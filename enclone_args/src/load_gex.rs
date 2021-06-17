// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Load gene expression and feature barcoding (antibody, antigen) data from
// Cell Ranger outputs.

use enclone_core::defs::*;
use enclone_core::slurp::*;
use hdf5::Dataset;
use io_utils::*;
use mirror_sparse_matrix::*;
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fs::{remove_file, File},
    io::BufRead,
    time::Instant,
};
use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn load_gex(
    ctl: &mut EncloneControl,
    gex_features: &mut Vec<Vec<String>>,
    gex_barcodes: &mut Vec<Vec<String>>,
    gex_matrices: &mut Vec<MirrorSparseMatrix>,
    fb_top_barcodes: &mut Vec<Vec<String>>,
    fb_top_matrices: &mut Vec<MirrorSparseMatrix>,
    cluster: &mut Vec<HashMap<String, usize>>,
    cell_type: &mut Vec<HashMap<String, String>>,
    cell_type_specified: &mut Vec<bool>,
    pca: &mut Vec<HashMap<String, Vec<f64>>>,
    gex_mults: &mut Vec<f64>,
    fb_mults: &mut Vec<f64>,
    gex_cell_barcodes: &mut Vec<Vec<String>>,
    have_gex: &mut bool,
    have_fb: &mut bool,
    h5_paths: &mut Vec<String>,
) -> Result<(), String> {
    let t = Instant::now();
    let mut results = Vec::<(
        usize,
        Vec<String>,
        Vec<String>,
        MirrorSparseMatrix,
        Option<f64>,
        Option<f64>,
        Vec<String>,
        HashMap<String, usize>,
        HashMap<String, String>,
        HashMap<String, Vec<f64>>,
        bool,
        String,
        String,
        MirrorSparseMatrix,
        Vec<String>,
    )>::new();
    for i in 0..ctl.origin_info.gex_path.len() {
        results.push((
            i,
            Vec::<String>::new(),
            Vec::<String>::new(),
            MirrorSparseMatrix::new(),
            None,
            None,
            Vec::<String>::new(),
            HashMap::<String, usize>::new(),
            HashMap::<String, String>::new(),
            HashMap::<String, Vec<f64>>::new(),
            false,
            String::new(),
            String::new(),
            MirrorSparseMatrix::new(),
            Vec::<String>::new(),
        ));
    }
    let gex_outs = &ctl.origin_info.gex_path;
    // Here and in other places, where an error message can be printed in a parallel loop, it
    // would be better if the thread could use a global lock to prevent multiple threads from
    // issuing an error message.
    //
    // A lot of time is spent in this parallel loop.  Some things are known about this:
    // 1. When running it over a large number of datasets, the observed load average is ~2, so
    //    somehow the parallelism is not working.
    // 2. We know where the time is spent in the loop, and this is marked below.
    results.par_iter_mut().for_each(|r| {
        let i = r.0;
        if gex_outs[i].len() > 0 {
            // First define the path where the GEX files should live, and make sure that the path
            // exists.

            let root = gex_outs[i].clone();
            let mut outs = root.clone();
            if root.ends_with("/outs") && path_exists(&root) {
                outs = root.clone();
            } else if root.ends_with("/outs") {
                outs = root.before("/outs").to_string();
                if !path_exists(&outs) {
                    r.11 = format!(
                        "\nThe directory\n{}\ndoes not exist.  Something must be amiss with \
                        the arguments to PRE and/or GEX and/or META.\n",
                        outs
                    );
                    return;
                }
            }

            // Define the file paths and test for their existence.

            let mut h5_path = String::new();
            let h5p = [
                "raw_feature_bc_matrix.h5",
                "raw_gene_bc_matrices_h5.h5",
                "multi/count/raw_feature_bc_matrix.h5",
            ];
            for x in h5p.iter() {
                let p = format!("{}/{}", outs, x);
                if path_exists(&p) {
                    h5_path = p;
                    break;
                }
            }
            if h5_path.len() == 0 {
                r.11 = format!(
                    "\nThe file raw_feature_bc_matrix.h5 is not in the directory\n{}\n\
                    and neither is the older-named version raw_gene_bc_matrices_h5.h5.  Perhaps \
                    something\nis amiss with the arguments to PRE and/or GEX and/or META.\n",
                    outs
                );
                return;
            }
            r.12 = h5_path.clone();
            let types_file = format!("{}/analysis_csv/celltypes/celltypes.csv", outs);

            // Define possible places for the analysis directory.

            let mut analysis = Vec::<String>::new();
            analysis.push(format!("{}/analysis_csv", outs));
            analysis.push(format!("{}/analysis", outs));
            analysis.push(format!("{}/count/analysis", outs));
            let pso = format!("{}/per_sample_outs", outs);
            if path_exists(&pso) {
                let samples = dir_list(&pso);
                if samples.solo() {
                    let a = format!("{}/{}/count/analysis", pso, samples[0]);
                    analysis.push(a);
                }
            }
            let pso2 = format!("{}/../per_sample_outs", outs);
            if path_exists(&pso2) {
                let samples = dir_list(&pso2);
                if samples.solo() {
                    let a = format!("{}/{}/count/analysis", pso2, samples[0]);
                    analysis.push(a);
                }
            }

            // Find the pca file.

            let mut pca_file = String::new();
            for x in analysis.iter() {
                pca_file = format!("{}/pca/10_components/projection.csv", x);
                if path_exists(&pca_file) {
                    break;
                }
            }

            // Find the cluster file.

            let mut cluster_file = String::new();
            for x in analysis.iter() {
                cluster_file = format!("{}/clustering/graphclust/clusters.csv", x);
                if path_exists(&cluster_file) {
                    break;
                }
            }

            // Proceed.

            let bin_file = format!("{}/feature_barcode_matrix.bin", outs);
            for f in [pca_file.clone(), cluster_file.clone()].iter() {
                if !path_exists(&f) {
                    r.11 = format!(
                        "\nThe file\n{}\ndoes not exist.  \
                        Perhaps one of your directories is missing some stuff.\n\n\
                        One possibility is that you ran \"cellranger count\" using only \
                        feature barcode (antibody) data,\nand you had less then ten antibodies.  \
                        Currently if you do this, cellranger will not run the\nsecondary \
                        analyses, so you'll be missing some files.  A workaround is to add \
                        some \"fake\" antibodies\nto pad out the total number to ten.\n\n\
                        Another possibility is that this is a multi run, and the path you \
                        provided\nis to a subdirectory of the outs folder.  In that case it may \
                        work to provide the path to outs\nor (equivalently) the parent \
                        directory.\n",
                        f
                    );
                    return;
                }
            }

            // Find metrics summary file.

            let mut csv = String::new();
            let mut csvs = Vec::<String>::new();
            csvs.push(format!("{}/metrics_summary.csv", outs));
            csvs.push(format!("{}/metrics_summary_csv.csv", outs));
            let pso = format!("{}/per_sample_outs", outs);
            if path_exists(&pso) {
                let samples = dir_list(&pso);
                if samples.solo() {
                    let a = format!("{}/{}/metrics_summary.csv", pso, samples[0]);
                    csvs.push(a);
                    let a = format!("{}/{}/metrics_summary_csv.csv", pso, samples[0]);
                    csvs.push(a);
                }
            }
            for i in 0..csvs.len() {
                let c = &csvs[i];
                if path_exists(&c) {
                    csv = c.clone();
                    break;
                }
            }
            if csv.len() == 0 {
                r.11 = format!(
                    "\nSomething wrong with GEX or META argument:\ncan't find the file \
                        metrics_summary.csv or metrics_summary_csv.csv in the directory\n\
                        {}",
                    outs
                );
                return;
            }

            // Determine the state of affairs of the bin file.  We end with one of three outcomes:
            //
            // 1. We're not using the bin file at all.
            // 2. We are reading the bin file.
            // 3. We are writing the bin file.

            let mut bin_file_state = 1;
            if !ctl.gen_opt.force_h5 {
                let bin_file_exists = path_exists(&bin_file);
                if !bin_file_exists {
                    if !ctl.gen_opt.h5 {
                        bin_file_state = 3;
                    }
                } else {
                    // THE FOLLOWING LINE HAS BEEN OBSERVED TO FAIL SPORADICALLY.  THIS HAS
                    // HAPPENED MULTIPLE TIMES.  THE FAIL WAS IN
                    // binary_read_to_ref::<u32>(&mut ff, &mut x[0], 11).unwrap();
                    // WHERE THE unwrap() FAILED ON
                    // UnexpectedEof, error: "failed to fill whole buffer".
                    //
                    // 2/15/21: this should now be fixed.

                    let v = get_code_version_from_file(&bin_file);
                    if v == 1 {
                        bin_file_state = 2;
                    } else {
                        bin_file_state = 3;
                    }
                }
            }

            // If we need to write feature_barcode_matrix.bin, make sure that's possible, before
            // spending a lot of time reading other stuff.

            if bin_file_state == 3 {
                let f = File::create(&bin_file);
                if !f.is_ok() {
                    r.11 = format!(
                        "\nenclone is trying to create the path\n{}\n\
                        but that path cannot be created.  This path is for the binary GEX \
                        matrix file that enclone can read\n\
                        faster than the hdf5 file.  Your options are:\n\
                        1. Make that location writable (or fix the path, if it's wrong).\n\
                        2. Find a new location where you can write.\n\
                        3. Don't specify NH5 (if you specified it).\n",
                        bin_file
                    );
                    return;
                }
                remove_file(&bin_file).unwrap();
            }

            // Read cell types.

            if path_exists(&types_file) {
                let f = open_userfile_for_read(&types_file);
                let mut count = 0;
                for line in f.lines() {
                    count += 1;
                    if count == 1 {
                        continue;
                    }
                    let s = line.unwrap();
                    let barcode = s.before(",");
                    let cell_type = s.after(",");
                    r.8.insert(barcode.to_string(), cell_type.to_string());
                    r.10 = true;
                }
            } else if ctl.gen_opt.mark_stats
                || ctl.gen_opt.mark_stats2
                || ctl.clono_filt_opt.marked_b
            {
                r.11 = format!(
                    "\nIf you use MARK_STATS or MARK_STATS2 or MARKED_B, celltypes.csv has to \
                    exist, and this file\n{}\ndoes not exist.\n",
                    types_file
                );
                return;
            }

            // Read PCA file.

            let f = open_userfile_for_read(&pca_file);
            let mut count = 0;
            for line in f.lines() {
                count += 1;
                if count == 1 {
                    continue;
                }
                let s = line.unwrap();
                let barcode = s.before(",");
                let x = s.after(",").split(',').collect::<Vec<&str>>();
                // This assert is turned off because in fact there are not always 10 components.
                // assert_eq!(x.len(), 10);
                let mut y = Vec::<f64>::new();
                for i in 0..x.len() {
                    y.push(x[i].force_f64());
                }
                r.9.insert(barcode.to_string(), y);
            }

            // Read graph clusters, and also get the cell barcodes from that.

            let f = open_userfile_for_read(&cluster_file);
            let mut count = 0;
            for line in f.lines() {
                count += 1;
                if count == 1 {
                    continue;
                }
                let s = line.unwrap();
                let (barcode, cluster) = (s.before(","), s.after(",").force_usize());
                r.7.insert(barcode.to_string(), cluster);
                r.6.push(barcode.to_string());
            }

            // Get the multipliers gene and feature barcode counts.

            let mut gene_mult = None;
            let mut fb_mult = None;
            let mut rpc = None;
            let mut fbrpc = None;
            let mut lines = Vec::<String>::new();
            {
                let f = open_userfile_for_read(&csv);
                for line in f.lines() {
                    let s = line.unwrap();
                    lines.push(s.to_string());
                }
            }
            if lines.is_empty() {
                r.11 = format!("\nThe file\n{}\nis empty.\n", csv);
                return;
            }
            let fields = parse_csv(&lines[0]);
            if fields.contains(&"Metric Name".to_string())
                && fields.contains(&"Metric Value".to_string())
                && fields.contains(&"Library Type".to_string())
            {
                let mut lib_field = 0;
                let mut name_field = 0;
                let mut value_field = 0;
                for i in 0..fields.len() {
                    if fields[i] == "Library Type" {
                        lib_field = i;
                    } else if fields[i] == "Metric Name" {
                        name_field = i;
                    } else if fields[i] == "Metric Value" {
                        value_field = i;
                    }
                }
                for j in 1..lines.len() {
                    let fields = parse_csv(&lines[j]);
                    if fields.len() < lib_field + 1
                        || fields.len() < name_field + 1
                        || fields.len() < value_field + 1
                    {
                        r.11 = format!(
                            "\nSomething appears to be wrong with the file\n{}:\n\
                            line {} doesn't have enough fields.\n",
                            csv,
                            j + 1,
                        );
                        return;
                    }
                    if fields[lib_field] == "Gene Expression"
                        && fields[name_field] == "Mean reads per cell"
                    {
                        let mut rpcx = fields[value_field].to_string();
                        rpcx = rpcx.replace(",", "");
                        rpcx = rpcx.replace("\"", "");
                        if !rpcx.parse::<usize>().is_ok() {
                            r.11 = format!(
                                "\nSomething appears to be wrong with the file\n{}:\n\
                                the Gene Expression Mean Reads per Cell value isn't an integer.\n",
                                csv
                            );
                            return;
                        }
                        rpc = Some(rpcx.force_usize() as isize);
                    // Note that where we have "Antibody Capture", we could hypothetically have
                    // "CRISPR Guide Capture" or "Custom Feature".
                    } else if fields[lib_field] == "Antibody Capture"
                        && fields[name_field] == "Mean reads per cell"
                    {
                        let mut fbrpcx = fields[value_field].to_string();
                        fbrpcx = fbrpcx.replace(",", "");
                        fbrpcx = fbrpcx.replace("\"", "");
                        if !fbrpcx.parse::<usize>().is_ok() {
                            r.11 = format!(
                                "\nSomething appears to be wrong with the file\n{}:\n\
                                the Antibody Capture Mean Reads per Cell value isn't an integer.\n",
                                csv
                            );
                            return;
                        }
                        fbrpc = Some(fbrpcx.force_usize() as isize);
                    }
                }
                if rpc.is_none() && fbrpc.is_none() {
                    r.11 = format!(
                        "\nGene expression or feature barcode data was expected, however the \
                        CSV file\n{}\n\
                        does not have values for Gene Expression Mean Reads per Cell or
                        Antibody Capture Mean Reads per Cell.\n\
                        This is puzzling.\n",
                        csv,
                    );
                    return;
                }
            } else {
                let mut rpc_field = None;
                let mut fbrpc_field = None;
                for line_no in 0..lines.len() {
                    let s = &lines[line_no];
                    let fields = parse_csv(&s);
                    if line_no == 0 {
                        for i in 0..fields.len() {
                            if fields[i] == "Mean Reads per Cell" {
                                rpc_field = Some(i);
                            } else if fields[i] == "Antibody: Mean Reads per Cell" {
                                fbrpc_field = Some(i);
                            }
                        }
                    } else if line_no == 1 {
                        if rpc_field.is_some() && rpc_field.unwrap() >= fields.len() {
                            r.11 = format!(
                                "\nSomething appears to be wrong with the file\n{}:\n\
                                the second line doesn't have enough fields.\n",
                                csv
                            );
                            return;
                        } else if rpc_field.is_some() {
                            let mut rpcx = fields[rpc_field.unwrap()].to_string();
                            rpcx = rpcx.replace(",", "");
                            rpcx = rpcx.replace("\"", "");
                            if !rpcx.parse::<usize>().is_ok() {
                                r.11 = format!(
                                    "\nSomething appears to be wrong with the file\n{}:\n\
                                    the Mean Reads per Cell field isn't an integer.\n",
                                    csv
                                );
                                return;
                            }
                            rpc = Some(rpcx.force_usize() as isize);
                        }
                        if fbrpc_field.is_some() && fbrpc_field.unwrap() >= fields.len() {
                            r.11 = format!(
                                "\nSomething appears to be wrong with the file\n{}:\n\
                                the second line doesn't have enough fields.\n",
                                csv
                            );
                            return;
                        } else if fbrpc_field.is_some() {
                            let mut fbrpcx = fields[fbrpc_field.unwrap()].to_string();
                            fbrpcx = fbrpcx.replace(",", "");
                            fbrpcx = fbrpcx.replace("\"", "");
                            if !fbrpcx.parse::<usize>().is_ok() {
                                r.11 = format!(
                                    "\nSomething appears to be wrong with the file\n{}:\n\
                                    the Antibody: Mean Reads per Cell field isn't an integer.\n",
                                    csv
                                );
                                return;
                            }
                            fbrpc = Some(fbrpcx.force_usize() as isize);
                        }
                    }
                }
                if rpc.is_none() && fbrpc.is_none() {
                    r.11 = format!(
                        "\nGene expression or feature barcode data was expected, however the \
                        CSV file\n{}\n\
                        does not have a field \"Mean Reads per Cell\" or \
                        \"Antibody: Mean Reads per Cell\".\n\
                        This is puzzling, and might be because a file within the Cell Ranger outs \
                        directory has been moved\n\
                        from its original location.\n",
                        csv,
                    );
                    return;
                }
            }
            if rpc.is_some() {
                const RPC_EXPECTED: f64 = 20_000.0;
                gene_mult = Some(RPC_EXPECTED / rpc.unwrap() as f64);
            }
            if fbrpc.is_some() {
                const FB_RPC_EXPECTED: f64 = 5_000.0;
                fb_mult = Some(FB_RPC_EXPECTED / fbrpc.unwrap() as f64);
            }
            r.4 = gene_mult;
            r.5 = fb_mult;

            // Read the top feature barcode matrix.

            let top_file = format!("{}/feature_barcode_matrix_top.bin", outs);
            if path_exists(&top_file) {
                read_from_file(&mut r.13, &top_file);
                let nrows = r.13.nrows();
                for i in 0..nrows {
                    r.14.push(r.13.row_label(i));
                }
            }

            // Read the binary matrix file if appropriate.

            if bin_file_state == 2 {
                read_from_file(&mut r.3, &bin_file);
                let (n, k) = (r.3.nrows(), r.3.ncols());
                for i in 0..n {
                    r.2.push(r.3.row_label(i));
                }
                for j in 0..k {
                    r.1.push(r.3.col_label(j));
                }

            // Otherwise we have to get stuff from the h5 file.
            } else {
                let mut matrix = Vec::<Vec<(i32, i32)>>::new();
                slurp_h5(
                    &h5_path,
                    bin_file_state == 3,
                    &mut r.2,
                    &mut r.1,
                    &mut matrix,
                );
                if bin_file_state == 3 {
                    r.3 = MirrorSparseMatrix::build_from_vec(&matrix, &r.2, &r.1);
                    write_to_file(&r.3, &bin_file);
                    // Note that if the dataset archive was complete, we would not need to do this.
                    if ctl.gen_opt.internal_run {
                        let earth = &ctl.gen_opt.config["earth"];
                        if !bin_file.starts_with(earth) {
                            let bin_file_alt =
                                format!("{}/current{}", earth, bin_file.after("current"));
                            write_to_file(&r.3, &bin_file_alt);
                        }
                    }
                }
            }
        }
        unique_sort(&mut r.6);
    });
    ctl.perf_stats(&t, "in load_gex main loop");

    // Test for error.

    let t = Instant::now();
    for i in 0..results.len() {
        if results[i].11.len() > 0 {
            return Err(results[i].11.clone());
        }
    }

    // Set have_gex and have_fb.

    for i in 0..results.len() {
        if results[i].4.is_some() {
            *have_gex = true;
        }
        if results[i].5.is_some() {
            *have_fb = true;
        }
    }
    for i in 0..results.len() {
        h5_paths.push(results[i].12.clone());
    }

    // Save results.  This avoids cloning, which saves a lot of time.

    let n = results.len();
    for (_i, (_x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, _x11, _x12, x13, x14)) in
        results.into_iter().take(n).enumerate()
    {
        gex_features.push(x1);
        gex_barcodes.push(x2);
        gex_matrices.push(x3);
        fb_top_matrices.push(x13);
        fb_top_barcodes.push(x14);
        let mut gex_mult = 1.0;
        if x4.is_some() {
            gex_mult = x4.unwrap();
        }
        gex_mults.push(gex_mult);
        let mut fb_mult = 1.0;
        if x5.is_some() {
            fb_mult = x5.unwrap();
        }
        fb_mults.push(fb_mult);
        gex_cell_barcodes.push(x6);
        cluster.push(x7);
        cell_type.push(x8);
        pca.push(x9);
        cell_type_specified.push(x10);
    }
    ctl.perf_stats(&t, "in load_gex tail");
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Get gene expression and feature barcoding counts.

pub fn get_gex_info(mut ctl: &mut EncloneControl) -> Result<GexInfo, String> {
    let mut gex_features = Vec::<Vec<String>>::new();
    let mut gex_barcodes = Vec::<Vec<String>>::new();
    let mut gex_matrices = Vec::<MirrorSparseMatrix>::new();
    let mut fb_top_barcodes = Vec::<Vec<String>>::new();
    let mut fb_top_matrices = Vec::<MirrorSparseMatrix>::new();
    let mut cluster = Vec::<HashMap<String, usize>>::new();
    let mut cell_type = Vec::<HashMap<String, String>>::new();
    let mut cell_type_specified = Vec::<bool>::new();
    let mut pca = Vec::<HashMap<String, Vec<f64>>>::new();
    let mut gex_mults = Vec::<f64>::new();
    let mut fb_mults = Vec::<f64>::new();
    let mut gex_cell_barcodes = Vec::<Vec<String>>::new();
    let mut have_gex = false;
    let mut have_fb = false;
    let mut h5_paths = Vec::<String>::new();
    load_gex(
        &mut ctl,
        &mut gex_features,
        &mut gex_barcodes,
        &mut gex_matrices,
        &mut fb_top_barcodes,
        &mut fb_top_matrices,
        &mut cluster,
        &mut cell_type,
        &mut cell_type_specified,
        &mut pca,
        &mut gex_mults,
        &mut fb_mults,
        &mut gex_cell_barcodes,
        &mut have_gex,
        &mut have_fb,
        &mut h5_paths,
    )?;
    let t = Instant::now();
    if ctl.gen_opt.gene_scan_test.is_some() && !ctl.gen_opt.accept_inconsistent {
        let mut allf = gex_features.clone();
        unique_sort(&mut allf);
        if allf.len() != 1 {
            let mut msg = format!(
                "\nCurrently, SCAN requires that all datasets have identical \
                 features, and they do not.\n\
                There are {} datasets and {} feature sets after removal of \
                 duplicates.\n",
                gex_features.len(),
                allf.len()
            );
            msg += &mut format!("Classification of features sets:\n\n");
            for i in 0..gex_features.len() {
                let p = bin_position(&allf, &gex_features[i]);
                msg += &mut format!("{} ==> {}\n", ctl.origin_info.dataset_id[i], p);
            }
            msg += "\n";
            return Err(msg);
        }
    }
    let mut h5_data = Vec::<Option<Dataset>>::new();
    let mut h5_indices = Vec::<Option<Dataset>>::new();
    let mut h5_indptr = Vec::<Vec<u32>>::new();
    if ctl.gen_opt.h5 {
        let gex_outs = &ctl.origin_info.gex_path;
        for i in 0..ctl.origin_info.dataset_path.len() {
            // let bin_file = format!("{}/feature_barcode_matrix.bin", gex_outs[i]);
            if gex_outs[i].len() > 0
            /* && !(path_exists(&bin_file) && !ctl.gen_opt.force_h5) */
            {
                let f = &h5_paths[i];
                let h = hdf5::File::open(&f).unwrap();
                h5_data.push(Some(h.dataset("matrix/data").unwrap()));
                h5_indices.push(Some(h.dataset("matrix/indices").unwrap()));
                let indptr = h.dataset("matrix/indptr").unwrap();
                let x: Vec<u32> = indptr.as_reader().read().unwrap().to_vec();
                h5_indptr.push(x);
            } else {
                h5_data.push(None);
                h5_indices.push(None);
                h5_indptr.push(Vec::<u32>::new());
            }
        }
    }
    fn compute_feature_id(gex_features: &Vec<String>) -> HashMap<String, usize> {
        let mut x = HashMap::<String, usize>::new();
        for j in 0..gex_features.len() {
            let f = &gex_features[j];
            let ff = f.split('\t').collect::<Vec<&str>>();
            for z in 0..2 {
                if ff[2].starts_with(&"Antibody") {
                    x.insert(format!("{}_ab", ff[z]), j);
                } else if ff[2].starts_with(&"CRISPR") {
                    x.insert(format!("{}_cr", ff[z]), j);
                } else if ff[2].starts_with(&"CUSTOM") {
                    x.insert(format!("{}_cu", ff[z]), j);
                } else if ff[2].starts_with(&"Gene") {
                    x.insert(format!("{}_g", ff[z]), j);
                }
            }
        }
        x
    }
    let n = gex_features.len();
    let pi = (0..n).into_par_iter();
    let mut feature_id = Vec::<HashMap<String, usize>>::new();
    pi.map(|i| compute_feature_id(&gex_features[i]))
        .collect_into_vec(&mut feature_id);
    let mut is_gex = Vec::<Vec<bool>>::new();
    for i in 0..gex_features.len() {
        is_gex.push(vec![false; gex_features[i].len()]);
        for j in 0..gex_features[i].len() {
            let f = &gex_features[i][j];
            let ff = f.split('\t').collect::<Vec<&str>>();
            if ff[2].starts_with(&"Gene") {
                is_gex[i][j] = true;
            }
        }
    }
    ctl.perf_stats(&t, "after load_gex");

    // Answer.

    Ok(GexInfo {
        gex_features: gex_features,
        gex_barcodes: gex_barcodes,
        gex_matrices: gex_matrices,
        fb_top_barcodes: fb_top_barcodes,
        fb_top_matrices: fb_top_matrices,
        cluster: cluster,
        cell_type: cell_type,
        cell_type_specified: cell_type_specified,
        pca: pca,
        gex_cell_barcodes: gex_cell_barcodes,
        gex_mults: gex_mults,
        fb_mults: fb_mults,
        h5_data: h5_data,
        h5_indices: h5_indices,
        h5_indptr: h5_indptr,
        is_gex: is_gex,
        feature_id: feature_id,
        have_gex: have_gex,
        have_fb: have_fb,
    })
}
