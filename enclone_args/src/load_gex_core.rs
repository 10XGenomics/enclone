// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Load gene expression and feature barcoding (antibody, antigen) data from
// Cell Ranger outputs.

use enclone_core::defs::*;
use enclone_core::slurp::*;
use io_utils::*;
use itertools::Itertools;
use mirror_sparse_matrix::*;
use rayon::prelude::*;
use serde_json::Value;
use std::{
    collections::HashMap,
    convert::TryInto,
    fs::{read_to_string, remove_file, File},
    io::{BufRead, BufReader, Read},
    time::Instant,
};
use string_utils::*;
use vector_utils::*;

// parse_csv_pure: same as parse_csv, but don't strip out quotes

pub fn parse_csv_pure(x: &str) -> Vec<String> {
    let mut y = Vec::<String>::new();
    let mut w = Vec::<char>::new();
    for c in x.chars() {
        w.push(c);
    }
    let (mut quotes, mut i) = (0, 0);
    while i < w.len() {
        let mut j = i;
        while j < w.len() {
            if quotes % 2 == 0 && w[j] == ',' {
                break;
            }
            if w[j] == '"' {
                quotes += 1;
            }
            j += 1;
        }
        let (start, stop) = (i, j);
        let mut s = String::new();
        for m in start..stop {
            s.push(w[m]);
        }
        y.push(s);
        i = j + 1;
    }
    if !w.is_empty() && *w.last().unwrap() == ',' {
        y.push(String::new());
    }
    y
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn load_gex(
    ctl: &mut EncloneControl,
    gex_features: &mut Vec<Vec<String>>,
    gex_barcodes: &mut Vec<Vec<String>>,
    gex_matrices: &mut Vec<MirrorSparseMatrix>,
    fb_top_barcodes: &mut Vec<Vec<String>>,
    fb_top_matrices: &mut Vec<MirrorSparseMatrix>,
    fb_total_umis: &mut Vec<u64>,
    fb_brn: &mut Vec<Vec<(String, u32, u32)>>,
    feature_refs: &mut Vec<String>,
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
    feature_metrics: &mut Vec<HashMap<(String, String), String>>,
    json_metrics: &mut Vec<HashMap<String, f64>>,
    metrics: &mut Vec<String>,
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
        Vec<String>,
        HashMap<(String, String), String>,
        HashMap<String, f64>,
        String,
        u64,
        Vec<(String, u32, u32)>,
        String,
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
            Vec::<String>::new(),
            HashMap::<(String, String), String>::new(),
            HashMap::<String, f64>::new(),
            String::new(),
            0,
            Vec::new(),
            String::new(),
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
        let pathlist = &mut r.15;
        let i = r.0;
        if !gex_outs[i].is_empty() {
            // First define the path where the GEX files should live, and make sure that the path
            // exists.

            let root = gex_outs[i].clone();
            let mut outs = root.clone();
            if root.ends_with("/outs") && path_exists(&root) {
                outs = root;
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
                    pathlist.push(p.clone());
                    h5_path = p;
                    break;
                }
            }
            if h5_path.is_empty() {
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
            analysis.push(outs.to_string());
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
                    pathlist.push(pca_file.clone());
                    break;
                }
            }

            // Find the json metrics file.

            let mut json_metrics_file = String::new();
            if !ctl.gen_opt.cellranger {
                for x in analysis.iter() {
                    let f = format!("{}/metrics_summary_json.json", x);
                    if path_exists(&f) {
                        json_metrics_file = f.clone();
                        pathlist.push(f);
                        break;
                    }
                }
            }

            // Find the feature metrics file.

            let mut feature_metrics_file = String::new();
            if !ctl.gen_opt.cellranger {
                for x in analysis.iter() {
                    let f = format!("{}/per_feature_metrics.csv", x);
                    if path_exists(&f) {
                        feature_metrics_file = f.clone();
                        pathlist.push(f);
                        break;
                    }
                }
            }

            // Find the metrics file.

            let mut metrics_file = String::new();
            if !ctl.gen_opt.cellranger {
                let summary_dir = format!("{}/../multi_web_summary_json/metrics_summary_csv", outs);
                if path_exists(&summary_dir) {
                    let list = dir_list(&summary_dir);
                    if list.solo() {
                        let path = format!("{}/{}", summary_dir, list[0]);
                        pathlist.push(path.clone());
                        metrics_file = path;
                    }
                }
            }

            // Find the cluster file.

            let mut cluster_file = String::new();
            for x in analysis.iter() {
                cluster_file = format!("{}/clustering/graphclust/clusters.csv", x);
                if path_exists(&cluster_file) {
                    pathlist.push(cluster_file.clone());
                    break;
                }
            }

            // Proceed.

            let bin_file = format!("{}/feature_barcode_matrix.bin", outs);
            for f in [pca_file.clone(), cluster_file.clone()].iter() {
                if !path_exists(f) {
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
                } else {
                    pathlist.push(f.to_string());
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
                if path_exists(c) {
                    csv = c.clone();
                    pathlist.push(c.to_string());
                    break;
                }
            }
            if csv.is_empty() {
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
                    pathlist.push(bin_file.clone());
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
                if f.is_err() {
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
                pathlist.push(types_file.clone());
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
                || ctl.clono_filt_opt_def.marked_b
            {
                r.11 = format!(
                    "\nIf you use MARK_STATS or MARK_STATS2 or MARKED_B, celltypes.csv has to \
                    exist, and this file\n{}\ndoes not exist.\n",
                    types_file
                );
                return;
            }

            // Read json metrics file.  Note that we do not enforce the requirement of this
            // file, so it may not be present.  Also it is not present in the outs folder of CS
            // pipelines, and a customer would have to rerun with --vdrmode=disable to avoid
            // deleting the file, and then move it to outs so enclone could find it.

            if !json_metrics_file.is_empty() {
                let m = std::fs::read_to_string(&json_metrics_file).unwrap();
                let v: Value = serde_json::from_str(&m).unwrap();
                let z = v.as_object().unwrap();
                for (var, value) in z.iter() {
                    if value.as_f64().is_some() {
                        let value = value.as_f64().unwrap();
                        r.17.insert(var.to_string(), value);
                    }
                }
            }

            // Read and parse metrics file.  Rewrite as metrics class, metric name, metric value.

            if !metrics_file.is_empty() {
                let m = std::fs::read_to_string(&metrics_file).unwrap();
                let fields = parse_csv_pure(m.before("\n"));
                let (mut class, mut name, mut value) = (None, None, None);
                for i in 0..fields.len() {
                    if fields[i] == "Library Type" {
                        class = Some(i);
                    } else if fields[i] == "Metric Name" {
                        name = Some(i);
                    } else if fields[i] == "Metric Value" {
                        value = Some(i);
                    }
                }
                let (class, name, value) = (class.unwrap(), name.unwrap(), value.unwrap());
                let mut lines = Vec::<String>::new();
                let mut first = true;
                for line in m.lines() {
                    if first {
                        first = false;
                    } else {
                        let fields = parse_csv_pure(line);
                        lines.push(format!(
                            "{},{},{}",
                            fields[class], fields[name], fields[value]
                        ));
                    }
                }
                r.18 = format!("{}\n", lines.iter().format("\n"));
            }

            // Read feature metrics file.  Note that we do not enforce the requirement of this
            // file, so it may not be present.

            if !feature_metrics_file.is_empty() {
                let mut count = 0;
                let f = open_for_read![&feature_metrics_file];
                let mut feature_pos = HashMap::<String, usize>::new();
                let mut xfields = Vec::<String>::new();
                for line in f.lines() {
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    if count == 0 {
                        for j in 0..fields.len() {
                            feature_pos.insert(fields[j].to_string(), j);
                        }
                        xfields = fields.clone();
                    } else {
                        let feature_type = &fields[feature_pos["feature_type"]];
                        let mut feature;
                        for pass in 1..=2 {
                            if pass == 1 {
                                feature = fields[feature_pos["feature_name"]].clone();
                            } else {
                                feature = fields[feature_pos["feature_id"]].clone();
                            }
                            if feature_type.starts_with(&"Antibody") {
                                feature += "_ab";
                            } else if feature_type.starts_with(&"CRISPR") {
                                feature += "_cr";
                            } else if feature_type.starts_with(&"CUSTOM") {
                                feature += "_cu";
                            } else if feature_type.starts_with(&"Gene") {
                                feature += "_g";
                            }
                            for j in 0..fields.len() {
                                if xfields[j] == "num_umis"
                                    || xfields[j] == "num_reads"
                                    || xfields[j] == "num_umis_cells"
                                    || xfields[j] == "num_reads_cells"
                                {
                                    r.16.insert(
                                        (feature.clone(), xfields[j].clone()),
                                        fields[j].clone(),
                                    );
                                }
                            }
                        }
                    }
                    count += 1;
                }
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
                        if rpcx.parse::<usize>().is_err() {
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
                        if fbrpcx.parse::<usize>().is_err() {
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
                    let fields = parse_csv(s);
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
                            if rpcx.parse::<usize>().is_err() {
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
                            if fbrpcx.parse::<usize>().is_err() {
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

            let mut top_file = format!("{}/../feature_barcode_matrix_top.bin", outs);
            if !path_exists(&top_file) {
                top_file = format!("{}//feature_barcode_matrix_top.bin", outs);
            }
            if path_exists(&top_file) {
                pathlist.push(top_file.clone());
                read_from_file(&mut r.13, &top_file);
                let nrows = r.13.nrows();
                for i in 0..nrows {
                    r.14.push(r.13.row_label(i));
                }
            }

            // Read the total UMIs.

            let mut top_file = format!("{}/../feature_barcode_matrix_top.total", outs);
            if !path_exists(&top_file) {
                top_file = format!("{}/feature_barcode_matrix_top.total", outs);
            }
            if path_exists(&top_file) {
                pathlist.push(top_file.clone());
                let mut f = open_for_read![&top_file];
                let mut bytes = Vec::<u8>::new();
                f.read_to_end(&mut bytes).unwrap();
                r.19 = u64::from_ne_bytes(bytes.try_into().unwrap());
            }

            // Read the barcode-ref-nonref UMI count file.

            let mut brn_file = format!("{}/../feature_barcode_matrix_top.brn", outs);
            if !path_exists(&brn_file) {
                brn_file = format!("{}//feature_barcode_matrix_top.brn", outs);
            }
            if path_exists(&brn_file) {
                pathlist.push(brn_file.clone());
                let f = open_for_read![&brn_file];
                for line in f.lines() {
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    r.20.push((
                        fields[0].to_string(),
                        fields[1].parse::<u32>().unwrap(),
                        fields[2].parse::<u32>().unwrap(),
                    ));
                }
            }

            // Read the feature reference file.

            let mut fref_file = format!("{}/../feature_reference.csv", outs);
            if !path_exists(&fref_file) {
                fref_file = format!("{}/feature_reference.csv", outs);
            }
            if path_exists(&fref_file) {
                pathlist.push(fref_file.clone());
                r.21 = read_to_string(&fref_file).unwrap();
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
    for i in 0..results.len() {
        ctl.pathlist.append(&mut results[i].15.clone());
    }
    ctl.perf_stats(&t, "in load_gex main loop");

    // Test for error.

    let t = Instant::now();
    for i in 0..results.len() {
        if !results[i].11.is_empty() {
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

    // Add some metrics.

    let extras = [
        (
            "ANTIBODY_G_perfect_homopolymer_frac",
            "Antibody Capture,G Homopolymer Frac",
        ),
        (
            "GRCh38_raw_rpc_20000_subsampled_filtered_bcs_median_unique_genes_detected",
            "Gene Expression,GRCh38 Median genes per cell (20k raw reads per cell)",
        ),
        (
            "GRCh38_raw_rpc_20000_subsampled_filtered_bcs_median_counts",
            "Gene Expression,GRCh38 Median UMI counts per cell (20k raw reads per cell)",
        ),
    ];
    for x in extras.iter() {
        let metric_name = x.0.to_string();
        let metric_display_name = x.1.to_string();
        let mut have = false;
        for i in 0..results.len() {
            if results[i].17.contains_key(&metric_name) {
                have = true;
            }
        }
        if have {
            for i in 0..results.len() {
                let mut value = String::new();
                if results[i].17.contains_key(&metric_name) {
                    value = format!("{:.3}", results[i].17[&metric_name]);
                }
                results[i].18 += &mut format!("{},{}\n", metric_display_name, value);
            }
        }
    }

    // Save results.  This avoids cloning, which saves a lot of time.

    let n = results.len();
    for (
        _i,
        (
            _x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            _x11,
            _x12,
            x13,
            x14,
            _x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
        ),
    ) in results.into_iter().take(n).enumerate()
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
        feature_metrics.push(x16);
        json_metrics.push(x17);
        metrics.push(x18);
        fb_total_umis.push(x19);
        fb_brn.push(x20);
        feature_refs.push(x21);
    }

    // Done.

    ctl.perf_stats(&t, "in load_gex tail");
    Ok(())
}
