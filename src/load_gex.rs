// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Load gene expression and antibody data.

use crate::defs::*;
use h5::Dataset;
use io_utils::*;
use load_feature_bc::*;
use mirror_sparse_matrix::*;
use perf_stats::*;
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fs::{copy, create_dir_all, remove_file, File},
    io::{BufRead, BufReader},
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
    gex_mults: &mut Vec<f64>,
    fb_mults: &mut Vec<f64>,
    gex_cell_barcodes: &mut Vec<Vec<String>>,
    have_gex: &mut bool,
    have_fb: &mut bool,
) {
    let pre = &ctl.gen_opt.pre;
    let comp = ctl.comp;
    let mut results = Vec::<(
        usize,
        Vec<String>,
        Vec<String>,
        MirrorSparseMatrix,
        Option<f64>,
        Option<f64>,
        Vec<String>,
    )>::new();
    for i in 0..ctl.sample_info.gex_path.len() {
        results.push((
            i,
            Vec::<String>::new(),
            Vec::<String>::new(),
            MirrorSparseMatrix::new(),
            None,
            None,
            Vec::<String>::new(),
        ));
    }
    let gex_outs = &ctl.sample_info.gex_path;
    results.par_iter_mut().for_each(|r| {
        let i = r.0;
        if gex_outs[i].len() > 0 {
            let mut root = gex_outs[i].clone();
            if root.ends_with("/outs") {
                root = root.rev_before("/outs").to_string();
            }
            if !path_exists(&root) {
                eprintln!(
                    "\nProbably something is wrong with the GEX argument:\nthe path{} \
                     does not exist.\n",
                    root
                );
                std::process::exit(1);
            }
            if !path_exists(&gex_outs[i]) {
                eprintln!(
                    "\nProbably something is wrong with the GEX argument:\nthe path{} \
                     does not exist.\n",
                    gex_outs[i]
                );
                std::process::exit(1);
            }
            let cb1 = format!("{}/filtered_feature_bc_matrix/barcodes.tsv.gz", gex_outs[i]);
            let cb2 = format!(
                "{}/filtered_gene_bc_matrices_mex/barcodes.tsv.gz",
                gex_outs[i]
            );
            let mut cb = cb1.clone();
            if !path_exists(&cb) {
                cb = cb2.clone();
                if !path_exists(&cb) {
                    eprintln!(
                        "\nSomething wrong with GEX argument:\nneither the path {}\nnor \
                         the path {} exists.\n",
                        cb1, cb2
                    );
                    std::process::exit(1);
                }
            }
            read_maybe_unzipped(&cb, &mut r.6);
            let csv1 = format!("{}/metrics_summary.csv", gex_outs[i]);
            let csv2 = format!("{}/metrics_summary_csv.csv", gex_outs[i]);
            if !path_exists(&csv1) && !path_exists(&csv2) {
                eprintln!(
                    "\nSomething wrong with GEX or META argument:\ncan't find the file \
                        metrics_summary.csv or metrics_summary_csv.csv in the directory\n\
                        {}",
                    gex_outs[i]
                );
                std::process::exit(1);
            }
            let tm = Instant::now();
            let mut gene_mult = None;

            let mut csv = csv1.clone();
            if !path_exists(&csv1) {
                csv = csv2.clone();
            }
            let f = open_for_read![&csv];
            let mut line_no = 0;
            let mut rpc_field = None;
            let mut rpc = None;
            let mut fbrpc_field = None;
            let mut fbrpc = None;
            for line in f.lines() {
                let s = line.unwrap();
                let fields = parse_csv(&s);
                line_no += 1;
                if line_no == 1 {
                    for i in 0..fields.len() {
                        if fields[i] == "Mean Reads per Cell" {
                            rpc_field = Some(i);
                        } else if fields[i] == "Antibody: Mean Reads per Cell" {
                            fbrpc_field = Some(i);
                        }
                    }
                } else if line_no == 2 {
                    if rpc_field.is_some() && rpc_field.unwrap() >= fields.len() {
                        eprintln!(
                            "\nSomething appears to be wrong with the file\n{}:\n\
                            the second line doesn't have enough fields.\n",
                            csv
                        );
                        std::process::exit(1);
                    } else if rpc_field.is_some() {
                        let mut rpcx = fields[rpc_field.unwrap()].to_string();
                        rpcx = rpcx.replace(",", "");
                        rpcx = rpcx.replace("\"", "");
                        if !rpcx.parse::<usize>().is_ok() {
                            eprintln!(
                                "\nSomething appears to be wrong with the file\n{}:\n\
                                the Mean Reads per Cell field isn't an integer.\n",
                                csv
                            );
                            std::process::exit(1);
                        }
                        rpc = Some(rpcx.force_usize() as isize);
                    }
                    if fbrpc_field.is_some() && fbrpc_field.unwrap() >= fields.len() {
                        eprintln!(
                            "\nSomething appears to be wrong with the file\n{}:\n\
                            the second line doesn't have enough fields.\n",
                            csv
                        );
                        std::process::exit(1);
                    } else if fbrpc_field.is_some() {
                        let mut fbrpcx = fields[fbrpc_field.unwrap()].to_string();
                        fbrpcx = fbrpcx.replace(",", "");
                        fbrpcx = fbrpcx.replace("\"", "");
                        if !fbrpcx.parse::<usize>().is_ok() {
                            eprintln!(
                                "\nSomething appears to be wrong with the file\n{}:\n\
                                the Antibody: Mean Reads per Cell field isn't an integer.\n",
                                csv
                            );
                            std::process::exit(1);
                        }
                        fbrpc = Some(fbrpcx.force_usize() as isize);
                    }
                }
            }
            if comp {
                println!(
                    "-- used {:.2} seconds getting reads per cell metric",
                    elapsed(&tm)
                );
            }
            if rpc.is_some() {
                const RPC_EXPECTED: f64 = 20_000.0;
                gene_mult = Some(RPC_EXPECTED / rpc.unwrap() as f64);
            }
            let mut fb_mult = None;
            if fbrpc.is_some() {
                const FB_RPC_EXPECTED: f64 = 5_000.0;
                fb_mult = Some(FB_RPC_EXPECTED / fbrpc.unwrap() as f64);
            }
            r.4 = gene_mult;
            r.5 = fb_mult;
            let (mut features, mut barcodes) = (Vec::<String>::new(), Vec::<String>::new());
            let bin_file = format!("{}/raw_feature_bc_matrix/matrix.bin", gex_outs[i]);
            let mut dir = format!("{}/raw_feature_bc_matrix", gex_outs[i]);
            if !path_exists(&dir) {
                dir = format!("{}/raw_gene_bc_matrices_mex", gex_outs[i]);
            }
            let tfb = Instant::now();
            if path_exists(&format!("{}/GRCh38", dir)) {
                read_maybe_unzipped(&format!("{}/GRCh38/barcodes.tsv.gz", dir), &mut barcodes);
            } else {
                read_maybe_unzipped(&format!("{}/barcodes.tsv.gz", dir), &mut barcodes);
            }
            if comp {
                println!(
                    "-- used {:.2} seconds reading features and barcodes",
                    elapsed(&tfb)
                );
            }
            let tfb = Instant::now();
            if path_exists(&format!("{}/GRCh38", dir)) {
                read_maybe_unzipped(&format!("{}/GRCh38/genes.tsv.gz", dir), &mut r.1);
            } else {
                read_maybe_unzipped(&format!("{}/features.tsv.gz", dir), &mut r.1);
            }
            if path_exists(&format!("{}/GRCh38", dir)) {
                read_maybe_unzipped(&format!("{}/GRCh38/barcodes.tsv.gz", dir), &mut r.2);
            } else {
                read_maybe_unzipped(&format!("{}/barcodes.tsv.gz", dir), &mut r.2);
            }
            if comp {
                println!(
                    "-- used {:.2} seconds reading features and barcodes",
                    elapsed(&tfb)
                );
            }
            if path_exists(&bin_file) && !ctl.gen_opt.force_h5 {
                let t = Instant::now();
                read_from_file(&mut r.3, &bin_file);
                if comp {
                    println!("-- used {:.2} seconds reading matrix.bin", elapsed(&t));
                }
            } else if !ctl.gen_opt.h5 {
                let mut matrix = Vec::<Vec<(i32, i32)>>::new();
                let mut dir = format!("{}/raw_feature_bc_matrix", gex_outs[i]);
                if !path_exists(&dir) {
                    dir = format!("{}/raw_gene_bc_matrices_mex", gex_outs[i]);
                }
                let mut matrix_file = format!("{}/matrix.mtx.gz", dir);
                if path_exists(&format!("{}/GRCh38", dir)) {
                    matrix_file = format!("{}/GRCh38/matrix.mtx.gz", dir);
                }
                if !path_exists(&matrix_file) {
                    eprintln!(
                        "\nYou've used the NH5 option, but a gene expression directory \
                         is incomplete for this purpose,\nbecause this file\n\
                         {}\ndoes not exist.\n",
                        matrix_file
                    );
                    std::process::exit(1);
                }
                load_feature_bc_matrix(&gex_outs[i], &mut features, &mut barcodes, &mut matrix);
                r.3 = MirrorSparseMatrix::build_from_vec(&matrix);
                if *pre != "".to_string() || !ctl.gen_opt.internal_run {
                    let mut go = ctl.sample_info.gex_path[i].clone();
                    if go.ends_with("/HEAD/outs") {
                        let id = go.rev_before("/HEAD/outs").rev_after("/");
                        go = format!("{}/{}/outs", pre, id);
                    }
                    let dir_new = format!("{}/raw_feature_bc_matrix", go);
                    let bin_file = format!("{}/raw_feature_bc_matrix/matrix.bin", go);
                    if !path_exists(&dir_new) {
                        create_dir_all(&dir_new).unwrap();
                    }
                    if path_exists(&bin_file) {
                        remove_file(&bin_file).unwrap();
                    }
                    write_to_file(&r.3, &bin_file);
                    if go != gex_outs[i] {
                        let old_json = format!("{}/metrics_summary_json.json", gex_outs[i]);
                        let new_json = format!("{}/metrics_summary_json.json", go);
                        copy(&old_json, &new_json).unwrap();
                        if path_exists(&format!("{}/GRCh38", dir)) {
                            copy(
                                &format!("{}/GRCh38/genes.tsv.gz", dir),
                                &format!("{}/GRCh38/genes.tsv.gz", dir_new),
                            )
                            .unwrap();
                        } else {
                            copy(
                                &format!("{}/features.tsv.gz", dir),
                                &format!("{}/features.tsv.gz", dir_new),
                            )
                            .unwrap();
                        }
                        if path_exists(&format!("{}/GRCh38", dir)) {
                            copy(
                                &format!("{}/GRCh38/barcodes.tsv.gz", dir),
                                &format!("{}/GRCh38/barcodes.tsv.gz", dir_new),
                            )
                            .unwrap();
                        } else {
                            copy(
                                &format!("{}/barcodes.tsv.gz", dir),
                                &format!("{}/barcodes.tsv.gz", dir_new),
                            )
                            .unwrap();
                        }
                        let mut fdir = format!("{}/filtered_feature_bc_matrix", gex_outs[i]);
                        if !path_exists(&fdir) {
                            fdir = format!("{}/filtered_gene_bc_matrices_mex", gex_outs[i]);
                        }
                        let fdir_new = format!("{}/filtered_feature_bc_matrix", go);
                        if !path_exists(&fdir_new) {
                            create_dir_all(&fdir_new).unwrap();
                        }
                        if path_exists(&format!("{}/GRCh38", fdir)) {
                            copy(
                                &format!("{}/GRCh38/barcodes.tsv.gz", fdir),
                                &format!("{}/GRCh38/barcodes.tsv.gz", fdir_new),
                            )
                            .unwrap();
                        } else {
                            copy(
                                &format!("{}/barcodes.tsv.gz", fdir),
                                &format!("{}/barcodes.tsv.gz", fdir_new),
                            )
                            .unwrap();
                        }
                    }
                }
            }
        }
        unique_sort(&mut r.6);
    });

    // Set have_gex and have_fb.

    for i in 0..results.len() {
        if results[i].4.is_some() {
            *have_gex = true;
        }
        if results[i].5.is_some() {
            *have_fb = true;
        }
    }

    // Save results.

    for i in 0..results.len() {
        // All this cloning seems pointless and inefficient.
        gex_features.push(results[i].1.clone());
        gex_barcodes.push(results[i].2.clone());
        gex_matrices.push(results[i].3.clone());
        let mut gex_mult = 1.0;
        if results[i].4.is_some() {
            gex_mult = results[i].4.unwrap();
        }
        gex_mults.push(gex_mult);
        let mut fb_mult = 1.0;
        if results[i].5.is_some() {
            fb_mult = results[i].5.unwrap();
        }
        fb_mults.push(fb_mult);
        gex_cell_barcodes.push(results[i].6.clone());
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Get gene expression and antibody counts.

pub fn get_gex_info(mut ctl: &mut EncloneControl) -> GexInfo {
    let mut gex_features = Vec::<Vec<String>>::new();
    let mut gex_barcodes = Vec::<Vec<String>>::new();
    let mut gex_matrices = Vec::<MirrorSparseMatrix>::new();
    let mut gex_mults = Vec::<f64>::new();
    let mut fb_mults = Vec::<f64>::new();
    let mut gex_cell_barcodes = Vec::<Vec<String>>::new();
    let mut have_gex = false;
    let mut have_fb = false;
    load_gex(
        &mut ctl,
        &mut gex_features,
        &mut gex_barcodes,
        &mut gex_matrices,
        &mut gex_mults,
        &mut fb_mults,
        &mut gex_cell_barcodes,
        &mut have_gex,
        &mut have_fb,
    );
    if ctl.gen_opt.gene_scan_test.is_some() && !ctl.gen_opt.accept_inconsistent {
        let mut allf = gex_features.clone();
        unique_sort(&mut allf);
        if allf.len() != 1 {
            eprintln!(
                "\nCurrently, SCAN requires that all datasets have identical \
                 features, and they do not."
            );
            eprintln!(
                "There are {} datasets and {} feature sets after removal of \
                 duplicates.",
                gex_features.len(),
                allf.len()
            );
            eprintln!("Classification of features sets:\n");
            for i in 0..gex_features.len() {
                let p = bin_position(&allf, &gex_features[i]);
                eprintln!("{} ==> {}", ctl.sample_info.dataset_id[i], p);
            }
            eprintln!("");
            std::process::exit(1);
        }
    }
    let mut h5_data = Vec::<Option<Dataset>>::new();
    let mut h5_indices = Vec::<Option<Dataset>>::new();
    let mut h5_indptr = Vec::<Vec<u32>>::new();
    if ctl.gen_opt.h5 {
        let gex_outs = &ctl.sample_info.gex_path;
        for i in 0..ctl.sample_info.dataset_path.len() {
            if gex_outs[i].len() > 0 {
                let mut f = format!("{}/raw_feature_bc_matrix.h5", gex_outs[i]);
                if !path_exists(&f) {
                    f = format!("{}/raw_gene_bc_matrices_h5.h5", gex_outs[i]);
                }
                if !path_exists(&f) {
                    eprintln!("\nThere's a missing input file:\n{}.\n", f);
                    std::process::exit(1);
                }
                let h = h5::File::open(&f, "r").unwrap();
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
    let mut feature_id = vec![HashMap::<String, usize>::new(); gex_features.len()];
    for i in 0..gex_features.len() {
        for j in 0..gex_features[i].len() {
            let f = &gex_features[i][j];
            let ff = f.split('\t').collect::<Vec<&str>>();
            for z in 0..2 {
                if ff[2].starts_with(&"Antibody") {
                    feature_id[i].insert(format!("{}_ab", ff[z]), j);
                } else if ff[2].starts_with(&"Antigen") {
                    feature_id[i].insert(format!("{}_ag", ff[z]), j);
                } else if ff[2].starts_with(&"CRISPR") {
                    feature_id[i].insert(format!("{}_cr", ff[z]), j);
                } else if ff[2].starts_with(&"CUSTOM") {
                    feature_id[i].insert(format!("{}_cu", ff[z]), j);
                } else if ff[2].starts_with(&"Gene") {
                    feature_id[i].insert(format!("{}_g", ff[z]), j);
                }
            }
        }
    }
    let mut is_gex = Vec::<Vec<bool>>::new();
    for i in 0..gex_features.len() {
        is_gex.push(vec![false; gex_features[i].len()]);
        for j in 0..gex_features[i].len() {
            let f = &gex_features[i][j];
            let ff = f.split('\t').collect::<Vec<&str>>();
            if !ff[2].starts_with(&"Antibody") {
                is_gex[i][j] = true;
            }
        }
    }

    // Answer.

    GexInfo {
        gex_features: gex_features,
        gex_barcodes: gex_barcodes,
        gex_matrices: gex_matrices,
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
    }
}
