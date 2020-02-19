// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Load gene expression and antibody data.

use crate::defs::*;
use crate::proc_args2::*;
use h5::Dataset;
use io_utils::*;
use load_feature_bc::*;
use mirror_sparse_matrix::*;
use perf_stats::*;
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fs::{copy, create_dir_all, remove_file},
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
            let root = gex_outs[i].rev_before("/outs");
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
            let json = format!("{}/metrics_summary_json.json", gex_outs[i]);
            if !path_exists(&json) {
                eprintln!(
                    "\nSomething wrong with GEX argument:\nthe path {} does not exist.\n",
                    json
                );
                std::process::exit(1);
            }
            let tm = Instant::now();
            let mut gene_mult = None;
            let rpc = get_metric_value(&json, &"reads_per_cell".to_string());
            if comp {
                println!(
                    "-- used {:.2} seconds getting reads per cell metric",
                    elapsed(&tm)
                );
            }
            if rpc.len() > 0 {
                if rpc.parse::<f64>().is_err() {
                    eprintln!(
                        "\nIn load_gex, failed to parse metric reads_per_cell from {}.\n",
                        json
                    );
                    eprintln!(
                        "This is unexpected behavior, if you can't figure out what \
                         happened, please let us know.\n"
                    );
                    std::process::exit(1);
                }
                let rpc = rpc.parse::<f64>().unwrap();
                const RPC_EXPECTED: f64 = 20_000.0;
                gene_mult = Some(RPC_EXPECTED / rpc);
            }
            let mut fb_mult = None;
            let fbrpc = get_metric_value(&json, &"ANTIBODY_reads_per_cell".to_string());
            if fbrpc.len() > 0 {
                if fbrpc.parse::<f64>().is_err() {
                    eprintln!(
                        "\nIn load_gex, failed to parse metric ANTIBODY_reads_per_cell from {}.\n",
                        json
                    );
                    eprintln!(
                        "This is unexpected behavior, if you can't figure out what \
                         happened, please let us know.\n"
                    );
                    std::process::exit(1);
                }
                let fbrpc = fbrpc.parse::<f64>().unwrap();
                const FB_RPC_EXPECTED: f64 = 5_000.0;
                fb_mult = Some(FB_RPC_EXPECTED / fbrpc);
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
            if !ctl.gen_opt.h5 {
                if path_exists(&bin_file) {
                    let t = Instant::now();
                    read_from_file(&mut r.3, &bin_file);
                    if comp {
                        println!("-- used {:.2} seconds reading matrix.bin", elapsed(&t));
                    }
                } else {
                    let mut matrix = Vec::<Vec<(i32, i32)>>::new();
                    load_feature_bc_matrix(&gex_outs[i], &mut features, &mut barcodes, &mut matrix);
                    r.3 = MirrorSparseMatrix::build_from_vec(&matrix);
                    if *pre != "".to_string() {
                        let go = ctl.sample_info.gex_path[i].clone();
                        let dir_new = format!("{}/raw_feature_bc_matrix", go);
                        let bin_file = format!("{}/raw_feature_bc_matrix/matrix.bin", go);
                        if !path_exists(&dir_new) {
                            create_dir_all(&dir_new).unwrap();
                        }
                        if path_exists(&bin_file) {
                            remove_file(&bin_file).unwrap();
                        }
                        write_to_file(&r.3, &bin_file);
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
                    }
                }
            }
        }
        unique_sort(&mut r.6);
    });

    // Check for lvars that don't make sense.

    let mut have_gex = false;
    let mut have_fb = false;
    for i in 0..results.len() {
        if results[i].4.is_some() {
            have_gex = true;
        }
        if results[i].5.is_some() {
            have_fb = true;
        }
    }
    for x in ctl.clono_print_opt.lvars.iter() {
        if *x == "gex_med".to_string()
            || *x == "gex_max".to_string()
            || x.ends_with("_g")
            || *x == "n_gex".to_string()
            || *x == "entropy".to_string()
        {
            if !have_gex {
                eprintln!(
                    "\nYou've supplied the lead column variable {},\nbut it would appear \
                     that you do not have gene expression data.\n",
                    *x
                );
                std::process::exit(1);
            }
        }
        if x.ends_with("_ab") || x.ends_with("_ag") || x.ends_with("_cr") || x.ends_with("_cu") {
            if !have_fb {
                eprintln!(
                    "\nYou've supplied the lead column variable {},\nbut it would appear \
                     that you do not have feature barcode data.\n",
                    *x
                );
                std::process::exit(1);
            }
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
    load_gex(
        &mut ctl,
        &mut gex_features,
        &mut gex_barcodes,
        &mut gex_matrices,
        &mut gex_mults,
        &mut fb_mults,
        &mut gex_cell_barcodes,
    );
    /*
    if ctl.gen_opt.gene_scan_test.is_some() {
        let mut allf = gex_features.clone();
        unique_sort(&mut allf);
        if allf.len() != 1 {
            eprintln!( "\nCurrently, SCAN requires that all datasets have identical \
                features, and they do not." );
            eprintln!( "There are {} datasets and {} feature sets after removal of \
                duplicates.", gex_features.len(), allf.len() );
            eprintln!( "Classification of features sets:\n" );
            for i in 0..gex_features.len() {
                let p = bin_position(&allf, &gex_features[i]);
                eprintln!( "{} ==> {}", ctl.sample_info.dataset_id[i], p );
            }
            eprintln!("");
            std::process::exit(1);
        }
    }
    */
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

    // Check lvars args.

    check_lvars(&mut ctl, &gex_features);

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
    }
}
