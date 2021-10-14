// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Copy a 10x pipestance, retaining only the files used by enclone, or which might otherwise
// be convenient to keep.

use enclone_core::slurp::slurp_h5;
use io_utils::{dir_list, path_exists};
use lz4::EncoderBuilder;
use mirror_sparse_matrix::{write_to_file, MirrorSparseMatrix};
use std::fs::{copy, remove_file, File};
use vector_utils::VecUtils;

fn copy_file(f: &str, dir1: &str, dir2: &str) {
    let x = copy(&format!("{}/{}", dir1, f), &format!("{}/{}", dir2, f));
    if x.is_err() {
        eprintln!("\nFailed to copy {} from\n{}\nto {}.\n", f, dir1, dir2);
    }
    x.unwrap();
}

fn mkdir(d: &str) {
    let x = std::fs::create_dir(&d);
    if x.is_err() {
        eprintln!("\nCan't make directory {}.\n", d);
    }
    x.unwrap();
}

fn mkdirp(d: &str) {
    std::fs::create_dir_all(&d).unwrap();
}

fn compress(source: &str, destination: &str) {
    let mut input_file = File::open(source).unwrap();
    let output_file = File::create(destination).unwrap();
    let mut encoder = EncoderBuilder::new().level(4).build(output_file).unwrap();
    std::io::copy(&mut input_file, &mut encoder).unwrap();
    let _ = encoder.finish();
}

fn lz4_file(f: &str) {
    compress(f, &format!("{}.lz4", f));
    remove_file(&f).unwrap();
}

pub fn copy_for_enclone(source: &str, target: &str) {
    let p = format!("{}/outs", source);

    // Define most of the gex files to be copied.

    let gex_files = [
        // "barcode_correction_csv.csv", // don't know why this was here
        // "barcode_summary.h5",         // don't know why this was here
        "filtered_barcodes.csv",
        "filtered_feature_bc_matrix.h5",
        "gene_properties.json",
        "metrics_summary_csv.csv",
        "metrics_summary_json.json",
        // "per_barcode_metrics.csv",    // don't know why this was here
        "per_feature_metrics.csv",
        "web_summary.html",
        "web_summary_pd.html",
    ]
    .to_vec();

    // Copy VDJ.

    if path_exists(&format!("{}/all_contig_annotations.json", p)) {
        mkdirp(&format!("{}/outs", target));

        // Copy and compress all_contig_annotations.json.

        let json_target = format!("{}/outs/all_contig_annotations.json", target);
        copy(&format!("{}/all_contig_annotations.json", p), &json_target).unwrap();
        lz4_file(&json_target);

        // Copy _invocation and _log.

        copy(
            &format!("{}/../_invocation", p),
            &format!("{}/_invocation", target),
        )
        .unwrap();
        copy(&format!("{}/../_log", p), &format!("{}/_log", target)).unwrap();

        // Copy these two files used for Immcantation.

        copy(
            &format!("{}/filtered_contig_annotations.csv", p),
            &format!("{}/outs/filtered_contig_annotations.csv", target),
        )
        .unwrap();
        copy(
            &format!("{}/filtered_contig.fasta", p),
            &format!("{}/outs/filtered_contig.fasta", target),
        )
        .unwrap();

        // Copy these files that are nice to have.

        copy(
            &format!("{}/cell_barcodes.json", p),
            &format!("{}/outs/cell_barcodes.json", target),
        )
        .unwrap();
        copy(
            &format!("{}/metrics_summary_json.json", p),
            &format!("{}/outs/metrics_summary_json.json", target),
        )
        .unwrap();

    // Copy GEX.
    } else if path_exists(&format!("{}/raw_feature_bc_matrix.h5", p))
        || path_exists(&format!("{}/raw_gene_bc_matrices_h5.h5", p))
    {
        mkdirp(&format!("{}/outs", target));

        // Copy _invocation.

        copy(
            &format!("{}/../_invocation", p),
            &format!("{}/_invocation", target),
        )
        .unwrap();

        // Copy h5 file, allowing for older-named version.

        let f;
        if path_exists(&format!("{}/raw_feature_bc_matrix.h5", p)) {
            f = "raw_feature_bc_matrix.h5";
        } else {
            f = "raw_gene_bc_matrices_h5.h5";
        }
        let h5_target = format!("{}/outs/{}", target, f);
        copy(&format!("{}/{}", p, f), &h5_target).unwrap();

        // Generate feature_barcode_matrix.bin.

        let mut barcodes = Vec::<String>::new();
        let mut features = Vec::<String>::new();
        let mut matrix = Vec::<Vec<(i32, i32)>>::new();
        slurp_h5(&h5_target, true, &mut barcodes, &mut features, &mut matrix);
        let msm = MirrorSparseMatrix::build_from_vec(&matrix, &barcodes, &features);
        let bin_file = format!("{}/outs/feature_barcode_matrix.bin", target);
        write_to_file(&msm, &bin_file);

        // Copy other files in outs.

        for f in gex_files.iter() {
            if path_exists(&format!("{}/{}", p, f)) {
                copy(&format!("{}/{}", p, f), &format!("{}/outs/{}", target, f)).unwrap();
            }
        }

        // Copy subdirectories of outs.

        {
            use fs_extra::dir::{copy, CopyOptions};
            let opt = CopyOptions::new();
            let mut d = "filtered_feature_bc_matrix";
            if !path_exists(&format!("{}/{}", p, d)) {
                d = "filtered_gene_bc_matrices_mex";
            }
            copy(&format!("{}/{}", p, d), &format!("{}/outs", target), &opt).unwrap();
            let dirs = vec!["analysis".to_string(), "analysis_csv".to_string()];
            for d in dirs.iter() {
                copy(&format!("{}/{}", p, d), &format!("{}/outs", target), &opt).unwrap();
            }
        }

    // Copy multi.
    } else if path_exists(&format!("{}/multi", p))
        && (path_exists(&format!("{}/count", p)) || path_exists(&format!("{}/count_pd", p)))
    {
        let count;
        if path_exists(&format!("{}/count_pd", p)) {
            count = "count_pd";
        } else {
            count = "count";
        }
        mkdirp(&format!("{}/outs/multi", target));
        mkdir(&format!("{}/outs/{}", target, count));
        std::fs::copy(
            &format!("{}/../_invocation", p),
            &format!("{}/_invocation", target),
        )
        .unwrap();
        std::fs::copy(&format!("{}/../_log", p), &format!("{}/_log", target)).unwrap();
        for vdj in ["vdj_b", "vdj_t"].iter() {
            let vdj_from = format!("{}/multi/{}", p, vdj);
            if path_exists(&vdj_from) {
                let vdj_to = format!("{}/outs/multi/{}", target, vdj);
                mkdir(&vdj_to);
                let json_target = format!("{}/all_contig_annotations.json", vdj_to);
                copy_file("all_contig_annotations.json", &vdj_from, &vdj_to);
                lz4_file(&json_target);
            }
        }
        if path_exists(&format!("{}/outs/multi/count", target)) {
            let count_to = format!("{}/outs/multi/count", target);
            mkdir(&count_to);
            let count_from = format!("{}/count", p);
            copy_file("raw_feature_bc_matrix.h5", &count_from, &count_to);
        }
        if path_exists(&format!("{}/vdj_reference", p)) {
            mkdirp(&format!("{}/outs/vdj_reference/fasta", target));
            let ref_from = format!("{}/vdj_reference", p);
            let ref_to = format!("{}/outs/vdj_reference", target);
            copy_file("reference.json", &ref_from, &ref_to);
            let fasta_from = format!("{}/vdj_reference/fasta", p);
            let fasta_to = format!("{}/outs/vdj_reference/fasta", target);
            copy_file("regions.fa", &fasta_from, &fasta_to);
        }

        let summary = format!("{}/multi_web_summary_json/metrics_summary_csv", p);
        if path_exists(&summary) {
            let summary_out = format!("{}/outs/multi_web_summary_json/metrics_summary_csv", target);
            mkdirp(&summary_out);
            let list = dir_list(&summary);
            for f in list.iter() {
                copy_file(f, &summary, &summary_out);
            }
        }

        let mut gex_files = gex_files.clone();
        gex_files.push("raw_feature_bc_matrix.h5");
        let count_pd_from = format!("{}/{}", p, count);
        let count_pd_to = format!("{}/outs/{}", target, count);
        for f in gex_files.iter() {
            if path_exists(&format!("{}/{}", count_pd_from, f)) {
                copy_file(f, &count_pd_from, &count_pd_to);
            }
        }

        // Generate feature_barcode_matrix.bin.

        let h5_target = format!("{}/raw_feature_bc_matrix.h5", count_pd_to);
        let mut barcodes = Vec::<String>::new();
        let mut features = Vec::<String>::new();
        let mut matrix = Vec::<Vec<(i32, i32)>>::new();
        slurp_h5(&h5_target, true, &mut barcodes, &mut features, &mut matrix);
        let msm = MirrorSparseMatrix::build_from_vec(&matrix, &barcodes, &features);
        let bin_file = format!("{}/feature_barcode_matrix.bin", count_pd_to);
        write_to_file(&msm, &bin_file);

        // Copy more.

        use fs_extra::dir::{copy, CopyOptions};
        let opt = CopyOptions::new();
        let dirs = ["filtered_feature_bc_matrix", "analysis_csv"];
        for d in dirs.iter() {
            copy(
                &format!("{}/{}/{}", p, count, d),
                &format!("{}/outs/{}", target, count),
                &opt,
            )
            .unwrap();
        }

        // Copy per_sample_outs.

        if path_exists(&format!("{}/per_sample_outs", p)) {
            let list = dir_list(&format!("{}/per_sample_outs", p));
            if list.solo() {
                let x = &list[0];
                let mut d = format!("per_sample_outs/{}/count/analysis", x);
                if !path_exists(&d) {
                    d = format!("per_sample_outs/{}/count/analysis_csv", x);
                }
                mkdirp(&format!("{}/outs/{}", target, d));
                copy(
                    &format!("{}/{}", p, d),
                    &format!("{}/outs/per_sample_outs/{}/count", target, x),
                    &opt,
                )
                .unwrap();
            }
        }

    // Copy new dir structure.
    } else if path_exists(&format!("{}/count", p)) || path_exists(&format!("{}/count_pd", p)) {
        let count;
        if path_exists(&format!("{}/count_pd", p)) {
            count = "count_pd";
        } else {
            count = "count";
        }
        mkdirp(&format!("{}/outs/{}", target, count));
        std::fs::copy(
            &format!("{}/../_invocation", p),
            &format!("{}/_invocation", target),
        )
        .unwrap();
        std::fs::copy(&format!("{}/../_log", p), &format!("{}/_log", target)).unwrap();
        for vdj in ["vdj_b", "vdj_t"].iter() {
            let vdj_from = format!("{}/{}", p, vdj);
            if path_exists(&vdj_from) {
                let vdj_to = format!("{}/outs/{}", target, vdj);
                mkdir(&vdj_to);
                let json_target = format!("{}/all_contig_annotations.json", vdj_to);
                copy_file("all_contig_annotations.json", &vdj_from, &vdj_to);
                lz4_file(&json_target);
            }
        }
        if path_exists(&format!("{}/vdj_reference", p)) {
            mkdirp(&format!("{}/outs/vdj_reference/fasta", target));
            let ref_from = format!("{}/vdj_reference", p);
            let ref_to = format!("{}/outs/vdj_reference", target);
            copy_file("reference.json", &ref_from, &ref_to);
            let fasta_from = format!("{}/vdj_reference/fasta", p);
            let fasta_to = format!("{}/outs/vdj_reference/fasta", target);
            copy_file("regions.fa", &fasta_from, &fasta_to);
        }
        let mut gex_files = gex_files.clone();
        gex_files.push("raw_feature_bc_matrix.h5");
        let count_pd_from = format!("{}/{}", p, count);
        let count_pd_to = format!("{}/outs/{}", target, count);
        for f in gex_files.iter() {
            if path_exists(&format!("{}/{}", count_pd_from, f)) {
                copy_file(f, &count_pd_from, &count_pd_to);
            }
        }

        // Generate feature_barcode_matrix.bin.

        let h5_target = format!("{}/raw_feature_bc_matrix.h5", count_pd_to);
        let mut barcodes = Vec::<String>::new();
        let mut features = Vec::<String>::new();
        let mut matrix = Vec::<Vec<(i32, i32)>>::new();
        slurp_h5(&h5_target, true, &mut barcodes, &mut features, &mut matrix);
        let msm = MirrorSparseMatrix::build_from_vec(&matrix, &barcodes, &features);
        let bin_file = format!("{}/feature_barcode_matrix.bin", count_pd_to);
        write_to_file(&msm, &bin_file);

        // Copy more.

        use fs_extra::dir::{copy, CopyOptions};
        let opt = CopyOptions::new();
        let dirs = ["filtered_feature_bc_matrix", "analysis_csv"];
        for d in dirs.iter() {
            copy(
                &format!("{}/{}/{}", p, count, d),
                &format!("{}/outs/{}", target, count),
                &opt,
            )
            .unwrap();
        }

    // Or something is wrong.
    } else {
        eprintln!(
            "\nThe directory {} is not a complete VDJ or GEX or MULTI directory.\n",
            p
        );
        std::process::exit(1);
    }
}
