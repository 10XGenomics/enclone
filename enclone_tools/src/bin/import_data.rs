// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Takes a single argument: a comma-separated list of ids, allowing hyphenated ranges.  Copy these
// to the internal collections.  If it has already been copied, it is moved aside, and at the end,
// deleted (assuming that the copy did not fail.
//
// You need to set quota first on the second internal location.
//
// This is probably not fully consistent with current pipeline structure.
//
// Run make_enclone_testlist_all after this to update the catalog.
//
// Optional second argument: FB_INFO: do nothing except attempt to create the
// feature barcode matrix.
//
// Optional second argument: FB_INFO_WRITE: do nothing except attempt to create the
// feature barcode matrix and write the corresponding files.
//
// For use at 10x Genomics.

use enclone_core::defs::get_config;
use enclone_core::packing::*;
use enclone_core::testlist::TEST_FILES_VERSION;
use enclone_tools::copy_for_enclone::copy_for_enclone;
use enclone_tools::feature_barcode_matrix::{
    feature_barcode_matrix, feature_barcode_matrix_seq_def, SequencingDef,
};
use io_utils::{fwriteln, open_for_read, open_for_write_new, path_exists};
use mirror_sparse_matrix::write_to_file;
use pretty_trace::PrettyTrace;
use serde_json::Value;
use std::collections::HashMap;
use std::env;
use std::fs::{copy, remove_dir_all, rename, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use string_utils::{parse_csv, TextUtils};

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let mut fb_info = false;
    let mut fb_info_write = false;
    let mut verbosity = 0;
    for i in 2..args.len() {
        if args[i] == "FB_INFO" {
            fb_info = true;
        } else if args[i] == "FB_INFO_WRITE" {
            fb_info_write = true;
        } else if args[i].starts_with("VERBOSITY=") {
            verbosity = args[i].after("VERBOSITY=").force_usize();
        } else {
            eprintln!("\nIllegal arg.\n");
            std::process::exit(1);
        }
    }
    let mut config = HashMap::<String, String>::new();
    let mut config_file = String::new();
    for (key, value) in env::vars() {
        if key == "ENCLONE_CONFIG" {
            config_file = value.to_string();
            if config_file.contains(',') {
                config_file = config_file.after(",").to_string();
            }
        }
    }
    let _ = get_config(&config_file, &mut config);
    let mut dests = Vec::<String>::new();
    dests.push(format!("{}/current{}", config["earth"], TEST_FILES_VERSION));
    dests.push(format!("{}/current{}", config["cloud"], TEST_FILES_VERSION));
    let ids0 = args[1].split(',').collect::<Vec<&str>>();
    let mut ids = Vec::<String>::new();
    for id in ids0.iter() {
        if id.contains('-') {
            let start = id.before("-").force_usize();
            let stop = id.after("-").force_usize();
            for p in start..=stop {
                ids.push(format!("{}", p));
            }
        } else {
            ids.push(id.to_string());
        }
    }
    let ones = &config["ones"];
    for id in ids.iter() {
        // Get path.

        let mut p = id.clone();
        assert!(p.parse::<usize>().is_ok());
        let url = format!("{}/{}", ones, p);
        let o = Command::new("curl")
            .arg(&url)
            .output()
            .expect("failed to execute http");
        let m = String::from_utf8(o.stdout).unwrap();
        if m.contains("502 Bad Gateway") {
            eprintln!(
                "\nWell this is sad.  The URL \
                {} yielded a 502 Bad Geteway \
                message.  Either try again later or ask someone for help.\n\n",
                url,
            );
            std::process::exit(1);
        }
        if m.contains("\"path\":\"") {
            let path = m.between("\"path\":\"", "\"").to_string();
            p = format!("{}/outs", path);
            if !path_exists(&p) {
                eprintln!(
                    "\nIt looks like you've provided an id {} for \
                    which\nthe pipeline outs folder has not yet been \
                    generated.\n\n",
                    p
                );
                std::process::exit(1);
            }
        } else {
            eprintln!(
                "\nIt looks like you've provided either an incorrect \
                id {} or else one for which\n\
                the pipeline outs folder has not yet been generated.\n\n",
                p
            );
            std::process::exit(1);
        }

        // Punt if the count pipeline appears to have succeeded but actually failed.

        if path_exists(&format!("{}/../SC_RNA_COUNTER_PD", p))
            && !path_exists(&format!("{}/analysis", p))
        {
            println!(
                "\nNo analysis directory found for\n{}\nso the outs directory is \
                incomplete, possibly because there were not enough reads.  Giving up.\n",
                p
            );
            continue;
        }

        // Move directories if they exist.

        let mut moved = false;
        if !fb_info && !fb_info_write {
            for dest in dests.iter() {
                if path_exists(&format!("{}/{}", dest, id)) {
                    if path_exists(&format!("{}/{}.aside", dest, id)) {
                        eprintln!("\nPlease remove {}/{}.aside.\n", dest, id);
                        std::process::exit(1);
                    }
                    rename(
                        &format!("{}/{}", dest, id),
                        &format!("{}/{}.aside", dest, id),
                    )
                    .unwrap();
                    moved = true;
                }
            }
        }

        // Start copy.

        if !fb_info && !fb_info_write {
            println!("copying {} using path = {}", id, p);
            for i in (0..dests.len()).rev() {
                let dest = &dests[i];
                let target = format!("{}/{}", dest, id);
                if path_exists(&target) {
                    eprintln!("\nPlease delete {}.\n", target);
                    std::process::exit(1);
                }
                copy_for_enclone(&format!("{}/..", p), &target);
            }
        }

        // Determine if the feature barcode matrix for the top feature barcodes should be
        // generated, and if so, what id to use.

        let mut seq_def = None;
        if path_exists(&format!("{}/../SC_RNA_COUNTER_PD", p)) {
            seq_def = feature_barcode_matrix_seq_def(id.force_usize());
        } else if path_exists(&format!("{}/../SC_MULTI_PD", p)) {
            let mut sample_indices = Vec::<String>::new();
            let mut lanes = Vec::<usize>::new();
            let read_path;
            {
                let mut antibody_seq_id = None;
                let inv = format!("{}/../_invocation", p);
                let f = open_for_read![inv];
                let mut lines = Vec::<String>::new();
                for line in f.lines() {
                    let s = line.unwrap();
                    lines.push(s);
                }
                for i in 0..lines.len() {
                    let fields = parse_csv(&lines[i]);
                    if fields.len() >= 5 && fields[4] == "Antibody Capture" {
                        antibody_seq_id = Some(fields[3].force_usize());
                        break;
                    }
                }
                if antibody_seq_id.is_some() {
                    let antibody_seq_id = antibody_seq_id.unwrap();
                    let pid = m.between("\"pipestance_id\":\"", "\"").to_string();
                    let meta = &config["meta"];
                    let url = format!("{}/{}", meta, pid);
                    println!("getting data for {} from {}", antibody_seq_id, url);
                    let o = Command::new("curl")
                        .arg(&url)
                        .output()
                        .expect("failed to execute curl for meta");
                    let mm = String::from_utf8(o.stdout).unwrap();
                    let v: Value = serde_json::from_str(&mm).unwrap();
                    let asi = format!("{}", antibody_seq_id);
                    let rrr = &v["sample_bag"]["sequencing_libraries"][&asi];
                    let lane = &rrr["metadata"]["lane"];
                    let lane = lane.to_string().between("\"", "\"").to_string();
                    let lane = lane.split(',').collect::<Vec<&str>>();
                    for x in lane.iter() {
                        lanes.push(x.force_usize());
                    }
                    let si_data = rrr["sample_indexes"].as_array().unwrap();
                    for j in 0..si_data.len() {
                        sample_indices.push(
                            si_data[j]["seq"]
                                .to_string()
                                .between("\"", "\"")
                                .to_string(),
                        );
                    }
                    let flowcell = rrr["sequencing_run"]["name"]
                        .to_string()
                        .between("\"", "\"")
                        .to_string();
                    read_path = v["fastq_paths"][&flowcell]
                        .to_string()
                        .between("\"", "\"")
                        .to_string();
                    seq_def = Some(SequencingDef {
                        read_path,
                        sample_indices,
                        lanes,
                    });
                }
            }
        }

        // Build feature barcode matrix for top feature barcodes.

        if seq_def.is_some() {
            // Get list of antibody reference feature barcodes.

            let mut ref_fb = Vec::<String>::new();
            {
                let mut f = format!("{}/multi/count/feature_reference.csv", p);
                if !path_exists(&f) {
                    let g = open_for_read![&format!("{}/../_invocation", p)];
                    for line in g.lines() {
                        let s = line.unwrap();
                        if s.contains("feature_reference") {
                            f = s.between("\"", "\"").to_string();
                        }
                    }
                }
                for i in (0..dests.len()).rev() {
                    let dest = &dests[i];
                    let target = format!("{}/{}", dest, id);
                    copy(&f, &format!("{}/outs/feature_reference.csv", target)).unwrap();
                }
                let f = open_for_read![&f];
                let mut seq_pos = 0;
                for (i, line) in f.lines().enumerate() {
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    if i == 0 {
                        for j in 0..fields.len() {
                            if fields[j] == "sequence" {
                                seq_pos = j;
                            }
                        }
                    } else {
                        ref_fb.push(fields[seq_pos].to_string());
                    }
                }
            }
            ref_fb.sort();

            // Keep going.

            let mut verb = verbosity;
            if verb == 0 && (fb_info || fb_info_write) {
                verb = 1;
            }
            let x = feature_barcode_matrix(&seq_def.unwrap(), id.force_usize(), verb, &ref_fb);
            if fb_info {
                std::process::exit(0);
            }
            if x.is_ok() {
                let (
                    m,
                    total,
                    brn,
                    common_gumi_freq,
                    common_gumi_content,
                    m_reads,
                    total_reads,
                    brnr,
                    bdcs,
                ) = x.unwrap();
                for i in (0..dests.len()).rev() {
                    let dest = &dests[i];
                    let target = format!("{}/{}", dest, id);
                    write_to_file(
                        &m,
                        &format!("{}/outs/feature_barcode_matrix_top.bin", target),
                    );
                    write_to_file(
                        &m_reads,
                        &format!("{}/outs/feature_barcode_matrix_top_reads.bin", target),
                    );
                    let mut f =
                        File::create(&format!("{}/outs/feature_barcode_matrix_top.total", target))
                            .unwrap();
                    f.write_all(&total.to_ne_bytes()).unwrap();
                    let mut f = File::create(&format!(
                        "{}/outs/feature_barcode_matrix_top.total_reads",
                        target
                    ))
                    .unwrap();
                    f.write_all(&total_reads.to_ne_bytes()).unwrap();
                    let mut f = open_for_write_new![&format!(
                        "{}/outs/feature_barcode_matrix_top.brn",
                        target
                    )];
                    for j in 0..brn.len() {
                        fwriteln!(f, "{},{},{}", brn[j].0, brn[j].1, brn[j].2);
                    }
                    let mut f = open_for_write_new![&format!(
                        "{}/outs/feature_barcode_matrix_top.brnr",
                        target
                    )];
                    for j in 0..brnr.len() {
                        fwriteln!(f, "{},{},{}", brnr[j].0, brnr[j].1, brnr[j].2);
                    }
                    let mut f = open_for_write_new![&format!(
                        "{}/outs/feature_barcode_matrix_top.bdcs",
                        target
                    )];
                    for j in 0..bdcs.len() {
                        fwriteln!(f, "{},{},{},{}", bdcs[j].0, bdcs[j].1, bdcs[j].2, bdcs[j].3);
                    }
                    let mut f = File::create(&format!(
                        "{}/outs/feature_barcode_matrix.common_gumis",
                        target
                    ))
                    .unwrap();
                    let mut bytes = Vec::<u8>::new();
                    bytes.append(&mut save_vec_f32(&common_gumi_freq));
                    bytes.append(&mut save_vec_vec_u8(&common_gumi_content));
                    f.write_all(&bytes).unwrap();
                }
            }
        }

        // Remove moved directories.

        if moved {
            for dest in dests.iter() {
                let movedir = format!("{}/{}.aside", dest, id);
                if path_exists(&movedir) {
                    remove_dir_all(&movedir).unwrap();
                }
            }
        }
    }
}
