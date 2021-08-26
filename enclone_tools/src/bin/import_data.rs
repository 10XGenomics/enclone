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
// For use at 10x Genomics.

use enclone_core::defs::*;
use enclone_core::testlist::*;
use enclone_tools::copy_for_enclone::*;
use enclone_tools::feature_barcode_matrix::*;
use io_utils::*;
use mirror_sparse_matrix::write_to_file;
use pretty_trace::*;
use serde_json::Value;
use std::collections::HashMap;
use std::env;
use std::fs::{remove_dir_all, rename, File};
use std::io::{BufRead, BufReader};
use std::process::Command;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
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
    let args: Vec<String> = env::args().collect();
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

        // Move directories if they exist.

        let mut moved = false;
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

        // Start copy.

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

        // Determine if the feature barcode matrix for the top feature barcodes should be
        // generated, and if so, what id to use.

        let mut seq_def = None;
        if path_exists(&format!("{}/../SC_RNA_COUNTER_PD", p)) {
            seq_def = Some(feature_barcode_matrix_seq_def(id.force_usize()));
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
                        read_path: read_path,
                        sample_indices: sample_indices,
                        lanes: lanes,
                    });
                }
            }
        }

        // Build feature barcode matrix for top feature barcodes.

        if seq_def.is_some() {
            let m = feature_barcode_matrix(&seq_def.unwrap(), id.force_usize(), false);
            if m.is_ok() {
                let m = m.unwrap();
                for i in (0..dests.len()).rev() {
                    let dest = &dests[i];
                    let target = format!("{}/{}", dest, id);
                    write_to_file(
                        &m,
                        &format!("{}/outs/feature_barcode_matrix_top.bin", target),
                    );
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
