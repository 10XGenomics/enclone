// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Build inputs for running other clonotyping software.
//
// Pass one argument, a list of numerical ids of internal datasets, allowing hyphenated ranges.
// Also, @test OK.
//
// Creates two files:
// filtered_contig.fasta
// filtered_contig_annotations.csv.

use enclone_core::defs::get_config;
use enclone_core::test_def::replace_at_test;
use enclone_core::testlist::*;
use io_utils::*;
use pretty_trace::*;
use std::collections::HashMap;
use std::env;
use std::io::{BufRead, Write};
use string_utils::*;
use vector_utils::VecUtils;

pub fn main() {
    PrettyTrace::new().on();

    // Get list of ids.

    let mut args: Vec<String> = env::args().collect();
    for i in 1..args.len() {
        replace_at_test(&mut args[i]);
        args[i] = args[i].replace(":", ",");
        args[i] = args[i].replace(";", ",");
    }
    let ids0 = args[1].split(',').collect::<Vec<&str>>();
    let mut ids = Vec::<usize>::new();
    for id in ids0.iter() {
        if id.contains('-') {
            let start = id.before("-").force_usize();
            let stop = id.after("-").force_usize();
            for p in start..=stop {
                ids.push(p);
            }
        } else {
            ids.push(id.force_usize());
        }
    }
    ids.sort();

    // Get configuration.

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

    // Location of files.

    let pre = format!("{}/current{}", config["cloud"], TEST_FILES_VERSION);

    // Concatenate the fasta files, prepending the id to the contig name.
    // Also concatenate the csv files, prepending the id to the contig name.

    let mut g1 = open_for_write_new!["filtered_contig.fasta"];
    let mut g = open_for_write_new!["filtered_contig_annotations.csv"];
    for (i, l) in ids.iter().enumerate() {
        // Find the filtered_contig input files.

        let outs = format!("{pre}/{l}/outs");
        let mut fdir = outs.clone();
        if path_exists(&format!("{outs}/per_sample_outs")) {
            let list = dir_list(&format!("{outs}/per_sample_outs"));
            if list.solo() {
                let x = &list[0];
                let vdj_b = format!("{outs}/per_sample_outs/{x}/vdj_b");
                if path_exists(&vdj_b) {
                    fdir = vdj_b;
                }
            }
        }

        // Edit and copy them.

        let mut count1 = 0;
        let f1 = open_for_read![&format!("{fdir}/filtered_contig.fasta")];
        for line in f1.lines() {
            let s = line.unwrap();
            if s.starts_with('>') {
                let bc = s.between(">", "-");
                let contig = s.rev_after("_");
                fwriteln!(g1, ">{}-{}_contig_{}", bc, l, contig);
                count1 += 1;
            } else {
                fwriteln!(g1, "{}", s);
            }
        }
        let f = open_for_read![&format!("{fdir}/filtered_contig_annotations.csv")];
        let mut count2 = 0;
        let mut first = true;
        for line in f.lines() {
            let s = line.unwrap();
            if i == 0 && first {
                fwriteln!(g, "{}", s);
            }
            if first {
                first = false;
                continue;
            }
            let fields = s.split(',').collect::<Vec<&str>>();
            assert!(fields[2].contains("_contig_"));
            for j in 0..fields.len() {
                if j > 0 {
                    fwrite!(g, ",");
                }
                if j == 0 {
                    fwrite!(g, "{}-{}", fields[j].before("-"), l);
                } else if j == 2 {
                    fwrite!(
                        g,
                        "{}-{}_contig_{}",
                        fields[j].before("-"),
                        l,
                        fields[j].rev_after("_")
                    );
                    count2 += 1;
                } else {
                    fwrite!(g, "{}", fields[j]);
                }
            }
            fwriteln!(g, "");
        }
        if count1 != count2 {
            eprintln!("\ninconsistency");
            eprintme!(l, count1, count2);
            eprintln!("");
            std::process::exit(1);
        }
    }
}
