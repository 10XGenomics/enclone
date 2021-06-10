// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Experimental code to find find unexpected feature barcodes.

use io_utils::*;
use itertools::Itertools;
use pretty_trace::*;
use std::env;
use std::io::{BufRead, BufReader};
use std::fs;
use std::fs::File;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let pipestance = &args[1];

    // Get read path and lanes and sample indices from the invocation file.

    let invocation = format!("{}/_invocation", pipestance);
    if !path_exists(&invocation) {
        println!("\n_invocation does not exist\n");
        std::process::exit(1);
    }
    let mut read_path = String::new();   // path to reads
    let mut si = Vec::<String>::new();   // sample indices
    let mut lanes = Vec::<usize>::new(); // lanes
    {
        let f = open_for_read![&invocation];
        let mut lines = Vec::<String>::new();
        for line in f.lines() {
            let s = line.unwrap();
            lines.push(s.to_string());
        }
        for j in 0..lines.len() {
            let s = &lines[j];
            if s.contains("\"read_path\": ") {
                let mut s = s.after("\"read_path\": ");
                if s.contains("\"") {
                    s = s.after("\"");
                    if s.contains("\"") {
                        read_path = s.before("\"").to_string();
                    }
                }
            } else if s.contains("\"lanes\": ") {
                let mut s = s.after("\"lanes\": ");
                if s.starts_with("[") && s.ends_with("],") {
                    s = s.between("[", "],");
                    if s.parse::<usize>().is_ok() {
                        lanes.push(s.force_usize());
                    }
                }
            }
        }
        if read_path.len() == 0 {
            eprintln!("\nfailed to find read path\n");
            std::process::exit(1);

        }
        if !path_exists(&read_path) {
            eprintln!("\nread path does not exist");
            std::process::exit(1);
        }
        let mut j = 0;
        let sample_indices_head = "\"sample_indices\": [";
        while j < lines.len() {
            if lines[j].contains(&sample_indices_head) {
                if lines[j].contains("]") {
                    si.push(lines[j].between("[\"", "\"]").to_string());
                    j += 1;
                } else {
                    j += 1;
                    while j < lines.len() && !lines[j].contains("]") {
                        si.push(lines[j].between("\"", "\"").to_string());
                        j += 1;
                    }
                }
            } else {
                j += 1;
            }
        }
    }
    if lanes.len() == 0 {
        eprintln!("\nfailed to find lanes\n");
        std::process::exit(1);
    }
    let attr = fs::metadata(&read_path).unwrap();
    if !attr.is_dir() {
        eprintln!("\nread path is not a directory\n");
        std::process::exit(1);
    }

    // Check to see if reads exist.

    let mut reads_exist = false;
    let x = dir_list(&read_path);
    for i in 0..x.len() {
        if !x[i].starts_with(".") {
            reads_exist = true;
            break;
        }
    }
    if !reads_exist {
        eprintln!("\nreads do not exist\n");
        std::process::exit(1);
    }

    // Report what we found.

    println!("\nread path = {}", read_path);
    println!("lanes = {}", lanes.iter().format(","));
    println!("sample indices = {}\n", si.iter().format(","));
}

// read-I1_si-TCACGTTGGG_lane-001-chunk-001.fastq.gz
// read-I2_si-TCACGTTGGG_lane-001-chunk-001.fastq.gz
// read-RA_si-TCACGTTGGG_lane-001-chunk-001.fastq.gz
