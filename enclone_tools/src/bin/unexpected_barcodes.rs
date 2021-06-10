// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Experimental code to find find unexpected feature barcodes.

use flate2::read::MultiGzDecoder;
use io_utils::*;
use itertools::Itertools;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::env;
use std::io::{BufRead, BufReader};
use std::fs;
use std::fs::File;
use std::time::Instant;
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    let t = Instant::now();
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

    // Find the read files.

    let mut read_files = Vec::<String>::new();
    let x = dir_list(&read_path);
    println!("");
    for f in x.iter() {
        for sample_index in si.iter() {
            for lane in lanes.iter() {
                if f.contains(&format!("read-RA_si-{}_lane-00{}-", sample_index, lane)) {
                    println!("{}", f);
                    read_files.push(f.clone());
                }
            }
        }
    }
    if read_files.is_empty() {
        eprintln!("\nreads do not exist\n");
        std::process::exit(1);
    }

    // Report what we found.

    println!("\nread path = {}", read_path);
    println!("lanes = {}", lanes.iter().format(","));
    println!("sample indices = {}", si.iter().format(","));
    println!("used {:.1} seconds\n", elapsed(&t));

    // Traverse the reads.

    let mut buf = Vec::<(Vec<u8>, Vec<u8>, Vec<u8>)>::new(); // {(barcode, umi, fb)}
    for rf in read_files.iter() {
        let f = format!("{}/{}", read_path, rf);
        let gz = MultiGzDecoder::new(File::open(&f).unwrap());
        let b = BufReader::new(gz);

        // Paired reads are in groups of eight lines.  Line 2 is the cell barcode-umi read, 
        // and line 6 is the read that contains the feature barcode.

        let mut count = 0;
        let mut barcode = Vec::<u8>::new();
        let mut umi = Vec::<u8>::new();
        let mut fb;
        for line in b.lines() {
            count += 1;
            if count % 8 == 2 || count % 8 == 6 {
                let s = line.unwrap();
                let s = s.as_bytes();
                if count % 8 == 2 {
                    assert!(s.len() >= 28);
                    barcode = s[0..16].to_vec();
                    umi = s[16..28].to_vec();
                } else {
                    fb = s[10..25].to_vec();
                    buf.push((barcode.clone(), umi.clone(), fb.clone()));
                }
            }
        }
    }
    println!("there are {} read pairs", buf.len());
    println!("\nused {:.1} seconds\n", elapsed(&t));

    // Unique sort.
    
    buf.par_sort();
    buf.dedup();
    println!("there are {} uniques", buf.len());
    println!("\nused {:.1} seconds\n", elapsed(&t));

    // Report common feature barcodes.

    println!("common feature barcodes\n");
    let mut fbx = Vec::<Vec<u8>>::new();
    for i in 0..buf.len() {
        fbx.push(buf[i].2.clone());
    }
    fbx.par_sort();
    let mut freq = Vec::<(u32, Vec<u8>)>::new();
    make_freq(&fbx, &mut freq);
    for i in 0..10 {
        println!("{} = {}", strme(&freq[i].1), freq[i].0);
    }
    println!("\nused {:.1} seconds\n", elapsed(&t));
}
