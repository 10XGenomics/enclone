// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Fetch the raw reads for a given barcode, yielding fasta as output.
//
// Usage:
// fetch_raw_reads DIR=... BC=...
// where DIR is the directory to look in and BC is the barcode (with suffix ignored).
//
// As implemented, this outputs one fasta record per read pair, and puts the R1 sequence
// in the header.  This makes some sense if R1 contains only the barcode and UMI.

extern crate io_utils;
extern crate pretty_trace;
extern crate string_utils;
extern crate vector_utils;

use io_utils::*;
use pretty_trace::*;
use std::env;
use std::fs::read_dir;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let mut dir = String::new();
    let mut bc = String::new();
    for i in 1..args.len() {
        if args[i].starts_with("DIR=") {
            dir = args[i].after("DIR=").to_string();
        }
        if args[i].starts_with("BC=") {
            bc = args[i].after("BC=").to_string();
        }
    }
    if bc.contains('-') {
        bc = bc.before("-").to_string();
    }
    let files = read_dir(&dir).unwrap();
    let mut fs = Vec::<String>::new();
    for f in files {
        let s: String = f.unwrap().file_name().into_string().unwrap();
        fs.push(s);
    }
    for f1 in fs.iter() {
        if f1.contains("_R1_") {
            for f2 in fs.iter() {
                if f2.contains("_R2_") {
                    let f21 = f2.replace("_R2_", "_R1_");
                    if *f1 == f21 {
                        let (g1, g2) = (format!("{}/{}", dir, f1), format!("{}/{}", dir, f2));
                        let (mut lines1, mut lines2) = (Vec::<String>::new(), Vec::<String>::new());
                        read_maybe_unzipped(&g1, &mut lines1);
                        read_maybe_unzipped(&g2, &mut lines2);
                        assert_eq!(lines1.len(), lines2.len());
                        for i in (0..lines1.len()).step_by(4) {
                            if lines1[i + 1].contains(&bc) {
                                // println!(">{}",lines1[i].between("@", " "));
                                // println!("{}", lines1[i+1]);
                                println!(">{}_{}", lines2[i].between("@", " "), lines1[i + 1]);
                                println!("{}", lines2[i + 1]);
                            }
                        }
                    }
                }
            }
        }
    }
}
