// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.
//
// Splice junctions from pgen.

use io_utils::*;
use pretty_trace::*;
use std::env;
use std::io::BufRead;
use string_utils::TextUtils;
use vdj_ann_ref::human_ref;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();

    // Get the standard human reference.

    let mut refnames = Vec::<String>::new();
    let mut refs = Vec::<Vec<u8>>::new();
    let mut utr = Vec::<bool>::new();
    {
        let href = human_ref();
        for (i, line) in href.lines().enumerate() {
            if i % 2 == 0 {
                let n = line.between("|", " ").to_string();
                utr.push(line.contains("5'UTR"));
                refnames.push(n);
            } else {
                refs.push(line.as_bytes().to_vec());
            }
        }
    }
    let nref = refs.len();

    // Load the input file.

    let mut jun = Vec::<Vec<u8>>::new();
    let mut hv = Vec::<String>::new();
    let mut hj = Vec::<String>::new();
    let f = open_for_read![&args[1]];
    for line in f.lines() {
        let s = line.unwrap();
        let fields = s.split(',').collect::<Vec<&str>>();
        jun.push(fields[0].as_bytes().to_vec());
        hv.push(fields[2].to_string());
        hj.push(fields[3].to_string());
    }
    let n = jun.len();

    // Process entries.

    // for i in 0..n {
    for i in 0..5 {
        if hv[i] == "IGHV3-NL1" {
            continue;
        }
        let mut vseq = Vec::<u8>::new();
        let mut jseq = Vec::<u8>::new();
        for j in 0..nref {
            if refnames[j] == hv[i] && !utr[j] {
                vseq = refs[j].clone();
                break;
            }
        }
        if vseq.is_empty() {
            println!("\nEntry {}, unable to find V gene named {}.\n", i + 1, hv[i]);
            std::process::exit(1);
        }
        for j in 0..nref {
            if refnames[j] == hj[i] {
                jseq = refs[j].clone();
                break;
            }
        }
        if jseq.is_empty() {
            println!("\nEntry {}, unable to find J gene named {}.\n", i + 1, hj[i]);
            std::process::exit(1);
        }
        let junv = &jun[i][0..6];
        use string_utils::strme;
        println!("\nV gene = {}", hv[i]);
        println!("\nlooking for {} in {}", strme(&junv), strme(&vseq));
        for k in (0..vseq.len() - junv.len()).rev() {
            if vseq[k..].starts_with(&junv) {
                println!("vstart for {} found at pos -{}", i + 1, vseq.len() - k);
            }
        }
    }
}
