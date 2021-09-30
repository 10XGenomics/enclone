// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Slow regression test code for IGHD constant region results.
//
// This will break when we change denovo.rs to not just run IGHD.

use io_utils::*;
use itertools::Itertools;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::io::Write;
use std::process::Command;
use std::time::Instant;
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    let t = Instant::now();
    let mut id_seq = Vec::<(String, String, String)>::new();
    let f = include_str!("../const_ighd_best_hits");
    for line in f.lines() {
        if line.starts_with('#') || line.contains("//") {
            continue;
        }
        if line.ends_with(" old") || line.ends_with(" new/good") || line.ends_with(" correct") {
            let id = line.before(":").to_string();
            let species = line.between(":", " ").to_string();
            let seq = line.rev_before(" ").rev_after(" ").to_string();
            id_seq.push((id, species, seq));
        }
    }
    println!();
    let mut results = Vec::<(usize, Vec<u8>, bool)>::new();
    for i in 0..id_seq.len() {
        results.push((i, Vec::<u8>::new(), false));
    }
    const NTHREADS: usize = 4;
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(NTHREADS)
        .build_global();
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let id = &id_seq[i].0;
        let species = &id_seq[i].1;
        let seq = &id_seq[i].2;
        let mut log = Vec::<u8>::new();
        fwriteln!(log, "running denovo {} = {}", id, species);
        let o = Command::new("denovo")
            .arg(&id)
            .output()
            .expect("failed to execute denovo");
        if o.status.code().unwrap() != 0 {
            fwriteln!(log, "\nOOPS FAILED!");
            res.2 = true;
        } else {
            let m = String::from_utf8(o.stdout).unwrap();
            fwriteln!(log, "{}", m);
            let mut seqs = Vec::<String>::new();
            for line in m.lines() {
                let fields = line.split(' ').collect::<Vec<&str>>();
                if fields[0] == "IGHD" {
                    let seqx = fields[8].rev_after("|").to_string();
                    seqs.push(seqx);
                }
            }
            unique_sort(&mut seqs);
            let s = format!("{}", seqs.iter().format(","));
            if s != *seq {
                fwriteln!(log, "old = {}, new = {}", seq, s);
                fwriteln!(log, "Answer is not correct.\n");
                res.2 = true;
            }
        }
        res.1 = log;
    });
    let mut fail = false;
    for i in 0..results.len() {
        if results[i].2 {
            fail = true;
        }
        print!("{}", strme(&results[i].1));
    }
    println!("used {:.2} minutes\n", elapsed(&t) / 60.0);
    if fail {
        println!("\nFAILED!\n");
        std::process::exit(1);
    }
}
