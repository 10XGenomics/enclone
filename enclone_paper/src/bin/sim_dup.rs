// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Predict the expected number of duplicates in the true data.
//
// Has some hardwired numbers.
//
// sim_dup hpgen_sim.csv
//
use io_utils::*;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::env;
use std::io::BufRead;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    // Data = pairs of heavy chain gene names and CDRH3 amino acid sequences.
    let mut data = Vec::<(String, String)>::new();
    for line in f.lines() {
        let s = line.unwrap();
        let fields = s.split(',').collect::<Vec<&str>>();
        data.push((fields[1].to_string(), fields[2].to_string()));
    }
    let mut counts = Vec::new();
    counts.push((364703 as f64 * 0.304).round() as usize);
    counts.push((293532 as f64 * 0.696).round() as usize);
    counts.push((352749 as f64 * 0.671).round() as usize);
    counts.push((397955 as f64 * 0.627).round() as usize);
    let mut dup = vec![false; data.len()];
    for m1 in 0..4 {
        let mut start1 = 0;
        for i in 0..m1 {
            start1 += counts[i];
        }
        let stop1 = start1 + counts[m1];
        // [start1, stop1) are the range of data entries for donor m1+1
        for m2 in m1 + 1..4 {
            let mut start2 = 0;
            for i in 0..m2 {
                start2 += counts[i];
            }
            let stop2 = start2 + counts[m2];
            // [start2, stop2) are the range of data entries for donor m2+1
            let mut results = Vec::<(usize, Vec<usize>)>::new();
            for k1 in start1..stop1 {
                results.push((k1, Vec::new()));
            }
            results.par_iter_mut().for_each(|res| {
                let k1 = res.0;
                for k2 in start2..stop2 {
                    if data[k1] == data[k2] {
                        // push the entries for both cells
                        res.1.push(k1);
                        res.1.push(k2);
                    }
                }
            });
            for i in 0..results.len() {
                for j in 0..results[i].1.len() {
                    dup[results[i].1[j]] = true;
                }
            }
        }
    }
    let mut dups = 0;
    for i in 0..dup.len() {
        if dup[i] {
            dups += 1;
        }
    }
    println!("dups = {}", dups);
}
