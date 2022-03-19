// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Optimize amino acid penalty matrix to increase joining given light chain coherence of at
// least 70%.

// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1,hcomp
//         > per_cell_stuff

use io_utils::*;
use perf_stats::elapsed;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use std::time::Instant;
use string_utils::TextUtils;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();
    let mut data = Vec::<(
        String,
        usize,
        Vec<u8>,
        String,
        String,
        usize,
    )>::new();
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("#") {
            continue;
        }
        let fields = s.split(',').collect::<Vec<&str>>();
        if first {
            for i in 0..fields.len() {
                tof.insert(fields[i].to_string(), i);
            }
            assert!(tof.contains_key("donors_cell"));
            assert!(tof.contains_key("v_name1"));
            assert!(tof.contains_key("v_name2"));
            assert!(tof.contains_key("dref"));
            assert!(tof.contains_key("cdr3_aa1"));
            first = false;
        } else {
            let dref = fields[tof["dref"]].force_usize();
            if dref > 0 {
                data.push((
                    fields[tof["v_name1"]].to_string(),
                    fields[tof["cdr3_aa1"]].len(),
                    fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                    fields[tof["donors_cell"]].to_string(),
                    fields[tof["v_name2"]].to_string(),
                    fields[tof["dref"]].force_usize(),
                ));
            }
        }
    }
    data.sort();

    // Replace paralogs.

    for i in 0..data.len() {
        data[i].4 = data[i].4.replace("D", "");
    }

    // Convert CDRH3 to [0,20).

    let aa = b"ACDEFGHIKLMNPQRSTVWY".to_vec();
    for i in 0..data.len() {
        let mut x = Vec::<u8>::new();
        for j in 0..data[i].2.len() {
            let c = data[i].2[j];
            let p = bin_position(&aa, &c) as u8;
            x.push(p);
        }
        data[i].2 = x;
    }
    
    // Define penalty matrix.

    let mut penalty = vec![vec![1.0; 20]; 20];
    for i in 0..20 {
        penalty[i][i] = 0.0;
    }

    // Define groups based on equal heavy chain gene names and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::<(usize, usize, (usize, usize))>::new();
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        bounds.push((i, j, (0, 0)));
        i = j;
    }

    // Results = for each percent identity, rounded down:
    // 1. count for equal light chain gene names and dref1 > 0 and dref2 > 0
    // 2. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.

    let t = Instant::now();
    bounds.par_iter_mut().for_each(|res| {
        let (i, j) = (res.0, res.1);
        for k1 in i..j {
            for k2 in k1 + 1..j {

                // Require different donors.

                if data[k1].3 == data[k2].3 {
                    continue;
                }

                // Compute stuff.

                let mut err = 0.0;
                for m in 0..data[k1].2.len() {
                    let c1 = data[k1].2[m] as usize;
                    let c2 = data[k2].2[m] as usize;
                    err += penalty[c1][c2];
                }
                err /= data[k1].2.len() as f64;
                err *= 100.0;
                let eq_light = data[k1].4 == data[k2].4;

                // Add to results.

                if err == 0.0 {
                    if eq_light {
                        res.2.0 += 1;
                    } else {
                        res.2.1 += 1;
                    }
                }
            }
        }
    });

    // Sum.

    let mut res = (0, 0);
    for i in 0..bounds.len() {
        res.0 += bounds[i].2.0;
        res.1 += bounds[i].2.1;
    }

    // Print.

    let n = res.0 + res.1;
    let nznz = 100.0 * res.0 as f64 / n as f64;
    println!("\n{nznz:.1}%");
    println!("used {:.1} seconds\n", elapsed(&t));
}
