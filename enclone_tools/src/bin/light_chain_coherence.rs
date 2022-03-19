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
        usize,
        String,
        usize,
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
            assert!(tof.contains_key("clonotype_ncells"));
            assert!(tof.contains_key("const1"));
            assert!(tof.contains_key("hcomp"));
            assert!(tof.contains_key("jun_ins"));
            first = false;
        } else {
            data.push((
                fields[tof["v_name1"]].to_string(),
                fields[tof["cdr3_aa1"]].len(),
                fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                fields[tof["donors_cell"]].to_string(),
                fields[tof["v_name2"]].to_string(),
                fields[tof["dref"]].force_usize(),
                fields[tof["clonotype_ncells"]].force_usize(),
                fields[tof["const1"]].to_string(),
                fields[tof["hcomp"]].force_usize(),
                fields[tof["jun_ins"]].force_usize(),
            ));
        }
    }
    data.sort();

    // Replace paralogs.

    for i in 0..data.len() {
        data[i].4 = data[i].4.replace("D", "");
    }

    // Define groups based on equal heavy chain gene names and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::<(usize, usize, Vec<Vec<(usize, usize, usize, usize)>>)>::new();
    let mut i = 0;
    while i < data.len() {
        // let j = next_diff12_9(&data, i as i32) as usize;
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        bounds.push((i, j, vec![vec![(0, 0, 0, 0); 11]; 7]));
        i = j;
    }

    // Results = for each percent identity, rounded down:
    // 1. count for equal light chain gene names and dref1 = 0 and dref2 = 0
    // 2. count for unequal light chain gene names and dref1 = 0 and dref2 = 0
    // 3. count for equal light chain gene names and dref1 > 0 and dref2 > 0
    // 4. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.
    //
    // Make one pass for all donors, and one pass each for each pair of donors.

    let t = Instant::now();
    bounds.par_iter_mut().for_each(|res| {
        let i = res.0;
        let j = res.1;
        for k1 in i..j {
            let dref1 = data[k1].5;
            if dref1 == 0 {
                continue;
            }
            for k2 in k1 + 1..j {
                let dref2 = data[k2].5;
                if dref2 == 0 {
                    continue;
                }

                // Require different donors.

                if data[k1].3 == data[k2].3 {
                    continue;
                }

                // Compute stuff.

                let mut same = 0;
                for m in 0..data[k1].2.len() {
                    if data[k1].2[m] == data[k2].2[m] {
                        same += 1;
                    }
                }
                let ident = 100.0 * same as f64 / data[k1].2.len() as f64;
                let ident = ident.floor() as usize;
                let ident = ident / 10;
                let eq_light = data[k1].4 == data[k2].4;

                // Go through passes.

                for pass in 0..1 {

                    // Add to results.

                    if dref1 > 0 && dref2 > 0 {
                        if eq_light {
                            res.2[pass][ident].2 += 1;
                        } else {
                            res.2[pass][ident].3 += 1;
                        }
                    }
                }
            }
        }
    });

    // Sum.

    let mut res = vec![vec![(0, 0, 0, 0); 11]; 7];
    for pass in 0..7 {
        for i in 0..bounds.len() {
            for j in 0..=10 {
                res[pass][j].0 += bounds[i].2[pass][j].0;
                res[pass][j].1 += bounds[i].2[pass][j].1;
                res[pass][j].2 += bounds[i].2[pass][j].2;
                res[pass][j].3 += bounds[i].2[pass][j].3;
            }
        }
    }

    // Print.

    let n = res[0][10].2 + res[0][10].3;
    let nznz = 100.0 * res[0][10].2 as f64 / n as f64;
    println!("\n{nznz:.1}%");
    println!("used {:.1} seconds", elapsed(&t));
}
