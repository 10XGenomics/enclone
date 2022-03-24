// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// See also private_light_chain_analysis.rs.
//
// Analyze light chains.  Supply a single file of data, with one line per cell, and fields
// including donors_cell,v_name1,v_name2,dref,cdr3_aa1.
//
// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1 > per_cell_stuff

use io_utils::*;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use string_utils::{stringme, TextUtils};

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();
    let mut data = Vec::<(String, usize, Vec<u8>, String, String, usize)>::new();
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
            first = false;
        } else {
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
    data.sort();

    // Replace paralogs.

    for i in 0..data.len() {
        data[i].4 = data[i].4.replace("D", "");
    }

    // Define penalty matrix.

    let f = include_bytes!["../../../enclone_paper/data/mat.575142"].to_vec();
    let f = stringme(&f);
    let mut m = Vec::<Vec<f64>>::new();
    for line in f.lines() {
        let mut s = line.to_string();
        let sb = s.replace(" ", "");
        if sb == "ACDEFGHIKLMNPQRSTVWY" {
            continue;
        }
        if s.len() > 2 && s.as_bytes()[0] >= b'A' {
            s = s[2..].to_string();
        }
        let fields = s.split(' ').collect::<Vec<&str>>();
        let mut row = Vec::<f64>::new();
        for i in 0..fields.len() {
            row.push(fields[i].force_f64());
        }
        m.push(row);
    }
    let mut penalty = vec![vec![0.0; 27]; 27];
    let aa = b"ACDEFGHIKLMNPQRSTVWY".to_vec();
    for i1 in 0..aa.len() {
        let c1 = (aa[i1] - b'A') as usize;
        for i2 in 0..aa.len() {
            let c2 = (aa[i2] - b'A') as usize;
            penalty[c1][c2] = m[i1][i2];
        }
    }
    let mut penalty_simple = vec![vec![1.0; 27]; 27];
    for i in 0..27 {
        penalty_simple[i][i] = 0.0;
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
        bounds.push((i, j, vec![vec![(0, 0, 0, 0); 21]; 7]));
        i = j;
    }

    // Results = for each percent identity, rounded down:
    // 3. count for equal light chain gene names and dref1 > 0 and dref2 > 0
    // 4. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.

    bounds.par_iter_mut().for_each(|res| {
        let i = res.0;
        let j = res.1;
        for k1 in i..j {
            for k2 in k1 + 1..j {
                if data[k1].3 == data[k2].3 {
                    continue;
                }
                for pass in 0..2 {
                    let penalty = if pass == 0 { &penalty_simple } else { &penalty };
                    let mut err = 0.0;
                    for m in 0..data[k1].2.len() {
                        let c1 = (data[k1].2[m] - b'A') as usize;
                        let c2 = (data[k2].2[m] - b'A') as usize;
                        err += penalty[c1][c2];
                    }
                    let same = data[k1].2.len() as f64 - err;
                    let ident = 100.0 * same / data[k1].2.len() as f64;
                    let ident = ident.floor() as usize;
                    let ident = ident / 5;
                    let (dref1, dref2) = (data[k1].5, data[k2].5);
                    let eq_light = data[k1].4 == data[k2].4;
                    if dref1 > 0 && dref2 > 0 {
                        for u in 0..=ident {
                            if eq_light {
                                res.2[pass][u].2 += 1;
                            } else {
                                res.2[pass][u].3 += 1;
                            }
                        }
                    }
                }
            }
        }
    });

    // Sum.

    let mut res = vec![vec![(0, 0, 0, 0); 21]; 7];
    for pass in 0..2 {
        for i in 0..bounds.len() {
            for j in 0..=20 {
                res[pass][j].0 += bounds[i].2[pass][j].0;
                res[pass][j].1 += bounds[i].2[pass][j].1;
                res[pass][j].2 += bounds[i].2[pass][j].2;
                res[pass][j].3 += bounds[i].2[pass][j].3;
            }
        }
    }

    // Print results.

    for pass in 0..2 {
        if pass == 0 {
            println!("\nsimple transition matrix\n");
        } else {
            println!("\nevolved transition matrix\n");
        }
        for j in 14..=20 {
            let n = res[pass][j].2 + res[pass][j].3;
            let nznz = 100.0 * res[pass][j].2 as f64 / n as f64;
            println!("{}% ==> {nznz:.1}% of {n}", 5 * j);
        }
    }
    println!("");
}
