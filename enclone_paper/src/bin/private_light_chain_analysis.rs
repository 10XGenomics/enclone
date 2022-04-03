// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// See also public_light_chain_analysis.rs.
//
// Analyze light chains.  Supply a single file of data, with one line per cell, and fields
// including donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1.
//
// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1,hcomp
//         > per_cell_stuff

use enclone_core::hcat;
use io_utils::*;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use string_utils::TextUtils;
use tables::print_tabular_vbox;

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
                fields[tof["donors_cell"]].to_string(),
                fields[tof["cdr3_aa1"]].len(),
                fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                fields[tof["v_name1"]].to_string(),
                fields[tof["v_name2"]].to_string(),
                fields[tof["dref"]].force_usize(),
            ));
        }
    }
    data.sort();

    // Replace paralogous gene names.

    for i in 0..data.len() {
        data[i].4 = data[i].4.replace("D", "");
    }

    // Define groups based on equal donor and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::<(usize, usize, Vec<Vec<(usize, usize, usize, usize)>>)>::new();
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        bounds.push((i, j, vec![vec![(0, 0, 0, 0); 11]; 5]));
        i = j;
    }

    // Results = for each percent identity, rounded down:
    // 1. count for equal light chain gene names and dref1 = 0 and dref2 = 0
    // 2. count for unequal light chain gene names and dref1 = 0 and dref2 = 0
    // 3. count for equal light chain gene names and dref1 > 0 and dref2 > 0
    // 4. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.
    //
    // Make one pass for each donor.

    bounds.par_iter_mut().for_each(|res| {
        let i = res.0;
        let j = res.1;
        let d = data[i].0.after("d").force_usize() - 1;
        for k1 in i..j {
            for k2 in k1 + 1..j {
                // Require different heavy chain V genes.

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
                let (dref1, dref2) = (data[k1].5, data[k2].5);
                let eq_light = data[k1].4 == data[k2].4;

                if dref1 == 0 && dref2 == 0 {
                    if eq_light {
                        res.2[d + 1][ident].0 += 1;
                        res.2[0][ident].0 += 1;
                    } else {
                        res.2[d + 1][ident].1 += 1;
                        res.2[0][ident].1 += 1;
                    }
                } else if dref1 > 0 && dref2 > 0 {
                    if eq_light {
                        res.2[d + 1][ident].2 += 1;
                        res.2[0][ident].2 += 1;
                    } else {
                        res.2[d + 1][ident].3 += 1;
                        res.2[0][ident].3 += 1;
                    }
                }
            }
        }
    });

    // Sum.

    let mut res = vec![vec![(0, 0, 0, 0); 11]; 5];
    for pass in 0..5 {
        for i in 0..bounds.len() {
            for j in 0..=10 {
                res[pass][j].0 += bounds[i].2[pass][j].0;
                res[pass][j].1 += bounds[i].2[pass][j].1;
                res[pass][j].2 += bounds[i].2[pass][j].2;
                res[pass][j].3 += bounds[i].2[pass][j].3;
            }
        }
    }

    // Print results.

    println!(
        "\nConsider two cells from the same donor that have the same CDR3H length,\n\
            and different heavy chain V genes."
    );
    println!("\nColumn 1: percent identity rounded down to nearest ten percent");
    println!("Column > 1: probability that light chain gene names are the same");
    let mut logs = Vec::<String>::new();
    for xpass in 1..=2 {
        let mut log = String::new();
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "CDR3H-AA".to_string(),
            "all".to_string(),
            "d1".to_string(),
            "d2".to_string(),
            "d3".to_string(),
            "d4".to_string(),
        ];
        rows.push(row);
        for j in 0..=10 {
            let row = vec!["\\hline".to_string(); 6];
            rows.push(row);
            let mut row = vec![format!("{}%", 10 * j)];
            for pass in 0..5 {
                if xpass == 1 {
                    let n = res[pass][j].2 + res[pass][j].3;
                    let nznz = 100.0 * res[pass][j].2 as f64 / n as f64;
                    row.push(format!("{nznz:.1}%"));
                } else {
                    let n = res[pass][j].0 + res[pass][j].1;
                    let nznz = 100.0 * res[pass][j].0 as f64 / n as f64;
                    row.push(format!("{nznz:.1}%"));
                }
            }
            rows.push(row);
        }
        print_tabular_vbox(&mut log, &rows, 0, &b"l|r|r|r|r|r".to_vec(), false, false);
        logs.push(log);
    }
    let mut logr = vec![Vec::<String>::new(); 2];
    for xpass in 0..2 {
        let r = logs[xpass].split('\n').map(str::to_owned).collect();
        logr[xpass] = r;
    }
    print!("\n both cells have dref > 0");
    print!("                               ");
    println!("both cells have dref = 0");
    let r = hcat(&logr[0], &logr[1], 3);
    for i in 0..r.len() {
        println!("{}", r[i]);
    }
}
