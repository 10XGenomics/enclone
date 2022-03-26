// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// See also private_light_chain_analysis.rs.
//
// Analyze light chains.  Supply a single file of data, with one line per cell, and fields
// including donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1.
//
// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1,hcomp,jun_ins
//         > per_cell_stuff
//
// enclone BIB=@training BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1,hcomp
//         > training_per_cell_stuff
//
// enclone BIB=1,2,3,29 BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1,hcomp
//         > training1_per_cell_stuff

use io_utils::*;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::cmp::max;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use std::mem::swap;
use string_utils::{stringme, strme, TextUtils};
use tables::*;
use vector_utils::unique_sort;

pub fn hcat(col1: &[String], col2: &[String], sep: usize) -> Vec<String> {
    let mut cat = Vec::<String>::new();
    let height = max(col1.len(), col2.len());
    let mut width1 = 0;
    for x in col1 {
        width1 = max(width1, x.chars().count() + sep);
    }
    for i in 0..height {
        let mut s = if i < col1.len() {
            col1[i].clone()
        } else {
            String::new()
        };
        while s.chars().count() < width1 {
            s += " ";
        }
        if i < col2.len() {
            s += &col2[i];
        }
        cat.push(s);
    }
    cat
}

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

    // Compute jun_ins frequency for memory and naive cells.

    let mut ins = vec![vec![0; 1000]; 2];
    let mut total = vec![0; 2];
    for pass in 0..2 {
        for i in 0..data.len() {
            let dref = data[i].5;
            let jun_ins = data[i].9;
            if (pass == 0 && dref > 0) || (pass == 1 && dref == 0) {
                ins[pass][jun_ins] += 1;
                total[pass] += 1;
            }
        }
    }
    println!("\njunction insertion length frequency for memory cells");
    for i in 0..=20 {
        if ins[0][i] > 0 {
            println!("{i}: {:.3}%", 100.0 * ins[0][i] as f64 / total[0] as f64);
        }
    }
    println!("\njunction insertion length frequency for naive cells");
    for i in 0..=20 {
        if ins[1][i] > 0 {
            println!("{i}: {:.3}%", 100.0 * ins[1][i] as f64 / total[1] as f64);
        }
    }

    // Compute mean jun_ins for various classes of cells.

    let mut n_memory = 0;
    let mut n_naive = 0;
    let mut ins_memory = 0;
    let mut ins_naive = 0;
    for i in 0..data.len() {
        let dref = data[i].5;
        let jun_ins = data[i].9;
        if dref == 0 {
            n_naive += 1;
            ins_naive += jun_ins
        } else {
            n_memory += 1;
            ins_memory += jun_ins
        }
    }
    println!(
        "\nmean junction insertion bases for memory = {:.1}",
        ins_memory as f64 / n_memory as f64
    );
    println!(
        "mean junction insertion bases for naive = {:.1}",
        ins_naive as f64 / n_naive as f64
    );
    let mut n_memory = 0;
    let mut n_naive = 0;
    let mut ins_memory = 0;
    let mut ins_naive = 0;
    let mut max_ins_memory = 0;
    let mut max_ins_memory_cdr3 = String::new();
    let mut max_ins_naive = 0;
    let mut i = 0;
    let mut last = String::new();
    let mut ins_mem = vec![0; 100];
    let mut ins = vec![vec![0; 1000]; 2];
    let mut total = vec![0; 2];
    println!("");
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].2 != data[i].2 {
                break;
            }
            j += 1;
        }
        let mut donors = Vec::<String>::new();
        for k in i..j {
            donors.push(data[k].3.clone());
        }
        unique_sort(&mut donors);
        if donors.len() > 1 {
            for k in i..j {
                let dref = data[k].5;
                let jun_ins = data[k].9;
                if dref == 0 {
                    n_naive += 1;
                    ins_naive += jun_ins;
                    max_ins_naive = max(max_ins_naive, jun_ins);
                    ins[1][jun_ins] += 1;
                    total[1] += 1;
                } else {
                    n_memory += 1;
                    ins_memory += jun_ins;
                    if k == i {
                        ins_mem[jun_ins] += 1;
                    }
                    if jun_ins >= 6 && last != strme(&data[k].2) {
                        println!(
                            "public memory with CDR3H = {} has jun_ins = {}",
                            stringme(&data[k].2),
                            jun_ins
                        );
                        last = stringme(&data[k].2);
                    }
                    if jun_ins > max_ins_memory {
                        max_ins_memory_cdr3 = stringme(&data[k].2);
                    }
                    max_ins_memory = max(max_ins_memory, jun_ins);
                    ins[0][jun_ins] += 1;
                    total[0] += 1;
                }
            }
        }
        i = j;
    }
    println!("\njunction insertion length frequency for public memory cells");
    for i in 0..=20 {
        if ins[0][i] > 0 {
            println!("{i}: {:.3}%", 100.0 * ins[0][i] as f64 / total[0] as f64);
        }
    }
    println!("\njunction insertion length frequency for public naive cells");
    for i in 0..=20 {
        if ins[1][i] > 0 {
            println!("{i}: {:.3}%", 100.0 * ins[1][i] as f64 / total[1] as f64);
        }
    }

    println!(
        "mean junction insertion bases for public memory = {:.1}",
        ins_memory as f64 / n_memory as f64
    );
    println!(
        "max = {} at CDR3H = {}",
        max_ins_memory, max_ins_memory_cdr3
    );
    println!(
        "mean junction insertion bases for public naive = {:.1}",
        ins_naive as f64 / n_naive as f64
    );
    println!("max = {}", max_ins_naive);
    println!("\npublic memory junction insertion length distribution, by cell group");
    for i in 0..ins_mem.len() {
        if ins_mem[i] > 0 {
            println!("{} ==> {}", i, ins_mem[i]);
        }
    }

    // Compute light chain coherence for memory cells as a function of insertion length.

    let mut n = vec![0; 100];
    let mut same = vec![0; 100];
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].2 != data[i].2 {
                break;
            }
            j += 1;
        }
        for k1 in i..j {
            for k2 in k1 + 1..j {
                if data[k1].5 > 0 && data[k2].5 > 0 {
                    if data[k1].3 != data[k2].3 {
                        let ins = data[k1].9;
                        if ins == data[k2].9 {
                            n[ins] += 1;
                            if data[k1].4 == data[k2].4 {
                                same[ins] += 1;
                            }
                        }
                    }
                }
            }
        }
        i = j;
    }
    println!(
        "\nlight chain coherence for public memory cells as function of junction insertion length"
    );
    for i in 0..=4 {
        if n[i] > 0 {
            println!(
                "ins = {}  ==> coherence = {:.1}% ({} of {})",
                i,
                100.0 * same[i] as f64 / n[i] as f64,
                same[i],
                n[i],
            );
        }
    }
    let mut n_big = 0;
    let mut same_big = 0;
    for i in 5..n.len() {
        n_big += n[i];
        same_big += same[i];
    }
    println!(
        "ins >= 5 ==> coherence = {:.1}% ({} of {})",
        100.0 * same_big as f64 / n_big as f64,
        same_big,
        n_big,
    );

    // Define groups based on equal heavy chain gene names and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::<(usize, usize, Vec<Vec<(usize, usize, usize, usize)>>)>::new();
    let mut bounds2 = Vec::<(usize, usize, Vec<usize>, Vec<usize>)>::new();
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
        bounds2.push((i, j, vec![0; 100], vec![0; 100]));
        i = j;
    }

    // Compute light chain coherence for memory cells as a function of insertion length, as
    // earlier but now using 90% CDRH3 identity.

    bounds2.par_iter_mut().for_each(|res| {
        let i = res.0;
        let j = res.1;
        for k1 in i..j {
            for k2 in k1 + 1..j {
                let mut samex = 0;
                for m in 0..data[k1].2.len() {
                    if data[k1].2[m] == data[k2].2[m] {
                        samex += 1;
                    }
                }
                let ident = 100.0 * samex as f64 / data[k1].2.len() as f64;
                let ident = ident.floor() as usize;
                if ident >= 90 {
                    if data[k1].5 > 0 && data[k2].5 > 0 {
                        if data[k1].3 != data[k2].3 {
                            let ins = data[k1].9;
                            if ins == data[k2].9 {
                                res.3[ins] += 1;
                                if data[k1].4 == data[k2].4 {
                                    res.2[ins] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    });
    let mut n = vec![0; 100];
    let mut same = vec![0; 100];
    for i in 0..bounds2.len() {
        for j in 0..100 {
            same[j] += bounds2[i].2[j];
            n[j] += bounds2[i].3[j];
        }
    }
    println!(
        "\nlight chain coherence for public (>= 90%) memory cells as function of junction \
        insertion length"
    );
    for i in 0..=9 {
        if n[i] > 0 {
            println!(
                "ins = {}  ==> coherence = {:.1}% ({} of {})",
                i,
                100.0 * same[i] as f64 / n[i] as f64,
                same[i],
                n[i],
            );
        }
    }
    let mut n_big = 0;
    let mut same_big = 0;
    for i in 10..n.len() {
        n_big += n[i];
        same_big += same[i];
    }
    println!(
        "ins >= 10 ==> coherence = {:.1}% ({} of {})",
        100.0 * same_big as f64 / n_big as f64,
        same_big,
        n_big,
    );

    // Results = for each percent identity, rounded down:
    // 1. count for equal light chain gene names and dref1 = 0 and dref2 = 0
    // 2. count for unequal light chain gene names and dref1 = 0 and dref2 = 0
    // 3. count for equal light chain gene names and dref1 > 0 and dref2 > 0
    // 4. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.
    //
    // Make one pass for all donors, and one pass each for each pair of donors.

    bounds.par_iter_mut().for_each(|res| {
        let i = res.0;
        let j = res.1;
        for k1 in i..j {
            for k2 in k1 + 1..j {
                // Require different donors.

                if data[k1].3 == data[k2].3 {
                    continue;
                }
                let (mut d1, mut d2) = (data[k1].3.clone(), data[k2].3.clone());
                if d1 > d2 {
                    swap(&mut d1, &mut d2);
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

                // Go through passes.

                for pass in 0..7 {
                    // Require specific donors.

                    if pass == 0 {
                    } else if pass == 1 {
                        if d1 != "d1" || d2 != "d2" {
                            continue;
                        }
                    } else if pass == 2 {
                        if d1 != "d1" || d2 != "d3" {
                            continue;
                        }
                    } else if pass == 3 {
                        if d1 != "d1" || d2 != "d4" {
                            continue;
                        }
                    } else if pass == 4 {
                        if d1 != "d2" || d2 != "d3" {
                            continue;
                        }
                    } else if pass == 5 {
                        if d1 != "d2" || d2 != "d4" {
                            continue;
                        }
                    } else {
                        if d1 != "d3" || d2 != "d4" {
                            continue;
                        }
                    }

                    // Add to results.

                    if dref1 == 0 && dref2 == 0 {
                        if eq_light {
                            res.2[pass][ident].0 += 1;
                        } else {
                            res.2[pass][ident].1 += 1;
                        }
                    } else if dref1 > 0 && dref2 > 0 {
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

    // Print results.

    println!(
        "\nConsider two cells from different donors that have the same heavy chain gene name \
        and CDR3H length."
    );
    println!("\nColumn 1: percent identity rounded down to nearest ten percent");
    println!("Column > 1: probability that light chain gene names are the same");
    let mut logs = Vec::<String>::new();
    for xpass in 1..=2 {
        let mut log = String::new();
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "CDR3H-AA".to_string(),
            "any".to_string(),
            "d1,d2".to_string(),
            "d1,d3".to_string(),
            "d1,d4".to_string(),
            "d2,d3".to_string(),
            "d2,d4".to_string(),
            "d3,d4".to_string(),
        ];
        rows.push(row);
        for j in 0..=10 {
            let row = vec!["\\hline".to_string(); 8];
            rows.push(row);
            let mut row = vec![format!("{}%", 10 * j)];
            for pass in 0..7 {
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
        print_tabular_vbox(
            &mut log,
            &rows,
            0,
            &b"l|r|r|r|r|r|r|r".to_vec(),
            false,
            false,
        );
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
