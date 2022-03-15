// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Analyze light chains.  Supply a single file of data, with one line per cell, and fields
// including donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1.
//
// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1,hcomp
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
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use std::mem::swap;
use string_utils::TextUtils;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();
    // data = {(v_name1, cdr3_len, cdr3_aa1, donor, v_name2, dref, clonotype_ncells, const1, hcomp)}
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
            ));
        }
    }
    data.sort();

    // Define groups based on equal heavy chain gene names and CDR3H length.
    // Plus placeholder for results, see next.

    let mut bounds = Vec::<(usize, usize, Vec<(usize, usize, usize, usize)>)>::new();
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
        bounds.push((i, j, vec![(0, 0, 0, 0); 11]));
        i = j;
    }

    // Results = for each percent identity, rounded down:
    // 1. count for equal light chain gene names and dref1 = 0 and dref2 = 0
    // 2. count for unequal light chain gene names and dref1 = 0 and dref2 = 0
    // 3. count for equal light chain gene names and dref1 > 0 and dref2 > 0
    // 4. count for unequal light chain gene names and dref1 > 0 and dref2 > 0.

    // Make one pass for all donors, and one pass each for each pair of donors.
    // This is seven times slower than it should be.

    for pass in 0..7 {
        for i in 0..bounds.len() {
            bounds[i].2 = vec![(0, 0, 0, 0); 11];
        }
        if pass == 0 {
            println!("\nany two donors");
        } else if pass == 1 {
            println!("donors 1 and 2");
        } else if pass == 2 {
            println!("donors 1 and 3");
        } else if pass == 3 {
            println!("donors 1 and 4");
        } else if pass == 4 {
            println!("donors 2 and 3");
        } else if pass == 5 {
            println!("donors 2 and 4");
        } else {
            println!("donors 3 and 4");
        }
        let mut res = vec![(0, 0, 0, 0); 11];
        bounds.par_iter_mut().for_each(|res| {
            let i = res.0;
            let j = res.1;
            for k1 in i..j {
                for k2 in k1 + 1..j {

                    // Require different donors.

                    let (mut d1, mut d2) = (data[k1].3.clone(), data[k2].3.clone());
                    if d1 > d2 {
                        swap(&mut d1, &mut d2);
                    }
                    if pass == 0 {
                        if d1 == d2 {
                            continue;
                        }
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
                            res.2[ident].0 += 1;
                        } else {
                            res.2[ident].1 += 1;
                        }
                    } else if dref1 > 0 && dref2 > 0 {
                        if eq_light {
                            res.2[ident].2 += 1;
                        } else {
                            res.2[ident].3 += 1;
                        }
                    }
                }
            }
        });
        for i in 0..bounds.len() {
            for j in 0..=10 {
                res[j].0 += bounds[i].2[j].0;
                res[j].1 += bounds[i].2[j].1;
                res[j].2 += bounds[i].2[j].2;
                res[j].3 += bounds[i].2[j].3;
            }
        }

        // Print results.

        println!(
            "\nConsider two cells from different donors that have the same heavy chain gene name\n\
            and CDR3H length."
        );
        println!(
            "\npercent identity on CDR3H-AA ==> \
            probability that light chain genes are the same"
        );
        println!("\npercent identity rounded down to nearest ten percent");
        println!("\nboth cells have dref > 0:\n");
        for j in 0..=10 {
            let n = res[j].2 + res[j].3;
            let nznz = 100.0 * res[j].2 as f64 / n as f64;
            let m = 10 * j;
            println!("{m}% ==> {nznz:.1}% (n = {n})");
        }
        println!("\nboth cells have dref == 0:\n");
        for j in 0..=10 {
            let n = res[j].0 + res[j].1;
            let nznz = 100.0 * res[j].0 as f64 / n as f64;
            let m = 10 * j;
            println!("{m}% ==> {nznz:.1}% (n = {n})");
        }
        println!("");
    }
}
