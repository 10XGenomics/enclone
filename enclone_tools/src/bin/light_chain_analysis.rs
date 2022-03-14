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
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use string_utils::TextUtils;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    let mut first = true;
    let mut tof = HashMap::<String, usize>::new();
    // data = {(v_name1, cdr3_aa1, donor, v_name2, dref, clonotype_ncells, const1, hcomp)}
    let mut data = Vec::<(String, String, String, String, usize, usize, String, usize)>::new();
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
                fields[tof["cdr3_aa1"]].to_string(),
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
    let mut all1 = 0;
    let mut same1 = 0;
    let mut all2 = 0;
    let mut same2 = 0;
    let mut mall1 = vec![vec![0; 4]; 4];
    let mut mall2 = vec![vec![0; 4]; 4];
    let mut msame1 = vec![vec![0; 4]; 4];
    let mut msame2 = vec![vec![0; 4]; 4];
    let mut total_hcomp1 = 0;
    let mut total_hcomp2 = 0;
    let mut i = 0;
    while i < data.len() {
        // let j = next_diff12_8(&data, i as i32) as usize;
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        for k1 in i..j {
            for k2 in i + 1..j {
                if data[k1].2 != data[k2].2 {
                    let lmatch = data[k1].3 == data[k2].3;
                    let (dref1, dref2) = (data[k1].4, data[k2].4);
                    // let (cncells1, cncells2) = (data[k1].5, data[k2].5);
                    // let (isotype1, isotype2) = (&data[k1].6, &data[k2].6);
                    let (hcomp1, hcomp2) = (data[k1].7, data[k2].7);
                    if dref1 == 0 || dref2 == 0 {
                        all1 += 1;
                        total_hcomp1 += hcomp1 + hcomp2;
                        if lmatch {
                            same1 += 1;
                        }
                    } else {
                        all2 += 1;
                        total_hcomp2 += hcomp1 + hcomp2;
                        if lmatch {
                            same2 += 1;
                        }
                    }
                    for z1 in 0..4 {
                        for z2 in z1 + 1..4 {
                            if data[k1].2 == format!("d{}", z1 + 1)
                                || data[k2].2 == format!("d{}", z1 + 1)
                            {
                                if data[k1].2 == format!("d{}", z2 + 1)
                                    || data[k2].2 == format!("d{}", z2 + 1)
                                {
                                    if dref1 == 0 || dref2 == 0 {
                                        mall1[z1][z2] += 1;
                                        if lmatch {
                                            msame1[z1][z2] += 1;
                                        }
                                    } else {
                                        mall2[z1][z2] += 1;
                                        if lmatch {
                                            msame2[z1][z2] += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        i = j;
    }
    println!(
        "interdonor light chain concordance for either naive = {:.1}% = {} of {}",
        100.0 * same1 as f64 / all1 as f64,
        same1,
        all1
    );
    println!(
        "for these cells, mean hcomp = {:.1}",
        total_hcomp1 as f64 / (all1 as f64 * 2.0)
    );
    println!(
        "interdonor light chain concordance for neither naive = {:.1}% = {} of {}",
        100.0 * same2 as f64 / all2 as f64,
        same2,
        all2
    );
    println!(
        "for these cells, mean hcomp = {:.1}",
        total_hcomp2 as f64 / (all2 as f64 * 2.0)
    );

    // The following doesn't make sense for the BIB=@training.

    for z1 in 0..4 {
        for z2 in z1 + 1..4 {
            println!("\nd{} versus d{}", z1 + 1, z2 + 1);
            println!(
                "interdonor light chain concordance for either naive = {:.1}% = {} of {}",
                100.0 * msame1[z1][z2] as f64 / mall1[z1][z2] as f64,
                msame1[z1][z2],
                mall1[z1][z2],
            );
            println!(
                "interdonor light chain concordance for neither naive = {:.1}% = {} of {}",
                100.0 * msame2[z1][z2] as f64 / mall2[z1][z2] as f64,
                msame2[z1][z2],
                mall2[z1][z2],
            );
        }
    }
}
