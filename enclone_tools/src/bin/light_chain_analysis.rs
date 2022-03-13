// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Analyze light chains.  Supply a single file of data, with one line per cell, and fields
// including donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1.
//
// Data from:
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1 
//         > per_cell_stuff
//
// enclone BIB=@training BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,const1 
//         > training_per_cell_stuff

use pretty_trace::PrettyTrace;
use io_utils::*;
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
    // data = {(v_name1, cdr3_aa1, donor, v_name2, dref, clonotype_ncells, const1)}
    let mut data = Vec::<(String, String, String, String, usize, usize, String)>::new();
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
            ));
        }
    }
    data.sort();
    let mut all1 = 0;
    let mut same1 = 0;
    let mut all2 = 0;
    let mut same2 = 0;
    let mut i = 0;
    while i < data.len() {
        // let j = next_diff12_7(&data, i as i32) as usize;
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        for k1 in i..j {
            for k2 in i+1..j {
                if data[k1].2 != data[k2].2 {
                    let lmatch = data[k1].3 == data[k2].3;
                    let (dref1, dref2) = (data[k1].4, data[k2].4);
                    // let (cncells1, cncells2) = (data[k1].5, data[k2].5);
                    // let (isotype1, isotype2) = (&data[k1].6, &data[k2].6);
                    if dref1 == 0 || dref2 == 0 {
                        all1 += 1;
                        if lmatch {
                            same1 += 1;
                        }
                    } else {
                        all2 += 1;
                        if lmatch {
                            same2 += 1;
                        }
                    }
                }
            }
        }
        i = j;
    }
    println!("interdonor light chain concordance for either naive = {:.1}% = {} of {}",
        100.0 * same1 as f64 / all1 as f64,
        same1,
        all1
    );
    println!("interdonor light chain concordance for neither naive = {:.1}% = {} of {}",
        100.0 * same2 as f64 / all2 as f64,
        same2,
        all2
    );
}
