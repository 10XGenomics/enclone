// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Analyze light chains.  Supply a single file of data, with one line per cell, and fields
// including donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells.

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
    // data = {(v_name1, cdr3_aa1, donor, v_name2, dref, clonotype_ncells)}
    let mut data = Vec::<(String, String, String, String, usize, usize)>::new();
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
            first = false;
        } else {
            data.push((
                fields[tof["v_name1"]].to_string(),
                fields[tof["cdr3_aa1"]].to_string(),
                fields[tof["donors_cell"]].to_string(),
                fields[tof["v_name2"]].to_string(),
                fields[tof["dref"]].force_usize(),
                fields[tof["clonotype_ncells"]].force_usize(),
            ));
        }
    }
    for i in 0..data.len() {
        if data[i].3.contains("H") {
            println!("{}, {}, {}, {}, {}, {}",
                data[i].0,
                data[i].1,
                data[i].2,
                data[i].3,
                data[i].4,
                data[i].5,
            );
        }
    }
    data.sort();
    let mut all = 0;
    let mut same = 0;
    let mut i = 0;
    while i < data.len() {
        // let j = next_diff12_6(&data, i as i32) as usize;
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
                    let (cncells1, cncells2) = (data[k1].5, data[k2].5);
                    // if dref1 >= 30 && dref2 >= 30 {
                    // if dref1 == 0 && dref2 == 0 {
                    // if cncells1 >= 10 && cncells2 >= 10 {
                    // if cncells1 == 1 && cncells2 == 1 && dref1 >= 20 && dref2 >= 20 { // 68.1%
                    // if cncells1 == 1 && cncells2 == 1 && dref1 >= 10 && dref2 >= 10 { // 71.6%
                    // if cncells1 == 1 && cncells2 == 1 && dref1 == 0 && dref2 == 0 { // 8.7%
                    println!("L = {}/{}, csize = {}/{}, dref = {}/{}",
                        data[k1].3, data[k2].3,
                        data[k1].5, data[k2].5,
                        data[k1].4, data[k2].4,
                    );
                    if true {
                        all += 1;
                        if lmatch {
                            same += 1;
                        }
                    }
                }
            }
        }
        i = j;
    }
    println!("interdonor light chain concordance = {:.1}% = {} of {}",
        100.0 * same as f64 / all as f64,
        same,
        all
    );
}
