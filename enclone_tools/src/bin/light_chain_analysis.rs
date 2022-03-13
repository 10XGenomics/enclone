// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Analyze light chains.  Supply a single file of data, with one line per cell, and fields
// donors,v_name1,v_name2,dref,cdr3_aa1.

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
    // data = {(v_name1, cdr3_aa1, donor, v_name2, dref)}
    let mut data = Vec::<(String, String, String, String, usize)>::new();
    for line in f.lines() {
        let s = line.unwrap();
        if s.len() == 0 || s.starts_with("enclone") {
            continue;
        }
        let fields = s.split(',').collect::<Vec<&str>>();
        if first {
            for i in 0..fields.len() {
                tof.insert(fields[i].to_string(), i);
            }
            assert!(tof.contains_key("donors"));
            assert!(tof.contains_key("v_name1"));
            assert!(tof.contains_key("v_name2"));
            assert!(tof.contains_key("dref"));
            assert!(tof.contains_key("cdr3_aa1"));
            first = false;
        } else {
            data.push((
                fields[tof["v_name1"]].to_string(),
                fields[tof["cdr3_aa1"]].to_string(),
                fields[tof["donors"]].to_string(),
                fields[tof["v_name2"]].to_string(),
                fields[tof["dref"]].force_usize(),
            ));
        }
    }
    data.sort();
    let mut all = 0;
    let mut same = 0;
    let mut i = 0;
    while i < data.len() {
        // let j = next_diff12_5(&data, i as i32) as usize;
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
                    if dref1 >= 40 && dref2 >= 40 {
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
