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
use string_utils::TextUtils;
use tables::*;

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
            first = false;
        } else {
            data.push((
                fields[tof["donors_cell"]].to_string(),
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
        data[i].3 = data[i].3.replace("D", "");
    }

    // Analyze the data.

    let mut i = 0;
    let mut all = 0;
    let mut same = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }

        // Now i..j define a group having the same donor and same CDR3H-AA.  Look at pairs of 
        // cells from this group.

        for k1 in i..j {
            for k2 in i + 1..j {
                if data[k1].2 == data[k2].2 { 
                    continue;
                }
                if data[k1].4 == 0 || data[k2].4 == 0 {
                    continue;
                }
                all += 1;
                if data[k1].3 == data[k2].3 {
                    println!("{}, {}", data[k1].2, data[k2].2);
                    same += 1;
                }
            }
        }
        i = j;
    }
    println!("{:.1}%", 100.0 * same as f64 / all as f64);
    println!("{} of {}", same, all);
}
