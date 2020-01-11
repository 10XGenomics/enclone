// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Given output of a clonotyping algorithm which took as inputs the pipeline outputs for the lenas 
// in enclone.testdata, find how many of the clonotypes cross donors.  The clonotyping algorithm
// should be a file with lines of the form:
// lena barcode clonotype_id.

/*
 
besp3% pwd
/mnt/home/wyatt.mcdonnell/immcantation_test/20191212_immcantation

besp3% cat clone | grep -v "^SEQ" | cut -f1,61 | tr '_' ' ' | tr '\t' ' ' | Col 1 2 5 | sort -u > /mnt/deck5/david.jaffe/imm/clone_summary

assess_clonotyping /mnt/deck5/david.jaffe/imm/clone_summary

669

The analogous number for enclone is 19.

enclone 123085 EXT=/mnt/deck5/david.jaffe/imm/clone_summary

---------------------------------------------------------------------------------------------------

cat /mnt/home/wyatt.mcdonnell/immcantation_test/123085/123085_ham/123085_ham_heavy_germ-pass.tab |  grep -v "^SEQ" | cut -f1,61 | sed 's/_/ /' | Col 1 3 | sed 's/^/123085 /' > /mnt/deck5/david.jaffe/imm/clone_123085_ham

enclone 123085 EXT=/mnt/deck5/david.jaffe/imm/clone_123085_ham

cat /mnt/home/wyatt.mcdonnell/immcantation_test/123085/123085_ham_aa/123085_ham_aa_heavy_germ-pass.tab |  grep -v "^SEQ" | cut -f1,61 | sed 's/_/ /' | Col 1 3 | sed 's/^/123085 /' > /mnt/deck5/david.jaffe/imm/clone_123085_ham_aa

cat /mnt/home/wyatt.mcdonnell/immcantation_test/123085/123085_hh_s1f/123085_hh_s1f_heavy_germ-pass.tab |  grep -v "^SEQ" | cut -f1,61 | sed 's/_/ /' | Col 1 3 | sed 's/^/123085 /' > /mnt/deck5/david.jaffe/imm/clone_123085_hh_s1f

cat /mnt/home/wyatt.mcdonnell/immcantation_test/123085/123085_hh_s5f/123085_hh_s5f_heavy_germ-pass.tab |  grep -v "^SEQ" | cut -f1,61 | sed 's/_/ /' | Col 1 3 | sed 's/^/123085 /' > /mnt/deck5/david.jaffe/imm/clone_123085_hh_s5f

*/

extern crate io_utils;
extern crate pretty_trace;
extern crate string_utils;
extern crate vector_utils;

use io_utils::*;
use pretty_trace::*;
use std::collections::HashMap;
use std::fs::File;
use std::env;
use std::io::{BufRead,BufReader};
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    let mut donor = HashMap::<String,usize>::new();
    let s = include_str!["../enclone.testdata"];
    let lines = s.split('\n').collect::<Vec<&str>>();
    for i in 0..lines.len() {
        let mut x = lines[i].to_string();
        if x.starts_with('#') {
            continue;
        }
        x = x.replace( " ", "" );
        let lenas = x.split(',').collect::<Vec<&str>>();
        for j in 0..lenas.len() {
            if lenas[j].contains('-') {
                let l1 = lenas[j].before("-").force_usize();
                let l2 = lenas[j].after("-").force_usize();
                for l in l1..=l2 {
                    donor.insert( format!( "{}", l ), i );
                }
            } else {
                donor.insert( lenas[j].to_string(), i );
            }
        }
    }
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    let mut lc = Vec::<(usize,String)>::new();
    for line in f.lines() {
        let s = line.unwrap();
        let fields = s.split(' ').collect::<Vec<&str>>();
        lc.push( ( fields[2].force_usize(), fields[0].to_string() ) );
    }
    lc.sort();
    let mut bads = 0;
    let mut i = 0;
    while i < lc.len() {
        let j = next_diff1_2( &lc, i as i32 ) as usize;
        let mut donors = Vec::<usize>::new();
        for k in i..j {
            donors.push( donor[&lc[k].1.to_string()] );
        }
        unique_sort(&mut donors);
        if donors.len() > 1 {
            bads += 1;
        }
        i = j;
    }
    println!( "bads = {}", bads );
}
