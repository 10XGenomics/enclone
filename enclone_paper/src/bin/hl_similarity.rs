// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Compute heavy chain similarity to light chain similarity, considering only memory cells from
// different donors.
//
// enclone BCR=@test BUILT_IN CHAINS_EXACT=2 CHAINS=2 NOPRINT POUT=stdout PCELL ECHOC
//         PCOLS=datasets_cell,donors_cell,v_name1,v_name2,dref,cdr3_aa1,clonotype_ncells,
//         const1,hcomp,jun_ins,d1_name1,vj_aa_nl1,vj_aa_nl2
//         > per_cell_stuff
//
// hl_similarity per_cell_stuff
//
// Generate a CSV file as output with fields as follows.
// class: 1 or 2; 2 has the additional restriction that the heavy chain gene names and
//        CDRH3 lengths are the same
// donor1: d1 or d2 or d3 or d4
// donor2: d1 or d2 or d3 or d4
// const1: const region name for first cell
// const2: const region name for first cell
// hd: heavy chain edit distance, excluding leader
// ld: light chain edit distance, excluding leader.

use io_utils::*;
use pretty_trace::PrettyTrace;
use rand_chacha;
use rand_chacha::rand_core::RngCore;
use rand_chacha::rand_core::SeedableRng;
use std::collections::{HashMap, HashSet};
use std::env;
use std::io::BufRead;
use string_utils::TextUtils;
use triple_accel::levenshtein;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Load data.

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
        String,
        String,
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
            assert!(tof.contains_key("datasets_cell"));
            assert!(tof.contains_key("donors_cell"));
            assert!(tof.contains_key("v_name1"));
            assert!(tof.contains_key("v_name2"));
            assert!(tof.contains_key("dref"));
            assert!(tof.contains_key("cdr3_aa1"));
            assert!(tof.contains_key("clonotype_ncells"));
            assert!(tof.contains_key("const1"));
            assert!(tof.contains_key("hcomp"));
            assert!(tof.contains_key("jun_ins"));
            assert!(tof.contains_key("d1_name1"));
            first = false;
        } else {
            data.push((
                /* 0 */ fields[tof["v_name1"]].to_string(),
                /* 1 */ fields[tof["cdr3_aa1"]].len(),
                /* 2 */ fields[tof["cdr3_aa1"]].to_string().as_bytes().to_vec(),
                /* 3 */ fields[tof["donors_cell"]].to_string(),
                /* 4 */ fields[tof["v_name2"]].to_string(),
                /* 5 */ fields[tof["dref"]].force_usize(),
                /* 6 */ fields[tof["clonotype_ncells"]].force_usize(),
                /* 7 */ fields[tof["const1"]].to_string(),
                /* 8 */ fields[tof["hcomp"]].force_usize(),
                /* 9 */ fields[tof["datasets_cell"]].force_usize(),
                /* 10 */ fields[tof["vj_aa_nl1"]].to_string(),
                /* 11 */ fields[tof["vj_aa_nl2"]].to_string(),
            ));
        }
    }
    data.sort();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Replace paralogs.

    for i in 0..data.len() {
        data[i].4 = data[i].4.replace("D", "");
    }

    // Print CSV header.

    println!("class,donor1,donor2,const1,const2,hd,ld");

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Select pairs of cells having different donors at random.  Compute their heavy and light
    // chain edit distances.

    const SAMPLE: usize = 5_000;
    let mut randme = rand_chacha::ChaCha8Rng::seed_from_u64(123456789);
    let mut seen = HashSet::<(usize, usize)>::new();
    while seen.len() < SAMPLE {
        let i1 = randme.next_u64() as usize % data.len();
        let i2 = randme.next_u64() as usize % data.len();
        if i1 < i2 && !seen.contains(&(i1, i2)) {
            let donor1 = &data[i1].3;
            let donor2 = &data[i2].3;
            let dref1 = data[i1].5;
            let dref2 = data[i2].5;
            if donor1 != donor2 && dref1 > 0 && dref2 > 0 {
                seen.insert((i1, i2));
                let h1 = &data[i1].10.as_bytes();
                let h2 = &data[i2].10.as_bytes();
                let hd = levenshtein(&h1, &h2) as usize;
                let l1 = &data[i1].11.as_bytes();
                let l2 = &data[i2].11.as_bytes();
                let ld = levenshtein(&l1, &l2) as usize;
                let class = 1;
                let const1 = &data[i1].7;
                let const2 = &data[i2].7;
                println!("{class},{donor1},{donor2},{const1},{const2},{hd},{ld}");
            }
        }
    }

    // Define groups based on equal heavy chain gene names and CDRH3 amino acid sequences.

    let mut bounds = Vec::<(usize, usize)>::new();
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        while j < data.len() {
            if data[j].0 != data[i].0 || data[j].1 != data[i].1 {
                break;
            }
            j += 1;
        }
        bounds.push((i, j));
        i = j;
    }

    // Find pairs.

    let mut bucket = Vec::<(usize, usize)>::new();
    for m in 0..bounds.len() {
        let i = bounds[m].0;
        let j = bounds[m].1;
        for k1 in i..j {
            for k2 in k1 + 1..j {
                let donor1 = &data[k1].3;
                let donor2 = &data[k2].3;
                let dref1 = data[k1].5;
                let dref2 = data[k2].5;
                if donor1 != donor2 && dref1 > 0 && dref2 > 0 {
                    let hname1 = &data[k1].0;
                    let hname2 = &data[k2].0;
                    if hname1 == hname2 && data[k1].2 == data[k2].2 {
                        bucket.push((k1, k2));
                    }
                }
            }
        }
    }

    // The same, but now assume that heavy chain gene names and CDRH3 lengths are the same.

    while seen.len() < 2 * SAMPLE {
        let m = randme.next_u64() as usize % bucket.len();
        let i1 = bucket[m].0;
        let i2 = bucket[m].1;
        if !seen.contains(&(i1, i2)) {
            seen.insert((i1, i2));
            let h1 = &data[i1].10.as_bytes();
            let h2 = &data[i2].10.as_bytes();
            let hd = levenshtein(&h1, &h2) as usize;
            let l1 = &data[i1].11.as_bytes();
            let l2 = &data[i2].11.as_bytes();
            let ld = levenshtein(&l1, &l2) as usize;
            let class = 2;
            let donor1 = &data[i1].3;
            let donor2 = &data[i2].3;
            let const1 = &data[i1].7;
            let const2 = &data[i2].7;
            println!("{class},{donor1},{donor2},{const1},{const2},{hd},{ld}");
        }
    }
}
