// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.
//
// Splice junctions from olga.

use debruijn::dna_string::DnaString;
use io_utils::*;
use pretty_trace::*;
use std::env;
use std::io::BufRead;
use string_utils::TextUtils;
use vdj_ann::annotate::{annotate_seq, get_cdr3_using_ann};
use vdj_ann::refx::make_vdj_ref_data_core;
use vdj_ann::refx::RefData;
use vdj_ann_ref::human_ref;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();

    // Get the human reference used by enclone.

    let mut refnames = Vec::<String>::new();
    let mut refs = Vec::<Vec<u8>>::new();
    let mut utr = Vec::<bool>::new();
    let href = human_ref();
    {
        for (i, line) in href.lines().enumerate() {
            if i % 2 == 0 {
                let n = line.between("|", " ").to_string();
                utr.push(line.contains("5'UTR"));
                refnames.push(n);
            } else {
                refs.push(line.as_bytes().to_vec());
            }
        }
    }

    // Add extra genes.  For IGHV3-NL1, we took the sequence in olga.seq, and found the leader
    // in GenBank AC245166.2:
    // ATGGAGAAATAGAGAGACTGAGTGTGAGTGAACAT
    // GAGTGAGAAAAACTGGATTTGTGTGGCATTTTCTGATAACGGTGTCCTTCTGTTTGCAGGTGTCCAGTGT.
    // However, this leader has stop codons.  Therefore we put in a FAKE leader, that for
    // IGHV3-11.

    refnames.push("IGHV3-43D".to_string());
    utr.push(false);
    refs.push(
        b"ATGGAGTTTGGACTGAGCTGGGTTTTCCTTGTTGCTATTTTAAAAGGTGTCCAGTGTGAAGTGCAGCTGGTGGAGTCT\
        GGGGGAGTCGTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTGATGATTATGCCATGCACTG\
        GGTCCGTCAAGCTCCGGGGAAGGGTCTGGAGTGGGTCTCTCTTATTAGTTGGGATGGTGGTAGCACCTACTATGCAGACTCTGTGA\
        AGGGTCGATTCACCATCTCCAGAGACAACAGCAAAAACTCCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACCGCCTTG\
        TATTACTGTGCAAAAGATA"
            .to_vec(),
    );
    refnames.push("IGHV3-30-3".to_string());
    utr.push(false);
    refs.push(
        b"ATGGAGTTTGGGCTGAGCTGGGTTTTCCTCGTTGCTCTTTTAAGAGGTGTCCAGTGTCAGGTGCAGCTGGTGGAGTCT\
        GGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCTATGCACTG\
        GGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGCAATAAATACTACGCAGACTCCGTGA\
        AGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTG\
        TATTACTGTGCGAGAGA"
            .to_vec(),
    );
    refnames.push("IGHV3-NL1".to_string());
    refs.push(
        b"ATGGAGTTTGGGCTGAGCTGGGTTTTCCTTGTTGCTATTATAAAAGGTG\
        TCCAGTGTCAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGGGGTCCCTGAGACT\
        CTCCTGTGCAGCGTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTCT\
        CAGTTATTTATAGCGGTGGTAGTAGCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAAC\
        ACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAAAGA"
            .to_vec(),
    );
    utr.push(false);
    refnames.push("IGHV4-30-2".to_string());
    refs.push(
        b"ATGAAACACCTGTGGTTCTTCCTCCTGCTGGTGGCAGCTCCCAGATGGGTCCTGTCCCAGCTGCAGCTGCAGGAGTCC\
        GGCTCAGGACTGGTGAAGCCTTCACAGACCCTGTCCCTCACCTGCGCTGTCTCTGGTGGCTCCATCAGCAGTGGTGGTTACTCCTG\
        GAGCTGGATCCGGCAGCCACCAGGGAAGGGCCTGGAGTGGATTGGGAGTATCTATTATAGTGGGAGCACCTACTACAACCCGTCCC\
        TCAAGAGTCGAGTCACCATATCCGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCTGCAGACACGGCT\
        GTGTATTACTGTGCGAGACA"
            .to_vec(),
    );
    utr.push(false);
    let nref = refs.len();
    let mut refdata = RefData::new();
    let ext_ref = String::new();
    make_vdj_ref_data_core(&mut refdata, &href, &ext_ref, false, true, None);

    // Load the input file.

    let mut jun = Vec::<Vec<u8>>::new();
    let mut hv = Vec::<String>::new();
    let mut hj = Vec::<String>::new();
    let mut cdr3 = Vec::<String>::new();
    let f = open_for_read![&args[1]];
    for line in f.lines() {
        let s = line.unwrap();
        let fields = s.split(',').collect::<Vec<&str>>();
        jun.push(fields[0].as_bytes().to_vec());
        hv.push(fields[2].to_string());
        hj.push(fields[3].to_string());
        cdr3.push(fields[1].to_string());
    }
    let n = jun.len();

    // Process entries.

    let mut fails = 0;
    // for i in 0..n {
    // for i in 0..20000 {
    for i in 0..10 {
        println!(
            "\n-------------------------------------------------------------------------\
            --------------------------"
        );
        println!("\nsequence {}", i + 1);
        let mut vseq = Vec::<u8>::new();
        let mut jseq = Vec::<u8>::new();
        for j in 0..nref {
            if refnames[j] == hv[i] && !utr[j] {
                vseq = refs[j].clone();
                break;
            }
        }
        if vseq.is_empty() {
            println!(
                "\nEntry {}, unable to find V gene named {}.\n",
                i + 1,
                hv[i]
            );
            std::process::exit(1);
        }
        for j in 0..nref {
            if refnames[j] == hj[i] {
                jseq = refs[j].clone();
                break;
            }
        }
        if jseq.is_empty() {
            println!(
                "\nEntry {}, unable to find J gene named {}.\n",
                i + 1,
                hj[i]
            );
            std::process::exit(1);
        }
        let junv = &jun[i][0..2];
        use string_utils::strme;
        println!("\nV gene = {}", hv[i]);
        println!("J gene = {}", hj[i]);
        println!("\nlooking for {} in {}", strme(&junv), strme(&vseq));
        let mut vstart = None;
        for k in (0..=vseq.len() - junv.len()).rev() {
            if vseq[k..].starts_with(&junv) {
                if k % 3 != 0 {
                    continue;
                }
                vstart = Some(k);
                println!("vstart for {} found at pos -{}", i + 1, vseq.len() - k);
                break;
            }
        }
        if vstart.is_none() {
            println!("\nfailed to find vstart for entry {}\n", i + 1);
            std::process::exit(1);
        }
        let mut jtail = jseq[jseq.len() - 31..].to_vec(); // concatenate to junction
        let mut full = vseq[0..vstart.unwrap()].to_vec();
        full.append(&mut jun[i].clone());
        full.append(&mut jtail);
        println!("\ncdr3[{}] = {}\n", i + 1, cdr3[i]);
        println!("full[{}] = {}\n", i + 1, strme(&full));
        let x = DnaString::from_dna_string(&strme(&full));
        let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
        annotate_seq(&x, &refdata, &mut ann, true, false, true);
        let mut cdr3x = Vec::<(usize, Vec<u8>, usize, usize)>::new();
        get_cdr3_using_ann(&x, &refdata, &ann, &mut cdr3x);
        if cdr3x.len() != 1 {
            println!("failed to find unique CDR3\n");
            println!("found {} CDR3s\n", cdr3x.len());
            fails += 1;
            continue;
        }
        println!("CDR3 = {}", strme(&cdr3x[0].1));
        if strme(&cdr3x[0].1) != cdr3[i] {
            println!("\nmismatch");
            std::process::exit(1);
        }
    }
    println!("\nThere were {} fails.\n", fails);
}
