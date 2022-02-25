// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// annotate_seqs fasta-file vdj-reference-file <number of reads to use>

use annotate::*;
use debruijn::dna_string::*;
use io_utils::*;
use pretty_trace::*;
use rayon::prelude::*;
use refx::*;
use std::env;
use std::io::{BufRead, Write};
use string_utils::*;
use vdj_ann::transcript::*;
use vdj_ann::*;

fn main() {
    // Set up and parse args.

    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let (fasta_file, ref_file) = (&args[1], &args[2]);
    let mut nreads = 0;
    if args.len() >= 4 {
        nreads = args[3].force_usize();
    }

    // Read sequences.

    let mut seqs = Vec::<DnaString>::new();
    let fin = open_for_read![&fasta_file];
    let mut last: String = String::new();
    let mut first = true;
    for line in fin.lines() {
        let s = line.unwrap();
        if first {
            if !s.starts_with(">") {
                panic!("fasta format failure");
            }
            first = false;
        } else {
            if s.starts_with(">") {
                seqs.push(DnaString::from_dna_string(&last));
                if seqs.len() == nreads {
                    break;
                }
                last.clear();
            } else {
                last += &s;
            }
        }
    }
    seqs.push(DnaString::from_dna_string(&last));
    if seqs.len() < nreads {
        eprintln!("\nOnly have {} sequences.\n", seqs.len());
        std::process::exit(1);
    }

    // Make reference data.

    let mut refdata = RefData::new();
    let refx = std::fs::read_to_string(&ref_file).unwrap();
    let ext_ref = String::new();
    make_vdj_ref_data_core(&mut refdata, &refx, &ext_ref, true, true, None);

    // Annotate the sequences.

    let mut results = Vec::<(usize, Vec<u8>, usize)>::new();
    for i in 0..seqs.len() {
        results.push((i, Vec::<u8>::new(), 0));
    }
    results.par_iter_mut().for_each(|r| {
        let i = r.0;
        let mut seq = seqs[i].clone();
        let len = seq.len();
        if seq.len() > 1000 {
            seq = seq.slice(0, 1000).to_owned();
        }
        let mut log = Vec::<u8>::new();
        fwriteln!(log, ", len = {}", len);
        print_annotations(&seq, &refdata, &mut log, false, true, false);
        let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
        annotate_seq(&seq, &refdata, &mut ann, true, false, true);
        let valid = is_valid(&seq, &refdata, &ann, false, &mut log, None);
        let mut shift = false;
        for i in 1..ann.len() {
            if ann[i - 1].2 == ann[i].2 {
                if ann[i].0 == ann[i - 1].0 + ann[i - 1].1 + 1
                    && ann[i].3 == ann[i - 1].3 + ann[i - 1].1
                {
                    shift = true;
                }
                if ann[i].0 == ann[i - 1].0 + ann[i - 1].1
                    && ann[i].3 == ann[i - 1].3 + ann[i - 1].1 + 1
                {
                    shift = true;
                }
            }
        }
        if shift {
            r.2 = 1;
        }
        let mut match_len = 0;
        for j in 0..ann.len() {
            match_len += ann[j].1;
        }
        let mut log0 = Vec::<u8>::new();
        print_cdr3_using_ann(&seq, &refdata, &ann, &mut log0);
        if (match_len >= 40 && !shift) || !log0.is_empty() {
            log.append(&mut log0);
            fwrite!(r.1, "{}", stringme(&log));
            if valid {
                fwrite!(r.1, "VALID!\n");
            }
            fwrite!(r.1, "\n");
        }
    });
    println!("");
    let (mut count, mut shift, mut valids) = (0, 0, 0);
    for i in 0..results.len() {
        shift += results[i].2;
        if results[i].1.len() > 0 {
            count += 1;
            let log = strme(&results[i].1);
            print!("sequence {} = {}{}", count, results[i].0 + 1, log);
            if log.contains("VALID") {
                valids += 1;
            }
        }
    }
    println!("{} valids", valids);
    println!("plus {} frameshifted transcripts\n", shift);
}
