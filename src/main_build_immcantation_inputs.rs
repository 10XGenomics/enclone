// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Build input files for Immcantation.  This creates two files:
// filtered_contig.fasta
// filtered_contig_annotations.csv.

use io_utils::*;
use pretty_trace::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use string_utils::*;

pub fn main_build_immcantation_inputs() {
    PrettyTrace::new().on();

    // Location of files.

    let pre = "/mnt/assembly/vdj/current12";

    // Get the list of lena ids.

    let testlines = include_str!["enclone.testdata"];
    let testlines = testlines.split('\n').collect::<Vec<&str>>();
    let mut lenas = Vec::<usize>::new();
    for x in testlines.iter() {
        if !x.starts_with('#') {
            let xx = x.replace(" ", "");
            let y = xx.split(',').collect::<Vec<&str>>();
            for z in y.iter() {
                if *z == "" {
                } else if !z.contains('-') {
                    lenas.push(z.force_usize());
                } else {
                    let l1 = z.before("-").force_usize();
                    let l2 = z.after("-").force_usize();
                    for l in l1..=l2 {
                        lenas.push(l);
                    }
                }
            }
        }
    }
    lenas.sort();

    // Concatenate the fasta files, prepending the lena id to the contig name.
    // Also concatenate the csv files, prepending the lena id to the contig name.

    let mut g1 = open_for_write_new!["filtered_contig.fasta"];
    let mut g = open_for_write_new!["filtered_contig_annotations.csv"];
    for (i, l) in lenas.iter().enumerate() {
        let mut count1 = 0;
        let f1 = open_for_read![&format!("{}//{}/outs/filtered_contig.fasta", pre, l)];
        for line in f1.lines() {
            let s = line.unwrap();
            if s.starts_with('>') {
                fwriteln!(g1, ">{}_{}", l, s.after(">"));
                count1 += 1;
            } else {
                fwriteln!(g1, "{}", s);
            }
        }
        let f = open_for_read![&format!(
            "{}//{}/outs/filtered_contig_annotations.csv",
            pre, l
        )];
        let mut count2 = 0;
        let mut first = true;
        for line in f.lines() {
            let s = line.unwrap();
            if i == 0 && first {
                fwriteln!(g, "{}", s);
            }
            if first {
                first = false;
                continue;
            }
            let fields = s.split(',').collect::<Vec<&str>>();
            assert!(fields[2].contains("_contig_"));
            for j in 0..fields.len() {
                if j > 0 {
                    fwrite!(g, ",");
                }
                if j == 2 {
                    fwrite!(g, "{}_", l);
                    count2 += 1;
                }
                fwrite!(g, "{}", fields[j]);
            }
            fwriteln!(g, "");
        }
        if count1 != count2 {
            eprintln!("\ninconsistency");
            eprintme!(l, count1, count2);
            eprintln!("");
            std::process::exit(1);
        }
    }
}
