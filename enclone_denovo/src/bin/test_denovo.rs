// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Test program to see if some contigs match de novo V genes better than our reference or IMGT.
//
// enclone BI=2 BUILT_IN NOPRINT POUT=stdout PCOLS=vj_seq1 > seqs
// or vj_seq2

use fasta_tools::*;
use io_utils::*;
use pretty_trace::PrettyTrace;
use std::cmp::min;
use std::io::BufRead;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let denovo = read_fasta_to_vec_vec_u8("xxx.fasta");
    let imgt = read_fasta_to_vec_vec_u8("imgt.fasta");
    let tenx = read_fasta_to_vec_vec_u8("regions.fa");
    let f = open_for_read!["seqs"];
    for (l, line) in f.lines().enumerate() {
        if l > 0 {
            let s = line.unwrap();
            let s = s.as_bytes();
            if s.len() == 0 {
                continue;
            }
            let mut best_mis = 1000000;
            let mut best_i = 0;
            for i in 0..denovo.len() / 2 {
                let i = 2 * i + 1;
                let mut mis = 0;
                for j in 0..min(denovo[i].len(), s.len()) {
                    if denovo[i][j] != s[j] {
                        mis += 1;
                    }
                }
                if mis < best_mis {
                    best_mis = mis;
                    best_i = i;
                }
            }
            let mut best_mis_imgt = 1000000;
            let mut best_i_imgt = 0;
            for i in 0..imgt.len() {
                if imgt[i].len() < 300 {
                    continue;
                }
                let mut mis = 0;
                for j in 0..min(imgt[i].len(), s.len()) {
                    if imgt[i][j] != s[j] {
                        mis += 1;
                    }
                }
                if mis < best_mis_imgt {
                    best_mis_imgt = mis;
                    best_i_imgt = i;
                }
            }
            let mut best_mis_tenx = 1000000;
            let mut best_i_tenx = 0;
            for i in 0..tenx.len() {
                if tenx[i].len() < 300 {
                    continue;
                }
                let mut mis = 0;
                for j in 0..min(tenx[i].len(), s.len()) {
                    if tenx[i][j] != s[j] {
                        mis += 1;
                    }
                }
                if mis < best_mis_tenx {
                    best_mis_tenx = mis;
                    best_i_tenx = i;
                }
            }
            if best_mis < best_mis_tenx && best_mis < best_mis_imgt {
                println!("");
                printme!(l, best_mis, best_mis_tenx, best_mis_imgt);
                println!("seq  = {}", strme(&s));
                println!("novo = {}", strme(&denovo[best_i]));
                println!("tenx = {}", strme(&tenx[best_i_tenx]));
                println!("imgt = {}", strme(&imgt[best_i_imgt]));
            }
        }
    }
}
