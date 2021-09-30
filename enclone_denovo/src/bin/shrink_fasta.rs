// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Read in a possibly gzipped fasta file and write it back out as a fasta file with runs
// of Ns truncated to length 1000.
//
// Usage shrink_fasta in-file out-file

use fasta_tools::*;
use io_utils::*;
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let (in_file, out_file) = (&args[1], &args[2]);
    let n = 1000;
    let mut x = read_fasta_to_vec_vec_u8(in_file);
    if n > 0 {
        for i in 0..x.len() {
            let mut j = 0;
            let mut y = Vec::<u8>::new();
            while j < x[i].len() {
                if x[i][j] != b'N' {
                    y.push(x[i][j]);
                    j += 1;
                } else {
                    let mut k = j + 1;
                    while k < x[i].len() && x[i][k] == b'N' {
                        k += 1;
                    }
                    let s = std::cmp::min(k - j, n);
                    for _ in 0..s {
                        y.push(b'N');
                    }
                    j = k;
                }
            }
            x[i] = y;
        }
    }
    let mut f = open_for_write_new![&out_file];
    for i in 0..x.len() {
        if i % 2 == 0 {
            fwriteln!(f, ">{}", strme(&x[i]));
        } else {
            for j in 0..x[i].len() {
                if j > 0 && j % 80 == 0 {
                    fwriteln!(f, "");
                }
                fwrite!(f, "{}", x[i][j] as char);
            }
            fwriteln!(f, "");
        }
    }
}
