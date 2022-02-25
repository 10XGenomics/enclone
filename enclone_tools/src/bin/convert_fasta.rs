// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Read in a possibly gzipped fasta file and write it back out as a binary Vec<Vec<u8>>.
//
// Usage convert_fasta in-file out-file

use binary_vec_io::*;
use fasta_tools::*;
use pretty_trace::*;
use std::env;
use std::fs::File;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let (in_file, out_file) = (&args[1], &args[2]);
    let x = read_fasta_to_vec_vec_u8(&in_file);
    binary_write_vec_vec::<u8>(&mut File::create(&out_file).unwrap(), &x).unwrap();
}
