// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Read in a possibly gzipped fasta file and write it back out as a bincoded Vec<Vec<u8>>.
//
// Usage convert_fasta in-file out-file

use fasta_tools::*;
use io_utils::write_obj;

use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let x = read_fasta_to_vec_vec_u8(&args[1]);
    write_obj(&x, &args[2]);
}
