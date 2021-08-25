// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Run feature_barcode_matrix on the given id.

use enclone_tools::feature_barcode_matrix::*;
use pretty_trace::*;
use std::env;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let _ = feature_barcode_matrix(args[1].force_usize(), true);
}
