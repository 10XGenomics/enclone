// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_visual::svg_to_geometry::*;
use pretty_trace::*;
use std::env;
use std::fs::read_to_string;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let s = read_to_string(&args[1]).unwrap();
    let g = svg_to_geometry(&s);
    if g.is_none() {
        eprintln!("failed");
    }
}
