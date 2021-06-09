// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_visual::svg_to_geometry::*;
use pretty_trace::*;
use std::fs::read_to_string;

fn main() {
    PrettyTrace::new().on();
    let s = read_to_string("~/p.svg").unwrap();
    let g = svg_to_geometry(&s);
    if g.is_none() {
        eprintln!("failed");
    }
}
