// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_visual::convert_svg_to_png::*;
use pretty_trace::*;

fn main() {
    PrettyTrace::new().on();
    let svg = std::fs::read("/Users/david.jaffe/plot.svg").unwrap();
    let png = convert_svg_to_png(&svg);
    std::fs::write("/Users/david.jaffe/plot.png", &png).unwrap();
}
