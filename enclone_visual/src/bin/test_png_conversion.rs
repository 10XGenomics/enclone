// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// usage: test_png_conversion whatever.svg
// ==> generates whatever.png

use enclone_core::convert_svg_to_png::*;
use pretty_trace::*;
use std::env;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let svg_file = &args[1];
    let png_file = format!("{}.png", svg_file.rev_before(".svg"));
    let svg = std::fs::read(&svg_file).unwrap();
    let png = convert_svg_to_png(&svg, 2000);
    std::fs::write(&png_file, &png).unwrap();
}
