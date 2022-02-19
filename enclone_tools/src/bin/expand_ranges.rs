// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// EXAMPLE
// expand_ranges 100-101,500,600-602
// 100,101,500,600,601,602
//
// Not robust to funny input.

use itertools::Itertools;
use pretty_trace::*;
use std::env;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let fields = args[1].split(',').collect::<Vec<&str>>();
    let mut x = Vec::<usize>::new();
    for i in 0..fields.len() {
        if !fields[i].contains("-") {
            x.push(fields[i].force_usize());
        } else {
            let a = fields[i].before("-").force_usize();
            let b = fields[i].after("-").force_usize();
            for m in a..=b {
                x.push(m);
            }
        }
    }
    println!("{}", x.iter().format(","));
}
