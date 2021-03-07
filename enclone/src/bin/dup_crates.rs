// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Count duplicated crates.

use io_utils::*;
use pretty_trace::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    let mut crates = Vec::<String>::new();
    let f = open_for_read!["Cargo.lock"];
    let mut lines = Vec::<String>::new();
    for line in f.lines() {
        let s = line.unwrap();
        lines.push(s.to_string());
    }
    for i in 0..lines.len() {
        if lines[i] == "[[package]]" {
            crates.push(lines[i + 1].clone());
        }
    }
    let n1 = crates.len();
    unique_sort(&mut crates);
    let n2 = crates.len();
    println!("{}", n1 - n2);
}
