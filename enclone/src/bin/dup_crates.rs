// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Find duplicated crates.

use io_utils::*;
use pretty_trace::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use string_utils::*;
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
    let mut dups = 0;
    let mut duplist = Vec::<(isize, String)>::new();
    let mut i = 0;
    while i < crates.len() {
        let j = next_diff(&crates, i);
        if j - i > 1 {
            dups += j - i - 1;
            duplist.push((
                i as isize - j as isize,
                crates[i].between("\"", "\"").to_string(),
            ));
        }
        i = j;
    }
    duplist.sort();
    println!("{}", dups);
    for x in duplist.iter() {
        println!("{} {}", -x.0, x.1);
    }
}
