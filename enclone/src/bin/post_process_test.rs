// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Post process cargo test results to remove stuff we don't need to see.
// This is done without buffering.

use std::io::{self, BufRead};

fn main() {
    let reject = [
        "running ",
        " Running ",
        " Doc-tests ",
        " filtered out",
        " Compiling ",
        " Finished ",
    ];
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let mut rejected = false;
        for r in reject.iter() {
            if line.contains(r) {
                rejected = true;
            }
        }
        if !rejected && line.len() > 0 {
            println!("{}", line);
        }
    }
}
