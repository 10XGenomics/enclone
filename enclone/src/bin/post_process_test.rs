// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Post process cargo test results to remove stuff we don't need to see.  This includes converting
// triple newlines to double newlines and removing leading newlines.  All this is done without
// buffering.

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
    let mut nulls = 0;
    let mut nonnull = false;
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let mut rejected = false;
        for r in reject.iter() {
            if line.contains(r) {
                rejected = true;
            }
        }
        if rejected {
            continue;
        }
        if line.len() == 0 {
            nulls += 1;
            if nulls > 1 || !nonnull {
                continue;
            }
        } else {
            nulls = 0;
            nonnull = true;
        }
        println!("{}", line);
    }
}
