// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Display a CSV file, given by the first argument.  Pipe to less.
//
// This can be used to examine the output of POUT.

use io_utils::open_for_read;
use pretty_trace::PrettyTrace;
use std::env;

use std::io::BufRead;
use string_utils::parse_csv;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let f = open_for_read![&args[1]];
    let mut fields = Vec::<String>::new();
    for (i, line) in f.lines().enumerate() {
        let s = line.unwrap();
        let values = parse_csv(&s);
        if i == 0 {
            fields = values;
        } else {
            for j in 0..fields.len() {
                println!(
                    "\nline {}, fields {} ({}) = {}",
                    i,
                    j + 1,
                    fields[j],
                    values[j]
                );
            }
        }
    }
}
