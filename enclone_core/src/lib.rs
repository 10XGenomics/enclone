// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod align_to_vdj_ref;
pub mod allowed_vars;
pub mod copy_for_enclone;
pub mod defs;
pub mod join_one;
pub mod linear_condition;
pub mod mammalian_fixed_len;
pub mod median;
pub mod opt_d;
pub mod print_tools;
pub mod run_test;
pub mod slurp;
pub mod testlist;
pub mod vdj_features;

use std::env;

const VERSION_STRING: &'static str = env!("VERSION_STRING");

pub fn version_string() -> String {
    VERSION_STRING.to_string()
}

// Parse a line, breaking at blanks, but not if they're in quotes.  And strip the quotes.
// Ridiculously similar to parse_csv, probably should refactor.

pub fn parse_bsv(x: &str) -> Vec<String> {
    let mut args = Vec::<String>::new();
    let mut w = Vec::<char>::new();
    for c in x.chars() {
        w.push(c);
    }
    let (mut quotes, mut i) = (0, 0);
    while i < w.len() {
        let mut j = i;
        while j < w.len() {
            if quotes % 2 == 0 && w[j] == ' ' {
                break;
            }
            if w[j] == '"' {
                quotes += 1;
            }
            j += 1;
        }
        let (mut start, mut stop) = (i, j);
        if stop - start >= 2 && w[start] == '"' && w[stop - 1] == '"' {
            start += 1;
            stop -= 1;
        }
        let mut s = String::new();
        for m in start..stop {
            s.push(w[m]);
        }
        args.push(s);
        i = j + 1;
    }
    args
}
