// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Find master.toml files repos in the given directory (e.g. ~/repos), and identify
// inconsistencies.

use io_utils::{dir_list, open_for_read, path_exists};
use pretty_trace::PrettyTrace;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use string_utils::TextUtils;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let repos = &args[1];
    let list = dir_list(repos);
    let mut origins = Vec::<String>::new();
    let mut masters = Vec::<Vec<String>>::new();
    for x in list.iter() {
        let master1 = format!("{}/{}/master.toml", repos, x);
        let master2 = format!("{}/{}/lib/rust/master.toml", repos, x);
        let mut master = master1.clone();
        if !path_exists(&master) {
            master = master2.clone();
        }
        if path_exists(&master) {
            origins.push(x.to_string());
            let f = open_for_read![&master];
            let mut y = Vec::<String>::new();
            for line in f.lines() {
                let s = line.unwrap();
                if !s.starts_with('#') && s.contains('=') {
                    y.push(s);
                }
            }
            masters.push(y);
        }
    }
    let mut count = 0;
    for i1 in 0..origins.len() {
        for i2 in i1 + 1..origins.len() {
            for m1 in 0..masters[i1].len() {
                for m2 in 0..masters[i2].len() {
                    if masters[i1][m1].before("=") == masters[i2][m2].before("=")
                        && masters[i1][m1] != masters[i2][m2]
                    {
                        println!(
                            "\ninconsistency between {} and {}",
                            origins[i1], origins[i2],
                        );
                        println!("{}\n{}", masters[i1][m1], masters[i2][m2]);
                        count += 1;
                    }
                }
            }
        }
    }
    println!("\ntotal inconsistencies = {}\n", count);
}
