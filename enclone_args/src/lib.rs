// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use io_utils::*;

pub mod load_gex;
pub mod load_gex_core;
pub mod proc_args;
pub mod proc_args2;
pub mod proc_args3;
pub mod proc_args_check;
pub mod proc_args_post;
pub mod process_special_arg;
pub mod read_json;

// parse_csv_pure: same as parse_csv, but don't strip out quotes

pub fn parse_csv_pure(x: &str) -> Vec<String> {
    let mut y = Vec::<String>::new();
    let mut w = Vec::<char>::new();
    for c in x.chars() {
        w.push(c);
    }
    let (mut quotes, mut i) = (0, 0);
    while i < w.len() {
        let mut j = i;
        while j < w.len() {
            if quotes % 2 == 0 && w[j] == ',' {
                break;
            }
            if w[j] == '"' {
                quotes += 1;
            }
            j += 1;
        }
        let (start, stop) = (i, j);
        let mut s = String::new();
        for m in start..stop {
            s.push(w[m]);
        }
        y.push(s);
        i = j + 1;
    }
    if !w.is_empty() && *w.last().unwrap() == ',' {
        y.push(String::new());
    }
    y
}

pub fn fnx(outs: &str, name: &str) -> String {
    let mut file = format!("{}/../{}", outs, name);
    if !path_exists(&file) {
        file = format!("{}/{}", outs, name);
    }
    file
}
