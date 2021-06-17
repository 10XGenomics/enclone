// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod align_to_vdj_ref;
pub mod allowed_vars;
pub mod defs;
pub mod join_one;
pub mod linear_condition;
pub mod mammalian_fixed_len;
pub mod median;
pub mod opt_d;
pub mod print_tools;
pub mod slurp;
pub mod testlist;
pub mod vdj_features;

use std::env;
use std::time::Duration;

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

pub fn fetch_url(url: &str) -> Result<String, String> {
    const TIMEOUT: u64 = 120; // timeout in seconds
    let req = attohttpc::get(url.clone()).read_timeout(Duration::new(TIMEOUT, 0));
    let response = req.send();
    if response.is_err() {
        panic!("Failed to access URL {}.", url);
    }
    let response = response.unwrap();
    if !response.is_success() {
        let msg = response.text().unwrap();
        if msg.contains("Not found") {
            return Err(format!(
                "\nAttempt to access the URL\n{}\nfailed with \"Not found\".  Could there \
                be something wrong with the id?\n",
                url
            ));
        }
        return Err(format!("Failed to access URL {}: {}.", url, msg));
    }
    Ok(response.text().unwrap())
}
