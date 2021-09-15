// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod align_to_vdj_ref;
pub mod allowed_vars;
pub mod cell_color;
pub mod combine_group_pics;
pub mod convert_svg_to_png;
pub mod defs;
pub mod join_one;
pub mod linear_condition;
pub mod logging;
pub mod main_testlist;
pub mod mammalian_fixed_len;
pub mod median;
pub mod opt_d;
pub mod prepare_for_apocalypse;
pub mod print_tools;
pub mod set_speakers;
pub mod slurp;
pub mod testlist;
pub mod var_reg;
pub mod vdj_features;

use lazy_static::*;
use std::env;
use std::io::BufRead;
use std::sync::Mutex;
use std::time::Duration;

lazy_static! {
    pub static ref BUG_REPORT_ADDRESS: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref REMOTE_HOST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
}

const VERSION_STRING: &'static str = env!("VERSION_STRING");

// WARNING: the version string will not be up to date unless enclone_core/build.rs is touched
// before compiling.  Use current_version_string() to always get the current version.

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

// Test to see if a line can be read from the given file f.  If not, return an error message
// the references arg, which is supposed to be the name of a command line argument from which
// f originated.

pub fn require_readable_file(f: &str, arg: &str) -> Result<(), String> {
    let x = std::fs::File::open(&f);
    if !x.is_ok() {
        return Err(format!(
            "\nThe file {} could not be opened because {}.  This came from \
            the command line argument {}.\n",
            f,
            x.err().unwrap(),
            arg,
        ));
    }
    let y = std::io::BufReader::new(x.unwrap());
    for line in y.lines() {
        if line.is_err() {
            let mut err = line.err().unwrap().to_string();
            if err.starts_with("Is a directory") {
                err = "it is a directory".to_string();
            }
            return Err(format!(
                "\nThe file {} could not be read because {}.\nThis came from \
                the command line argument {}.\n",
                f, err, arg,
            ));
        }
        break;
    }
    Ok(())
}
