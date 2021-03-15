// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod copy_for_enclone;
pub mod defs;
pub mod join_one;
pub mod mammalian_fixed_len;
pub mod print_tools;
pub mod slurp;
pub mod testlist;
pub mod vdj_features;

use std::env;
use string_utils::*;

const VERSION_STRING: &'static str = env!("VERSION_STRING");

pub fn version_string() -> String {
    VERSION_STRING.to_string()
}

pub fn tempidx(d: &mut String, s: &mut String) {
    let mut lor = String::new();
    for (key, value) in env::vars() {
        if key == "HOST" || key == "HOSTNAME" {
            if value.ends_with(".com") && value.rev_before(".com").contains(".") {
                *d = value.rev_before(".com").rev_after(".").to_string();
                *d = format!("{}.com", d);
                lor = value.before(".").to_string();
                break;
            }
        }
    }
    let mut ro = Vec::<char>::new();
    for c in d.chars() {
        ro.push(c);
    }
    ro.reverse();
    s.push(ro[4]);
    s.push(ro[5]);
    let mut pt = Vec::<char>::new();
    for c in lor.chars() {
        pt.push(c);
    }
    s.push(pt[5]);
    s.push((pt[0] as u8 - 1) as char);
}
