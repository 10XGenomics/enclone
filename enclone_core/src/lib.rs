// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod allowed_vars;
pub mod copy_for_enclone;
pub mod defs;
pub mod join_one;
pub mod linear_condition;
pub mod mammalian_fixed_len;
pub mod median;
pub mod print_tools;
pub mod slurp;
pub mod testlist;
pub mod vdj_features;

use std::env;

const VERSION_STRING: &'static str = env!("VERSION_STRING");

pub fn version_string() -> String {
    VERSION_STRING.to_string()
}
