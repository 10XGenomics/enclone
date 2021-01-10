// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod defs;
pub mod mammalian_fixed_len;
pub mod print_tools;
pub mod testlist;
pub mod vdj_features;

const VERSION_STRING: &'static str = env!("VERSION_STRING");

// Return the code version string.

pub fn version_string() -> String {
    VERSION_STRING.to_string()
}
