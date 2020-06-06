// Copyright (c) 2020 10x Genomics, Inc. All rights reserved.

pub mod defs;
pub mod print_tools;
pub mod testlist;
pub mod types;

const VERSION_STRING: &'static str = env!("VERSION_STRING");

pub fn version_string() -> String {
    VERSION_STRING.to_string()
}
