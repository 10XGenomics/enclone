// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Parse the vars file to test if it's valid.

use enclone_vars::var::*;
use pretty_trace::*;

fn main() {
    PrettyTrace::new().on();
    let old = std::fs::read_to_string("enclone_vars/src/vars").unwrap();
    let _ = parse_variables(&old);
}
