// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Actually export code.

// Read the vars file and export code.  This is a partial implementation.

use enclone_vars::export_code::*;
use io_utils::*;
use pretty_trace::PrettyTrace;

use std::io::Write;

fn main() {
    PrettyTrace::new().on();
    let outs = export_code(0);
    for i in 0..outs.len() {
        let mut f = open_for_write_new![&outs[i].0];
        fwrite!(f, "{}", outs[i].1);
    }
}
