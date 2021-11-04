// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Actually export code.

// Read the vars file and export code.  This is a partial implementation.

use enclone_vars::export_code::*;
use pretty_trace::PrettyTrace;

fn main() {
    PrettyTrace::new().on();
    export_code();
}
