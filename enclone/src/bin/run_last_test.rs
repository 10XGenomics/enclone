// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Run just the last test in the TESTS category in testlist.rs.

use enclone::run_test::*;
use enclone_core::testlist::*;
use pretty_trace::*;

fn main() {
    PrettyTrace::new().on();
    let mut out = String::new();
    let mut ok = false;
    let mut logx = String::new();
    run_test(
        "enclone",
        TESTS.len() - 1,
        "",
        &TESTS[TESTS.len() - 1],
        "test",
        &mut ok,
        &mut logx,
        &mut out,
    );
    print!("{}", logx);
}
