// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Update the output of all main tests.
//
// Do not do this unless you're highly confident that it's safe to do so, without
// manually examining the outputs.
//
// NOTE: you have to run this from the enclone_main directory.  Otherwise it won't work.

use enclone_core::run_test::*;
use enclone_core::testlist::*;
use io_utils::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

fn main() {
    PrettyTrace::new().on();
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..TESTS.len() {
        if !TESTS[i].contains("EXPECT_FAIL") && !TESTS[i].contains("EXPECT_OK") {
            results.push((i, false, String::new()));
        }
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let testname = &TESTS[it];
        let mut out = String::new();
        run_test(
            "enclone", it, "", &testname, "test", &mut res.1, &mut res.2, &mut out,
        );
        if !res.1 {
            let out_file = format!("testx/inputs/outputs/enclone_test{}_output", it + 1);
            let mut f = open_for_write_new![&out_file];
            fwrite!(f, "{}", out);
        }
    });
}
