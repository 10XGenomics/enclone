// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Walk through the main tests, show which have changed, and give option to update results.
// (First version, just prints output and pipes to less.)
//
// NOTE: you have to run this from the enclone_main directory.  Otherwise it won't work.

use enclone::misc1::setup_pager;
use enclone_core::main_testlist::*;
use enclone_tools::run_test::*;
use pretty_trace::*;
use rayon::prelude::*;

fn main() {
    PrettyTrace::new().on();
    setup_pager(true);
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..TESTS.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let mut out = String::new();
        run_test(
            "enclone",
            res.0,
            "",
            &TESTS[res.0],
            "test",
            &mut res.1,
            &mut res.2,
            &mut out,
        );
    });
    for i in 0..TESTS.len() {
        if !results[i].1 {
            print!("{}", results[i].2);
        }
    }
}
