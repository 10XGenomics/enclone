// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Update the output of all main tests.
//
// Do not do this unless you're highly confident that it's safe to do so, without
// manually examining the outputs.
//
// Note duplication of some code with enclone/tests/enclone_test.rs.

use enclone::testlist::*;
use io_utils::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::Command;
use string_utils::*;

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
        let mut test = TESTS[it].to_string();
        test = test.replace("\n", "");
        for _ in 0..3 {
            test = test.replace("  ", " ");
        }
        let out_file = format!("test/inputs/outputs/enclone_test{}_output", it + 1);

        // Get arguments, by parsing command, breaking at blanks, but not if they're in quotes.
        // This is identical to parse_csv, except for the splitting character.
        // Should refactor.

        let mut args = Vec::<String>::new();
        let mut w = Vec::<char>::new();
        for c in test.chars() {
            w.push(c);
        }
        let (mut quotes, mut i) = (0, 0);
        while i < w.len() {
            let mut j = i;
            while j < w.len() {
                if quotes % 2 == 0 && w[j] == ' ' {
                    break;
                }
                if w[j] == '"' {
                    quotes += 1;
                }
                j += 1;
            }
            let (mut start, mut stop) = (i, j);
            if stop - start >= 2 && w[start] == '"' && w[stop - 1] == '"' {
                start += 1;
                stop -= 1;
            }
            let mut s = String::new();
            for m in start..stop {
                s.push(w[m]);
            }
            args.push(s);
            i = j + 1;
        }

        // Form the command and execute it.

        let mut new = Command::new("target/release/enclone");
        let mut new = new.arg(format!("PRE=test/inputs/version{}", TEST_FILES_VERSION));
        for i in 0..args.len() {
            new = new.arg(&args[i]);
        }
        let new = new
            .arg("FORCE_EXTERNAL")
            .arg("MAX_CORES=24")
            .output()
            .expect(&format!("failed to execute enclone for test{}", it + 1));
        let new2 = stringme(&new.stdout);
        let mut f = open_for_write_new![&out_file];
        fwrite!(f, "{}", new2);
    });
}
