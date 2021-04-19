// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

use ansi_escape::*;
use enclone::html::*;
use enclone_core::defs::*;
use enclone_core::run_test::*;
use enclone_core::testlist::*;
use enclone_core::*;
use enclone_proto::proto_io::{read_proto, ClonotypeIter};
use enclone_proto::types::EncloneOutputs;
use failure::Error;
use flate2::read::GzDecoder;
use io_utils::*;
use itertools::Itertools;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use serde_json::Value;
use sha2::{Digest, Sha256};
use stats_utils::*;
use std::cmp::min;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::{metadata, read_dir, read_to_string, remove_dir_all, remove_file, File};
use std::io;
use std::io::prelude::*;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::process::{Command, Stdio};
use std::thread;
use std::time;
use std::time::{Duration, Instant};
use string_utils::*;
use vector_utils::*;

const LOUPE_OUT_FILENAME: &str = "testx/__test_proto";

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test peak memory.  This is designed for one server at 10x Genomics
//
// This is in a separate file because doing so makes it possible to run just this one test,
// and have it run by ./test, and have it not run by the CI.
//
// You can run this test alone by typing
// cargo test --test enclone_test_peak_mem -- --nocapture
// from the root of the repo.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_peak_memory() {
    PrettyTrace::new().on();

    // Specify mem requirements.

    let dataset = "BCR=123085";
    let expected_mb = 367.1;

    let max_percent_dev = 0.5;

    // Only run internally.

    let mut internal_run = false;
    for (key, value) in env::vars() {
        if key.contains("TELEPORT") && value.contains("10xgenomics.com") {
            internal_run = true;
        }
    }
    if !internal_run {
        return;
    }

    // Define function to get mem.

    fn mem_val(dataset: &str) -> f64 {
        let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
            .arg(&dataset)
            .args(&["NOPRINT", "COMP"])
            .output()
            .expect(&format!("failed to execute peak mem test"));
        if new.status.code() != Some(0) {
            eprint!(
                "\nenclone_site_examples: nonzero status code for peak mem test, stderr =\n{}",
                strme(&new.stderr),
            );
            std::process::exit(1);
        }
        let out = strme(&new.stdout);
        let mut mem = None;
        for line in out.lines() {
            if line.starts_with("peak mem usage = ") && line.ends_with(" MB") {
                mem = Some(line.between("peak mem usage = ", " MB").force_f64());
            }
        }
        if mem.is_none() {
            eprintln!("\nfailed to find mem usage for peak mem\n");
            std::process::exit(1);
        }
        mem.unwrap()
    }

    // Run the test.

    let mut ok = false;
    let mut mems = Vec::<f64>::new();
    for _tries in 0..3 {
        let mem = mem_val(&dataset);
        mems.push(mem);
        let dev = 100.0 * (mem - expected_mb).abs() / expected_mb;
        if dev <= max_percent_dev {
            ok = true;
            break;
        }
    }
    if !ok {
        for _tries in 3..10 {
            mems.push(mem_val(&dataset));
        }
    }
    let mut mean = 0.0;
    for i in 0..mems.len() {
        mean += mems[i];
    }
    mean /= mems.len() as f64;
    let dev = 100.0 * (mean - expected_mb).abs() / expected_mb;
    let msg = format!(
        "Peak mem of {:.1} MB observed, which differs from expected value of {} by {:.2}%.",
        mean, expected_mb, dev
    );
    if !ok && dev > max_percent_dev {
        eprintln!("\n{}\n", msg);
        eprintln!("observed values:\n");
        for i in 0..mems.len() {
            eprintln!("{:.1}", mems[i]);
        }
        eprintln!(
            "Please note that this test was designed to work correctly from a single server\n\
            at 10x Genomics.  \
            If you're running from a different server, the expected memory value\n\
            may need to be changed.  This might also be the case if that server was changed.\n\
            Otherwise, your options are:\n\
            1. Change the value of expected_mb in enclone_test.rs to {:.1}.\n\
            2. Optimize to reduce mem usage (if value exceeds expected).\n\
            3. Repeat the test, but this is very unlikely to succeed.\n",
            mean
        );
        std::process::exit(1);
    } else {
        println!("{}", msg);
    }
}
