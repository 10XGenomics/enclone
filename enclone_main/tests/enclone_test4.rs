// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

use enclone_core::defs::*;
use io_utils::*;
use pretty_trace::*;
use stats_utils::*;
use std::fs::{metadata, File};
use std::io::{BufRead, BufReader};
use std::process::Command;
use string_utils::*;
use vector_utils::*;

const LOUPE_OUT_FILENAME: &str = "testx/__test_proto";

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 27. Test Linux executable size.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_executable_size() {
    PrettyTrace::new().on();
    const ENCLONE_SIZE: usize = 80298672;
    const ENCLONE_SIZE_MAX_PER_DIFF: f64 = 1.0;
    let f = format!("../target/debug/enclone");
    let n = metadata(&f).unwrap().len() as usize;
    let delta = 100.0 * abs_diff(ENCLONE_SIZE, n) as f64 / ENCLONE_SIZE as f64;
    if delta > ENCLONE_SIZE_MAX_PER_DIFF {
        eprintln!(
            "\nenclone executable is only allowed to change by {}%, but it has changed \
            by {:.1}%.\n",
            ENCLONE_SIZE_MAX_PER_DIFF, delta,
        );
        eprintln!("old size = {}, new size = {}\n", ENCLONE_SIZE, n);
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 28. Test cpu usage.  This is designed for one server at 10x Genomics.  It runs
// single-threaded and measures total instructions used.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_cpu_usage() {
    PrettyTrace::new().on();
    let args = [
        "stat",
        "-e",
        "instructions:u",
        "enclone",
        "BCR=123085",
        "NOPRINT",
        "MAX_CORES=1",
    ];
    let new = Command::new("perf")
        .args(&args)
        .output()
        .expect(&format!("failed to execute test_cpu_usage"));
    if new.status.code() != Some(0) {
        eprint!(
            "\netest_cpu_usage: failed to execute, stderr =\n{}",
            strme(&new.stderr),
        );
        std::process::exit(1);
    }
    let mut gi = 0.0;
    let out = strme(&new.stderr);
    for line in out.lines() {
        let mut line = line.to_string();
        if line.contains("instructions:u") {
            for _ in 0..3 {
                line = line.replace("  ", " ");
            }
            line = line.replace(",", "");
            line = line.between(" ", " ").to_string();
            gi = line.force_f64() / 1_000_000_000.0;
        }
    }
    const REQUIRED_GI: f64 = 18.7520;
    let err = ((gi - REQUIRED_GI) / REQUIRED_GI).abs();
    let report = format!(
        "Observed GI = {:.4}, versus required GI = {:.4}, err = {:.2}%, versus max \
        allowed err = 0.10%.",
        gi,
        REQUIRED_GI,
        100.0 * err
    );
    if err > 0.001 {
        eprintln!("\n{}\n", report);
        eprintln!(
            "Possible causes of failure:\n\
            1. You running on a server other than the one this is designed for.  Won't work.\n\
            2. The server was changed.\n\
            3. A code change altered its performance.\n\
            4. You got very unlucky."
        );
        eprintln!("\nIf it makes sense, you can change REQUIRED_GI.\n");
        std::process::exit(1);
    } else {
        println!("{}", report);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 29. Test source code file length.  Cap the length in lines of the longest .rs file.  We do this
// because long files tend to increase compilation time.  They should be split up where possible.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_source_code_file_length() {
    PrettyTrace::new().on();
    const MAX_RS_LINES: usize = 1000;
    let top = dir_list("..");
    let mut dirs = Vec::<String>::new();
    for d in top.iter() {
        for x in ["src", "src/bin", "tests"].iter() {
            let d = format!("../{}/{}", d, x);
            if path_exists(&d) {
                dirs.push(d.clone());
            }
        }
    }
    let mut fail = false;
    for d in dirs.iter() {
        let fs = dir_list(d);
        for x in fs.iter() {
            if x.ends_with(".rs") {
                let y = format!("{}/{}", d, x);
                let f = open_for_read![&y];
                let mut n = 0;
                for _ in f.lines() {
                    n += 1;
                }
                if n > MAX_RS_LINES {
                    eprintln!(
                        "\nSource code file {} has {} lines, which exceeds the allowed max \
                        of {}.\n",
                        x, n, MAX_RS_LINES,
                    );
                    fail = true;
                }
            }
        }
    }
    if fail {
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 30. Cap number of duplicated crates.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_dupped_crates() {
    PrettyTrace::new().on();
    let mut crates = Vec::<String>::new();
    let f = open_for_read!["../Cargo.lock"];
    let mut lines = Vec::<String>::new();
    for line in f.lines() {
        let s = line.unwrap();
        lines.push(s.to_string());
    }
    for i in 0..lines.len() {
        if lines[i] == "[[package]]" {
            crates.push(lines[i + 1].clone());
        }
    }
    let n1 = crates.len();
    unique_sort(&mut crates);
    let n2 = crates.len();
    let d = n1 - n2;
    const DUPPED_CRATES: usize = 2;
    if d != DUPPED_CRATES {
        eprintln!(
            "\nThe number of duplicated crates is {}, but the required number is {}.\n",
            d, DUPPED_CRATES,
        );
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 30. Make sure that help page list is correct.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_help_page_list() {
    let (mut help1, mut help2) = (Vec::<String>::new(), Vec::<String>::new());
    let autox = dir_list("../pages/auto");
    for i in 0..autox.len() {
        if autox[i].starts_with("help") {
            help1.push(autox[i].between(".", ".").to_string());
        }
    }
    for x in HELP_PAGES.iter() {
        help2.push(x.to_string());
    }
    help1.sort();
    help2.sort();
    for x in help2.iter() {
        if !bin_member(&help1, &x) {
            eprintln!(
                "\nHelp page for {} is in HELP_PAGES but not in enclone/pages/auto.\n",
                x
            );
            std::process::exit(1);
        }
    }
    for x in help1.iter() {
        if !bin_member(&help2, &x) {
            eprintln!(
                "\nHelp page for {} is in enclone/pages/auto but not in HELP_PAGES.\n",
                x
            );
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 31. Test workspace dependency structure.  These restrictions are there to reduce compile time.
// To be expanded over time.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_dependency_structure() {
    // Don't allow enclone_core to reach to any other enclone crate.

    let f = include_str!["../../enclone_core/Cargo.toml"];
    for line in f.lines() {
        if !line.starts_with("name =")
            && line.starts_with("enclone")
            && !line.starts_with("enclone_proto")
        {
            eprintln!(
                "\nenclone_core should not depend on any other enclone crate\n\
                except enclone_proto.  This restriction is there to reduce compile time.\n"
            );
            std::process::exit(1);
        }
    }

    // Don't allow any crate except enclone_main to reach the enclone crate.

    let top = dir_list("..");
    for d in top.iter() {
        if d.starts_with("enclone") && d != "enclone_main" {
            let toml = format!("../{}/Cargo.toml", d);
            if path_exists(&toml) {
                let f = open_for_read![&toml];
                for line in f.lines() {
                    let s = line.unwrap();
                    if s.starts_with("enclone =") {
                        eprintln!(
                            "\nThe crate {} has the crate enclone as a dependency.  In an \
                            attempt to reduce\ncompile time, we only allow this for the crate \
                            enclone_main.\n",
                            d
                        );
                        std::process::exit(1);
                    }
                }
            }
        }
    }
}
