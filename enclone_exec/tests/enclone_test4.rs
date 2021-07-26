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
    const ENCLONE_SIZE: usize = 93925704;
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
    const REQUIRED_GI: f64 = 18.4775;
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
    const DUPPED_CRATES: usize = 29;
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

    // Don't allow any crate except enclone_main and enclone_tools to reach the enclone crate.

    let top = dir_list("..");
    for d in top.iter() {
        if d.starts_with("enclone") && d != "enclone_main" && d != "enclone_tools" {
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 32. Test that the Rust version used by GitHub Actions is the same as the Rust version
// that's in use.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_rust_version() {
    PrettyTrace::new().on();
    let new = Command::new("rustc")
        .arg("--version")
        .output()
        .expect(&format!("failed to execute test_cpu_usage"));
    if new.status.code() != Some(0) {
        eprint!(
            "\netest_rust_version: failed to execute, stderr =\n{}",
            strme(&new.stderr),
        );
        std::process::exit(1);
    }
    let version = strme(&new.stdout).between(" ", " ");
    let test_yaml = format!("../.github/workflows/test.yaml");
    let f = open_for_read![&test_yaml];
    let mut version_found = false;
    for line in f.lines() {
        let s = line.unwrap();
        if s.contains("rustup default ") {
            version_found = true;
            let ga_version = s.after("rustup default ");
            if ga_version != version {
                eprintln!(
                    "\nThe version of Rust you are running is {}, but GitHub Actions is \
                    using version {} in at\nleast one place.  \
                    Please update the file .github/workflows/test.yaml.\n",
                    version, ga_version,
                );
                std::process::exit(1);
            }
        }
    }
    if !version_found {
        eprintln!("\nFailed to find Rust version in .github/workflows/test.yaml.  Weird.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 33. Don't allow exit in code in some crates.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_exit() {
    PrettyTrace::new().on();
    let exit_free_crates = [
        "enclone",
        "enclone_args",
        "enclone_core",
        "enclone_help",
        "enclone_main",
        "enclone_print",
        "enclone_proto",
        "enclone_tail",
    ];
    for cname in exit_free_crates.iter() {
        let files = dir_list(&format!("../{}/src", cname));
        for f in files.iter() {
            if f.ends_with(".rs") {
                let g = open_for_read![&format!("../{}/src/{}", cname, f)];
                for line in g.lines() {
                    let s = line.unwrap();
                    if s.contains("process::exit") {
                        eprintln!("exit not allowed in crate {}, but is in file {}", cname, f);
                        std::process::exit(1);
                    }
                }
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 34. Require that authors for all crates are the same.  We do not track authors on a
// crate-by-crate basis, and always just list them all.  Also some people have contributed
// to enclone, but not to specific crates.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_authors() {
    PrettyTrace::new().on();
    let dirs = dir_list(&format!("../.."));
    let mut aud = Vec::<(String, Vec<String>)>::new();
    for d in dirs.iter() {
        let toml = format!("../../{}/Cargo.toml", d);
        if path_exists(&toml) {
            let f = open_for_read![&toml];
            let mut au = Vec::<String>::new();
            let mut started = false;
            for line in f.lines() {
                let s = line.unwrap();
                if s.starts_with("authors =") {
                    started = true;
                }
                if started {
                    au.push(s.to_string());
                }
                if s.contains("]") {
                    break;
                }
            }
            aud.push((d.clone(), au));
        }
    }
    for i in 1..aud.len() {
        if aud[i].1 != aud[0].1 {
            eprintln!(
                "\nThe author lists for crates {} and {} are different.\n\
                They are required to be the same.\n",
                aud[0].0, aud[i].0,
            );
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 35. Test HONEY_OUT and HONEY_IN.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_honey() {
    let cmd1 = "BCR=123085:123089 PLOT=\"gui_stdout,s1->blue,s2->red\" \
        HONEY_OUT=testx/outputs/honey NOPRINT";
    let cmd2 = "BCR=123085:123089 PLOT_BY_ISOTYPE=gui_stdout HONEY_IN=testx/outputs/honey NOPRINT";
    let args1 = cmd1.split(' ').collect::<Vec<&str>>();
    let args2 = cmd2.split(' ').collect::<Vec<&str>>();
    let new1 = Command::new(env!("CARGO_BIN_EXE_enclone"))
        .args(&args1)
        .output()
        .expect(&format!("failed to execute test_honey 1"));
    if new1.status.code() != Some(0) {
        eprintln!(
            "\ntest_honey 1: failed to execute, stderr =\n{}",
            strme(&new1.stderr),
        );
        std::process::exit(1);
    }
    let new2 = Command::new(env!("CARGO_BIN_EXE_enclone"))
        .args(&args2)
        .output()
        .expect(&format!("failed to execute test_honey 2"));
    if new2.status.code() != Some(0) {
        eprintln!(
            "\ntest_honey 2: failed to execute, stderr =\n{}",
            strme(&new2.stderr),
        );
        std::process::exit(1);
    }
    let out = strme(&new2.stdout);
    let expected = include_str!["../testx/inputs/outputs/test_honey.svg"];
    if out != expected {
        eprintln!("\ntest honey yielded changed answer\n");
        std::process::exit(1);
    }
}
