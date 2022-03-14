// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

use enclone_core::defs::*;
use enclone_main::main_enclone::main_enclone;
use enclone_ranger::main_enclone::main_enclone_ranger;
use enclone_vars::export_code::*;
use enclone_vars::var::*;
use enclone_vars::*;
use io_utils::*;
use pretty_trace::*;
use stats_utils::*;
use std::fs::{metadata, read_to_string, remove_file, rename, File};
use std::io::{BufRead, BufReader, Read};
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
    const ENCLONE_SIZE: usize = 100777216;
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
//
// Using BUILT_IN because it increases code coverage.

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
        "BUILT_IN",
        "TREE",
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
    const REQUIRED_GI: f64 = 69.8962;
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
// Exception ..._auto.rs, for now.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_source_code_file_length() {
    PrettyTrace::new().on();
    const MAX_RS_LINES: usize = 1000;
    let top = dir_list("..");
    let mut dirs = Vec::<String>::new();
    for d in top.iter() {
        if d != "enclone_denovo" {
            for x in ["src", "src/bin", "tests"].iter() {
                let d = format!("../{}/{}", d, x);
                if path_exists(&d) {
                    dirs.push(d.clone());
                }
            }
        }
    }
    let mut fail = false;
    for d in dirs.iter() {
        let fs = dir_list(d);
        for x in fs.iter() {
            if x.ends_with(".rs") && !x.ends_with("_auto.rs") {
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
    const DUPPED_CRATES: usize = 39;
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

#[cfg(not(feature = "cpu"))]
#[test]
fn test_dependency_structure() {
    // Restrict crates reached by enclone_core.

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

    // Restrict crates reaching enclone.

    let top = dir_list("..");
    for d in top.iter() {
        if d.starts_with("enclone")
            && d != "enclone_main"
            && d != "enclone_tools"
            && d != "enclone_denovo"
            && d != "enclone_ranger"
            && d != "enclone_stuff"
        {
            let toml = format!("../{}/Cargo.toml", d);
            if path_exists(&toml) {
                let f = open_for_read![&toml];
                for line in f.lines() {
                    let s = line.unwrap();
                    if s.starts_with("enclone =") {
                        eprintln!(
                            "\nThe crate {} has the crate enclone as a dependency.  In an \
                            attempt to reduce\ncompile time, we only allow this for certain \
                            crates.\n",
                            d
                        );
                        std::process::exit(1);
                    }
                }
            }
        }
    }

    // Restrict crates reaching enclone_tail.

    let top = dir_list("..");
    for d in top.iter() {
        if d.starts_with("enclone")
            && d != "enclone_main"
            && d != "enclone_visual"
            && d != "enclone_paper"
        {
            let toml = format!("../{}/Cargo.toml", d);
            if path_exists(&toml) {
                let f = open_for_read![&toml];
                for line in f.lines() {
                    let s = line.unwrap();
                    if s.starts_with("enclone_tail =") {
                        eprintln!(
                            "\nThe crate {} has the crate enclone_tail as a dependency.  In an \
                            attempt to reduce\ncompile time, we only allow this for certain \
                            crates.\n",
                            d
                        );
                        std::process::exit(1);
                    }
                }
            }
        }
    }

    // Restrict crates reaching enclone_help.

    let top = dir_list("..");
    for d in top.iter() {
        if d.starts_with("enclone") && d != "enclone_main" {
            let toml = format!("../{}/Cargo.toml", d);
            if path_exists(&toml) {
                let f = open_for_read![&toml];
                for line in f.lines() {
                    let s = line.unwrap();
                    if s.starts_with("enclone_help =") {
                        eprintln!(
                            "\nThe crate {} has the crate enclone_help as a dependency.  In an \
                            attempt to reduce\ncompile time, we only allow this for certain \
                            crates.\n",
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
                        let mut chars = Vec::<char>::new();
                        for c in s.chars() {
                            chars.push(c);
                        }
                        let mut comment = false;
                        for i in 0..chars.len() - 1 {
                            if chars[i] == ' ' {
                            } else if chars[i] == '/' && chars[i + 1] == '/' {
                                comment = true;
                            } else {
                                break;
                            }
                        }
                        if !comment {
                            eprintln!("exit not allowed in crate {}, but is in file {}", cname, f);
                            std::process::exit(1);
                        }
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
//
// Do this to rebuild if needed:
// enclone BCR=123085:123089 PLOT="gui_stdout,s1->blue,s2->red" HONEY_OUT=enclone_exec/testx/outputs/honey NOPRINT > /dev/null
// enclone BCR=123085:123089 PLOT_BY_ISOTYPE=gui_stdout HONEY_IN=enclone_exec/testx/outputs/honey NOPRINT > enclone_exec/testx/inputs/outputs/test_honey.svg

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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 36. Test to see if GitHub Actions scripts are valid YAML.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_yaml() {
    PrettyTrace::new().on();
    let content = include_str!["../../.github/workflows/release.yaml"];
    let yaml = yaml_rust::YamlLoader::load_from_str(&content);
    if yaml.is_err() {
        eprintln!("\nrelease.yaml is not valid YAML.\n");
        std::process::exit(1);
    }
    let content = include_str!["../../.github/workflows/test.yaml"];
    let yaml = yaml_rust::YamlLoader::load_from_str(&content);
    if yaml.is_err() {
        eprintln!("\ntest.yaml is not valid YAML.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 37. Test to see if vars file is sorted and formatted correctly.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_vars() {
    let old = std::fs::read_to_string("../enclone_vars/src/vars").unwrap();
    let new = sort_vars(&old);
    if new != old {
        eprintln!("\nPlease run var_sort to sort the variables file.\n");
        std::process::exit(1);
    }
    let _ = parse_variables(&old);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 38. Test to see that main_enclone_ranger and main_enclone do the same thing on one test case.
// An attempt to make this basic failed.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_ranger() {
    if path_exists("testx/outputs/tiny_multi_CS_6.1") {
        std::fs::remove_dir_all("testx/outputs/tiny_multi_CS_6.1").unwrap();
    }
    let options = fs_extra::dir::CopyOptions::new();
    fs_extra::dir::copy(
        "../enclone-data/big_inputs/version15/tiny_multi_CS_6.1",
        "testx/outputs",
        &options,
    )
    .unwrap();
    rename(
        "testx/outputs/tiny_multi_CS_6.1/outs/multi/vdj_b/all_contig_annotations.json.lz4",
        "testx/outputs/tiny_multi_CS_6.1/outs/multi/vdj_b/contig_annotations.json.lz4",
    )
    .unwrap();
    let proto_out = "testx/outputs/test1.proto";
    let donor_ref_out = "testx/outputs/test1.donor_ref.fasta";
    let mut files = vec![vec![Vec::<u8>::new(); 2]; 2];
    for pass in 1..=2 {
        for f in [&proto_out, &donor_ref_out].iter() {
            if path_exists(f) {
                remove_file(f).unwrap();
            }
        }
        let mut args = Vec::<String>::new();
        args.push("enclone".to_string());
        args.push("CELLRANGER".to_string());
        args.push("PRE=".to_string());
        args.push("FORCE_EXTERNAL".to_string());
        args.push("NOPAGER".to_string());
        args.push("NOPRINT".to_string());
        args.push("MAX_CORES=8".to_string());
        args.push("BCR=testx/outputs/tiny_multi_CS_6.1".to_string());
        args.push(
            "REF=../enclone-data/big_inputs/version15/tiny_multi_CS_6.1/outs/vdj_reference/\
            fasta/regions.fa"
                .to_string(),
        );
        args.push(format!("PROTO={}", proto_out));
        args.push(format!("DONOR_REF_FILE={}", donor_ref_out));
        if pass == 1 {
            main_enclone_ranger(&args).unwrap();
        } else {
            main_enclone(&args).unwrap();
        }
        if !path_exists(&proto_out) {
            eprintln!("pass = {}, can't find proto_out", pass);
            std::process::exit(1);
        }
        if !path_exists(&donor_ref_out) {
            eprintln!("pass = {}, can't find donor_ref_out", pass);
            std::process::exit(1);
        }
        let mut f = File::open(&proto_out).unwrap();
        f.read_to_end(&mut files[0][pass - 1]).unwrap();
        let mut f = File::open(&donor_ref_out).unwrap();
        f.read_to_end(&mut files[1][pass - 1]).unwrap();
    }
    for i in 0..2 {
        if files[i][0] != files[i][1] {
            eprintln!("files are different");
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 39. Test to see if there are unpushed commits.  The reason we do this is that one can
// accidentally merge a PR for a branch that has unpushed commits, and that causes grief.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_unpushed() {
    let new = Command::new("git")
        .arg("status")
        .output()
        .expect(&format!("failed to execute git status"));
    if new.status.code() != Some(0) {
        eprint!(
            "\netest_unpushed: failed to execute, stderr =\n{}",
            strme(&new.stderr),
        );
        std::process::exit(1);
    }
    let out = strme(&new.stdout);
    if out.contains("Your branch is ahead of") {
        eprintln!("\nYour branch has unpushed commits.  Please push them.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 40. Test to see that export_code has been run.
// An attempt to make this basic failed.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_export_code() {
    PrettyTrace::new().on();
    let outs = export_code(1);
    for i in 0..outs.len() {
        let f = format!("../{}", outs[i].0);
        let current = std::fs::read_to_string(&f).unwrap();
        if outs[i].1 != current {
            eprintln!("\nexport_code output {} has changed.\n", outs[i].0);
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 41. Test a funny command with a newline that asserted at one point.
// An attempt to make this basic failed.

// NOT BASIC

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_newline_bad() {
    let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
        .arg("METAX=\"bcr;86237;\n     null:83808\"")
        .output()
        .expect(&format!("failed to execute test_newline_bad"));
    if new.status.code() != Some(0) {
        eprintln!("\ntest_newline_bad: failed\n");
        eprintln!("stderr = {}\n", strme(&new.stderr));
        std::process::exit(1);
    }
}
