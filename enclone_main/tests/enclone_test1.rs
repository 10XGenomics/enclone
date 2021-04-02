// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

// There are three categories of tests here:
// 1. basic tests (feature = basic), runs without additional data requirements
// 2. nonbasic tests, requires extended dataset distributed with enclone
// 3. speed test (feature = cpu), requires non-public datasets.

use ansi_escape::*;
use enclone::html::*;
use enclone::misc3::parse_bsv;
use enclone::run_test::*;
use enclone_core::testlist::*;
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

fn valid_link(link: &str) -> bool {
    use attohttpc::*;
    let req = attohttpc::get(link.clone()).read_timeout(Duration::new(10, 0));
    let response = req.send();
    if response.is_err() {
        return false;
    } else {
        let response = response.unwrap();
        if response.is_success() {
            return true;
        }
        return false;
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 1. Test for redundancy of columns in parseable output.  While we're at it, test that if POUT
// is set to a file name, then variables of the form abbr:name use the field label abbr.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_for_parseable_redundancy() {
    let test = r###"BCR=123085 GEX=123217 LVARSP="IG%:IG.*_g_%" MIN_CHAINS_EXACT=2 CDR3=CAREGGVGVVTATDWYFDLW POUT=testx/outputs/redundancy_out"###;
    let pre_arg = format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    );
    let args = parse_bsv(&test);
    let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
        .arg(&pre_arg)
        .args(&args)
        .output()
        .expect(&format!("failed to execute test_for_parseable_redundancy"));
    if new.status.code() != Some(0) {
        eprint!(
            "\neparseable redundancy test: failed to execute, stderr =\n{}",
            strme(&new.stderr),
        );
        std::process::exit(1);
    }
    let f = open_for_read!["testx/outputs/redundancy_out"];
    let mut found = false;
    for line in f.lines() {
        let s = line.unwrap();
        let mut fields = s.split(',').collect::<Vec<&str>>();
        fields.sort();
        assert!(fields.len() > 0);
        let mut i = 0;
        while i < fields.len() {
            let j = next_diff(&fields, i);
            if j - i > 1 {
                eprintln!(
                    "\nParseable output field {} appears more than once.\n",
                    fields[i]
                );
                std::process::exit(1);
            }
            i = j;
        }
        for i in 0..fields.len() {
            if fields[i] == "IG%" {
                found = true;
            }
        }
        if !found {
            eprintln!("\nParseable output abbreviated field IG% not found.\n");
        }
        break;
    }
    let _ = remove_file("testx/outputs/redundancy_out");
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 2.

// NOT BASIC

// Make sure all help pages have been edited.
// Also make sure that they do not contain backslashes, which are likely errors.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_help_pages_edited() {
    let all = read_dir("../pages/auto").unwrap();
    let mut fail = false;
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        if f.contains(".help.") {
            let mut edited = false;
            let h = open_for_read![&format!("{}", f)];
            for line in h.lines() {
                let s = line.unwrap();
                if s.contains("googletag") {
                    edited = true;
                }
            }
            if !edited {
                eprintln!(
                    "\nThe page {} has not been edited.  Please run ./build.\n",
                    f
                );
                std::process::exit(1);
            }
        }
        let h = open_for_read![&format!("{}", f)];
        for line in h.lines() {
            let s = line.unwrap();
            if s.contains("\\") {
                eprintln!("\nIllegal backslash in {}, line is:\n{}\n", f, s);
                fail = true;
            }
        }
    }
    if fail {
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 3.

// NOT BASIC

// Make sure that two css files still exist.  These can never be deleted because they are
// accessed by certain html output of enclone.  These files could be out in the wild and we
// don't want to break them.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_css_existence() {
    let _ = include_str!["../../pages/enclone.css"];
    let _ = include_str!["../../pages/enclone_css_v2.css"];
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 4. Make sure that if sync_master was run, nothing would change.
//
// A bit ugly because of duplicated code.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_sync_master() {
    let mut version = HashMap::<String, String>::new();
    let f = open_for_read!["../master.toml"];
    for line in f.lines() {
        let s = line.unwrap();
        if !s.starts_with('#') && s.contains("=") {
            version.insert(s.before(" = ").to_string(), s.after(" = ").to_string());
        }
    }
    let all = read_dir("..").unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        let toml = format!("{}/Cargo.toml", f);
        if path_exists(&toml) {
            let g = open_for_read![&toml];
            for line in g.lines() {
                let s = line.unwrap();
                if s.contains(" =") {
                    let cratex = s.before(" =").to_string();
                    if version.contains_key(&cratex) {
                        let t = format!("{} = {}", cratex, version[&cratex]);
                        if t != s {
                            eprintln!("\nFound change in {}.\nold: {}\nnew: {}", toml, s, t);
                            eprintln!("You probably need to run sync_to_master\n");
                            std::process::exit(1);
                        }
                    }
                }
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 5.

// Run the download command on the landing page and make sure it works.
//
// Runs with "small", and passes second argument so we can put outputs in a defined place.
//
// There are three passes.  The first pass tests the copy of install.sh that is one master, and
// the second pass tests the local version; the third tests with wget forced.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_curl_command() {
    let mut internal_run = false;
    for (key, value) in env::vars() {
        if key.contains("TELEPORT") && value.contains("10xgenomics.com") {
            internal_run = true;
        }
    }
    if internal_run {
        if !path_exists("testx/outputs") {
            eprintln!(
                "\ntest_curl_command:\n\
                You need to create the directory enclone_main/testx/outputs.\n\
                If you run \"./build\" this will be done for you.\n"
            );
            std::process::exit(1);
        }
        for pass in 1..=3 {
            for f in ["enclone", "bin", ".profile", ".subversion"].iter() {
                let g = format!("testx/outputs/{}", f);
                if path_exists(&g) {
                    if !metadata(&g).unwrap().is_dir() {
                        remove_file(&g).unwrap();
                    } else {
                        remove_dir_all(&g).unwrap();
                    }
                }
            }
            let command;
            let version;
            if pass == 1 {
                command = "curl -sSf -L bit.ly/enclone_install | bash -s small testx/outputs";
                version = "master";
            } else if pass == 2 {
                command = "cat ../install.sh | bash -s small testx/outputs";
                version = "local";
            } else {
                command = "cat ../install.sh | bash -s small testx/outputs force_wget";
                version = "local";
            }
            let o = Command::new("sh").arg("-c").arg(&command).output().unwrap();
            if o.status.code().unwrap() != 0 {
                eprintln!(
                    "\nAttempt to run enclone install command using {} version of \
                    install.sh failed.\n",
                    version
                );
                eprint!("stdout:\n{}", strme(&o.stdout));
                eprint!("stderr:\n{}", strme(&o.stderr));
                std::process::exit(1);
            }
            // The following test is there because a bash script can fail without setting
            // a nonzero exit code.
            if strme(&o.stderr).starts_with("bash:") {
                eprintln!(
                    "\nOn pass {}, it would appear that the install script failed, because \
                    \"bash:\" appears in the stderr.\n\
                    Here is stderr:\n\n{}",
                    pass,
                    strme(&o.stderr)
                );
                std::process::exit(1);
            }
            let req = [
                "bin/enclone",
                "enclone/datasets/123085/outs/all_contig_annotations.json.lz4",
                "enclone/datasets_small_checksum",
                "enclone/version",
            ];
            for (jf, f) in req.iter().enumerate() {
                if !path_exists(&format!("testx/outputs/{}", f)) {
                    eprintln!(
                        "\nAttempt to run enclone install command using {} version of \
                        install.sh failed to fetch {}.\n",
                        version, f
                    );
                    eprintln!("results of install command:");
                    eprint!("stdout:\n{}", strme(&o.stdout));
                    eprint!("stderr:\n{}", strme(&o.stderr));
                    std::process::exit(1);
                }
                // Make sure that the all contigs file is not "essentially empty".  This happened.

                if pass >= 2 && jf == 1 {
                    let p = format!("testx/outputs/{}", f);
                    let len = metadata(&p).unwrap().len();
                    if len < 10_000_000 {
                        eprintln!(
                            "\nDownload yielded truncated all_contig_annotations.json.lz4 \
                            file, size = {} bytes.\n",
                            len
                        );
                        eprintln!("The command was\n{}\n", command);
                        if len < 100 {
                            eprintln!("file contents = ");
                            let f = open_for_read![&p];
                            for line in f.lines() {
                                let s = line.unwrap();
                                eprintln!("{}", s);
                            }
                            eprintln!("");
                        }
                        std::process::exit(1);
                    }
                }
            }
            if path_exists("testx/outputs/.subversion") {
                eprintln!(
                    "\nAttempt to run enclone install command using {} version of \
                    install.sh created .subversion.\n",
                    version
                );
                std::process::exit(1);
            }
            for f in ["enclone", "bin", ".profile", ".subversion"].iter() {
                let g = format!("testx/outputs/{}", f);
                if path_exists(&g) {
                    if *f == ".profile" {
                        remove_file(&g).unwrap();
                    } else {
                        remove_dir_all(&g).unwrap();
                    }
                }
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 6.

// NOT BASIC

// Make sure that the dataset checksum files are current.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_datasets_sha256() {
    let sha_command1 = format!(
        "git -C ../enclone-data write-tree --prefix=big_inputs/version{}",
        TEST_FILES_VERSION
    );
    let sha_command2 = "cat ../datasets_medium_checksum";
    let sha1 = Command::new("sh")
        .arg("-c")
        .arg(&sha_command1)
        .output()
        .unwrap();
    let sha1_status = sha1.status.code().unwrap();
    if sha1_status != 0 {
        eprintln!(
            "\nsha_command1 = {}\nfailed for datasets_medium_checksum\n",
            sha_command1
        );
        std::process::exit(1);
    }
    let sha1 = sha1.stdout;
    let sha2 = Command::new("sh")
        .arg("-c")
        .arg(&sha_command2)
        .output()
        .unwrap()
        .stdout;
    if sha1 != sha2 {
        eprintln!(
            "\nThe file datasets_medium_checksum is not current.  You can update it by typing\n\
            ./build\ndatasets_medium_checksum = {}\ncomputed sha             = {}",
            strme(&sha2),
            strme(&sha1),
        );
        std::process::exit(1);
    }
    let sha_command1 = format!(
        "git -C ../enclone-data write-tree --prefix=big_inputs/version{}/123085",
        TEST_FILES_VERSION
    );
    let sha_command2 = "cat ../datasets_small_checksum";
    let sha1 = Command::new("sh")
        .arg("-c")
        .arg(&sha_command1)
        .output()
        .unwrap()
        .stdout;
    let sha2 = Command::new("sh")
        .arg("-c")
        .arg(&sha_command2)
        .output()
        .unwrap()
        .stdout;
    if sha1 != sha2 {
        eprintln!(
            "\nThe file datasets_small_checksum is not current.  You can update it by typing\n\
            ./build\ndatasets_small_checksum = {}\ncomputed sha            = {}",
            strme(&sha2),
            strme(&sha1),
        );
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 7.

// SPEED (AND NOT BASIC)
// calibrated for one server at 10x Genomics, and requires linux
// cargo test --test enclone_test --features cpu -- --nocapture
// from enclone_main directory
// or just ./speed from root directory

#[cfg(not(feature = "basic"))]
#[cfg(feature = "cpu")]
#[test]
fn test_cpu() {
    PrettyTrace::new().on();

    // Introductory comments.

    println!(
        "\nSPEED TESTS\n\n\
        • These are calibrated for a particular server at \
        10x Genomics.  If this code is run\nusing a different server, or if that server is \
        changed, the tests will need to be recalibrated.\n\
        • These tests also use 10x Genomics datasets that are not distributed publicly\n\
        (although perhaps could be).\n\
        • Finally note that the datasets \
        themselves could be changed without changing this code,\nand that could affect results."
    );
    println!(
        "\nThese tests are expected to fail intermittently simply because of stochastic variation\n\
        in computational performance (or competing load on the server).  Note also that they\n\
        depend on files being in a cached state."
    );

    // Speed test 1.

    let it = 1;
    let test =
        "BI=10 NCROSS NGEX NOPRINT PRINT_CPU NCORES BUILT_IN EXPECT_OK EXPECT_NULL NO_PRE NFORCE";
    let expect = 16638;
    let percent_dev = 6.0;
    println!("\nSpeed test 1");
    println!(
        "\nThis tests cpu cycles.  If the code is parallelized better, this test may get \n\
        slower.  Such changes should be accepted if they reduce wallclock."
    );
    let mut out = String::new();
    let mut ok = false;
    let mut log = String::new();
    let mut cpu_all_start = 0;
    {
        let f = open_for_read!["/proc/stat"];
        for line in f.lines() {
            let s = line.unwrap();
            let mut t = s.after("cpu");
            while t.starts_with(' ') {
                t = t.after(" ");
            }
            cpu_all_start = t.before(" ").force_usize();
            break;
        }
    }
    run_test(
        env!("CARGO_BIN_EXE_enclone"),
        it,
        "",
        &test,
        "cpu",
        &mut ok,
        &mut log,
        &mut out,
    );
    let this_used = out.before("\n").force_usize();
    let mut cpu_all_stop = 0;
    {
        let f = open_for_read!["/proc/stat"];
        for line in f.lines() {
            let s = line.unwrap();
            let mut t = s.after("cpu");
            while t.starts_with(' ') {
                t = t.after(" ");
            }
            cpu_all_stop = t.before(" ").force_usize();
            break;
        }
    }
    let all_used = cpu_all_stop - cpu_all_start;
    let dev = 100.0 * (this_used as f64 - expect as f64) / (expect as f64);
    println!(
        "\nused cpu = {} = {:.1}% of total, dev = {:.1}%\n",
        this_used,
        percent_ratio(this_used, all_used),
        dev
    );
    if dev.abs() > percent_dev {
        eprintln!("cpu deviation exceeded max of {}%\n", percent_dev);
        std::process::exit(1);
    }

    // Speed test 2.

    let it = 2;
    let test =
        "BI=1-2,5-12 MIX_DONORS NOPRINT PRINT_CPU NCORES BUILT_IN EXPECT_OK EXPECT_NULL NO_PRE NFORCE";
    let expect = 137.6;
    let percent_dev = 6.0;
    println!("Speed test 2");
    println!(
        "\nThis tests wall clock.  It is thus particularly susceptible to competing load \
        on the server.\nIt will also be very slow unless it has been run recently so files \
        are in cache.\nThis test takes about two minutes and will trigger a warning from cargo \
        after 60 seconds.\n"
    );
    let t = Instant::now();
    let mut out = String::new();
    let mut ok = false;
    let mut log = String::new();
    run_test(
        env!("CARGO_BIN_EXE_enclone"),
        it,
        "",
        &test,
        "cpu",
        &mut ok,
        &mut log,
        &mut out,
    );
    let this_used = elapsed(&t);
    let dev = 100.0 * (this_used as f64 - expect as f64) / (expect as f64);
    println!(
        "\nused wallclock = {:.1} seconds, dev = {:.1}%\n",
        this_used, dev
    );
    if dev.abs() > percent_dev {
        eprintln!("cpu deviation exceeded max of {}%\n", percent_dev);
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 8.

// NOT BASIC

// Test licenses of included packages and their dependencies.
//
// The following rules are applied:
// 1. If the license field in Cargo.toml is set to MIT or ISC or Zlib or WTFPL or MPL-2.0
//    or CC0-1.0, or is a logical expression for which one of those is sufficient, then there is
//    no problem.  Note that for MPL-2.0, we inform people how to get the source code for
//    dependent crates.
// 2. If both license and license_field are null, then there is no problem.
// 3. If the license field is Apache-2.0, or a logical expression for which that is sufficient,
//    and there is no NOTICE file, then there is no problem.  Note that we include the
//    Apache-2.0 license as part of this repo in third_party.
// 4. If the package is owned by 10x, then there is no problem.
// 5. arrayref and cloudabi OK because we've included the license for it.
// 6. fuchsia-cprng OK because Cargo.toml refers to a BSD-style license, in a file LICENSE,
//    at https://fuchsia.googlesource.com/fuchsia/+/master/LICENSE, which we include in
//    third_party under fuchsia.
// 7. ring OK because we acknowledge OpenSSL in the file acknowledgements and because we include
//    the ring license.
// 8. webpki OK because we include the webpki license and also that for chromium.
// 9. instant OK because we include the license.
// 10. nalgebra OK because currently licensed under Apache-2.0 and there is not NOTICE file.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_licenses() {
    const ACCEPTABLE_LICENSE_TYPES: [&str; 6] =
        ["MIT", "ISC", "Zlib", "WTFPL", "MPL-2.0", "CC0-1.0"];
    const A2: &str = "Apache-2.0";
    const ACCEPTABLE_10X_PACKAGES: [&str; 2] = ["exons", "vdj_ann"];
    const ACCEPTABLE_OTHER_PACKAGES: [&str; 7] = [
        "arrayref",
        "cloudabi",
        "fuchsia-cprng",
        "instant",
        "nalgebra",
        "ring",
        "webpki",
    ];
    let new = Command::new("cargo-license").arg("-d").arg("-j").output();
    if new.is_err() {
        eprintln!(
            "\nFailed to execute cargo-license.  This means that either you have not \
            installed cargo-license,\nor that you have not added it to your PATH.  \
            To install it, type:\n\
            cargo install cargo-license\n\
            When it is done installing, it will tell you where it put the binary, and you\n\
            should add that path to your PATH.\n\
            You can also avoid this test entirely by running instead \
            \"cd enclone; cargo test basic -- --nocapture\".\n"
        );
        std::process::exit(1);
    }
    let lic = &new.unwrap().stdout;
    let mut f = &lic[..];
    let mut fails = Vec::<String>::new();
    loop {
        match read_vector_entry_from_json(&mut f) {
            None => break,
            Some(x) => {
                let v: Value = serde_json::from_str(strme(&x)).unwrap();
                let package = v["name"].to_string().between("\"", "\"").to_string();
                let version = v["version"].to_string().between("\"", "\"").to_string();
                let mut license = String::new();
                if v.get("license").is_some() {
                    license = v["license"].to_string();
                    if license.contains('"') {
                        license = license.between("\"", "\"").to_string();
                    }
                }
                let mut license_file = String::new();
                if v.get("license_file").is_some() {
                    license_file = v["license_file"].to_string();
                    if license_file.contains('"') {
                        license_file = license_file.between("\"", "\"").to_string();
                    }
                }
                if license == "null" && license_file == "null" {
                    continue;
                }
                let mut repo = String::new();
                if v.get("repository").is_some() {
                    repo = v["repository"].to_string();
                    if repo.contains('"') {
                        repo = repo.between("\"", "\"").to_string();
                    }
                }
                let mut ok = false;
                if package.starts_with("enclone") {
                    ok = true;
                }
                for y in ACCEPTABLE_10X_PACKAGES.iter() {
                    if package == *y {
                        ok = true;
                    }
                }
                for y in ACCEPTABLE_OTHER_PACKAGES.iter() {
                    if package == *y {
                        ok = true;
                    }
                }
                for y in ACCEPTABLE_LICENSE_TYPES.iter() {
                    if license == *y {
                        ok = true;
                    }
                    if license.ends_with(&format!(" OR {}", y)) {
                        ok = true;
                    }
                    if license.starts_with(&format!("{} OR ", y)) {
                        ok = true;
                    }
                }
                if !ok {
                    let (mut x1, mut x2) = (false, false);
                    if repo.starts_with("https://github.com") {
                        let f1 = format!("{}/blob/master/Cargo.toml", repo);
                        if valid_link(&f1) {
                            x1 = true;
                        }
                        let f2 = format!("{}/blob/master/NOTICE", repo);
                        if valid_link(&f2) {
                            x2 = true;
                        }
                    }
                    let a2 = license == A2
                        || license.ends_with(&format!(" OR {}", A2))
                        || license.starts_with(&format!("{} OR ", A2));
                    if a2 && x1 && !x2 {
                        continue;
                    }
                    fails.push(format!("{}, {}, {}, {}", package, version, license, repo));
                }
            }
        }
    }
    if fails.len() > 0 {
        fails.sort();
        let mut msg = format!("\nLicense check failed.  The following packages had problems:\n");
        for i in 0..fails.len() {
            msg += &format!("{}. {}\n", i + 1, fails[i]);
        }
        eprintln!(
            "{}\nYou may want to retry the test, since the license checks fails sporadically \
            at a low rate.\n",
            msg
        );
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 9. Test that files are rustfmt'ed.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_formatting() {
    let new = Command::new("cargo-fmt")
        .arg("--all")
        .arg("--")
        .arg("--check")
        .output()
        .expect(&format!("failed to execute test_formatting"));
    if new.status.code().unwrap() != 0 {
        eprintln!("\nYou need to run rustfmt.\n");
        eprintln!("{}\n", strme(&new.stdout));
        std::process::exit(1);
    }
}
