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
#[cfg(not(feature = "mem"))]
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
#[cfg(not(feature = "mem"))]
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
#[cfg(not(feature = "mem"))]
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
#[cfg(not(feature = "mem"))]
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
#[cfg(not(feature = "mem"))]
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
                    std::process::exit(1);
                }
                // Make sure that the all contigs file is not "essentially empty".  This happened.

                if pass >= 2 && jf == 1 {
                    let len = metadata(&format!("testx/outputs/{}", f)).unwrap().len();
                    if len < 10_000_000 {
                        eprintln!(
                            "\nDownload yielded truncated all_contig_annotations.json.lz4 \
                            file, size = {} bytes.\n",
                            len
                        );
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
#[cfg(not(feature = "mem"))]
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
// calibrated for bespin1, and requires linux
// cargo test --test enclone_test --features cpu -- --nocapture
// from enclone_main directory
// or just ./speed from root directory

#[cfg(not(feature = "basic"))]
#[cfg(feature = "cpu")]
#[cfg(not(feature = "mem"))]
#[test]
fn test_cpu() {
    // Introductory comments.

    println!(
        "\nSPEED TESTS\n\n\
        • These are calibrated for a particular server, bespin1 at \
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
    let test = "BI=10 NCROSS NGEX NOPRINT PRINT_CPU NCORES EXPECT_OK EXPECT_NULL NO_PRE NFORCE";
    let expect = 7700;
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
        "BI=1-2,5-12 MIX_DONORS NOPRINT PRINT_CPU NCORES EXPECT_OK EXPECT_NULL NO_PRE NFORCE";
    let expect = 59.0;
    let percent_dev = 6.0;
    println!("Speed test 2");
    println!(
        "\nThis tests wall clock.  It is thus particularly susceptible to competing load \
        on the server.\nIt will also be very slow unless it has been run recently so files \
        are in cache.\nThis test takes about a minute and may trigger a warning from cargo \
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
#[cfg(not(feature = "mem"))]
#[test]
fn test_licenses() {
    const ACCEPTABLE_LICENSE_TYPES: [&str; 6] =
        ["MIT", "ISC", "Zlib", "WTFPL", "MPL-2.0", "CC0-1.0"];
    const A2: &str = "Apache-2.0";
    const ACCEPTABLE_10X_PACKAGES: [&str; 6] = [
        "enclone",
        "enclone_print",
        "enclone_tail",
        "enclone_versions",
        "exons",
        "vdj_ann",
    ];
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
        eprintln!("{}", msg);
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 9. Test that files are rustfmt'ed.

#[cfg(not(feature = "cpu"))]
#[cfg(not(feature = "mem"))]
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 10. The following is a single test, containing many subtests, each of which is a regression test
// for a given enclone command line.
//
// If you ever need to change the output of all tests, use the main program
// update_all_main_tests.rs in enclone/src/bin.  Note that there is some duplicated code there.

#[cfg(not(feature = "cpu"))]
#[cfg(not(feature = "mem"))]
#[test]
fn test_enclone() {
    PrettyTrace::new().on();
    let t = Instant::now();
    println!("running tests using {}", env!("CARGO_BIN_EXE_enclone"));
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..TESTS.len() {
        results.push((i, false, String::new()));
    }
    let this = include_str!("../../enclone_core/src/testlist.rs");
    let mut tracking = false;
    let mut comments = Vec::<String>::new();
    let mut lines = Vec::<String>::new();
    for line in this.lines() {
        if line.starts_with("pub const TESTS: ") {
            tracking = true;
            continue;
        }
        if tracking {
            if line == "];" {
                break;
            }
            lines.push(line.to_string());
        }
    }
    let mut i = 0;
    while i < lines.len() {
        let mut j = i + 1;
        while j < lines.len() {
            if lines[j].starts_with("    // ") {
                let c = lines[j].after("    // ").as_bytes()[0];
                if c >= b'0' && c <= b'9' {
                    break;
                }
            }
            j += 1;
        }
        comments.push(lines[i..j].iter().format("\n").to_string());
        i = j;
    }

    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = TESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            &comments[it],
            &test,
            "test",
            &mut res.1,
            &mut res.2,
            &mut out,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ntotal time for {} enclone subtests = {:.2} seconds\n",
        TESTS.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 11.

// NOT BASIC

// Regression tests using the extended public dataset collection.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[cfg(not(feature = "mem"))]
#[test]
fn test_extended() {
    PrettyTrace::new().on();
    let t = Instant::now();
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..EXTENDED_TESTS.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = EXTENDED_TESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            "",
            &test,
            "ext_test",
            &mut res.1,
            &mut res.2,
            &mut out,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\nextended tests total time for {} enclone subtests = {:.2} seconds\n",
        EXTENDED_TESTS.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 12.

// NOT BASIC

// Crash tests.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[cfg(not(feature = "mem"))]
#[test]
fn test_crash() {
    PrettyTrace::new().on();
    let t = Instant::now();
    let mut crash_tests = Vec::<String>::new();
    for i in 0..CRASH_SETS.len() {
        crash_tests.push(format!("{} {} {}", CRASH_DATA, CRASH_SETS[i], CRASH_OPTS));
    }
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..crash_tests.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = &crash_tests[it];
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            "",
            &test,
            "crash_test",
            &mut res.1,
            &mut res.2,
            &mut out,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ncrash tests total time for {} enclone subtests = {:.2} seconds\n",
        crash_tests.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 13.

// NOT BASIC

// Regression tests for internal features.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[cfg(not(feature = "mem"))]
#[test]
fn test_internal() {
    PrettyTrace::new().on();
    let t = Instant::now();
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..INTERNAL_TESTS.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = INTERNAL_TESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
            "",
            &test,
            "internal_test",
            &mut res.1,
            &mut res.2,
            &mut out,
        );
    });
    for i in 0..results.len() {
        print!("{}", results[i].2);
        if !results[i].1 {
            std::process::exit(1);
        }
    }
    println!(
        "\ninternal tests total time for {} enclone subtests = {:.2} seconds\n",
        INTERNAL_TESTS.len(),
        elapsed(&t)
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 14.

// NOT BASIC

// Test site for broken links and spellcheck.
//
// Two approaches for checking broken links left in place for now, to delete one, and the
// corresponding crate from Cargo.toml.
//
// This looks for
// ▓<a href="..."▓
// ▓http:...[, '")}<#\n]▓
// ▓https:...[, '")}<#\n]▓
// ▓<img src="..."▓.
// (These also test termination by ". ".)
// SHOULD also look for at least:
// ▓ href="..."▓
// ▓ href='...'▓.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[cfg(not(feature = "mem"))]
#[test]
fn test_for_broken_links_and_spellcheck() {
    extern crate attohttpc;
    use std::time::Duration;

    // Set up link exceptions.  These are links that have been observed to break periodically.
    // The web.archive.org one is probably just too slow, and we should allow for that rather
    // than have it on the unreliable list.  The "period" version is because of a parsing bug.
    // Also, obviously, these links should be retested to determine if they are permanently broken
    // rather than just unreliable.

    let unreliable_links = [
        "http://www.abybank.org/abdb",
        "http://www.abybank.org/abdb/Data/NR_LH_Combined_Martin.tar.bz2",
        "http://www.bioinf.org.uk/abs/info.html",
        "http://www.bioinf.org.uk/abs/info.html#martinnum",
        "https://web.archive.org/web/20200803185732/http://www.bioinf.org.uk/abs/info.html",
        "https://web.archive.org/web/20200803185732/http://www.bioinf.org.uk/abs/info.html.",
    ];

    // Set up dictionary exceptions.  We should rewrite the code to avoid looking in certain
    // places and reduce the dictionary exceptions accordingly.

    let extra_words =
        "abybank adefghiklmnpqrstvwy amazonaws anarci barcode barcodes barcoding bcn \
        bioinf cdiff cellranger chmod clonotype clonotypes \
        clonotyping codebase colorn contig contigs cqvwdsssdhpyvf cred crispr \
        csv ctrlc cvar cvars datalayer dejavusansmono dref dyiid enclone executables false fcell \
        fixedtextbox foursie foursies frameshifted frameshifts fwr fwyh ganesh \
        genomics germline github githubusercontent google googletagmanager grok gz html \
        hypermutation hypermutations igh igk igl ighm igkc imgt \
        indel indels inkt jsdelivr json krh levenshtein linux loh lvars \
        macbook mait metadata mkdir \
        moresies multiomic nall ncbi ncross ndoublet newick nimproper \
        nopager noprint nqual nwhitef oligos onesie onesies parseable pbmc \
        pcell pdb pgas phad phylip \
        plasmablast preinstalled prepends pwm pwms redownloads \
        researchsquare samtools screenshot segn \
        sloooooooow spacebar stackexchange standalone stdout sthnqedkr subclonotype \
        subclonotypes svg thresholding timepoint tracebacks trb tsv twosie ubuntu \
        umi umis underperforming unicode untarring vdj vilella vilfwym vilm website wget wikimedia \
        wikipedia workaround workflow xf xhtml xkcd xxxxxxxxxxx xxxxxxxxxxxxxxxxxxxxxxx zenodo zx";
    let extra_words = extra_words.split(' ').collect::<Vec<&str>>();

    // Set up dictionary.

    let dictionary0 = read_to_string("../enclone-data/english_wordlist").unwrap();
    let dictionary0 = dictionary0.split('\n').collect::<Vec<&str>>();
    let mut dictionary = Vec::<String>::new();
    for w in dictionary0.iter() {
        let mut x = w.to_string();
        x.make_ascii_lowercase();
        dictionary.push(x);
    }
    for w in extra_words {
        dictionary.push(w.to_string());
    }
    unique_sort(&mut dictionary);

    // Find html pages on site.

    let mut htmls = vec!["../index.html".to_string()];
    let pages = read_dir("../pages").unwrap();
    for page in pages {
        let page = page.unwrap().path();
        let page = page.to_str().unwrap();
        if page.ends_with(".html") {
            htmls.push(format!("{}", page));
        }
    }
    let auto = read_dir("../pages/auto").unwrap();
    for page in auto {
        let page = page.unwrap().path();
        let page = page.to_str().unwrap();
        if page.ends_with(".html") {
            htmls.push(format!("{}", page));
        }
    }

    // Hardcoded exceptions to link testing, because of slowness.

    let mut tested = HashSet::<String>::new();
    tested.insert("https://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd".to_string());
    tested.insert("http://www.w3.org/1999/xhtml".to_string());

    // Hardcode exception for funny svn URL.

    tested.insert("https://github.com/10XGenomics/enclone/trunk".to_string());

    // Test each html.

    let mut dict_fail = false;
    for x in htmls {
        let mut bads = HashSet::<String>::new();
        let f = open_for_read![x];
        let depth = x.matches('/').count();
        for line in f.lines() {
            let mut s = line.unwrap();

            // Test spelling.  Case insensitive.

            let mut s0 = s.replace(',', " ");
            s0 = s0.replace('.', " ");
            s0 = s0.replace(';', " ");
            let words = s0.split(' ').collect::<Vec<&str>>();
            for i in 0..words.len() {
                let mut ok = true;
                let w = words[i].to_string();
                for c in w.chars() {
                    if !c.is_ascii_alphabetic() {
                        ok = false;
                    }
                }
                if w.is_empty() || !ok {
                    continue;
                }
                let mut wl = w.clone();
                wl.make_ascii_lowercase();
                if !bin_member(&dictionary, &wl.to_string()) && !bads.contains(&wl.to_string()) {
                    let mut wu = w.clone();
                    wu.make_ascii_uppercase();
                    if w != wu || w.len() < 20 {
                        // arbitrary long uppercase strings allowed
                        bads.insert(wl.to_string());
                        eprintln!(
                            "\nthe word \"{}\" in file {} isn't in the dictionary",
                            wl, x
                        );
                        dict_fail = true;
                    }
                }
            }

            // Check links.

            let mut links = Vec::<String>::new();
            let mut chars = Vec::<char>::new();
            for c in s.chars() {
                chars.push(c);
            }
            let mut i = 0;
            let terminators = vec![',', ' ', '\'', '"', ')', '}', '<', '#', '\n'];
            while i < chars.len() {
                let http = chars[i..].starts_with(&vec!['h', 't', 't', 'p', ':']);
                let https = chars[i..].starts_with(&vec!['h', 't', 't', 'p', 's', ':']);
                if http || https {
                    for j in i + 5..chars.len() {
                        if terminators.contains(&chars[j])
                            || (chars[j] == '.' && j < chars.len() - 1 && chars[j + 1] == ' ')
                        {
                            let mut link = String::new();
                            for k in i..j {
                                link.push(chars[k]);
                            }
                            if !tested.contains(&link.to_string()) {
                                links.push(link.clone());
                                tested.insert(link.to_string());
                            }
                            i = j - 1;
                            break;
                        }
                    }
                }
                i += 1;
            }
            let s2 = s.clone();
            while s.contains("<a href=\"") {
                let link = s.between("<a href=\"", "\"");
                if tested.contains(&link.to_string()) {
                    s = s.after("<a href=\"").to_string();
                    continue;
                }
                tested.insert(link.to_string());

                // Allow mailto to enclone.

                if link == "mailto:enclone@10xgenomics.com" {
                    s = s.after("<a href=\"").to_string();
                    continue;
                }

                // Otherwise if not http..., assume it's a file path.

                if !link.starts_with("http") {
                    let mut link = link.to_string();
                    if link.contains('#') {
                        link = link.before("#").to_string();
                    }
                    let mut z = link.clone();
                    for _ in 0..depth - 1 {
                        if !z.starts_with("../") {
                            eprintln!("something wrong with file {} on page {}", link, x);
                            std::process::exit(1);
                        }
                        z = z.after("../").to_string();
                    }
                    z = format!("../{}", z);
                    if !path_exists(&z) {
                        eprintln!("failed to find file {} on page {}", link, x);
                        std::process::exit(1);
                    }
                    s = s.after("<a href=\"").to_string();
                    continue;
                }

                // And finally do http....

                links.push(link.to_string());
                s = s.after("<a href=\"").to_string();
            }
            s = s2;
            while s.contains("<img src=\"") {
                let path = s.between("<img src=\"", "\"");
                if tested.contains(&path.to_string()) {
                    s = s.after("<img src=\"").to_string();
                    continue;
                }
                tested.insert(path.to_string());
                let path = path.to_string();
                let mut z = path.clone();
                for _ in 0..depth - 1 {
                    if !path.starts_with("../") {
                        eprintln!("something wrong with file {} on page {}", path, x);
                        std::process::exit(1);
                    }
                    z = z.after("../").to_string();
                }
                z = format!("../{}", z);
                if !path_exists(&z) {
                    eprintln!("failed to find file {} on page {}", path, x);
                    std::process::exit(1);
                }
                s = s.after("<img src=\"").to_string();
            }
            for link in links {
                // Temporary workaround.

                if link == "https://10xgenomics.github.io/enclone/install.sh" {
                    continue;
                }

                // Test for known unreliable links.

                let mut unreliable = false;
                for l in unreliable_links.iter() {
                    if *l == link {
                        unreliable = true;
                    }
                }
                if unreliable {
                    continue;
                }

                // eprintln!("checking link \"{}\"", link);

                // Approach 1 to testing if link works.  This seemed to hang once in spite of
                // the timeout.

                use attohttpc::*;
                const LINK_RETRIES: usize = 5;
                for i in 0..LINK_RETRIES {
                    if i > 0 {
                        thread::sleep(time::Duration::from_millis(100));
                        eprintln!("retrying link {}, attempt {}", link, i);
                    }
                    let req = attohttpc::get(link.clone()).read_timeout(Duration::new(10, 0));
                    let response = req.send();
                    if response.is_err() {
                        eprintln!("\ncould not read link {} on page {}\n", link, x);
                        if i == LINK_RETRIES - 1 {
                            std::process::exit(1);
                        }
                    } else {
                        let response = response.unwrap();
                        if response.is_success() {
                            break;
                        }
                        eprintln!("\ncould not read link {} on page {}\n", link, x);
                        if i == LINK_RETRIES - 1 {
                            std::process::exit(1);
                        }
                    }
                }

                // Approach 2 to testing if link works.  This may not have a timeout and does
                // not auto retry like approach 1.  Also may not compile anymore.

                /*
                use reqwest::StatusCode;
                let req = reqwest::blocking::get(link);
                if req.is_err() {
                    eprintln!("\ncould not read link {} on page {}\n", link, x);
                    std::process::exit(1);
                }
                if req.unwrap().status() == StatusCode::NOT_FOUND {
                    eprintln!("\ncould not read link {} on page {}\n", link, x);
                    std::process::exit(1);
                }
                */
            }
        }
    }
    if dict_fail {
        eprintln!("");
        std::process::exit(1);
    }
}
