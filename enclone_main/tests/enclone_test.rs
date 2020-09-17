// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

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

// 1. Test for redundancy of columns in parseable output.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_for_parseable_redundancy() {
    let test = "BCR=85333 CDR3=CARDLRVEGFDYW POUT=testx/outputs/redundancy_out";
    let args = parse_bsv(&test);
    let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
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

// NOT BASIC

// Run the download command on the landing page and make sure it works.
//
// Runs with "small", and passes second argument so we can put outputs in a defined place.
//
// There are two passes.  The first pass tests the copy of install.sh that is one master, and
// the second pass tests the local version.
//
// Not sure if this needs to be internal-only.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_curl_command() {
    let mut internal_run = false;
    for (key, value) in env::vars() {
        if (key == "HOST" || key == "HOSTNAME") && value.ends_with(".fuzzplex.com") {
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
        for pass in 1..=2 {
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
                command = "curl -sSf -L bit.ly/enclone_install | sh -s small testx/outputs";
                version = "master";
            } else {
                command = "cat ../install.sh | sh -s small testx/outputs";
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
            let req = [
                "bin/enclone",
                "enclone/datasets/123085/outs/all_contig_annotations.json.lz4",
                "enclone/datasets_small_checksum",
                "enclone/version",
            ];
            for f in req.iter() {
                if !path_exists(&format!("testx/outputs/{}", f)) {
                    eprintln!(
                        "\nAttempt to run enclone install command using {} version of \
                        install.sh failed to fetch {}.\n",
                        version, f
                    );
                    std::process::exit(1);
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
// calibrated for bespin1, and requires linux
// cargo test --test enclone_test --features cpu -- --nocapture
// from enclone_main directory
// or just ./speed from root directory

#[cfg(not(feature = "basic"))]
#[cfg(feature = "cpu")]
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

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
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
    const ACCEPTABLE_OTHER_PACKAGES: [&str; 5] =
        ["arrayref", "cloudabi", "fuchsia-cprng", "ring", "webpki"];
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
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let test = TESTS[it].to_string();
        let mut out = String::new();
        run_test(
            env!("CARGO_BIN_EXE_enclone"),
            it,
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

// Regression tests for internal features.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
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

// 13.

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
#[test]
fn test_for_broken_links_and_spellcheck() {
    extern crate attohttpc;
    use std::time::Duration;

    // Set up dictionary exceptions.  We should rewrite the code to avoid looking in certain
    // places and reduce the dictionary exceptions accordingly.

    let extra_words = "abybank amazonaws anarci barcode barcodes barcoding bcn \
        bioinf cdiff cellranger chmod clonotype clonotypes \
        clonotyping codebase colorn contig contigs cqvwdsssdhpyvf cred crispr \
        csv ctrlc cvar cvars datalayer dejavusansmono dref dyiid enclone executables false fcell \
        fixedtextbox foursie foursies frameshifts fwr \
        genomics germline github githubusercontent google googletagmanager grok gz html \
        hypermutation hypermutations igh igk igl ighm igkc imgt \
        indel indels inkt jsdelivr json levenshtein linux loh lvars macbook mait metadata mkdir \
        moresies multiomic nall ncbi ncross newick nimproper \
        nopager noprint nqual nwhitef oligos onesie onesies parseable pbmc pcell pdb phylip \
        plasmablast preinstalled prepends pwm pwms redownloads samtools screenshot segn \
        sloooooooow spacebar stackexchange standalone stdout subclonotype \
        subclonotypes svg thresholding timepoint tracebacks trb twosie ubuntu \
        umi umis underperforming unicode untarring vdj website wget wikimedia \
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 14.

// NOT BASIC

// Test site examples to make sure they are what they claim to be, and that the
// merged html files are correct.

#[cfg(not(feature = "basic"))]
#[cfg(not(feature = "cpu"))]
#[test]
fn test_site_examples() {
    for i in 0..SITE_EXAMPLES.len() {
        let example_name = SITE_EXAMPLES[i].0;
        let test = SITE_EXAMPLES[i].1;
        let in_file = format!("../{}", example_name);
        let in_stuff = read_to_string(&in_file).expect(&format!("couldn't find {}", in_file));
        let args = parse_bsv(&test);
        let new = Command::new(env!("CARGO_BIN_EXE_enclone"))
            .args(&args)
            .output()
            .expect(&format!("failed to execute test_site_examples"));
        if new.status.code() != Some(0) {
            eprint!(
                "\nenclone_site_examples: example {} failed to execute, stderr =\n{}",
                i + 1,
                strme(&new.stderr),
            );
            std::process::exit(1);
        }
        let out_stuff = strme(&new.stdout);
        if in_stuff != out_stuff {
            eprintln!("\nThe output for site example {} has changed.\n", i + 1);
            eprintln!("stderr:\n{}", strme(&new.stderr));
            let old_lines = in_stuff.split('\n').collect::<Vec<&str>>();
            let new_lines = out_stuff.split('\n').collect::<Vec<&str>>();
            eprintln!(
                "old stdout has {} lines; new stdout has {} lines",
                old_lines.len(),
                new_lines.len(),
            );
            for i in 0..min(old_lines.len(), new_lines.len()) {
                if old_lines[i] != new_lines[i] {
                    eprintln!(
                        "first different stdout line is line {}\nold = {}\nnew = {}",
                        i + 1,
                        old_lines[i],
                        new_lines[i],
                    );
                    break;
                }
            }
            let save = format!("testx/outputs/{}", example_name.rev_after("/"));
            {
                let mut f = open_for_write_new![&save];
                fwrite!(f, "{}", out_stuff);
            }
            let mut in_filex = in_file.clone();
            if in_filex.starts_with("../") {
                in_filex = in_filex.after("../").to_string();
            }
            eprintln!("\nPlease diff {} enclone_main/{}.", in_filex, save);
            eprintln!(
                "\nPossibly this could be because you're running \"cargo t\" in an \
                environment without the\n\
                extended dataset collection.  Possibly you should run \
                \"cd enclone; cargo test basic -- --nocapture\" instead.\n\n\
                Otherwise, if you're satisfied with the new output, you can update using\n\n\
                enclone {} > {}.\n",
                args.iter().format(" "),
                example_name
            );
            std::process::exit(1);
        }
    }

    insert_html(
        "../pages/index.html.src",
        "testx/outputs/index.html",
        true,
        0,
    );
    insert_html(
        "../pages/expanded.html.src",
        "testx/outputs/expanded.html",
        true,
        2,
    );
    let new_index = read_to_string("testx/outputs/index.html").unwrap();
    if read_to_string("../index.html").unwrap() != new_index {
        eprintln!("\nContent of index.html has changed.");
        {
            let mut f = open_for_write_new!["testx/outputs/index.html.new"];
            fwrite!(f, "{}", new_index);
        }
        eprintln!("Please diff index.html enclone_main/testx/outputs/index.html.new.\n");
        std::process::exit(1);
    }
    /*
    if read_to_string("../pages/auto/expanded.html").unwrap()
        != edit_html(&read_to_string("testx/outputs/expanded.html").unwrap())
    */
    if read_to_string("../pages/auto/expanded.html").unwrap()
        != read_to_string("testx/outputs/expanded.html").unwrap()
    {
        eprintln!("\nContent of expanded.html has changed.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 15. Test that examples are what we claim they are.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_enclone_examples() {
    PrettyTrace::new().on();
    for t in 0..EXAMPLES.len() {
        let testn = format!("{}", EXAMPLES[t]);
        let out_file = format!("../enclone_help/src/example{}", t + 1);
        let old = read_to_string(&out_file).unwrap();
        let args = testn.split(' ').collect::<Vec<&str>>();
        let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
        let mut new = new.arg(format!(
            "PRE=../enclone-data/big_inputs/version{}",
            TEST_FILES_VERSION
        ));
        for i in 0..args.len() {
            new = new.arg(&args[i]);
        }
        let new = new
            .arg("FORCE_EXTERNAL")
            .output()
            .expect(&format!("failed to execute test_enclone_examples"));
        let new2 = stringme(&new.stdout);
        if new.status.code() != Some(0) {
            eprint!(
                "\nenclone_test_examples: example{} failed to execute, stderr =\n{}",
                t + 1,
                strme(&new.stderr),
            );
            std::process::exit(1);
        }
        if old != new2 {
            eprintln!(
                "\nenclone_test_examples: the file example{} is not up to date\n",
                t + 1
            );
            eprintln!("old output =\n{}", old);
            eprintln!("new output =\n{}\n", new2);
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 16. Test that references to the dataset version in README.md are current.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_version_number_in_readme() {
    PrettyTrace::new().on();
    let readme = read_to_string("../README.md").unwrap();
    let fields = readme.split('/').collect::<Vec<&str>>();
    for x in fields {
        if x.starts_with("version") {
            let y = x.after("version");
            if y.parse::<usize>().is_ok() {
                let v = y.force_usize();
                assert_eq!(v, TEST_FILES_VERSION as usize);
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 17. Test that the DejaVuSansMono definition in enclone_css_v2.css has not changed.  We put this
// here because that definition has to be manually tested, and we don't want it accidentally
// changed and broken.  This is really gross, but it's not clear how to do it better.
//
// Absolutely hideous implementation to verify that
// cat ../pages/enclone_css_v2.css | head -36 = "2474276863 1467".
//
// Only works with high probability.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_dejavu() {
    PrettyTrace::new().on();
    let mut cat_output_child = Command::new("cat")
        .arg("../pages/enclone_css_v2.css")
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    if let Some(cat_output) = cat_output_child.stdout.take() {
        let mut head_output_child = Command::new("head")
            .arg("-36")
            .stdin(cat_output)
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        cat_output_child.wait().unwrap();
        if let Some(head_output) = head_output_child.stdout.take() {
            let cksum_output_child = Command::new("cksum")
                .stdin(head_output)
                .stdout(Stdio::piped())
                .spawn()
                .unwrap();
            let cksum_stdout = cksum_output_child.wait_with_output().unwrap();
            head_output_child.wait().unwrap();
            let cksum = String::from_utf8(cksum_stdout.stdout).unwrap();
            // println!("cksum = {}", cksum);
            assert!(cksum == "2474276863 1467\n".to_string());
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 18. Test that help output hasn't changed.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_help_output() {
    PrettyTrace::new().on();
    let pages = vec![
        "setup",
        "main",
        "quick",
        "how",
        "command",
        "glossary",
        "example1",
        "example2",
        "input",
        "input_tech",
        "parseable",
        "filter",
        "special",
        "lvars",
        "cvars",
        "amino",
        "display",
        "indels",
        "color",
        "faq",
        "developer",
        "all",
    ];
    for p in pages {
        let mut command = format!("enclone help {}", p);
        if p == "setup" {
            command = "enclone help".to_string();
        } else if p == "main" {
            command = "enclone".to_string();
        }
        let out_file = format!("../pages/auto/help.{}.html", p);
        let old = read_to_string(&out_file).unwrap();
        let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
        let mut new = new.arg("HTML");
        if p == "setup" {
            new = new.arg("help");
        } else if p == "main" {
        } else {
            new = new.arg("help");
            new = new.arg(p);
        }
        new = new.arg("STABLE_DOC");
        new = new.arg("NOPAGER");
        let new = new
            .arg("FORCE_EXTERNAL")
            .output()
            .expect(&format!("failed to execute test_help_output"));
        if new.status.code() != Some(0) {
            eprintln!("Attempt to run {} failed.\n", command);
            std::process::exit(1);
        }
        let new2 = edit_html(&stringme(&new.stdout));
        if old != new2 {
            eprintme!(old.len(), new2.len());
            eprintln!(
                "\nHelp test failed on {}.\n\
                 You need to update help output by typing \"./build\", \
                    assuming that the change is expected.\n",
                p
            );
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 19. Test that enclone help all HTML works (without STABLE_DOC).

#[cfg(not(feature = "cpu"))]
#[test]
fn test_help_no_stable() {
    PrettyTrace::new().on();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg("help");
    new = new.arg("all");
    new = new.arg("HTML");
    new = new.arg("NOPAGER");
    let new = new
        .arg("FORCE_EXTERNAL")
        .output()
        .expect(&format!("failed to execute test_help_output"));
    if new.status.code() != Some(0) {
        eprintln!("Attempt to run enclone help all without STABLE_DOC failed.\n");
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 20. Test that PREBUILD works.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_enclone_prebuild() {
    PrettyTrace::new().on();
    let t = Instant::now();
    let mb = format!(
        "../enclone-data/big_inputs/version{}/123749/outs/raw_feature_bc_matrix/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    );
    if path_exists(&mb) {
        remove_file(&mb).unwrap();
    }

    // First pass: run with NH5.

    let test_id = 48;
    let it = test_id - 1;
    let testn = format!("{} NH5", TESTS[it]);
    let out_file = format!("testx/inputs/outputs/enclone_test{}_output", test_id);
    let old = read_to_string(&out_file).unwrap();
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg(format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    ));
    for i in 0..args.len() {
        new = new.arg(&args[i]);
    }
    // dubious use of expect:
    let new = new
        .arg("FORCE_EXTERNAL")
        .output()
        .expect(&format!("failed to execute test_enclone_prebuild"));
    // let new_err = strme(&new.stderr).split('\n').collect::<Vec<&str>>();
    let new2 = stringme(&new.stdout);
    if old != new2 {
        eprintln!(
            "\nenclone_test_prebuild: first pass output has changed.\n\
            If you are OK with the new output, it should work to type:\n\
            enclone {} > enclone_main/{}\n",
            testn, out_file
        );
        eprintln!("old output =\n{}\n", old);
        eprintln!("new output =\n{}\n", new2);
        std::process::exit(1);
    }
    if !path_exists(&format!(
        "../enclone-data/big_inputs/version{}/123749/outs/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    )) {
        panic!("\nenclone_test_prebuild: did not create feature_barcode_matrix.bin.");
    }

    // Second pass: run without PREBUILD but using the feature_barcode_matrix.bin that the first
    // pass created.

    let testn = TESTS[it];
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new(env!("CARGO_BIN_EXE_enclone"));
    let mut new = new.arg(format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    ));
    for i in 0..args.len() {
        new = new.arg(&args[i]);
    }
    // dubious use of expect:
    let new = new
        .arg("FORCE_EXTERNAL")
        .output()
        .expect(&format!("failed to execute enclone_test_prebuild"));
    // let new_err = strme(&new.stderr).split('\n').collect::<Vec<&str>>();
    let new2 = stringme(&new.stdout);
    if old != new2 {
        eprintln!(
            "\nenclone_test_prebuild: second pass output has changed.\n\
             You may want to add more info to this failure message.\n\
             And don't forget to remove feature_barcode_matrix.bin.\n"
        );
        eprintln!("new output =\n{}\n", new2);
        std::process::exit(1);
    }

    // Clean up: delete feature_barcode_matrix.bin.

    std::fs::remove_file(&format!(
        "../enclone-data/big_inputs/version{}/123749/outs/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    ))
    .unwrap();
    println!("\nused {:.2} seconds in enclone_test_prebuild", elapsed(&t));
}

fn check_enclone_outs_consistency(enclone_outs: &EncloneOutputs) {
    let uref_items = &enclone_outs.universal_reference.items;
    for cl in &enclone_outs.clonotypes {
        for cl_chain in &cl.chains {
            assert!(uref_items.len() > cl_chain.v_idx as usize);
            assert!(uref_items.len() > cl_chain.j_idx as usize);
            assert_eq!(uref_items[cl_chain.v_idx as usize].region, 1);
            assert_eq!(uref_items[cl_chain.j_idx as usize].region, 3);
            if let Some(d_idx) = cl_chain.d_idx {
                assert!(uref_items.len() > d_idx as usize);
                assert_eq!(uref_items[d_idx as usize].region, 2);
            }
            if let Some(u_idx) = cl_chain.u_idx {
                assert!(uref_items.len() > u_idx as usize);
                assert_eq!(uref_items[u_idx as usize].region, 0);
            }
            for ex_cl in &cl.exact_clonotypes {
                for ex_cl_chain in ex_cl.chains.iter().map(|info| &info.chain) {
                    if let Some(c_region_idx) = ex_cl_chain.c_region_idx {
                        assert!(uref_items.len() > c_region_idx as usize);
                        assert_eq!(uref_items[c_region_idx as usize].region, 4);
                    }
                }
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 21. This test runs enclone for two test inputs, with LOUPE output
// turned on. It will then read both the bincode and proto file created
// and asserts that we get the same data structure either way.
//
// It also tests to make sure that the LOUPE output is unchanged.  If it changed for a good
// reason, update the output file.  Otherwise perhaps something has gone wrong!

#[cfg(not(feature = "cpu"))]
#[test]
fn test_proto_write() -> Result<(), Error> {
    let tests = vec!["BCR=123085", "TCR=101287"];
    let pre_arg = format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    );

    let bin_file = format!("{}.binary", LOUPE_OUT_FILENAME);
    let proto_file = format!("{}.proto", LOUPE_OUT_FILENAME);

    let binary_arg = format!("BINARY={}", bin_file);
    let proto_arg = format!("PROTO={}", proto_file);
    for t in tests.iter() {
        // FIXME: It would be nicer to use the enclone API here
        std::process::Command::new(env!("CARGO_BIN_EXE_enclone"))
            .args(&[&pre_arg, *t, &binary_arg, &proto_arg])
            .output()
            .expect(&format!("failed to execute enclone for test_proto_write"));

        // Test to make sure proto and bin are consistent.

        let outputs_proto = read_proto(&proto_file)?;
        let outputs_bin: EncloneOutputs = io_utils::read_obj(&bin_file);
        if outputs_proto != outputs_bin {
            eprintln!("\noutputs_proto is not equal to outputs_bin\n");
            std::process::exit(1);
        }

        // Test to make sure that the clonotype iterator works

        let clonotypes: Vec<_> = ClonotypeIter::from_file(&proto_file).unwrap().collect();
        assert!(clonotypes == outputs_proto.clonotypes);

        // Check consistency

        check_enclone_outs_consistency(&outputs_proto);

        // Test to make sure output is unchanged.

        let oldx = format!("testx/inputs/{}.binary.sha256", t.after("="));
        let mut fold = std::fs::File::open(&oldx)?;
        let mut cksum_old = Vec::<u8>::new();
        fold.read_to_end(&mut cksum_old)?;
        let mut fnew = std::fs::File::open(&bin_file)?;
        let mut cksum_new = Vec::<u8>::new();
        fnew.read_to_end(&mut cksum_new)?;
        let cksum_new = format!("{:x}", sha2::Sha256::digest(&cksum_new));
        std::fs::remove_file(&proto_file)?;
        std::fs::remove_file(&bin_file)?;
        if strme(&cksum_old) != cksum_new {
            eprintln!(
                "\nThe binary output of enclone on {} has changed.  If this is expected,\n\
                please run the command\n\
                echo -n {} > enclone_main/testx/inputs/{}.binary.sha256",
                t,
                &cksum_new,
                t.after("=")
            );
            std::process::exit(1);
        }
    }
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 22. Test for no change to the output of the command that appears in the enclone annotated
// example on the landing page.

#[cfg(not(feature = "cpu"))]
#[test]
fn test_annotated_example() {
    PrettyTrace::new().on();
    let it = 0;
    let test = "BCR=123085 CDR3=CTRDRDLRGATDAFDIW";
    let mut out = String::new();
    let mut log = String::new();
    let mut ok = false;
    run_test(
        env!("CARGO_BIN_EXE_enclone"),
        it,
        &test,
        "annotated_example_test",
        &mut ok,
        &mut log,
        &mut out,
    );
    print!("{}", log);
    if !ok {
        let mut log = Vec::<u8>::new();
        emit_red_escape(&mut log);
        fwriteln!(
            log,
            "Oh no: the results for the annotated example on the landing page have \
            changed.  Assuming that\nthe change is intentional, to fix this you \
            need to do two things:\n\
            1. Update the test output.\n\
            2. Manually update the annotated example output.\n\
            Because the second item is such a big pain, we stopped doing it, but you should\n\
            check to make sure that it has not changed too much from what is shown."
        );
        emit_end_escape(&mut log);
        eprintln!("{}", strme(&log));
        std::process::exit(1);
    }
}

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
