// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

#![allow(unused_imports, dead_code)]

// somewhat misplaced comments:
//
// Run enclone on a test case and verify that the output is identical to what
// was gotten before.  If the output is different, look at it
// and decide if the change is justified, and if so update the output file.
//
// This test only runs if you use cargo test --release.  The test for
// not debug_assertions is a proxy for that.
//
// To test just this test, use:
//
// cargo test --release -p enclone enclone -- --nocapture

use ansi_escape::*;
use enclone::proto_io::read_proto;
use enclone::testlist::*;
use enclone::types::EncloneOutputs;
use failure::Error;
use flate2::read::GzDecoder;
use io_utils::*;
use perf_stats::*;
use pretty_trace::*;
use rayon::prelude::*;
use std::cmp::min;
use std::fs::{read_to_string, remove_file, File};
use std::io::{Read, Write};
use std::process::Command;
use std::time::Instant;
use string_utils::*;

const LOUPE_OUT_FILENAME: &str = "test/__test_proto";

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[cfg(debug_assertions)]
#[test]
fn test_enclone_fail() {
    println!("\n\"cargo test\" deliberately fails here because without running in release mode,");
    println!("the test in enclone would be too slow.  We could simply elide the test, but");
    println!("then you wouldn't know that you're missing an important test.  If you really want");
    println!("to run all tests except this test in debug(dev) mode, please use");
    println!("\"cargo test --all --exclude enclone\"");
    println!("However please also note that even with the extra test, \"cargo test --release\"");
    println!("will be faster then the above.\n");
    assert!(0 == 1);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// The following is a single test, containing many subtests, each of which is a regression test
// for a given enclone command line.
//
// If you ever need to change the output of all tests, use the main program
// update_all_main_tests.rs in enclone/src/bin.  Note that there is some duplicated code there.

#[cfg(not(debug_assertions))]
#[test]
fn test_enclone() {
    PrettyTrace::new().on();
    let t = Instant::now();
    //                       id     ok    output
    let mut results = Vec::<(usize, bool, String)>::new();
    for i in 0..TESTS.len() {
        results.push((i, false, String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let it = res.0;
        let mut test = TESTS[it].to_string();
        test = test.replace("\n", "");
        for _ in 0..3 {
            test = test.replace("  ", " ");
        }
        let mut expect_null = false;
        let mut expect_fail = false;
        let mut expect_ok = false;
        if test.contains(" EXPECT_NULL") {
            test = test.replace(" EXPECT_NULL", "");
            expect_null = true;
        }
        if test.contains(" EXPECT_FAIL") {
            test = test.replace(" EXPECT_FAIL", "");
            expect_fail = true;
        }
        if test.contains(" EXPECT_OK") {
            test = test.replace(" EXPECT_OK", "");
            expect_ok = true;
        }
        let mut log = Vec::<u8>::new();
        let out_file = format!("test/inputs/outputs/enclone_test{}_output", it + 1);
        if !path_exists(&out_file) && !expect_fail && !expect_ok {
            fwriteln!(log, "\nYou need to create the output file {}.\n", out_file);
            fwriteln!(
                log,
                "Do this by executing the following command from \
                 cellranger/lib/rust/enclone:\n"
            );
            emit_bold_escape(&mut log);
            fwriteln!(
                log,
                "enclone PRE=test/inputs/version{} {} \
                 > test/inputs/outputs/enclone_test{}_output\n",
                TEST_FILES_VERSION,
                test,
                it + 1
            );
            emit_end_escape(&mut log);
            fwriteln!(log, "and then adding/committing the new file.");
            res.2 = stringme(&log);
        } else {
            let mut old = String::new();
            if !expect_fail && !expect_ok {
                old = read_to_string(&out_file).unwrap();
            }

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
            // dubious use of expect:
            let new = new
                .arg("FORCE_EXTERNAL")
                // Cap number of cores at 24.  Surprisingly, for testing on a 64-core
                // server, this significantly reduces wallclock.  And substituting either
                // 16 or 64 is slower.  Slower at the time of testing!  As we add tests or
                // change the algorithms, this may change.
                .arg("MAX_CORES=24")
                .output()
                .expect(&format!("failed to execute enclone for test{}", it + 1));
            let new_err = strme(&new.stderr).split('\n').collect::<Vec<&str>>();
            let new2 = stringme(&new.stdout);

            // Process tests that were supposed to fail or supposed to succeed.

            if expect_fail || expect_ok {
                res.1 = false;
                if new.status.code().is_none() {
                    fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                    fwriteln!(
                        log,
                        "Something really funky happened, status code unavailable.\n"
                    );
                } else {
                    let status = new.status.code().unwrap();
                    if expect_fail {
                        if status == 0 {
                            fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                            fwriteln!(
                                log,
                                "That test was supposed to have failed, but instead \
                                 succeeded.\n"
                            );
                        } else if status != 1 {
                            fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                            fwriteln!(
                                log,
                                "That test was supposed to have failed with exit status 1,\n\
                                 but instead failed with exit status {}.\n",
                                status
                            );
                        } else {
                            res.1 = true;
                        }
                    } else {
                        if status != 0 {
                            fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                            fwriteln!(
                                log,
                                "That test was supposed to have succeeded, but instead \
                                 failed.\n"
                            );
                        } else {
                            res.1 = true;
                        }
                    }
                }
                res.2 = stringme(&log);

            // Process tests that yield the expected stdout.
            } else if old == new2 {
                res.1 = true;
                if old.len() <= 1 && !expect_null {
                    fwriteln!(
                        log,
                        "\nWarning: old output for subtest {} has {} bytes.\n",
                        it + 1,
                        old.len()
                    );
                }
                if new.stderr.len() > 0 {
                    fwriteln!(log, "Command for subtest {} failed.\n", it + 1);
                    fwriteln!(log, "stderr has {} bytes:", new.stderr.len());
                    fwrite!(log, "{}", strme(&new.stderr));
                    res.1 = false;
                }
                res.2 = stringme(&log);

            // Process tests that yield unexpected stdout.
            } else {
                fwriteln!(log, "\nSubtest {}: old and new differ", it + 1);
                fwriteln!(
                    log,
                    "old has u8 length {} and new has u8 length {}",
                    old.len(),
                    new2.len()
                );
                let mut oldc = Vec::<char>::new();
                let mut newc = Vec::<char>::new();
                for c in old.chars() {
                    oldc.push(c);
                }
                for c in new2.chars() {
                    newc.push(c);
                }
                fwriteln!(
                    log,
                    "old has char length {} and new has char length {}",
                    oldc.len(),
                    newc.len()
                );
                for i in 0..min(oldc.len(), newc.len()) {
                    if oldc[i] != newc[i] {
                        fwriteln!(
                            log,
                            "the first difference is at character {}: old = \"{}\", \
                             new = \"{}\"\n",
                            i,
                            oldc[i],
                            newc[i]
                        );
                        break;
                    }
                }
                fwrite!(log, "old:\n{}", old);
                fwrite!(log, "new:\n{}", new2);
                if new_err.len() != 1 || new_err[0].len() != 0 {
                    fwriteln!(log, "stderr has {} lines:", new_err.len());
                    for i in 0..new_err.len() {
                        fwriteln!(log, "{}", new_err[i]);
                    }
                }
                // let f = format!(
                //     "test/inputs/version{}/{}/outs/all_contig_annotations.json.lz4",
                //         version, args[0].after("=") );
                // if !path_exists(&f) {
                //     println!( "Perhaps you forgot to lz4 compress the json file.\n" );
                //     std::process::exit(1);
                // }
                // println!( "The size of {} is {} bytes.", f, fs::metadata(&f).unwrap().len() );

                fwriteln!(
                    log,
                    "enclone subtest {} failed.  If you are happy with the new output, \
                     you can replace the\noutput by executing the following command from \
                     cellranger/lib/rust/enclone (essential!):\n",
                    it + 1
                );
                emit_bold_escape(&mut log);
                fwriteln!(
                    log,
                    "enclone PRE=test/inputs/version{} {} \
                     > test/inputs/outputs/enclone_test{}_output\n",
                    TEST_FILES_VERSION,
                    test,
                    it + 1
                );
                emit_end_escape(&mut log);
                fwrite!(log, "and then committing the changed file.  ");
                fwriteln!(
                    log,
                    "You can then retest using:\n\n\
                     cargo test --release -p enclone enclone  -- --nocapture"
                );
                if new2.len() > 0 {
                    fwriteln!(log, "");
                    res.2 = stringme(&log);
                } else if old != new2 {
                    fwriteln!(log, "old != new");
                    res.2 = stringme(&log);
                }
            }
        }
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

// Test that examples are what we claim they are.

#[cfg(not(debug_assertions))]
#[test]
fn test_enclone_examples() {
    PrettyTrace::new().on();
    for t in 0..EXAMPLES.len() {
        let testn = format!("{}", EXAMPLES[t]);
        let out_file = format!("src/example{}", t + 1);
        let old = read_to_string(&out_file).unwrap();
        let args = testn.split(' ').collect::<Vec<&str>>();
        let mut new = Command::new("target/release/enclone");
        let mut new = new.arg(format!("PRE=test/inputs/version{}", TEST_FILES_VERSION));
        for i in 0..args.len() {
            new = new.arg(&args[i]);
        }
        let new = new
            .arg("FORCE_EXTERNAL")
            .output()
            .expect(&format!("failed to execute test_enclone_examples"));
        let new2 = stringme(&new.stdout);
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

// Test that references to the dataset version in README.md are current.

#[cfg(not(debug_assertions))]
#[test]
fn test_version_number_in_readme() {
    PrettyTrace::new().on();
    let readme = read_to_string("README.md").unwrap();
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

// Test that enclone help all HTML works (without STABLE_DOC).

#[cfg(not(debug_assertions))]
#[test]
fn test_help_no_stable() {
    PrettyTrace::new().on();
    let mut new = Command::new("target/release/enclone");
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

// Test that help output hasn't changed.

#[cfg(not(debug_assertions))]
#[test]
fn test_help_output() {
    PrettyTrace::new().on();
    let out_file = format!("src/help.all.html");
    let old = read_to_string(&out_file).unwrap();
    let mut new = Command::new("target/release/enclone");
    let mut new = new.arg("help");
    new = new.arg("all");
    new = new.arg("HTML");
    new = new.arg("STABLE_DOC");
    new = new.arg("NOPAGER");
    let new = new
        .arg("FORCE_EXTERNAL")
        .output()
        .expect(&format!("failed to execute test_help_output"));
    if new.status.code() != Some(0) {
        eprintln!("Attempt to run enclone help all failed.\n");
        std::process::exit(1);
    }
    let new2 = stringme(&new.stdout);
    if old != new2 {
        eprintln!(
            "\nYou need to update help output by typing \
                \"enclone help all HTML STABLE_DOC > src/help.all.html\", \
                assuming that the change is expected.\n"
        );
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test that PREBUILD works.  

#[cfg(not(debug_assertions))]
#[test]
fn test_enclone_prebuild() {
    PrettyTrace::new().on();
    let t = Instant::now();
    let mb = format!(
        "test/inputs/version{}/123749/outs/raw_feature_bc_matrix/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    );
    if path_exists(&mb) {
        remove_file(&mb).unwrap();
    }

    // First pass: run with NH5.

    let test_id = 48;
    let it = test_id - 1;
    let testn = format!("{} NH5", TESTS[it]);
    let out_file = format!("test/inputs/outputs/enclone_test{}_output", test_id);
    let old = read_to_string(&out_file).unwrap();
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new("target/release/enclone");
    let mut new = new.arg(format!("PRE=test/inputs/version{}", TEST_FILES_VERSION));
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
             You may want to add more info to this failure message.\n"
        );
        eprintln!("old output =\n{}\n", old);
        eprintln!("new output =\n{}\n", new2);
        std::process::exit(1);
    }
    if !path_exists(&format!(
        "test/inputs/version{}/123749/outs/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    )) {
        panic!("\nenclone_test_prebuild: did not create feature_barcode_matrix.bin.");
    }

    // Second pass: run without PREBUILD but using the feature_barcode_matrix.bin that the first
    // pass created.

    let testn = TESTS[it];
    let args = testn.split(' ').collect::<Vec<&str>>();
    let mut new = Command::new("target/release/enclone");
    let mut new = new.arg(format!("PRE=test/inputs/version{}", TEST_FILES_VERSION));
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
        "test/inputs/version{}/123749/outs/feature_barcode_matrix.bin",
        TEST_FILES_VERSION
    ))
    .unwrap();
    println!("\nused {:.2} seconds in enclone_test_prebuild", elapsed(&t));
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This test runs enclone for two test inputs, with LOUPE output
// turned on. It will then read both the bincode and proto file created
// and asserts that we get the same data structure either way.
//
// It also tests to make sure that the LOUPE output is unchanged.  If it changed for a good
// reason, update the output file.  Otherwise perhaps something has gone wrong!

#[cfg(not(debug_assertions))]
#[test]
fn test_proto_write() -> Result<(), Error> {
    let tests = vec!["BCR=123085", "TCR=101287"];
    let pre_arg = format!("PRE=test/inputs/version{}", TEST_FILES_VERSION);
    let binary_arg = format!("BINARY={}.bin", LOUPE_OUT_FILENAME);
    let proto_arg = format!("PROTO={}.proto", LOUPE_OUT_FILENAME);
    for t in tests.iter() {
        // FIXME: It would be nicer to use the enclone API here
        std::process::Command::new("target/release/enclone")
            .args(&[&pre_arg, *t, &binary_arg, &proto_arg])
            .output()
            .expect(&format!("failed to execute enclone for test_proto_write"));

        // Test to make sure output is unchanged.

        let oldx = format!("test/inputs/{}.binary.gz", t.after("="));
        let newx = format!("{}.bin", LOUPE_OUT_FILENAME);
        let mut f = File::open(&oldx)?;
        let mut oldbufgz = Vec::<u8>::new();
        f.read_to_end(&mut oldbufgz)?;
        let mut gz = GzDecoder::new(&oldbufgz[..]);
        let mut oldbuf = Vec::<u8>::new();
        gz.read_to_end(&mut oldbuf)?;
        let mut f = File::open(&newx)?;
        let mut newbuf = Vec::<u8>::new();
        f.read_to_end(&mut newbuf)?;
        if oldbuf != newbuf {
            eprintln!(
                "\nThe binary output of enclone on {} has changed.  If this is expected,\n\
                please regenerate the file {}\n",
                t, oldx
            );
            std::process::exit(1);
        }

        // Test to make sure proto and bin are consistent.

        let outputs_proto = read_proto(format!("{}.proto", LOUPE_OUT_FILENAME))?;
        let outputs_bin: EncloneOutputs = io_utils::read_obj(format!("{}.bin", LOUPE_OUT_FILENAME));
        std::fs::remove_file(format!("{}.proto", LOUPE_OUT_FILENAME))?;
        std::fs::remove_file(format!("{}.bin", LOUPE_OUT_FILENAME))?;
        assert!(outputs_proto == outputs_bin);
    }

    Ok(())
}
