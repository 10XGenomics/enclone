// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use ansi_escape::*;
use enclone_core::parse_bsv;
use enclone_core::testlist::*;
use io_utils::*;
use itertools::Itertools;
use std::cmp::min;
use std::fs::read_to_string;
use std::io::Write;
use std::process::Command;
use string_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Run an enclone test.

pub fn run_test(
    enclone: &str,     // name of the enclone executable
    it: usize,         // test number
    comments: &str,    // info about the test
    test: &str,        // arguments for the test
    testname: &str,    // test category e.g. "test" or "ext_test"
    ok: &mut bool,     // true if test passes
    logx: &mut String, // logging from test
    out: &mut String,  // stdout of test
) {
    let orig_test = test.to_string();
    let mut test = test.replace("\n", "");
    for _ in 0..3 {
        test = test.replace("  ", " ");
    }
    let mut expect_null = false;
    let mut expect_fail = false;
    let mut expect_ok = false;
    let mut set_in_stone = false;
    let mut no_pre = false;
    let mut nforce = false;
    let mut ncores = false;
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
    if test.contains(" SET_IN_STONE") {
        test = test.replace(" SET_IN_STONE", "");
        set_in_stone = true;
    }
    if test.contains(" NO_PRE") {
        test = test.replace(" NO_PRE", "");
        no_pre = true;
    }
    if test.contains(" NFORCE") {
        test = test.replace(" NFORCE", "");
        nforce = true;
    }
    if test.contains(" NCORES") {
        test = test.replace(" NCORES", "");
        ncores = true;
    }
    test = test.replace("{TEST_FILES_VERSION}", &format!("{}", TEST_FILES_VERSION));
    while test.contains("  ") {
        test = test.replace("  ", " ");
    }
    let mut log = Vec::<u8>::new();
    let out_file = format!("testx/inputs/outputs/enclone_{}{}_output", testname, it + 1);
    let mut pre_arg = format!(
        "PRE=../enclone-data/big_inputs/version{}",
        TEST_FILES_VERSION
    );
    let mut local_pre_arg = format!(
        "PRE=enclone-data/big_inputs/version{},enclone_exec",
        TEST_FILES_VERSION
    );
    if no_pre {
        pre_arg = String::new();
        local_pre_arg = String::new();
    }
    if !path_exists(&out_file) && !expect_fail && !expect_ok {
        fwriteln!(log, "\nYou need to create the output file {}.\n", out_file);
        fwriteln!(
            log,
            "Do this by executing the following command from \
             the top level of the enclone repo:\n"
        );
        emit_bold_escape(&mut log);
        fwriteln!(
            log,
            "enclone {} {} > enclone_exec/testx/inputs/outputs/enclone_{}{}_output; \
             git add enclone_exec/testx/inputs/outputs/enclone_{}{}_output\n",
            local_pre_arg,
            test,
            testname,
            it + 1,
            testname,
            it + 1
        );
        fwriteln!(
            log,
            "If you just added a test to the TESTS group in testlist.rs, it would have been \
            faster if you had\nsimply typed run_last_test to get this information.  But you first \
            need to run ./build.\n",
        );
        emit_end_escape(&mut log);
        *logx = stringme(&log);
    } else {
        let mut old = String::new();
        if !expect_fail && !expect_ok {
            old = read_to_string(&out_file).unwrap();
        }
        let args = parse_bsv(&test);

        // Try to prevent a sporadic failure mode.

        let (mut p1, mut p2) = (false, false);
        for i in 0..args.len() {
            if args[i] == "GEX=123217" {
                p1 = true;
            } else if args[i] == "H5" {
                p2 = true;
            }
        }
        if p1 && !p2 {
            eprintln!(
                "\nFound GEX=123217 without H5.  Because of the PREBUILD test, this \
                can cause sporadic failures.\nHere is the test:\n{}\n",
                orig_test
            );
            std::process::exit(1);
        }

        // Form the command and execute it.

        let mut new = Command::new(&enclone);
        let mut new = new.arg(&args[0]);
        if !no_pre {
            new = new.arg(&pre_arg);
        }
        for i in 1..args.len() {
            new = new.arg(&args[i]);
        }
        if !nforce {
            new = new.arg("FORCE_EXTERNAL")
        }
        if !ncores {
            // Cap number of cores at 24.  Surprisingly, for testing on a 64-core
            // server, this significantly reduces wallclock.  And substituting either
            // 16 or 64 is slower.  Slower at the time of testing!  As we add tests or
            // change the algorithms, this may change.
            new = new.arg("MAX_CORES=24")
        }

        // Turn off process killing, which is needed to get enclone to emit the correct exit
        // code, but doesn't play nicely with the testing machinery.

        new = new.arg("NO_KILL");

        // dubious use of expect:
        let new = new
            .output()
            .unwrap_or_else(|_| panic!("failed to execute enclone for test{}", it + 1));
        let new_err = strme(&new.stderr).split('\n').collect::<Vec<&str>>();
        let new2 = stringme(&new.stdout);
        *out = new2.clone();

        // Process tests that were supposed to fail or supposed to succeed.

        if expect_fail || expect_ok {
            *ok = false;
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
                        *ok = true;
                    }
                } else if status != 0 {
                    fwriteln!(log, "\nCommand for subtest {} failed.", it + 1);
                    fwrite!(
                        log,
                        "That test was supposed to have succeeded, but instead \
                         failed, with stderr = {}",
                        new_err.iter().format("\n")
                    );
                    fwriteln!(
                        log,
                        "The command was\n\nenclone {} {}\n",
                        test,
                        local_pre_arg
                    );
                    fwriteln!(log, "stdout = {}", strme(&new.stdout));
                    fwriteln!(log, "stderr = {}", strme(&new.stderr));
                } else {
                    *ok = true;
                }
            }
            *logx = stringme(&log);

        // Process tests that yield the expected stdout.
        } else if old == new2 {
            *ok = true;
            if old.len() <= 1 && !expect_null {
                fwriteln!(
                    log,
                    "\nWarning: old output for subtest {} has {} bytes.\n",
                    it + 1,
                    old.len()
                );
            }
            if !new.stderr.is_empty() {
                fwriteln!(log, "Command for subtest {} failed.\n", it + 1);
                fwriteln!(log, "stderr has {} bytes:", new.stderr.len());
                fwrite!(log, "{}", strme(&new.stderr));
                *ok = false;
            }
            *logx = stringme(&log);

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
                    fwriteln!(
                        log,
                        "Note that this can be confusing if the differing character is in an \
                            ANSI escape sequence.\n"
                    );
                    let mut jold = i;
                    while jold > 0 && oldc[jold] != '\n' {
                        jold -= 1;
                    }
                    let mut jnew = i;
                    while jnew > 0 && newc[jnew] != '\n' {
                        jnew -= 1;
                    }
                    let mut kold = i + 1;
                    while kold < oldc.len() && oldc[kold] != '\n' {
                        kold += 1;
                    }
                    let mut knew = i + 1;
                    while knew < newc.len() && newc[knew] != '\n' {
                        knew += 1;
                    }
                    let mut old_line = String::new();
                    for z in jold + 1..kold {
                        old_line.push(oldc[z]);
                    }
                    let mut new_line = String::new();
                    for z in jnew + 1..knew {
                        new_line.push(newc[z]);
                    }
                    fwriteln!(
                        log,
                        "first differing line:\nold = {}\nnew = {}\n",
                        old_line,
                        new_line,
                    );
                    break;
                }
            }
            fwrite!(log, "old:\n{}", old);
            fwrite!(log, "new:\n{}", new2);
            if new_err.len() != 1 || !new_err[0].is_empty() {
                fwriteln!(log, "stderr has {} lines:", new_err.len());
                for i in 0..new_err.len() {
                    fwriteln!(log, "{}", new_err[i]);
                }
            }
            // let f = format!(
            //     "testx/inputs/version{}/{}/outs/all_contig_annotations.json.lz4",
            //         version, args[0].after("=") );
            // if !path_exists(&f) {
            //     println!( "Perhaps you forgot to lz4 compress the json file.\n" );
            //     std::process::exit(1);
            // }
            // println!( "The size of {} is {} bytes.", f, fs::metadata(&f).unwrap().len() );

            if !comments.is_empty() {
                fwriteln!(
                    log,
                    "enclone subtest {} failed.  Here are the comments for that test:\n\n{}\n\n\
                     If you are happy with the new output, \
                     you can replace the\noutput by executing the following command from \
                     the top level of the enclone repo (essential):\n",
                    it + 1,
                    comments
                );
            } else {
                fwriteln!(
                    log,
                    "enclone subtest {} failed.  If you are happy with the new output, \
                     you can replace the\noutput by executing the following command from \
                     the top level of the enclone repo (essential):\n",
                    it + 1
                );
            }
            if set_in_stone {
                fwriteln!(
                    log,
                    "🔴 However, the output of this test was not supposed to have changed.\n\
                     🔴 Please be extremely careful if you change it.\n",
                );
            }
            emit_bold_escape(&mut log);
            fwriteln!(
                log,
                "enclone {} {} \
                 > enclone_exec/testx/inputs/outputs/enclone_{}{}_output\n",
                local_pre_arg,
                test,
                testname,
                it + 1
            );
            emit_end_escape(&mut log);
            fwrite!(log, "and then committing the changed file.  ");
            fwriteln!(
                log,
                "You can then retest using:\n\n\
                 cargo test -p enclone enclone  -- --nocapture"
            );
            if !new2.is_empty() {
                fwriteln!(log, "");
                *logx = stringme(&log);
            } else if old != new2 {
                fwriteln!(log, "old != new");
                *logx = stringme(&log);
            }
        }
    }
}
