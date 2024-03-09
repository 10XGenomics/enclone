// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::parse_bsv;
use enclone_testlist::TEST_FILES_VERSION;
use io_utils::{fwrite, fwriteln, path_exists};
use rayon::prelude::*;
use std::cmp::min;
use std::fmt::Write;
use std::fs::read_to_string;
use std::process::Command;
use string_utils::strme;

/// A test case is declared with a fixed test index, a comment, and the code to run.
pub type TestCase = (i32, &str, &str);

/// Run a collection of test cases in parallel.
/// TODO: replace this by declaring test cases individually.
pub fn run_tests(enclone: &str, test_name: &str, max_cores: usize, tests: &[TestCase]) {
    println!("running tests using {enclone}");
    let failures: Vec<_> = tests
        .par_iter()
        .map(|(number, comments, args)| {
            run_test(enclone, *number, comments, args, test_name, max_cores)
        })
        .filter_map(Result::err)
        .collect();
    for f in &failures {
        print!("{}", f.log);
    }
    assert!(failures.is_empty());
}

#[derive(Default)]
pub struct TestResult {
    /// Logging.
    pub log: String,
    /// Contents of stdout.
    pub stdout: String,
}

/// Run a single test case.
pub fn run_test(
    enclone: &str,    // name of the enclone executable
    testname: &str,   // test category e.g. "test" or "ext_test"
    max_cores: usize, // max cores, or if 0, default value determined here
    num: usize,       // test number
    comments: &str,   // info about the test
    test: &str,       // arguments for the test
) -> Result<TestResult, TestResult> {
    let mut test = test.replace('\n', "");
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
    let out_file = format!("testx/inputs/outputs/enclone_{}{}_output", testname, num);
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
        return Err(TestResult {
            log: format!(
                "\nYou need to create the output file {out_file}.\n\
                 Do this by executing the following command from \
                 the top level of the enclone repo:\n\n\
                 enclone {local_pre_arg} {test} > enclone_exec/testx/inputs/outputs/enclone_{testname}{test_index_from_1}_output; \
                 git add enclone_exec/testx/inputs/outputs/enclone_{testname}{test_index_from_1}_output\n\n\
                 If you just added a test to the TESTS group in testlist.rs, it would have been \
                 faster if you had\nsimply typed run_last_test to get this information.  But you first \
                 need to run ./build.",
            ),
            ..Default::default()
        });
    }
    let mut old = String::new();
    if !expect_fail && !expect_ok {
        old = read_to_string(&out_file).unwrap();
    }
    let args = parse_bsv(&test);

    // Form the command and execute it.

    let mut new = Command::new(enclone);
    let mut new = new.arg(args[0]);
    if !no_pre {
        new = new.arg(&pre_arg);
    }
    for arg in args.iter().skip(1) {
        new = new.arg(arg);
    }
    if !nforce {
        new = new.arg("FORCE_EXTERNAL")
    }
    if !ncores && max_cores == 0 {
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
        .unwrap_or_else(|_| panic!("failed to execute enclone for test {}", num));
    let mut res = TestResult {
        stdout: String::from_utf8_lossy(&new.stdout).to_string(),
        ..Default::default()
    };

    let Some(status) = new.status.code() else {
        res.log = format!("Failed to get status code for command for subtest {num}.");
        return Err(res);
    };

    // Process tests that were supposed to fail or supposed to succeed.
    if expect_fail {
        return match status {
            1 => Ok(res),
            0 => {
                res.log = format!(
                    "Command for subtest {test_index_from_1} failed.\n\
            The test command was expected to fail but it succeeded unexpectedly."
                );
                Err(res)
            }
            _ => {
                res.log = format!(
                    "Command for subtest {test_index_from_1} failed.\n
            The test was expected to fail with with exit status 1,\n\
            but instead failed with exit status {status}."
                );
                Err(res)
            }
        };
    }
    if expect_ok {
        return if status == 0 {
            Ok(res)
        } else {
            res.log = format!(
                "Command for subtest {test_index_from_1} failed with exit status {status}.\n\
                The command was\n\nenclone {test} {local_pre_arg}\n\n
                stdout = {}\n\n\
                stderr = {}",
                strme(&new.stdout),
                strme(&new.stderr)
            );
            Err(res)
        };
    }

    // Process tests that yield the expected stdout.
    if old == res.stdout {
        if old.len() <= 1 && !expect_null {
            fwriteln!(
                &mut res.log,
                "Warning: old output for subtest {test_index_from_1} has {} bytes.\n",
                old.len()
            );
        }
        if !new.stderr.is_empty() {
            fwriteln!(
                &mut res.log,
                "Command for subtest {test_index_from_1} failed.\n"
            );
            fwriteln!(&mut res.log, "stderr has {} bytes:", new.stderr.len());
            fwrite!(&mut res.log, "{}", strme(&new.stderr));
            return Err(res);
        }
        return Ok(res);
    }

    // Process tests that yield unexpected stdout.
    fwriteln!(&mut res.log, "\nSubtest {}: old and new differ", num);
    fwriteln!(
        &mut res.log,
        "old has u8 length {} and new has u8 length {}",
        old.len(),
        res.stdout.len()
    );
    let oldc: Vec<_> = old.chars().collect();
    let newc: Vec<_> = res.stdout.chars().collect();
    fwriteln!(
        &mut res.log,
        "old has char length {} and new has char length {}",
        oldc.len(),
        newc.len()
    );
    for i in 0..min(oldc.len(), newc.len()) {
        if oldc[i] != newc[i] {
            fwriteln!(
                &mut res.log,
                "the first difference is at character {}: old = \"{}\", \
                         new = \"{}\"\n",
                i,
                oldc[i],
                newc[i]
            );
            fwriteln!(
                &mut res.log,
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
                &mut res.log,
                "first differing line:\nold = {}\nnew = {}\n",
                old_line,
                new_line,
            );
            break;
        }
    }
    fwrite!(&mut res.log, "old:\n{}", old);
    fwrite!(&mut res.log, "new:\n{}", res.stdout);
    let new_err = String::from_utf8_lossy(&new.stderr).to_string();
    if !new_err.is_empty() {
        let lines = new_err.chars().filter(|c| *c == '\n').count() + 1;
        fwriteln!(&mut res.log, "stderr has {lines} lines:\n{new_err}");
    }

    if !comments.is_empty() {
        fwriteln!(
            &mut res.log,
            "enclone subtest {} failed.  Here are the comments for that test:\n\n{}\n\n\
                     If you are happy with the new output, \
                     you can replace the\noutput by executing the following command from \
                     the top level of the enclone repo (essential):\n",
            num,
            comments
        );
    } else {
        fwriteln!(
            &mut res.log,
            "enclone subtest {} failed.  If you are happy with the new output, \
                     you can replace the\noutput by executing the following command from \
                     the top level of the enclone repo (essential):\n",
            num
        );
    }
    if set_in_stone {
        fwriteln!(
            &mut res.log,
            "🔴 However, the output of this test was not supposed to have changed.\n\
                     🔴 Please be extremely careful if you change it.\n",
        );
    }
    fwriteln!(
        &mut res.log,
        "enclone {} {} \
                 > enclone_exec/testx/inputs/outputs/enclone_{}{}_output\n",
        local_pre_arg,
        test,
        testname,
        num
    );
    fwrite!(&mut res.log, "and then committing the changed file.  ");
    fwriteln!(
        &mut res.log,
        "You can then retest using:\n\n\
                 cargo test -p enclone enclone  -- --nocapture"
    );
    Err(res)
}
