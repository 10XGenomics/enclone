// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Check that an out-of-range reference within a rayon parallel loop yields a correct traceback.

use enclone_build::set_panic_handler;
use rayon::prelude::*;

fn main() {
    set_panic_handler(&[]);
    let z = vec![0; 100];
    let mut x = vec![0; 100];
    x.par_iter_mut().for_each(|r| {
        let _ = z[100 + *r];
    });
}

#[test]
fn test_traceback1() {
    use assert_cmd::prelude::*;
    use std::process::Command;
    let mut cmd = Command::cargo_bin("traceback1").expect(
        "\nAttempt to run traceback1 failed.  The most likely explanation for this is that\n\
        somehow you did not run \"cargo b\".  Please try that now, and be sure you are doing\n\
        it from the top-level enclone directory.\n",
    );
    let cmd = cmd
        .output()
        .unwrap_or_else(|_| panic!("{}", "very strange, failed to execute test_traceback1"));

    let err = std::str::from_utf8(&cmd.stderr).unwrap();
    let source = "panicked at enclone_build/src/bin/traceback1.rs:13";
    let count = err.matches(source).count();
    assert_eq!(
        1, count,
        "expected to find exactly one instance of traceback source \"{source}\" but found {count}"
    );
}
