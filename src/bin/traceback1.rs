// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Check that an out-of-range reference within a rayon parallel loop yields a correct traceback.
// Make sure that line 14 stays as line 14.  Otherwise change the reference to traceback1.rs:14.

extern crate rayon;

use rayon::prelude::*;

fn main() {
    let z = vec![0; 100];
    let mut x = vec![0; 100];
    x.par_iter_mut().for_each(|r| {
        let _ = z[100 + *r]; // line 14!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    });
}

#[test]
fn test_traceback1() {
    extern crate assert_cmd;
    use assert_cmd::prelude::*;
    use std::{env, process::Command};
    env::set_var("RUST_BACKTRACE", "1");
    let mut cmd = Command::cargo_bin("traceback1").unwrap();
    let cmd = cmd
        .output()
        .expect(&format!("very strange, failed to execute test_traceback1"));
    let morsel = "traceback1.rs:14";
    if !std::str::from_utf8(&cmd.stderr).unwrap().contains(&morsel) {
        panic!("test_traceback1 failed because did not find {} as expected", morsel);
    }
}
