// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// The purpose of this file is to make some version information available so that it can be
// printed out at appropriate points by enclone.  This files is a slightly modified version
// of https://vallentin.dev/2019/06/06/versioning.

extern crate chrono;
extern crate string_utils;

use chrono::prelude::*;
use std::env::consts::{ARCH, OS};
use std::process::Command;
use string_utils::*;

#[cfg(debug_assertions)]
const BUILD_TYPE: &'static str = "debug";
#[cfg(not(debug_assertions))]
const BUILD_TYPE: &'static str = "release";

fn main() {
    let version_string = format!(
        "{} : {}{} : {} : {} : {} : {}",
        get_branch_name(),
        get_commit_hash(),
        if is_working_tree_clean() { "" } else { "+" },
        get_commit_date(),
        BUILD_TYPE,
        OS,
        ARCH
    );
    println!("cargo:rustc-env=VERSION_STRING={}", version_string);
}

fn get_commit_hash() -> String {
    match std::env::var("GITHUB_SHA") {
        Ok(v) => return v[0..7].to_string(),
        _ => (),
    }

    let output = Command::new("git")
        .arg("log")
        .arg("-1")
        .arg("--pretty=format:%h") // Abbreviated commit hash
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .unwrap();
    assert!(output.status.success());
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn is_github() -> bool {
    match std::env::var("GITHUB_SHA") {
        Ok(_) => true,
        _ => false,
    }
}

// We used to have the commit date here but this is easier and serves the same purpose for
// the version string.

fn get_commit_date() -> String {
    Local::now().to_string().before(" ").to_string()
}

fn get_branch_name() -> String {
    if is_github() {
        match std::env::var("GITHUB_REF") {
            Ok(v) => return v,
            _ => return "master".to_string(),
        }
    }

    let output = Command::new("git")
        .arg("rev-parse")
        .arg("--abbrev-ref")
        .arg("HEAD")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .unwrap();
    assert!(output.status.success());
    String::from_utf8_lossy(&output.stdout)
        .trim_end()
        .to_string()
}

fn is_working_tree_clean() -> bool {
    match std::env::var("GITHUB_SHA") {
        Ok(_) => return true,
        _ => (),
    }

    let status = Command::new("git")
        .arg("diff")
        .arg("--quiet")
        .arg("--exit-code")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .status()
        .unwrap();
    status.code().unwrap() == 0
}
