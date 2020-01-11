// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// slightly modified version of https://vallentin.dev/2019/06/06/versioning
extern crate prost_build;

use prost_build::Config;
use std::env::consts::{ARCH, OS};
use std::process::Command;

#[cfg(debug_assertions)]
const BUILD_TYPE: &'static str = "debug";
#[cfg(not(debug_assertions))]
const BUILD_TYPE: &'static str = "release";

fn main() {
    let version_string = format!(
        "{}:{}{}, {}, {} [{}]",
        get_branch_name(),
        get_commit_hash(),
        if is_working_tree_clean() { "" } else { "+" },
        BUILD_TYPE,
        OS,
        ARCH
    );
    println!("cargo:rustc-env=VERSION_STRING={}", version_string);
    let mut config = Config::new();
    config.type_attribute(".", "#[derive(::serde::Serialize, ::serde::Deserialize)]");
    config.compile_protos(&["types.proto"], &["."]).unwrap();
}

fn get_commit_hash() -> String {
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

fn get_branch_name() -> String {
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
    let status = Command::new("git")
        .arg("diff")
        .arg("--quiet")
        .arg("--exit-code")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .status()
        .unwrap();
    status.code().unwrap() == 0
}
