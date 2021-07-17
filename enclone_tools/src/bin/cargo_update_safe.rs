// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// For each dependent crate, cargo update it and run cargo b.  If that succeeds,
// commit the change.  Otherwise, undo.
//
// Note that this could in principle be updated.
//
// Takes a long time.

use io_utils::*;
use pretty_trace::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::Command;
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    let f = open_for_read!["Cargo.lock"];
    let mut crates = Vec::<String>::new();
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("name = \"") {
            let cratex = s.between("name = \"", "\"");
            crates.push(cratex.to_string());
        }
    }
    unique_sort(&mut crates);
    for cratex in crates.iter() {
        let cratex = &*cratex;
        println!("updating {}", cratex);
        let new = Command::new("cargo")
            .arg("update")
            .arg("-p")
            .arg(&cratex)
            .output()
            .expect(&format!("failed to execute cargo update"));
        if new.status.code() != Some(0) {
            println!("update failed");
            continue;
        }

        let new = Command::new("git")
            .arg("status")
            .output()
            .expect(&format!("failed to execute git status"));
        if new.status.code() != Some(0) {
            println!("git status failed, something is wrong");
            std::process::exit(1);
        }
        if !strme(&new.stdout).contains("Cargo.lock") {
            println!("nothing changed");
            continue;
        }
        println!("compiling");
        let new = Command::new("cargo")
            .arg("b")
            .output()
            .expect(&format!("failed to execute cargo b"));
        if new.status.code() != Some(0) {
            println!("compilation failed, resetting");
            let new = Command::new("git")
                .arg("checkout")
                .arg("Cargo.lock")
                .output()
                .expect(&format!("failed to execute git checkout"));
            if new.status.code() != Some(0) {
                println!("git checkout failed, something is wrong");
                std::process::exit(1);
            }
            continue;
        }
        let new = Command::new("git")
            .arg("commit")
            .arg("-a")
            .arg("-m")
            .arg(&format!("update cratex {}", cratex))
            .output()
            .expect(&format!("failed to execute git commit"));
        if new.status.code() != Some(0) {
            println!("git commit failed, something is wrong");
            std::process::exit(1);
        }
    }
}
