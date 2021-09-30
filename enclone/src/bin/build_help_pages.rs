// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Build help pages, called by ./build.

use enclone_core::defs::*;
use io_utils::*;
use pretty_trace::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::Command;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let new = Command::new("target/debug/enclone")
        .arg("HTML")
        .arg("STABLE_DOC")
        .output()
        .unwrap_or_else(|_| panic!("{}", "failed to execute enclone".to_string()));
    if new.status.code() != Some(0) {
        panic!();
    }
    let mut f = open_for_write_new!["pages/auto/help.main.html"];
    fwrite!(f, "{}", strme(&new.stdout));
    let new = Command::new("target/debug/enclone")
        .arg("help")
        .arg("HTML")
        .arg("STABLE_DOC")
        .output()
        .unwrap_or_else(|_| panic!("{}", "failed to execute enclone".to_string()));
    if new.status.code() != Some(0) {
        eprintln!("\nbuild_help_pages failed");
        eprintln!("stdout = {}", strme(&new.stdout));
        eprintln!("stderr = {}\n", strme(&new.stderr));
        std::process::exit(1);
    }
    let mut f = open_for_write_new!["pages/auto/help.setup.html"];
    fwrite!(f, "{}", strme(&new.stdout));
    for x in HELP_PAGES.iter() {
        let new;
        if *x == "setup" {
            new = Command::new("target/debug/enclone")
                .arg("help")
                .arg("HTML")
                .arg("STABLE_DOC")
                .output()
                .unwrap_or_else(|_| panic!("{}", "failed to execute enclone".to_string()));
        } else {
            new = Command::new("target/debug/enclone")
                .arg("help")
                .arg(x)
                .arg("HTML")
                .arg("STABLE_DOC")
                .output()
                .unwrap_or_else(|_| panic!("{}", "failed to execute enclone".to_string()));
        }
        if new.status.code() != Some(0) {
            eprintln!("\nbuild_help_pages failed on {}", x);
            eprintln!("stdout:\n{}", strme(&new.stdout));
            eprintln!("stderr:\n{}", strme(&new.stderr));
            std::process::exit(1);
        }
        let mut f = open_for_write_new![format!("pages/auto/help.{}.html", x)];
        fwrite!(f, "{}", strme(&new.stdout));
    }
}
