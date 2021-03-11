// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

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
        .expect(&format!("failed to execute enclone"));
    if new.status.code() != Some(0) {
        panic!();
    }
    let mut f = open_for_write_new!["pages/auto/help.main.html"];
    fwrite!(f, "{}", strme(&new.stdout));
    let new = Command::new("target/debug/enclone")
        .arg("HTML")
        .arg("STABLE_DOC")
        .arg("help")
        .output()
        .expect(&format!("failed to execute enclone"));
    if new.status.code() != Some(0) {
        panic!();
    }
    let mut f = open_for_write_new!["pages/auto/help.setup.html"];
    fwrite!(f, "{}", strme(&new.stdout));
    for x in HELP_PAGES.iter() {
        let new = Command::new("target/debug/enclone")
            .arg("help")
            .arg(x)
            .arg("HTML")
            .arg("STABLE_DOC")
            .output()
            .expect(&format!("failed to execute enclone"));
        if new.status.code() != Some(0) {
            panic!();
        }
        let mut f = open_for_write_new![format!("pages/auto/help.{}.html", x)];
        fwrite!(f, "{}", strme(&new.stdout));
    }
}
