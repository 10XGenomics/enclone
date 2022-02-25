// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Build help pages, called by ./build.

use enclone_core::defs::HELP_PAGES;
use io_utils::{fwrite, open_for_write_new};
use pretty_trace::PrettyTrace;

use std::io::Write;
use std::process::Command;
use string_utils::strme;

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
        let mut f = open_for_write_new![format!("pages/auto/help.{}.html", x).as_str()];
        fwrite!(f, "{}", strme(&new.stdout));
    }
}
