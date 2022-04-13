// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Check executable versions.  For now, only checks git.

use pretty_trace::PrettyTrace;
use std::process::Command;
use string_utils::strme;
use string_utils::TextUtils;

fn main() {
    PrettyTrace::new().on();

    let new = Command::new("git")
        .arg("--version")
        .output()
        .expect("git version test failed");
    if new.status.code() != Some(0) {
        eprintln!("\ngit version test didn't execute\n");
        std::process::exit(1);
    }
    let s = strme(&new.stdout);
    let mut v = s.after("git version ").to_string();
    if v.contains(" ") {
        v = v.before(" ").to_string();
    }
    if v.contains("\n") {
        v = v.before("\n").to_string();
    }
    let x = v.split('.').collect::<Vec<&str>>();
    if x.len() < 3
        || !x[0].parse::<usize>().is_ok()
        || !x[1].parse::<usize>().is_ok()
        || !x[2].parse::<usize>().is_ok()
    {
        eprintln!("\nCould not parse git version {v}.\n");
        std::process::exit(1);
    }
    if x[0].force_usize() < 2 || (x[0].force_usize() == 2 && x[1].force_usize() < 24) {
        eprintln!("\ngit version must be at least 2.24.z\n");
        std::process::exit(1);
    }
}
