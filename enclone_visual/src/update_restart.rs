// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use std::env;
use std::process::Command;
use string_utils::*;

pub fn restart_enclone() {
    let args: Vec<String> = env::args().collect();
    let mut args1 = Vec::<String>::new();
    for i in 1..args.len() {
        args1.push(args[i].clone());
    }
    let mut o = Command::new("enclone")
        .args(&args1)
        .spawn()
        .expect("failed to execute enclone restart");
    let _ = o.wait().expect("failed to wait on child");
    std::process::exit(0);
}

pub fn update_enclone() {
    println!("Automatically updating enclone, following the instructions at bit.ly/enclone.\n");
    let mut home = String::new();
    for (key, value) in env::vars() {
        if key == "HOME" {
            home = value.clone();
        }
    }
    if home.len() == 0 {
        eprintln!("Weird, unable to determine your home directory.\n");
        std::process::exit(1);
    }
    let o = Command::new("curl")
        .arg("-s")
        .arg("-L")
        .arg(
            "https://github.com/10XGenomics/enclone/\
            releases/latest/download/enclone_macos",
        )
        .arg("--output")
        .arg(&format!("{}/bin/enclone", home))
        .output()
        .expect("failed to execute curl");
    if o.status.code() != Some(0) {
        eprintln!(
            "Update failed with the following error message:\n{}",
            strme(&o.stderr)
        );
        std::process::exit(1);
    }
}
