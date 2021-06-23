// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use std::env;
use std::process::Command;

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
