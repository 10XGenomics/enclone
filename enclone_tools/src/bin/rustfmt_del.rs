// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Run rustfmt on all source files.

use io_utils::*;
use pretty_trace::PrettyTrace;
use rayon::prelude::*;
use std::process::Command;

fn main() {
    PrettyTrace::new().on();
    let mut paths = Vec::<String>::new();
    let all = dir_list(".");
    for d in all.iter() {
        let subs = ["src", "src/bin", "tests"];
        for sub in subs.iter() {
            let s = format!("{}/{}", d, sub);
            if path_exists(&s) {
                let all = dir_list(&s);
                for f in all.iter() {
                    if f.ends_with(".rs") {
                        let path = format!("{}/{}", s, f);
                        paths.push(path);
                    }
                }
            }
        }
    }
    let mut groups = Vec::<Vec<String>>::new();
    let batch = 10;
    for i in (0..paths.len()).step_by(batch) {
        let j = std::cmp::min(i + batch, paths.len());
        let mut batch = Vec::<String>::new();
        for k in i..j {
            batch.push(paths[k].clone());
        }
        groups.push(batch);
    }
    groups.par_iter_mut().for_each(|group| {
        let new = Command::new("rustfmt")
            .arg("--edition")
            .arg("2018")
            .args(&*group)
            .output()
            .unwrap_or_else(|_| panic!("{}", "failed to execute rustfmt"));
        if new.status.code() != Some(0) {
            eprintln!("\nrustfmt failed");
            std::process::exit(1);
        }
    });
}
