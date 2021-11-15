// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Sync all the crate versions in the workspace to the versions defined in the file master.toml
// in the top-level of the workspace.

use io_utils::{fwrite, open_for_read, open_for_write_new, path_exists};
use itertools::Itertools;
use pretty_trace::PrettyTrace;
use std::collections::HashMap;
use std::fs::read_dir;
use std::io::{BufRead, Write};
use string_utils::TextUtils;

fn main() {
    PrettyTrace::new().on();
    let mut version = HashMap::<String, String>::new();
    let f = open_for_read!["master.toml"];
    for line in f.lines() {
        let s = line.unwrap();
        if !s.starts_with('#') && s.contains('=') {
            version.insert(s.before(" = ").to_string(), s.after(" = ").to_string());
        }
    }
    let all = read_dir(".").unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        let toml = format!("{}/Cargo.toml", f);
        if path_exists(&toml) {
            let mut newg = Vec::<String>::new();
            {
                let g = open_for_read![&toml];
                for line in g.lines() {
                    let s = line.unwrap();
                    let mut t = s.clone();
                    if s.contains(" =") && !s.after(" =").starts_with(" [") {
                        let cratex = s.before(" =").to_string();
                        if version.contains_key(&cratex) {
                            t = format!("{} = {}", cratex, version[&cratex]);
                        }
                    }
                    newg.push(t);
                }
            }
            let mut g = open_for_write_new![&toml];
            fwrite!(g, "{}\n", newg.iter().format("\n"));
        }
    }
}
