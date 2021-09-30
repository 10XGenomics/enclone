// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// sync_to_genome_list dir

// Find all the filenames in dir that contain a colon, and match them to entries in genome_list
// based on their accession, then move them to the location in genome_list.
//
// This provides a vehicle for managing the names of the genome files.

use pretty_trace::PrettyTrace;
use std::env;
use std::fs::{read_dir, rename};
use string_utils::TextUtils;
use vector_utils::{bin_member, sort_sync2};

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let dir = &args[1];
    let listx = include_str!("../genome_list");
    let (mut acc1, mut fns1) = (Vec::<String>::new(), Vec::<String>::new());
    for line in listx.lines() {
        if line.contains(':') && !line.starts_with('#') {
            acc1.push(line.before(":").to_string());
            fns1.push(line.to_string());
        }
    }
    let (mut acc2, mut fns2) = (Vec::<String>::new(), Vec::<String>::new());
    let all = read_dir(&dir).unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        let f = f.rev_after("/");
        if f.contains(':') {
            acc2.push(f.before(":").to_string());
            fns2.push(f.to_string());
        }
    }
    sort_sync2(&mut acc1, &mut fns1);
    sort_sync2(&mut acc2, &mut fns2);
    for i in 0..acc1.len() {
        if !bin_member(&acc2, &acc1[i]) {
            eprintln!("\n{} in genome_list but not in dir\n", acc1[i]);
            std::process::exit(1);
        }
    }
    for i in 0..acc2.len() {
        if !bin_member(&acc1, &acc2[i]) {
            eprintln!("\n{} in dir but not in genome_list\n", acc2[i]);
            std::process::exit(1);
        }
    }
    assert_eq!(acc1, acc2);
    for i in 0..acc1.len() {
        if fns1[i] != fns2[i] {
            let f2 = format!("{}/{}", dir, fns2[i]);
            let f1 = format!("{}/{}", dir, fns1[i]);
            println!("renaming {} to {}", f2, f1);
            let res = rename(&f1, &f2);
            if res.is_err() {
                eprintln!("rename failed");
                std::process::exit(1);
            }
        }
    }
}
