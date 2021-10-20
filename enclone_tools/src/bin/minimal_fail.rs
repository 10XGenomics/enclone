// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Given a BCR  all_contig_annotations.json file, and other arguments for enclone, such that 
// the corresponding enclone command crashes, iteratively attempt to make the json file smaller,
// while maintaining the crash.
//
// Argument 1: path to the all_contig_annotations.json file.
// Argument 2: working directory.
// Argument 3-: additional arguments.

use io_utils::*;
use pretty_trace::PrettyTrace;
use serde_json::Value;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::Command;
use vector_utils::*;

fn panics(
    entries: &Vec<Value>, 
    using: &Vec<bool>, 
    work: &str, 
    working_json: &str, 
    extra: &Vec<String>,
) -> bool {
    let n = entries.len();
    let mut entriesx = entries.clone();
    let mut to_delete = vec![false; n];
    for i in 0..n {
        if !using[i] {
            to_delete[i] = true;
        }
    }
    erase_if(&mut entriesx, &to_delete);
    {
        let mut f = open_for_write_new![&working_json];
        fwriteln!(f, "[");
        for i in 0..entriesx.len() {
            fwrite!(f, "{}", entriesx[i]);
            if i < entriesx.len() - 1 {
                fwrite!(f, ",");
            }
            fwriteln!(f, "");
        }
        fwriteln!(f, "]");
    }
    let o = Command::new("enclone")
        .arg(&format!("BCR={}", work))
        .args(&*extra)
        .output()
        .expect("failed to execute enclone");
    let status = o.status.code().unwrap();
    let panicked = status == 101;
    panicked
}

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let json = &args[1];
    let work = &args[2];
    let mut extra = Vec::<String>::new();
    for i in 3..args.len() {
        extra.push(args[i].clone());
    }
    let json = std::fs::read_to_string(&json).unwrap();
    let v: Value = serde_json::from_str(&json).unwrap();
    let entries = v.as_array().unwrap();
    let n = entries.len();
    println!("len = {}", n);
    let using = vec![true; n];

    let outs = format!("{}/outs", work);
    if !path_exists(&outs) {
        std::fs::create_dir(&outs).unwrap();
    }
    let working_json = format!("{}/all_contig_annotations.json", work);

    // First try deleting blocks of 1000.

    loop {
        let panicked = panics(&entries, &using, &work, &working_json, &extra);
        println!("{}", panicked);
        break;
    }
}
