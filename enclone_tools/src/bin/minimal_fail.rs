// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Given a BCR  all_contig_annotations.json file, and other arguments for enclone, such that
// the corresponding enclone command fails, iteratively attempt to make the json file smaller,
// while maintaining the fails.
//
// We require that the enclone passes, whereas enclone.old fails.
//
// Fail may be defined by (a) output containing specified text or (b) panic.
//
// This can be used to make tests.
//
// Argument 1: text (possibly in quotes) or "#PANIC".
// Argument 2: path to the all_contig_annotations.json file.
// Argument 3: working directory.
// Argument 4-: additional arguments.
//
// Some extensions you might make by temporarily editing this file:
// 1. Change the name of the enclone executable.
// 2. Add a second test using a different enclone executable.

use io_utils::*;
use pretty_trace::PrettyTrace;
use serde_json::Value;
use std::env;
use std::fs::{copy, File};
use std::io::{BufWriter, Write};
use std::process::Command;
use string_utils::*;
use vector_utils::*;

fn fails(
    fail_condition: &str,
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

    // Execute new enclone.  Require OK.

    let o = Command::new("enclone")
        .arg(&format!("BCR={}", work))
        .args(&*extra)
        .output()
        .expect("failed to execute new enclone");
    if o.status.code() != Some(0) {
        return false;
    }

    // Execute old enclone.

    let o = Command::new("enclone.old")
        .arg(&format!("BCR={}", work))
        .args(&*extra)
        .output()
        .expect("failed to execute old enclone");

    // Test for lack of status code.  It is not clear if this can happen.  It is possible that
    // it only happened on old code run with new rust, and reflects some sort of incompatibility
    // therein.

    if o.status.code().is_none() {
        return false;
    }

    // Check for fail.

    let status = o.status.code().unwrap();
    let panicked = status == 101;
    let out = strme(&o.stdout);

    let mut failed = false;
    if fail_condition == "#PANIC" {
        failed = panicked;
    } else if out.contains(&fail_condition) {
        failed = true;
    }
    if failed {
        copy(&working_json, &format!("{}.good", working_json)).unwrap();
    }
    failed
}

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let fail_condition = &args[1];
    let json = &args[2];
    let work = &args[3];
    let mut extra = Vec::<String>::new();
    for i in 4..args.len() {
        extra.push(args[i].clone());
    }
    let json = std::fs::read_to_string(&json).unwrap();
    let v: Value = serde_json::from_str(&json).unwrap();
    let entries = v.as_array().unwrap();
    let n = entries.len();
    println!("len = {}", n);
    let mut using = vec![true; n];

    let outs = format!("{}/outs", work);
    if !path_exists(&outs) {
        std::fs::create_dir(&outs).unwrap();
    }
    let working_json = format!("{}/all_contig_annotations.json", outs);

    let mut total = 0;
    loop {
        let mut progress = false;
        for k in [10000, 1000, 100, 10, 1].iter() {
            let k = *k;
            if total == 0 || total > k {
                if n >= k {
                    for i in (0..=n - k).step_by(k) {
                        let mut usingx = using.clone();
                        for j in 0..k {
                            usingx[i + j] = false;
                        }
                        if usingx != using {
                            println!("trying to delete {}", k);
                            let panicked = fails(
                                &fail_condition,
                                &entries,
                                &usingx,
                                &work,
                                &working_json,
                                &extra,
                            );
                            if panicked {
                                using = usingx;
                                total = 0;
                                for x in using.iter() {
                                    if *x {
                                        total += 1;
                                    }
                                }
                                println!("deleted {}, leaving {}", k, total);
                                progress = true;
                            }
                        }
                    }
                }
            }
        }
        if !progress {
            break;
        }
    }
    println!("\nfinal = {}\n", total);
}
