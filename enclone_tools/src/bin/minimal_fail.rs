// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Given a BCR  all_contig_annotations.json file, and other arguments for enclone, such that
// the corresponding enclone command panics, iteratively attempt to make the json file smaller,
// while maintaining the panic.
//
// This can be used to make tests.
//
// Argument 1: path to the all_contig_annotations.json file.
// Argument 2: working directory.
// Argument 3-: additional arguments.

use io_utils::*;
use itertools::Itertools;
use pretty_trace::PrettyTrace;
use serde_json::Value;
use std::env;
use std::fs::{copy, File};
use std::io::{BufWriter, Write};
use std::process::Command;
use string_utils::*;
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
    // not clear how this can happen
    if o.status.code().is_none() {
        // eprintln!("weird, failed to get status code");
        return false;
    }
    if o.status.code().is_none() {
        eprint!("\nfailed to get status code, stdout =\n{}\nstderr =\n{}", strme(&o.stdout), strme(&o.stderr));
        eprintln!("\nThe command was enclone BCR={} {}.\n", work, extra.iter().format(" "));
        std::process::exit(1);
    }
    let status = o.status.code().unwrap();
    let panicked = status == 101;
    if panicked {
        // println!("\nvery good, panicked with stderr =\n{}", strme(&o.stderr));
        copy(&working_json, &format!("{}.good", working_json)).unwrap();
    /*
    } else if status != 0 {
        println!("\noh dear, status = {} and stderr =\n{}", status, strme(&o.stderr));
        std::process::exit(1);
    */
    }
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
    let mut using = vec![true; n];

    let outs = format!("{}/outs", work);
    if !path_exists(&outs) {
        std::fs::create_dir(&outs).unwrap();
    }
    let working_json = format!("{}/all_contig_annotations.json", work);

    loop {
        let mut progress = false;
        for k in [1000, 100, 10, 1].iter() {
            let k = *k;
            if n >= k {
                for i in (0..=n - k).step_by(k) {
                    let mut usingx = using.clone();
                    for j in 0..k {
                        usingx[i + j] = false;
                    }
                    if usingx != using {
                        println!("trying to delete {}", k);
                        let panicked = panics(&entries, &usingx, &work, &working_json, &extra);
                        if panicked {
                            using = usingx;
                            let mut total = 0;
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
        if !progress {
            break;
        }
    }
}
