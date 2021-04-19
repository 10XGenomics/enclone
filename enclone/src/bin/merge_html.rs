// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Build html files by inserting other html files.
//
// If supplied the single argument BUILD, also rebuild from source.

use enclone::html::*;
use enclone_core::testlist::SITE_EXAMPLES;
use enclone_core::*;
use io_utils::*;
use itertools::Itertools;
use pretty_trace::*;
use rayon::prelude::*;
use std::env;
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();

    // Build from source.

    let args: Vec<String> = env::args().collect();
    if args.len() == 2 && args[1] == "BUILD" {
        let mut results = Vec::<(usize, String)>::new();
        for i in 0..SITE_EXAMPLES.len() {
            results.push((i, String::new()));
        }
        results.par_iter_mut().for_each(|r| {
            let i = r.0;
            let test = SITE_EXAMPLES[i].1;
            let args = parse_bsv(&test);
            let new = Command::new("target/debug/enclone")
                .args(&args)
                .arg("MAX_CORES=24")
                .output()
                .expect(&format!("failed to execute enclone"));
            if new.status.code() != Some(0) {
                eprint!(
                    "\nenclone_site_examples: example {} failed to execute, stderr =\n{}",
                    i + 1,
                    strme(&new.stderr),
                );
                std::process::exit(1);
            }
            r.1 = stringme(&new.stdout);
        });
        for i in 0..SITE_EXAMPLES.len() {
            let example_name = SITE_EXAMPLES[i].0;
            let out_file = format!("{}", example_name);
            let mut f = open_for_write_new![&out_file];
            fwrite!(&mut f, "{}", results[i].1);
        }
    }
    let mut site_ex = Vec::<String>::new();
    for i in 0..SITE_EXAMPLES.len() {
        let example_name = SITE_EXAMPLES[i].0;
        let out_file = format!("{}", example_name);
        site_ex.push(out_file);
    }
    unique_sort(&mut site_ex);

    // Apply insert_html.

    let all = read_dir("pages").unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        if f.ends_with(".html.src") {
            let mut level = 2;
            let mut target = format!("pages/auto/{}.html", f.between("/", "."));
            if target == "pages/auto/index.html".to_string() {
                target = "index.html".to_string();
                level = 0;
            }
            insert_html(&f, &target, false, level);
        }
    }

    // This is ugly.  Edit the html pages.

    let all = read_dir("pages/auto").unwrap();
    let mut allx = vec!["index.html".to_string()];
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap().to_string();
        if f.ends_with(".html") && !bin_member(&site_ex, &f) {
            allx.push(f);
        }
    }
    for f in allx {
        let mut edited = false;
        let mut lines = Vec::<String>::new();
        {
            let h = open_for_read![&format!("{}", f)];
            for line in h.lines() {
                let s = line.unwrap();
                lines.push(s.clone());
                if s.contains("googletag") {
                    edited = true;
                }
            }
        }
        if !edited {
            let x = format!("{}\n", lines.iter().format("\n"));
            let mut h = open_for_write_new![&format!("{}", f)];
            fwrite!(h, "{}", edit_html(&x));
        }
    }
}
