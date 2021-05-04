// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Takes a single argument: a comma-separated list of ids, allowing hyphenated ranges.  Copy these
// to the internal collections.  If it has already been copied, it is moved aside, and at the end,
// deleted (assuming that the copy did not fail.
//
// You need to set quota first on the second internal location.
//
// This is probably not fully consistent with current pipeline structure.
//
// Run make_enclone_testlist_all after this to update the catalog.
//
// For use at 10x Genomics.

use enclone_core::copy_for_enclone::*;
use enclone_core::defs::*;
use enclone_core::testlist::*;
use io_utils::*;
use pretty_trace::*;
use std::collections::HashMap;
use std::env;
use std::fs::{remove_dir_all, rename};
use std::process::Command;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let mut config = HashMap::<String, String>::new();
    let _ = get_config(&mut config);
    let mut dests = Vec::<String>::new();
    dests.push(format!("{}/current{}", config["earth"], TEST_FILES_VERSION));
    dests.push(format!("{}/current{}", config["cloud"], TEST_FILES_VERSION));
    let args: Vec<String> = env::args().collect();
    let ids0 = args[1].split(',').collect::<Vec<&str>>();
    let mut ids = Vec::<String>::new();
    for id in ids0.iter() {
        if id.contains('-') {
            let start = id.before("-").force_usize();
            let stop = id.after("-").force_usize();
            for p in start..=stop {
                ids.push(format!("{}", p));
            }
        } else {
            ids.push(id.to_string());
        }
    }
    let ones = &config["ones"];
    for id in ids.iter() {
        // Get path.

        let mut p = id.clone();
        assert!(p.parse::<usize>().is_ok());
        let url = format!("{}/{}", ones, p);
        let o = Command::new("curl")
            .arg(&url)
            .output()
            .expect("failed to execute http");
        let m = String::from_utf8(o.stdout).unwrap();
        if m.contains("502 Bad Gateway") {
            eprintln!(
                "\nWell this is sad.  The URL \
                {} yielded a 502 Bad Geteway \
                message.  Either try again later or ask someone for help.\n\n",
                url,
            );
            std::process::exit(1);
        }
        if m.contains("\"path\":\"") {
            let path = m.between("\"path\":\"", "\"").to_string();
            p = format!("{}/outs", path);
            if !path_exists(&p) {
                eprintln!(
                    "\nIt looks like you've provided an id {} for \
                    which\nthe pipeline outs folder has not yet been \
                    generated.\n\n",
                    p
                );
                std::process::exit(1);
            }
        } else {
            eprintln!(
                "\nIt looks like you've provided either an incorrect \
                id {} or else one for which\n\
                the pipeline outs folder has not yet been generated.\n\n",
                p
            );
            std::process::exit(1);
        }

        // Move directories if they exist.

        let mut moved = false;
        for dest in dests.iter() {
            if path_exists(&format!("{}/{}", dest, id)) {
                if path_exists(&format!("{}/{}.aside", dest, id)) {
                    eprintln!("\nPlease remove {}/{}.aside.\n", dest, id);
                    std::process::exit(1);
                }
                rename(
                    &format!("{}/{}", dest, id),
                    &format!("{}/{}.aside", dest, id),
                )
                .unwrap();
                moved = true;
            }
        }

        // Start copy.

        println!("copying {}", id);
        for i in (0..dests.len()).rev() {
            let dest = &dests[i];
            let target = format!("{}/{}", dest, id);
            if path_exists(&target) {
                eprintln!("\nPlease delete {}.\n", target);
                std::process::exit(1);
            }
            copy_for_enclone(&format!("{}/..", p), &target);
        }

        // Remove moved directories.

        if moved {
            for dest in dests.iter() {
                remove_dir_all(&format!("{}/{}.aside", dest, id)).unwrap();
            }
        }
    }
}
