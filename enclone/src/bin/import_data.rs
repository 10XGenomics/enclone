// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Takes a single argument: a comma-separated list of ids, allowing hyphenated ranges.  Copy these
// to the internal collections.  If it has already been copied, you'll need to delete it first.
//
// For use at 10x Genomics.

use enclone_core::copy_for_enclone::*;
use enclone_core::testlist::*;
use io_utils::*;
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::Command;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let mut dests = vec![format!("/mnt/assembly/vdj/current{}", TEST_FILES_VERSION)];
    let mut domain = String::new();
    for (key, value) in env::vars() {
        if key == "HOST" || key == "HOSTNAME" {
            if value.ends_with(".com") && value.rev_before(".com").contains(".") {
                let d = value.rev_before(".com").rev_after(".");
                domain = format!("{}.com", d);
            }
        }
    }
    let cloud_loc = "/mnt/assembly/vdj/cloud";
    if !path_exists(&cloud_loc) {
        eprintln!("\nCan't find the cloud (1).\n");
        std::process::exit(1);
    }
    let f = open_for_read![&cloud_loc];
    for line in f.lines() {
        let cloud = line.unwrap();
        let cloud_path = format!("{}/current{}", cloud, TEST_FILES_VERSION);
        if !path_exists(&cloud_path) {
            eprintln!("\nCan't find the cloud (2).\n");
            std::process::exit(1);
        }
        dests.push(cloud_path);
    }

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
    for id in ids.iter() {
        // Get path.

        let mut p = id.clone();
        assert!(p.parse::<usize>().is_ok());
        let url = format!("https://xena.{}/api/analyses/{}", domain, p);
        let o = Command::new("curl")
            .arg(url)
            .output()
            .expect("failed to execute xena http");
        let m = String::from_utf8(o.stdout).unwrap();
        if m.contains("502 Bad Gateway") {
            eprintln!(
                "\nWell this is sad.  The URL \
                http://xena.{}/api/analyses/{} yielded a 502 Bad Geteway \
                message.  Either try again later or ask someone for help.\n\n",
                domain, p,
            );
            std::process::exit(1);
        }
        if m.contains("\"path\":\"") {
            let path = m.between("\"path\":\"", "\"").to_string();
            p = format!("{}/outs", path);
            if !path_exists(&p) {
                eprintln!(
                    "\nIt looks like you've provided a xena id {} for \
                    which\nthe pipeline outs folder has not yet been \
                    generated.\n\n",
                    p
                );
                std::process::exit(1);
            }
        } else {
            eprintln!(
                "\nIt looks like you've provided either an incorrect \
                xena id {} or else one for which\n\
                the pipeline outs folder has not yet been generated.\n\n",
                p
            );
            std::process::exit(1);
        }

        // Start copy.

        println!("copying {}", id);
        for dest in dests.iter() {
            let target = format!("{}/{}", dest, id);
            if path_exists(&target) {
                eprintln!("\nPlease delete {}.\n", target);
                std::process::exit(1);
            }
            copy_for_enclone(&format!("{}/..", p), &target);
        }
    }
}
