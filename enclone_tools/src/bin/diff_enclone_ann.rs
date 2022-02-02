// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// This program is designed for the assessment of changes to the annotation algorithm
// (in crate vdj_ann) or the reference sequence.
//
// This only makes sense for BCR.
//
// Run enclone twice on the same dataset, where the two runs differ either by
// • an extra command-line argument (specified by OLD_PARAM=...)
// • or a code version (specified by OLD_EXEC=..., as comparsed to whatever PATH gives).
// Then diff the outputs.
//
// usage diff_enclone_ann <argument list>
// where exactly one of OLD_PARAM=... or OLD_EXEC... must appear in the list.
//
// This measures changes using the dref variable.  Higher dref is assumed to be worse.  Use of dref
// is sort of a poor man's proxy for alignment score (e.g. deducible from ALIGN1 or ALIGN2), but
// has the advantage of being more reflective of enclone's view of the data.  For example, ALIGN1
// will correctly reflect two indels, but dref will reflect enclone's messed up view of them.
//
// Typically one would run this with just the dataset specification and OLD_EXEC.

use pretty_trace::PrettyTrace;
use std::env;
use std::process::Command;
use string_utils::{strme, TextUtils};

fn main() {
    PrettyTrace::new().on();
    let mut args: Vec<String> = env::args().collect();
    let (mut old_param, mut old_exec) = (None, None);
    let mut args2 = Vec::<String>::new();
    let blacklist = [
        "PCELL",
        "POUT",
        "PCOLS",
        "BARCODE",
        "NOPAGER",
        "NOPRINT",
        "BUILT_IN",
        "CHAINS_EXACT",
    ];
    for i in 1..args.len() {
        let mut a = args[i].clone();
        if a.contains('=') {
            a = a.before("=").to_string();
        }
        for j in 0..blacklist.len() {
            if a == blacklist[j] {
                eprintln!("\nPlease do not use {} as an argument.\n", a);
                std::process::exit(1);
            }
        }
        if args[i].starts_with("OLD_PARAM=") {
            old_param = Some(args[i].after("OLD_PARAM=").to_string());
        } else if args[i].starts_with("OLD_EXEC=") {
            old_exec = Some(args[i].after("OLD_EXEC=").to_string());
        } else {
            args2.push(args[i].clone());
        }
    }
    args = args2;
    if !(old_param.is_some() ^ old_exec.is_some()) {
        eprintln!("\nExactly one of OLD_PARAM=... or OLD_EXEC=... must be specified.\n");
        std::process::exit(1);
    }
    let pcols = "datasets,barcode,dref";

    // Run enclone twice, once for "old" and once for "new".

    let mut lines = Vec::<String>::new();
    for pass in 1..=2 {
        let mut argsp = args.clone();
        if pass == 1 && old_param.is_some() {
            argsp.push(old_param.clone().unwrap());
        }
        let mut new;
        if pass == 2 || old_param.is_some() {
            new = Command::new("enclone");
        } else {
            new = Command::new(old_exec.clone().unwrap());
        }
        let new = new
            .args(&argsp)
            .arg("BUILT_IN")
            .arg("CHAINS_EXACT=2")
            .arg("PCELL")
            .arg(&format!("PCOLS={}", pcols))
            .arg("POUT=stdout")
            .arg("NOPRINT")
            .output()
            .unwrap_or_else(|_| panic!("{}", "failed to execute enclone".to_string()));
        if new.status.code() != Some(0) {
            eprint!(
                "\nenclone call returned nonzero status code, stderr =\n{}",
                strme(&new.stderr),
            );
            std::process::exit(1);
        }
        let out = strme(&new.stdout);
        for (i, x) in out.lines().enumerate() {
            if i > 0 {
                let fields = x.split(',').collect::<Vec<&str>>();
                let dataset = &fields[0];
                let barcode = &fields[1];
                let dref = &fields[2];
                let state = if pass == 1 { "old" } else { "new" };
                lines.push(format!("{}:{},{},{}", dataset, barcode, dref, state));
            }
        }
    }

    // Identify the differences.

    lines.sort();
    let mut i = 0;
    let (mut old_dref, mut new_dref) = (0, 0);
    let (mut lost_old, mut lost_new) = (0, 0);
    let mut changed_cells = 0;
    while i < lines.len() {
        if i < lines.len() - 1 && lines[i].before(",") == lines[i + 1].before(",") {
            if lines[i].rev_before(",") != lines[i + 1].rev_before(",") {
                let mut l1 = lines[i].clone();
                let mut l2 = lines[i + 1].clone();
                let state1 = l1.rev_after(",");
                if state1 == "new" {
                    std::mem::swap(&mut l1, &mut l2);
                }
                let dref1 = l1.between(",", ",").force_usize();
                let dref2 = l2.between(",", ",").force_usize();
                if dref1 > dref2 {
                    // println!("\n{}\n{}  BETTER", l1, l2);
                } else {
                    if dref2 - dref1 >= 5 {
                        println!("\n{}\n{}  WORSE", l1, l2);
                    }
                }
                old_dref += dref1;
                new_dref += dref2;
                changed_cells += 1;
            }
            i += 2;
        } else {
            let state = lines[i].rev_after(",");
            if state == "new" {
                lost_old += 1;
                // println!("\n{}  BETTER", lines[i]);
            } else {
                lost_new += 1;
                println!("\n{}  WORSE", lines[i]);
            }
            i += 1;
        }
    }
    println!("\nchanged cells = {}", changed_cells);
    println!("\nold dref sum = {}", old_dref);
    println!("new dref sum = {}", new_dref);
    println!("improvement = {}", old_dref - new_dref);
    println!("\ncells lost by old = {}", lost_old);
    println!("cells lost by new = {}\n", lost_new);
}
