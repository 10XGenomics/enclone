// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Run enclone twice on the same dataset, where the two runs differ either by
// • an extra command-line argument (specified by OLD_PARAM=...)
// • or a code version (specified by OLD_EXEC=..., as comparsed to whatever PATH gives).
// Then diff the outputs, piping an organized version of visual output to less.
//
// usage diff_clonotypes <argument list>
// where exactly one of OLD_PARAM=... or OLD_EXEC... must appear in the list.
//
// As implemented, this only detects differences where a clonotype changes either its
// barcodes or its number of chains.
//
// example:
// diff_enclone BCR=123085 OLD_PARAM=NDOUBLET
//
// Do not use any of the grouping arguments.

use enclone::misc1::*;
use equiv::EquivRel;
use itertools::Itertools;
use pretty_trace::*;
use std::collections::HashMap;
use std::env;
use std::process::Command;
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    setup_pager(true);
    let mut args: Vec<String> = env::args().collect();
    let (mut old_param, mut old_exec) = (None, None);
    let mut args2 = Vec::<String>::new();
    let blacklist = ["PCELL", "POUT", "PCOLS", "BARCODE", "NOPAGER", "NOPRINT"];
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
    let pcols = "barcode,group_id,exact_subclonotype_id,nchains";

    // Run enclone twice, once for "old" and once for "new".  Return the parseable output lines
    // as outs.

    let mut outs = Vec::<Vec<String>>::new();
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
        let mut lines = Vec::<String>::new();
        for x in out.lines() {
            lines.push(x.to_string());
        }
        outs.push(lines);
    }
    let outs_orig = outs.clone();

    // Map barcodes to groups and then remove the groups.  And sort.

    let mut to_group = vec![HashMap::<String, String>::new(); 2];
    let mut to_groupx = HashMap::<String, (String, usize)>::new();
    for v in 0..2 {
        for i in 0..outs[v].len() {
            let barcode = outs[v][i].before(",").to_string();
            let group = outs[v][i].after(",").before(",").to_string();
            to_group[v].insert(barcode.clone(), group.clone());
            to_groupx.insert(barcode, (group, v));
            outs[v][i] = format!(
                "{},{}",
                outs[v][i].before(","),
                outs[v][i].after(",").after(",")
            );
        }
        outs[v].sort();
    }

    // Find barcodes that have different entries in the two outs.

    let mut bc = Vec::<String>::new();
    let mut group = vec![Vec::<(String, String)>::new(); 2];
    for v1 in 0..2 {
        let v2 = (v1 + 1) % 2;
        let (lines1, lines2) = (&outs[v1], &outs[v2]);
        for j in 0..lines1.len() {
            if !bin_member(lines2, &lines1[j]) {
                let b = lines1[j].before(",").to_string();
                bc.push(b.clone());
                group[v1].push((to_group[v1][&b].clone(), b));
            }
        }
        group[v1].sort();
    }
    unique_sort(&mut bc);

    // Group the barcodes.  If two barcodes are in the same clonotype for either run, then
    // they get grouped together.

    let mut e = EquivRel::new(bc.len() as i32);
    for v in 0..2 {
        let mut i = 0;
        while i < group[v].len() {
            let j = next_diff1_2(&group[v], i as i32) as usize;
            for k in i + 1..j {
                e.join(
                    bin_position(&bc, &group[v][i].1) as i32,
                    bin_position(&bc, &group[v][k].1) as i32,
                );
            }
            i = j;
        }
    }

    // Define the groups.

    let mut groups = Vec::<Vec<String>>::new();
    let mut reps = Vec::<i32>::new();
    e.orbit_reps(&mut reps);
    for i in 0..reps.len() {
        let mut o = Vec::<i32>::new();
        e.orbit(reps[i], &mut o);
        let mut b = Vec::<String>::new();
        for j in 0..o.len() {
            b.push(bc[o[j] as usize].clone());
        }
        groups.push(b);
    }

    // Expand the groups to include the full clonotypes.

    for j in 0..groups.len() {
        let mut groupsx = Vec::<(String, usize)>::new();
        for i in 0..groups[j].len() {
            groupsx.push(to_groupx[&groups[j][i]].clone());
        }
        unique_sort(&mut groupsx);
        for v in 0..2 {
            for i in 0..outs_orig[v].len() {
                let barcode = outs_orig[v][i].before(",").to_string();
                let group = outs_orig[v][i].after(",").before(",").to_string();
                if bin_member(&groupsx, &(group, v)) {
                    groups[j].push(barcode.clone());
                }
            }
        }
        unique_sort(&mut groups[j]);
    }

    // For each group of barcodes, run enclone using the group, for both old and new.

    println!();
    for i in 0..groups.len() {
        let bc_arg = format!("{}", groups[i].iter().format(","));
        for pass in 1..=2 {
            if pass == 1 {
                println!("OLD {}", i + 1);
            } else {
                println!("NEW {}", i + 1);
            }
            let mut argsp = args.clone();
            if pass == 1 && old_param.clone().is_some() {
                argsp.push(old_param.clone().unwrap());
            }
            let mut new;
            if pass == 2 || old_param.is_some() {
                new = Command::new("enclone")
            } else {
                new = Command::new(old_exec.clone().unwrap());
            }
            let new = new
                .args(&argsp)
                .arg(&format!("BARCODE={}", bc_arg))
                .arg("NOPAGER")
                .output()
                .unwrap_or_else(|_| panic!("{}", "failed to execute enclone".to_string()));
            if new.status.code() != Some(0) {
                eprint!(
                    "\nenclone call returned nonzero status code, stderr =\n{}",
                    strme(&new.stderr),
                );
                std::process::exit(1);
            }
            print!("{}", strme(&new.stdout));
        }
    }
}
