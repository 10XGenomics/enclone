// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Assess immcantation clonotyping.  Pass a single argument, which is to be an immcantation
// output file.
//
// 1. First run build_immcantation_inputs.
// 2. Then run immcantation.
// 3. Then run this code.
//
// We assume that the argument to build_immcantation_inputs is a list of ids from
// enclone.testdata.bcr.gex.
//
// We ignore lines for which clone_id is unspecified.

use io_utils::*;
use pretty_trace::*;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use string_utils::*;
use vector_utils::*;

pub fn main() {
    PrettyTrace::new().on();

    // Get dataset classification by sample.

    let mut to_donor = HashMap::<usize, usize>::new();
    {
        let f = include_str!["../../../enclone/src/enclone.testdata.bcr.gex"];
        let mut donor = 0;
        for s in f.lines() {
            if s.starts_with("DONOR=") {
                if s.after("DONOR=").parse::<usize>().is_err() {
                    break;
                }
                donor = s.after("DONOR=").force_usize();
            }
            if s.starts_with("BCR=") {
                let ids0 = s.after("BCR=").split(',').collect::<Vec<&str>>();
                let mut ids = Vec::<usize>::new();
                for id in ids0.iter() {
                    if id.contains('-') {
                        let start = id.before("-").force_usize();
                        let stop = id.after("-").force_usize();
                        for p in start..=stop {
                            ids.push(p);
                        }
                    } else {
                        ids.push(id.force_usize());
                    }
                }
                for id in ids.iter() {
                    to_donor.insert(*id, donor);
                }
            }
        }
    }

    // Parse the immcantation output file.

    let args: Vec<String> = env::args().collect();
    if args.len() != 2 || (!args[1].ends_with(".csv") && !args[1].ends_with(".tsv")) {
        eprintln!("\nPlease run with a single argument which is a path to a CSV or TSV file.\n");
        std::process::exit(1);
    }
    let f = open_for_read![&args[1]];
    let mut n_seq = None;
    let mut n_clone = None;
    let mut assignments = Vec::<(usize, String, String)>::new();
    let mut datasets = Vec::<String>::new();
    let mut unassigned = 0;
    for (i, line) in f.lines().enumerate() {
        let s = line.unwrap();
        let fields;
        if args[1].ends_with(".csv") {
            fields = parse_csv(&s);
        } else {
            fields = s.split('\t').map(str::to_owned).collect();
        }
        if i == 0 {
            for j in 0..fields.len() {
                if fields[j] == "sequence_id" {
                    n_seq = Some(j);
                } else if fields[j] == "clone_id" {
                    n_clone = Some(j);
                }
            }
        } else {
            let seq_id = &fields[n_seq.unwrap()];
            let clone_id = &fields[n_clone.unwrap()];
            if clone_id == "" || clone_id == "NA" {
                unassigned += 1;
                continue;
            }
            if !clone_id.parse::<usize>().is_ok() {
                eprintln!("\nProblem parsing line {}, clone_id = {}.", i + 1, clone_id);
            }
            let clone_id = clone_id.force_usize();
            let dataset = seq_id.before("_").to_string();
            datasets.push(dataset.clone());
            let barcode = seq_id.between("_", "_").to_string();
            assignments.push((clone_id, dataset, barcode));
        }
    }
    unique_sort(&mut datasets);
    use itertools::Itertools;
    println!("\ndatasets: {}", datasets.iter().format(","));
    unique_sort(&mut assignments);
    println!("\n{} cells assigned to clonotypes", assignments.len());
    println!("{} lines with clonotype unspecified\n", unassigned);
    println!("top clonotype sizes:\n");
    let mut sizes = Vec::<usize>::new();
    let mut i = 0;
    while i < assignments.len() {
        let j = next_diff1_3(&assignments, i as i32) as usize;
        sizes.push(j - i);
        i = j;
    }
    reverse_sort(&mut sizes);
    for i in 0..10 {
        println!("{}", sizes[i]);
    }
    println!("");
}
