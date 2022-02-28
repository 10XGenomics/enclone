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
use itertools::Itertools;
use pretty_trace::*;
use std::cmp::max;
use std::collections::HashMap;
use std::env;
use std::io::BufRead;
use string_utils::*;
use vector_utils::*;

pub fn main() {
    PrettyTrace::new().on();

    // Get dataset classification by sample.

    let mut to_donor = HashMap::<String, usize>::new();
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
                let mut ids = Vec::<String>::new();
                for id in ids0.iter() {
                    if id.contains('-') {
                        let start = id.before("-").force_usize();
                        let stop = id.after("-").force_usize();
                        for p in start..=stop {
                            ids.push(p.to_string());
                        }
                    } else {
                        ids.push(id.to_string());
                    }
                }
                for id in ids.iter() {
                    to_donor.insert(id.to_string(), donor);
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
            let dataset = seq_id.between("-", "_").to_string();
            datasets.push(dataset.clone());
            let barcode = format!("{}-1", seq_id.before("-"));
            assignments.push((clone_id, dataset, barcode));
        }
    }
    unique_sort(&mut datasets);
    unique_sort(&mut assignments);

    // Create clonotypes.  Note that we allow a cell to be assigned to more than one clonotype.

    let mut max_id = 0;
    for i in 0..assignments.len() {
        max_id = max(max_id, assignments[i].0);
    }
    let mut clono = vec![Vec::<(usize, String, String)>::new(); max_id + 1];
    for i in 0..assignments.len() {
        let dataset = &assignments[i].1;
        if !to_donor.contains_key(dataset) {
            eprintln!("\nCan't find donor for dataset {}.\n", dataset);
        }
        let donor = to_donor[dataset];
        let barcode = assignments[i].2.clone();
        clono[assignments[i].0].push((donor, dataset.to_string(), barcode));
    }
    for i in 0..clono.len() {
        clono[i].sort();
    }

    let mut mc = 0;
    let mut mci = 0;
    for i in 0..clono.len() {
        if clono[i].len() > mc {
            mc = clono[i].len();
            mci = i;
        }
    }
    let mut bcs = Vec::<String>::new();
    let mut datasets_top = Vec::<String>::new();
    for x in clono[mci].iter() {
        bcs.push(x.2.clone());
        datasets_top.push(x.1.clone());
    }
    unique_sort(&mut datasets_top);
    println!("\nbarcodes in top clonotype = {}", bcs.iter().format(","));
    println!(
        "\ndatasets in top clonotype = {}",
        datasets_top.iter().format(",")
    );

    // Generate some stats.

    println!("\ndatasets: {}", datasets.iter().format(","));
    println!("\n{} cells assigned to clonotypes", assignments.len());
    println!("{} lines with clonotype unspecified\n", unassigned);
    let mut sizes = Vec::<usize>::new();
    for i in 0..clono.len() {
        if clono[i].len() > 0 {
            sizes.push(clono[i].len());
        }
    }
    reverse_sort(&mut sizes);
    if sizes.len() > 10 {
        sizes.truncate(10);
    }
    println!("top clonotype sizes: {}\n", sizes.iter().format(", "));

    // Compute sensitivity/specificity stats.

    let mut max_donor = 0;
    for i in 0..clono.len() {
        for j in 0..clono[i].len() {
            max_donor = max(max_donor, clono[i][j].0);
        }
    }
    let mut cells_by_donor = vec![0 as usize; max_donor + 1];
    let mut merges2 = 0;
    let mut mixes = 0;
    let mut wrongotypes = 0;
    let mut clono2 = 0;
    for i in 0..clono.len() {
        let mut wrong = false;
        let mut cells_by_donor_this = vec![0; max_donor + 1];
        for c in clono[i].iter() {
            cells_by_donor[c.0] += 1;
            cells_by_donor_this[c.0] += 1;
        }
        for j1 in 0..clono[i].len() {
            for j2 in j1 + 1..clono[i].len() {
                if clono[i][j1].0 != clono[i][j2].0 {
                    mixes += 1;
                    wrong = true;
                }
            }
        }
        if wrong {
            wrongotypes += 1;
        }
        if clono[i].len() > 1 {
            clono2 += 1;
        }
        for j in 0..cells_by_donor_this.len() {
            let n = cells_by_donor_this[j];
            if n > 1 {
                merges2 += (n * (n - 1)) / 2;
            }
        }
    }
    let mut cross = 0;
    let mut intra = 0;
    for i1 in 0..cells_by_donor.len() {
        if cells_by_donor[i1] > 1 {
            intra += cells_by_donor[i1] * (cells_by_donor[i1] - 1) / 2;
        }
        for i2 in i1 + 1..cells_by_donor.len() {
            cross += cells_by_donor[i1] * cells_by_donor[i2];
        }
    }
    println!("number of intradonor comparisons = {}", add_commas(intra));
    println!(
        "number of intradonor cell-cell merges (quadratic) = {}",
        add_commas(merges2)
    );
    println!("number of cross-donor comparisons = {}", add_commas(cross));
    println!(
        "number of cross-donor comparisons that mix donors = {}",
        add_commas(mixes)
    );
    println!("number of mixed clonotypes = {}", wrongotypes);
    println!(
        "number of clonotypes having at least two cells = {}",
        clono2
    );
    let rate = (mixes as f64) * 1_000_000_000.0 / (cross as f64);
    println!("rate of cross donor mixing = {:.2} x 10^-9", rate);
    let bogus = (intra as f64) * (mixes as f64) / (cross as f64);
    println!(
        "estimated number of false intradonor merges = {}\n",
        add_commas(bogus.round() as usize)
    );
}
