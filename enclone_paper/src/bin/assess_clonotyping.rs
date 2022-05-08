// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Assess alternative clonotyping algorithm.
//
// Usage is something like this:
//
// 1. Let ids by the list of ids.  These need to be in enclone.testdata.bcr.gex.
//
// 2. build_clonotyping_inputs ids
//
// 3. Run alternative clonotyping software to yield clones.tsv.
//
// 4. enclone BCR=ids MIX_DONORS MIN_CHAINS_EXACT=2 NOPRINT POUT=stdout PCELL
//            PCOLS=datasets_cell,barcode PCOLS_SHOW=dataset,barcode > post_filter.csv
//
// 5. enclone BCR=@test MIX_DONORS MIN_CHAINS_EXACT=2 NOPRINT SUMMARY
//
// to get cross = number of cross-donor comparisons.
//
// 6. assess_clonotyping clones.tsv post_filter.csv cross
//
// Optional argument: VERBOSE.  Print the mixed clonotypes and exit.

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

    // Parse arguments.

    let args: Vec<String> = env::args().collect();
    if (args.len() != 4 && args.len() != 5)
        || (!args[1].ends_with(".csv") && !args[1].ends_with(".tsv"))
        || !args[2].ends_with(".csv")
        || !args[3].parse::<usize>().is_ok()
        || (args.len() == 5 && args[4] != "VERBOSE")
    {
        eprintln!("\nPlease read the usage in the source file.\n");
        std::process::exit(1);
    }
    let clone_file = &args[1];
    let filter_file = &args[2];
    let cross = args[3].force_usize();
    let verbose = args.len() == 5;

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

    // Load post filter data.

    let mut post_filter = Vec::<(String, String)>::new();
    let f = open_for_read![&filter_file];
    for (i, line) in f.lines().enumerate() {
        let s = line.unwrap();
        if i == 0 {
            assert_eq!(s, "dataset,barcode");
        } else {
            post_filter.push((s.before(",").to_string(), s.after(",").to_string()));
        }
    }
    post_filter.sort();

    // Parse the clonotyping output file.

    let mut assignments = Vec::<(String, String, String)>::new();
    let mut clone_ids = Vec::<String>::new();
    let mut datasets = Vec::<String>::new();
    let mut unassigned = 0;
    {
        let mut n_seq = None;
        let mut n_clone = None;
        let f = open_for_read![&clone_file];
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
                let mut seq_id = fields[n_seq.unwrap()].clone();
                if seq_id.contains("-1-") {
                    seq_id = seq_id.replace("-1-", "-");
                }
                let clone_id = fields[n_clone.unwrap()].clone();
                if clone_id == "" || clone_id == "NA" {
                    unassigned += 1;
                    continue;
                }
                let dataset = seq_id.between("-", "_").to_string();

                datasets.push(dataset.clone());
                let barcode = format!("{}-1", seq_id.before("-"));
                if bin_member(&post_filter, &(dataset.clone(), barcode.clone())) {
                    assignments.push((clone_id.clone(), dataset, barcode));
                    clone_ids.push(clone_id);
                }
            }
        }
        unique_sort(&mut datasets);
        unique_sort(&mut assignments);
    }
    unique_sort(&mut clone_ids);
    let mut assignments2 = Vec::<(usize, String, String)>::new();
    for i in 0..assignments.len() {
        let id = bin_position(&clone_ids, &assignments[i].0) as usize;
        assignments2.push((id, assignments[i].1.clone(), assignments[i].2.clone()));
    }
    unique_sort(&mut assignments2);
    let assignments = assignments2;

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

    // Print mixed clonotypes.

    if verbose {
        let mut mixes = 0;
        for i in 0..clono.len() {
            let c = &clono[i];
            let mut donors = Vec::<usize>::new();
            for j in 0..c.len() {
                donors.push(c[j].0);
            }
            unique_sort(&mut donors);
            if donors.len() > 1 {
                mixes += 1;
                println!("\nMIXED CLONOTYPE {mixes}");
                for j in 0..c.len() {
                    println!("[{}] {} {} {}", j + 1, c[j].0, c[j].1, c[j].2,);
                }
            }
        }
        std::process::exit(0);
    }

    // Print barcodes in top clonotype.

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
}
