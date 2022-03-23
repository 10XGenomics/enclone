// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Analyze donor reference and exit.
//
// This displays tables, which are somewhat mangled unless there are at least four donors.

use debruijn::dna_string::DnaString;
use enclone_core::defs::EncloneControl;
use enclone_core::version_string;
use itertools::Itertools;
use std::cmp::min;
use std::env;
use string_utils::*;
use tables::*;
use vdj_ann::refx::{make_vdj_ref_data_core, RefData};
use vector_utils::*;

pub fn analyze_donor_ref(
    refdata: &RefData,
    ctl: &EncloneControl,
    alt_refs: &Vec<(usize, usize, DnaString, usize, bool)>,
) {
    // Analyze donor reference.

    if ctl.gen_opt.external_ref.len() > 0 {
        if ctl.gen_opt.echo {
            let args: Vec<String> = env::args().collect();
            println!("\n{} : {}", env!("CARGO_PKG_VERSION"), version_string());
            println!("{}", args.iter().format(" "));
        }
        let mut erefdata = RefData::new();
        let f = std::fs::read_to_string(&ctl.gen_opt.external_ref).unwrap();
        make_vdj_ref_data_core(&mut erefdata, &f, "", true, true, None);
        let mut refs = Vec::<(String, String, Vec<u8>)>::new(); // {(gene, allele, seq)}

        // Store the external (IMGT) alleles.

        for i in 0..erefdata.refs.len() {
            if erefdata.is_v(i) {
                let allele = erefdata.rheaders_orig[i].between("*", " ");
                refs.push((
                    erefdata.name[i].clone(),
                    allele.to_string(),
                    erefdata.refs[i].to_ascii_vec(),
                ));
            }
        }

        // Store the enclone reference alleles.

        for i in 0..refdata.refs.len() {
            if refdata.is_v(i) {
                refs.push((
                    refdata.name[i].clone(),
                    format!("uref{}", i),
                    refdata.refs[i].to_ascii_vec(),
                ));
            }
        }

        // Store the donor reference alleles;

        for i in 0..alt_refs.len() {
            let donor = alt_refs[i].0;
            let ref_id = alt_refs[i].1;
            let name = &refdata.name[ref_id];
            let alt_seq = &alt_refs[i].2;
            refs.push((
                name.clone(),
                format!("dref{}_{}", i, donor),
                alt_seq.to_ascii_vec(),
            ));
        }

        // Sort the alleles and group by gene.

        refs.sort();
        let mut i = 0;
        while i < refs.len() {
            let j = next_diff1_3(&refs, i as i32) as usize;
            let gene = &refs[i].0;
            let mut alleles = Vec::<(Vec<u8>, String)>::new(); // (sequence, name)
            let mut have_alt = false;
            for k in i..j {
                if refs[k].1.starts_with("dref") {
                    have_alt = true;
                }
                alleles.push((refs[k].2.clone(), refs[k].1.clone()));
            }

            // Delete reference alleles having very low count relative to others.

            let mut to_delete = vec![false; alleles.len()];
            let mut mm = 0;
            for k in 0..alleles.len() {
                if alleles[k].1.starts_with("dref") {
                    let ii = alleles[k].1.between("dref", "_").force_usize();
                    mm = std::cmp::max(mm, alt_refs[ii].3);
                }
            }
            for k in 0..alleles.len() {
                if alleles[k].1.starts_with("dref") {
                    let ii = alleles[k].1.between("dref", "_").force_usize();
                    if alt_refs[ii].4 && alt_refs[ii].3 * 10 < mm {
                        to_delete[k] = true;
                    }
                }
            }
            erase_if(&mut alleles, &to_delete);

            // Proceed.

            if have_alt {
                // Truncate alleles so that they all have the same length.

                let mut m = 1000000;
                for r in 0..alleles.len() {
                    m = min(m, alleles[r].0.len());
                }
                for r in 0..alleles.len() {
                    alleles[r].0.truncate(m);
                }

                // Now alleles = all the alleles for one gene, and there is at least one
                // donor reference allele.  Combine identical alleles, and reorder.

                alleles.sort();
                let mut allelesg = Vec::<(Vec<String>, Vec<u8>)>::new();
                let mut r = 0;
                while r < alleles.len() {
                    let s = next_diff1_2(&alleles, r as i32) as usize;
                    let mut names = Vec::<String>::new();
                    for t in r..s {
                        names.push(alleles[t].1.clone());
                    }
                    allelesg.push((names, alleles[r].0.clone()));
                    r = s;
                }

                // Find the positions at which the alleles differ, and make variant matrix.

                let mut dp = Vec::<usize>::new();
                for p in 0..m {
                    let mut bases = Vec::<u8>::new();
                    for r in 0..allelesg.len() {
                        bases.push(allelesg[r].1[p]);
                    }
                    unique_sort(&mut bases);
                    if bases.len() > 1 {
                        dp.push(p);
                    }
                }
                let mut dm = vec![vec![0; dp.len()]; allelesg.len()];
                for u in 0..dp.len() {
                    for r in 0..allelesg.len() {
                        dm[r][u] = allelesg[r].1[dp[u]];
                    }
                }

                // Make donor matrix.

                let ndonors = ctl.origin_info.donor_list.len();
                let mut dd = vec![vec![false; ndonors]; allelesg.len()];
                for r in 0..allelesg.len() {
                    for k in 0..allelesg[r].0.len() {
                        let n = &allelesg[r].0[k];
                        if n.starts_with("dref") {
                            let d = n.after("_").force_usize();
                            dd[r][d] = true;
                        }
                    }
                }

                // Make IMGT matrix.

                let mut imgts = Vec::<String>::new();
                for r in 0..allelesg.len() {
                    for n in allelesg[r].0.iter() {
                        if !n.starts_with("d") && !n.starts_with("u") {
                            imgts.push(n.to_string());
                        }
                    }
                }
                unique_sort(&mut imgts);
                let nimgt = imgts.len();
                let mut im = vec![vec![false; nimgt]; allelesg.len()];
                for r in 0..allelesg.len() {
                    for k in 0..allelesg[r].0.len() {
                        let n = &allelesg[r].0[k];
                        let p = bin_position(&imgts, n);
                        if p >= 0 {
                            im[r][p as usize] = true;
                        }
                    }
                }

                // Make table, if it won't be too wide.

                let mut log = String::new();
                if dp.len() <= 20 {
                    let mut rows = Vec::<Vec<String>>::new();
                    let mut row = Vec::<String>::new();
                    row.push("allele".to_string());
                    row.push("donor".to_string());
                    for _ in 0..ndonors - 1 {
                        row.push("\\ext".to_string());
                    }
                    if nimgt > 0 {
                        row.push("IMGT".to_string());
                        for _ in 0..nimgt - 1 {
                            row.push("\\ext".to_string());
                        }
                    }
                    if dp.len() > 0 {
                        row.push("position".to_string());
                        for _ in 0..dp.len() - 1 {
                            row.push("\\ext".to_string());
                        }
                    }
                    rows.push(row);
                    let mut row = vec!["".to_string()];
                    row.append(&mut vec!["\\hline".to_string(); ndonors + nimgt + dp.len()]);
                    rows.push(row);
                    let mut row = Vec::<String>::new();
                    row.push("".to_string());
                    for d in 0..ndonors {
                        row.push(format!("{}", d + 1));
                    }
                    for k in 0..nimgt {
                        row.push(imgts[k].clone());
                    }
                    for u in 0..dp.len() {
                        row.push(dp[u].to_string());
                    }
                    rows.push(row);
                    let row = vec!["\\hline".to_string(); ndonors + nimgt + dp.len() + 1];
                    rows.push(row);
                    for r in 0..allelesg.len() {
                        let mut row = Vec::<String>::new();
                        let allele_name = (b'A' + r as u8) as char;
                        let mut an = String::new();
                        an.push(allele_name);
                        for n in allelesg[r].0.iter() {
                            if n.starts_with("uref") {
                                an.push('*');
                                break;
                            }
                        }
                        row.push(an);
                        for d in 0..ndonors {
                            if dd[r][d] {
                                row.push("▓".to_string());
                            } else {
                                row.push(" ".to_string());
                            }
                        }
                        for k in 0..nimgt {
                            if im[r][k] {
                                row.push("▓".to_string());
                            } else {
                                row.push(" ".to_string());
                            }
                        }
                        for u in 0..dp.len() {
                            row.push((allelesg[r].1[dp[u]] as char).to_string());
                        }
                        rows.push(row);
                    }
                    let mut just = b"l|".to_vec();
                    just.append(&mut vec![b'l'; ndonors]);
                    if nimgt > 0 {
                        just.push(b'|');
                        just.append(&mut vec![b'l'; nimgt]);
                    }
                    if dp.len() > 0 {
                        just.push(b'|');
                        just.append(&mut vec![b'l'; dp.len()]);
                    }
                    print_tabular_vbox(&mut log, &rows, 1, &just, false, false);
                }

                // Print.

                println!("\nworking on {}, have {} seqs", gene, alleles.len());
                println!(
                    "alleles differ at {} positions = {}",
                    dp.len(),
                    dp.iter().format(",")
                );
                if log.len() > 0 {
                    log.truncate(log.len() - 1);
                    println!("\n{}", log);
                    println!("* = a universal reference\n");
                }
                for m1 in 0..alleles.len() {
                    for m2 in m1 + 1..alleles.len() {
                        let a1 = &alleles[m1];
                        let a2 = &alleles[m2];
                        let mut diffs = 0;
                        for p in 0..min(a1.0.len(), a2.0.len()) {
                            if a1.0[p] != a2.0[p] {
                                diffs += 1;
                            }
                        }
                        println!(
                            "{} = {} vs {} = {} ==> {} diffs",
                            m1 + 1,
                            a1.1,
                            m2 + 1,
                            a2.1,
                            diffs
                        );
                    }
                }
                for m1 in 0..alleles.len() {
                    let mut best = 1_000_000;
                    let a1 = &alleles[m1];
                    if !a1.1.starts_with("dref") {
                        continue;
                    }
                    for m2 in 0..alleles.len() {
                        let a2 = &alleles[m2];
                        if a2.1.starts_with("dref") {
                            continue;
                        }
                        let mut diffs = 0;
                        for p in 0..min(a1.0.len(), a2.0.len()) {
                            if a1.0[p] != a2.0[p] {
                                diffs += 1;
                            }
                        }
                        best = min(best, diffs);
                    }
                    println!("{} is distance {} from a reference", a1.1, best);
                }
            }
            i = j;
        }
        std::process::exit(0);
    }
}
