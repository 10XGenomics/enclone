// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

use debruijn::dna_string::DnaString;
use enclone_core::defs::EncloneControl;
use itertools::Itertools;
use std::cmp::min;
use string_utils::*;
use tables::*;
use vdj_ann::refx::{make_vdj_ref_data_core, RefData};
use vector_utils::*;

pub fn analyze_donor_ref(
    refdata: &RefData,
    ctl: &EncloneControl,
    alt_refs: &Vec<(usize,usize,DnaString)>,
) {

    // Analyze donor reference.

    if ctl.gen_opt.external_ref.len() > 0 {
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
            let _donor = alt_refs[i].0;
            let ref_id = alt_refs[i].1;
            let name = &refdata.name[ref_id];
            let alt_seq = &alt_refs[i].2;
            refs.push((name.clone(), format!("dref{}", i), alt_seq.to_ascii_vec()));
        }

        // Sort the alleles and group by gene.

        refs.sort();
        let mut i = 0;
        while i < refs.len() {
            let j = next_diff1_3(&refs, i as i32) as usize;
            let gene = &refs[i].0;
            let mut alleles = Vec::<(Vec<u8>, String)>::new(); // (name, sequence)
            let mut have_alt = false;
            for k in i..j {
                if refs[k].1.starts_with("dref") {
                    have_alt = true;
                }
                alleles.push((refs[k].2.clone(), refs[k].1.clone()));
            }

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

                // Find the positions at which the alleles differ, and make a matrix.

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

                // Make table, if it won't be too wide.

                let mut log = String::new();
                if dp.len() <= 20 {
                    let mut rows = Vec::<Vec<String>>::new();
                    let mut row = Vec::<String>::new();
                    row.push("allele".to_string());
                    for u in 0..dp.len() {
                        row.push(dp[u].to_string());
                    }
                    rows.push(row);
                    for r in 0..allelesg.len() {
                        let mut row = Vec::<String>::new();
                        row.push(format!("{}", allelesg[r].0.iter().format(",")));
                        for u in 0..dp.len() {
                            row.push((allelesg[r].1[dp[u]] as char).to_string());
                        }
                        rows.push(row);
                    }
                    let mut just = b"l|".to_vec();
                    just.append(&mut vec![b'l'; dp.len()]);
                    print_tabular_vbox(&mut log, &rows, 2, &just, false, false);
                }

                // Print.

                println!("\nworking on {}, have {} seqs", gene, alleles.len());
                println!("alleles differ at {} positions = {}", dp.len(), dp.iter().format(","));
                println!("{}", log);
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
