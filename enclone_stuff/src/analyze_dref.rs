// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

use debruijn::dna_string::DnaString;
use enclone_core::defs::EncloneControl;
use std::cmp::min;
use string_utils::*;
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
            let mut alleles = Vec::<(String, Vec<u8>)>::new();
            let mut have_alt = false;
            for k in i..j {
                if refs[k].1.starts_with("dref") {
                    have_alt = true;
                }
                alleles.push((refs[k].1.clone(), refs[k].2.clone()));
            }

            if have_alt {

                // Now alleles = all the alleles for one gene, and there is at least one
                // donor reference allele.

                println!("working on {}, have {} seqs", gene, alleles.len());
                for m1 in 0..alleles.len() {
                    for m2 in m1 + 1..alleles.len() {
                        let a1 = &alleles[m1];
                        let a2 = &alleles[m2];
                        let mut diffs = 0;
                        for p in 0..min(a1.1.len(), a2.1.len()) {
                            if a1.1[p] != a2.1[p] {
                                diffs += 1;
                            }
                        }
                        println!(
                            "{} = {} vs {} = {} ==> {} diffs",
                            m1 + 1,
                            a1.0,
                            m2 + 1,
                            a2.0,
                            diffs
                        );
                    }
                }
                for m1 in 0..alleles.len() {
                    let mut best = 1_000_000;
                    let a1 = &alleles[m1];
                    if !a1.0.starts_with("dref") {
                        continue;
                    }
                    for m2 in 0..alleles.len() {
                        let a2 = &alleles[m2];
                        if a2.0.starts_with("dref") {
                            continue;
                        }
                        let mut diffs = 0;
                        for p in 0..min(a1.1.len(), a2.1.len()) {
                            if a1.1[p] != a2.1[p] {
                                diffs += 1;
                            }
                        }
                        best = min(best, diffs);
                    }
                    println!("{} is distance {} from a reference", a1.0, best);
                }
            }
            i = j;
        }
        std::process::exit(0);
    }
}
