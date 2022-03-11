// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Compute sensitivity and specificity statistics.

use enclone_core::defs::ExactClonotype;
use io_utils::fwriteln;
use std::io::Write;
use string_utils::{add_commas, stringme, TextUtils};
use vdj_ann::refx::RefData;
use vector_utils::*;

pub fn sens_spec(
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    refdata: &RefData,
) -> String {
    // Compute intra_mergeable and cross_mergeable, the number of pairs of cells that share
    // a heavy chain with the same V/J gene name and CDR3 length, and likewise for a light chain,
    // and are either from the same donor (in the intra_mergeable case), or from different donors
    // (in the cross_mergeable) case.  We only consider cells having two or three chains.
    //
    // Also compute intra_merged and cross_merged, the number of pairs of cells as above that
    // are placed in the same clonotype.

    let mut intra_mergeable = 0;
    let mut cross_mergeable = 0;
    let mut intra_merged = 0;
    let mut cross_merged = 0;
    // (heavy V+J name, light V+J name, CDR3H_len, CDR3L_len, donor, CDR3H_dna)
    let mut nnlld = Vec::<(String, String, usize, usize, usize, Vec<u8>)>::new();
    for i in 0..exacts.len() {
        let mut donors = Vec::<usize>::new();
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            for i1 in 0..ex.share.len() {
                if ex.share[i1].left {
                    let mut vname = refdata.name[ex.share[i1].v_ref_id].clone();
                    if vname.contains("*") {
                        vname = vname.before("*").to_string();
                    }
                    let mut jname = refdata.name[ex.share[i1].j_ref_id].clone();
                    if jname.contains("*") {
                        jname = jname.before("*").to_string();
                    }
                    let vjh = format!("{}+{}", vname, jname);
                    for i2 in 0..ex.share.len() {
                        if !ex.share[i2].left {
                            let mut vname = refdata.name[ex.share[i2].v_ref_id].clone();
                            if vname.contains("*") {
                                vname = vname.before("*").to_string();
                            }
                            let mut jname = refdata.name[ex.share[i2].j_ref_id].clone();
                            if jname.contains("*") {
                                jname = jname.before("*").to_string();
                            }
                            let vjl = format!("{}+{}", vname, jname);
                            for k in 0..ex.clones.len() {
                                let donor = ex.clones[k][0].donor_index;
                                if donor.is_some() {
                                    nnlld.push((
                                        vjh.clone(),
                                        vjl.clone(),
                                        ex.share[i1].cdr3_aa.len(),
                                        ex.share[i2].cdr3_aa.len(),
                                        donor.unwrap(),
                                        ex.share[i1].cdr3_dna.clone().as_bytes().to_vec(),
                                    ));
                                    donors.push(donor.unwrap());
                                }
                            }
                        }
                    }
                }
            }
        }
        donors.sort();
        let mut freq = Vec::<(u32, usize)>::new();
        make_freq(&donors, &mut freq);
        for m in 0..freq.len() {
            let n = freq[m].0 as usize;
            intra_merged += (n * (n - 1)) / 2;
        }
        for m1 in 0..freq.len() {
            for m2 in m1 + 1..freq.len() {
                let n1 = freq[m1].0 as usize;
                let n2 = freq[m2].0 as usize;
                cross_merged += n1 * n2;
            }
        }
    }
    nnlld.sort();
    let mut i = 0;
    while i < nnlld.len() {
        let mut j = i + 1;
        while j < nnlld.len() {
            if nnlld[j].0 != nnlld[i].0
                || nnlld[j].1 != nnlld[i].1
                || nnlld[j].2 != nnlld[i].2
                || nnlld[j].3 != nnlld[i].3
            {
                break;
            }
            j += 1;
        }
        for k1 in i..j {
            for k2 in k1 + 1..j {
                let mut matches = 0;
                for u in 0..nnlld[k1].5.len() {
                    if nnlld[k1].5[u] == nnlld[k2].5[u] {
                        matches += 1;
                    }
                }
                if matches as f64 / nnlld[k1].5.len() as f64 >= 0.8 {
                    if nnlld[k1].4 == nnlld[k2].4 {
                        intra_mergeable += 1;
                    } else {
                        cross_mergeable += 1;
                    }
                }
            }
        }
        i = j;
    }

    // Report stats.

    let mut log = Vec::<u8>::new();
    fwriteln!(
        log,
        "   • mergeable pairs of +cells from same donor = {}",
        add_commas(intra_mergeable)
    );
    fwriteln!(
        log,
        "   • mergeable pairs of +cells from different donors = {}",
        add_commas(cross_mergeable)
    );
    fwriteln!(
        log,
        "   • merged pairs of +cells from same donor = {}",
        add_commas(intra_merged)
    );
    fwriteln!(
        log,
        "   • merged pairs of +cells from different donors = {}",
        add_commas(cross_merged)
    );
    let error_rate = cross_merged as f64 / cross_mergeable as f64;
    fwriteln!(log, "   • error rate = {:.2} x 10^-3", error_rate * 1_000.0);
    let inferred_errors = (error_rate * intra_mergeable as f64).round() as usize;
    fwriteln!(
        log,
        "   • inferred errors = {}",
        add_commas(inferred_errors)
    );
    let adjusted = intra_merged as isize - inferred_errors as isize;
    if adjusted >= 0 {
        fwriteln!(
            log,
            "   • adjusted true merges = {}",
            add_commas(adjusted as usize)
        );
    } else {
        fwriteln!(
            log,
            "   • adjusted true merges = -{}",
            add_commas(-adjusted as usize)
        );
    }
    stringme(&log)
}
