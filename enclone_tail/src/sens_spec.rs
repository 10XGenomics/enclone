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
    let mut nnlld = Vec::<(String, String, usize, usize, usize)>::new();
    for i in 0..exacts.len() {
        let mut donors = Vec::<usize>::new();
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            let mut heavy = Vec::<(String, usize)>::new();
            let mut light = Vec::<(String, usize)>::new();
            for i in 0..ex.share.len() {
                let mut vname = refdata.name[ex.share[i].v_ref_id].clone();
                if vname.contains("*") {
                    vname = vname.before("*").to_string();
                }
                let mut jname = refdata.name[ex.share[i].j_ref_id].clone();
                if jname.contains("*") {
                    jname = jname.before("*").to_string();
                }
                let vj = format!("{}+{}", vname, jname);
                if ex.share[i].left {
                    heavy.push((vj, ex.share[i].cdr3_aa.len()));
                } else {
                    light.push((vj, ex.share[i].cdr3_aa.len()));
                }
            }
            unique_sort(&mut heavy);
            unique_sort(&mut light);
            if heavy.len() > 0 || light.len() > 0 {
                for j in 0..ex.clones.len() {
                    if ex.clones[j][0].donor_index.is_some() {
                        donors.push(ex.clones[j][0].donor_index.unwrap());
                    }
                }
            }
            for i1 in 0..heavy.len() {
                for i2 in 0..light.len() {
                    for j in 0..ex.clones.len() {
                        if ex.clones[j][0].donor_index.is_some() {
                            nnlld.push((
                                heavy[i1].0.clone(),
                                light[i2].0.clone(),
                                heavy[i1].1,
                                light[i2].1,
                                ex.clones[j][0].donor_index.unwrap(),
                            ));
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
        let mut donors = Vec::<usize>::new();
        for k in i..j {
            donors.push(nnlld[k].4);
        }
        donors.sort();
        let mut freq = Vec::<(u32, usize)>::new();
        make_freq(&donors, &mut freq);
        for m in 0..freq.len() {
            let n = freq[m].0 as usize;
            intra_mergeable += (n * (n - 1)) / 2;
        }
        for m1 in 0..freq.len() {
            for m2 in m1 + 1..freq.len() {
                let n1 = freq[m1].0 as usize;
                let n2 = freq[m2].0 as usize;
                cross_mergeable += n1 * n2;
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
    fwriteln!(
        log,
        "   • error rate = {:.2} x 10^-6",
        error_rate * 1_000_000.0
    );
    let inferred_errors = (error_rate * intra_mergeable as f64).round() as usize;
    fwriteln!(
        log,
        "   • inferred errors = {}",
        add_commas(inferred_errors)
    );
    fwriteln!(
        log,
        "   • adjusted true merges = {}",
        add_commas(intra_merged - inferred_errors)
    );
    stringme(&log)
}
