// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Find the optimal D segment, the runner up, and the delta between the scores.  This uses
// the donor V and J segments that are assigned to the clonotype.  Note that the optimal D
// segment may be null.  This is obvious from looking at data.

use crate::align_to_vdj_ref::*;
use crate::defs::*;
use enclone_proto::types::*;
use std::cmp::min;
use vdj_ann::refx::*;

pub fn opt_d(
    ex: &ExactClonotype,
    col: usize,
    u: usize,
    rsi: &ColInfo,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    opt: &mut Vec<usize>,
    opt2: &mut Vec<usize>,
    delta: &mut f64,
) {
    let mid = rsi.mat[col][u].unwrap();
    assert!(ex.share[mid].left);
    let mut comp = 1000000.0;
    let td = &ex.share[mid];
    let tig = &td.seq_del;

    // Go through passes.  If IGH/TRB, we go through every D segment.  Otherwise
    // there is just one pass.

    let z = refdata.ds.len();
    let mut ds = Vec::<Option<usize>>::new();
    let mut counts = Vec::<f64>::new();
    for di in 0..=z {
        let mut d = 0;
        if di < z {
            d = refdata.ds[di];
        }
        const LFLANK: usize = 15;
        const RFLANK: usize = 35;

        // Start to build reference concatenation.  First append the V segment.

        let mut concat = Vec::<u8>::new();
        let mut vref = refdata.refs[rsi.vids[col]].to_ascii_vec();
        if rsi.vpids[col].is_none() {
        } else {
            vref = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
        }
        let mut vstart = 0;
        if vref.len() >= LFLANK {
            vstart = vref.len() - LFLANK;
        }
        vref = vref[vstart..vref.len()].to_vec();
        concat.append(&mut vref.clone());

        // Append the D segment if IGH/TRB.

        let mut dref = Vec::<u8>::new();
        let mut drefname = String::new();
        if ex.share[mid].left && di < z {
            dref = refdata.refs[d].to_ascii_vec();
            drefname = refdata.name[d].clone();
            concat.append(&mut dref.clone());
        }

        // Append the J segment.

        let jref = refdata.refs[rsi.jids[col]].to_ascii_vec();
        let jend = min(RFLANK, jref.len());

        // Align the V..J sequence on the contig to the reference concatenation.

        let mut seq_start = vstart as isize;
        // probably not exactly right
        if ex.share[mid].annv.len() > 1 {
            let q1 = ex.share[mid].annv[0].0 + ex.share[mid].annv[0].1;
            let q2 = ex.share[mid].annv[1].0;

            seq_start += q1 as isize - q2 as isize;
        }
        let seq_end = tig.len() - (jref.len() - jend);
        let seq = tig[seq_start as usize..seq_end].to_vec();
        let jref = jref[0..jend].to_vec();
        concat.append(&mut jref.clone());

        let (_ops, count) = align_to_vdj_ref(&seq, &vref, &dref, &jref, &drefname, true);

        counts.push(count);
        if di < z {
            ds.push(Some(d));
        } else {
            ds.push(None);
        }
        if count > comp {
            comp = count;
        }
    }

    // Reverse sort sync (counts, ds).

    let mut counts_ds = Vec::new();
    for i in 0..counts.len() {
        counts_ds.push((counts[i], ds[i]));
    }
    counts_ds.sort_by(|a, b| b.partial_cmp(a).unwrap()); // reverse sort
    counts.clear();
    ds.clear();
    for i in 0..counts_ds.len() {
        counts.push(counts_ds[i].0);
        ds.push(counts_ds[i].1);
    }

    // Proceed.

    let mut comp = 0.0;
    let mut best_d = None;
    if counts.len() > 0 {
        comp = counts[0];
        best_d = ds[0];
    }
    let mut second_comp = 0.0;
    let mut best_d2 = None;
    if counts.len() > 1 {
        second_comp = counts[1];
        best_d2 = ds[1];
    }
    if best_d.is_some() {
        *opt = vec![best_d.unwrap()];
    } else {
        *opt = vec![];
    }
    if best_d2.is_some() {
        *opt2 = vec![best_d2.unwrap()];
    } else {
        *opt2 = vec![];
    }
    *delta = comp - second_comp;
}
