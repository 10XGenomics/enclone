// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Find the optimal D segment, the runner up, and the delta between the scores.  This uses
// the donor V and J segments that are assigned to the clonotype.  Note that the optimal D
// segment may be null.  This is obvious from looking at data.

use crate::align_to_vdj_ref::*;
use crate::defs::*;
use enclone_proto::types::*;
use std::cmp::min;
use vdj_ann::refx::*;
use vector_utils::*;

pub fn opt_d(
    ex: &ExactClonotype,
    col: usize,
    u: usize,
    rsi: &ColInfo,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    opt: &mut Option<usize>,
    opt2: &mut Option<usize>,
    delta: &mut usize,
) {
    let mid = rsi.mat[col][u].unwrap();
    assert!(ex.share[mid].left);
    let mut comp = 1000000;
    let td = &ex.share[mid];
    let tig = &td.seq;

    // Go through passes.  If IGH/TRB, we go through every D segment.  Otherwise
    // there is just one pass.

    let z = refdata.ds.len();
    let mut ds = Vec::<Option<usize>>::new();
    let mut counts = Vec::<i32>::new();
    for di in 0..=z {
        let mut d = 0;
        if di < z {
            d = refdata.ds[di];
        }
        const FLANK: usize = 35;

        // Start to build reference concatenation.  First append the V segment.

        let mut concat = Vec::<u8>::new();
        let mut vref = refdata.refs[rsi.vids[col]].to_ascii_vec();
        if rsi.vpids[col].is_none() {
        } else {
            vref = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
        }
        let mut vstart = 0;
        if vref.len() >= FLANK {
            vstart = vref.len() - FLANK;
        }
        vref = vref[vstart..vref.len()].to_vec();
        let mut _refsplits = Vec::<usize>::new();
        _refsplits.push(0);
        concat.append(&mut vref.clone());
        _refsplits.push(concat.len());

        // Append the D segment if IGH/TRB.

        let mut dref = Vec::<u8>::new();
        if ex.share[mid].left && di < z {
            dref = refdata.refs[d].to_ascii_vec();
            concat.append(&mut dref.clone());
            _refsplits.push(concat.len());
        }

        // Append the J segment.

        let jref = refdata.refs[rsi.jids[col]].to_ascii_vec();
        let jend = min(FLANK, jref.len());
        let jref = jref[0..jend].to_vec();
        concat.append(&mut jref.clone());
        _refsplits.push(concat.len());

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

        let al = align_to_vdj_ref(&seq, &vref, &dref, &jref);

        let count = al.score;
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
    sort_sync2(&mut counts, &mut ds);
    counts.reverse();
    ds.reverse();
    let mut comp = 0 as i32;
    let mut best_d = None;
    if counts.len() > 0 {
        comp = counts[0];
        best_d = ds[0];
    }
    let mut second_comp = 0 as i32;
    let mut best_d2 = None;
    if counts.len() > 1 {
        second_comp = counts[1];
        best_d2 = ds[1];
    }
    *opt = best_d;
    *opt2 = best_d2;
    *delta = (comp - second_comp) as usize;
}
