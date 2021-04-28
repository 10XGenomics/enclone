// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Find the optimal D segment, the runner up, and the delta between the scores.  This uses
// the donor V and J segments that are assigned to the clonotype.  Note that the optimal D
// segment may be null.  This is obvious from looking at data.

use crate::align_to_vdj_ref::*;
use crate::defs::*;
use enclone_proto::types::*;
use vdj_ann::refx::*;

pub fn vflank(_seq: &[u8], vref: &[u8]) -> usize {
    let mut flank = 15;
    if flank > vref.len() {
        flank = vref.len();
    }
    flank
}

pub fn jflank(seq: &[u8], jref: &[u8]) -> usize {
    let flank = 15;
    if flank > jref.len() {
        return jref.len();
    }
    const MATCHLEN: usize = 5;
    let mut start = 0;
    for i in 0..=jref.len() - MATCHLEN {
        let mut matchlen = 0;
        for j in 0..MATCHLEN {
            if seq[seq.len() - jref.len() + i + j] != jref[i + j] {
                break;
            }
            matchlen += 1;
        }
        if matchlen == MATCHLEN {
            start = i;
            break;
        }
    }
    flank + start
}

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

    // Go through every D segment, or possibly every concatenation of D segments.

    let mut todo = Vec::<Vec<usize>>::new();
    todo.push(vec![]);
    for i in refdata.ds.iter() {
        todo.push(vec![*i]);
    }
    if ex.share[mid].cdr3_aa.len() >= 20 {
        for i1 in refdata.ds.iter() {
            for i2 in refdata.ds.iter() {
                todo.push(vec![*i1, *i2]);
            }
        }
    }
    let mut ds = Vec::<Vec<usize>>::new();
    let mut counts = Vec::<f64>::new();
    for di in 0..todo.len() {
        // Start to build reference concatenation.  First append the V segment.

        let mut concat = Vec::<u8>::new();
        let mut vref = refdata.refs[rsi.vids[col]].to_ascii_vec();
        if rsi.vpids[col].is_none() {
        } else {
            vref = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
        }
        let vstart = vref.len() - vflank(&tig, &vref);
        vref = vref[vstart..vref.len()].to_vec();
        concat.append(&mut vref.clone());

        // Append the D segment or segments.

        let mut dref = Vec::<u8>::new();
        let mut d2ref = Vec::<u8>::new();
        let mut drefname = String::new();
        for j in 0..todo[di].len() {
            let d = todo[di][j];
            if j == 0 {
                dref = refdata.refs[d].to_ascii_vec();
            } else if j == 1 {
                d2ref = refdata.refs[d].to_ascii_vec();
            }
            if j > 0 {
                drefname += ":";
            }
            drefname += &mut refdata.name[d].clone();
        }
        concat.append(&mut dref.clone());
        concat.append(&mut d2ref.clone());

        // Append the J segment.

        let jref = refdata.refs[rsi.jids[col]].to_ascii_vec();
        let jend = jflank(&tig, &jref);

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

        let (_ops, count) = align_to_vdj_ref(&seq, &vref, &dref, &d2ref, &jref, &drefname, true);

        counts.push(count);
        ds.push(todo[di].clone());
        if count > comp {
            comp = count;
        }
    }

    // Reverse sort sync (counts, ds).

    let mut counts_ds = Vec::new();
    for i in 0..counts.len() {
        counts_ds.push((counts[i], ds[i].clone()));
    }
    counts_ds.sort_by(|a, b| b.partial_cmp(a).unwrap()); // reverse sort
    counts.clear();
    ds.clear();
    for i in 0..counts_ds.len() {
        counts.push(counts_ds[i].0);
        ds.push(counts_ds[i].1.clone());
    }

    // Proceed.

    let mut comp = 0.0;
    let mut best_d = Vec::new();
    if counts.len() > 0 {
        comp = counts[0];
        best_d = ds[0].clone();
    }
    let mut second_comp = 0.0;
    let mut best_d2 = Vec::new();
    if counts.len() > 1 {
        second_comp = counts[1];
        best_d2 = ds[1].clone();
    }
    *opt = best_d;
    *opt2 = best_d2;
    *delta = comp - second_comp;
}
