// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Find the optimal D segment, the runner up, and the delta between the scores.  This uses
// The donor V and J segments that are assigned to the clonotype.

use crate::defs::*;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
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
    opt: &mut usize,
    opt2: &mut usize,
    delta: &mut usize,
) {
    let mid = rsi.mat[col][u].unwrap();
    assert!(ex.share[mid].left);
    let mut comp = 1000000;
    let td = &ex.share[mid];
    let tig = &td.seq;
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let mut aligner = Aligner::new(-6, -1, &score);

    // Go through passes.  If IGH/TRB, we go through every D segment.  Otherwise
    // there is just one pass.

    let z = refdata.ds.len();
    let mut ds = Vec::<usize>::new();
    let mut counts = Vec::<usize>::new();
    for di in 0..z {
        let d = refdata.ds[di];
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
        concat.append(&mut vref.clone());

        // Append the D segment if IGH/TRB.

        if ex.share[mid].left {
            let mut x = refdata.refs[d].to_ascii_vec();
            concat.append(&mut x);
        }

        // Append the J segment.

        let jref = refdata.refs[rsi.jids[col]].to_ascii_vec();
        let jend = min(FLANK, jref.len());
        let mut x = jref[0..jend].to_vec();
        concat.append(&mut x);

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
        let al = aligner.semiglobal(&seq, &concat);
        let mut m = 0;
        let mut pos = al.xstart;
        let mut count = 0;
        let start = td.cdr3_start - td.ins_len() - seq_start as usize;
        let stop = td.j_stop - td.v_start - seq_start as usize;
        while m < al.operations.len() {
            let n = next_diff(&al.operations, m);
            match al.operations[m] {
                Match => {
                    pos += 1;
                }
                Subst => {
                    if pos >= start && pos < stop {
                        count += 1;
                    }
                    pos += 1;
                }
                Del => {
                    if pos >= start && pos < stop {
                        count += n - m;
                    }
                    pos += n - m;
                    m = n - 1;
                }
                Ins => {
                    if pos >= start && pos < stop {
                        count += n - m;
                    }
                    m = n - 1;
                }
                _ => {}
            };
            m += 1;
        }
        counts.push(count);
        ds.push(d);
        if count < comp {
            comp = count;
        }
    }
    sort_sync2(&mut counts, &mut ds);
    let mut comp = 0;
    let mut best_d = 0;
    if counts.len() > 0 {
        comp = counts[0];
        best_d = ds[0];
    }
    let mut second_comp = 0;
    let mut best_d2 = 0;
    if counts.len() > 1 {
        second_comp = counts[1];
        best_d2 = ds[1];
    }
    *opt = best_d;
    *opt2 = best_d2;
    *delta = second_comp - comp;
}
