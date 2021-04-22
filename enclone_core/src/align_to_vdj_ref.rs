// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Align a sequence to a concatenated V(D)J reference, encouraging gaps at the junction points.
//
// results of current version:
//
// dataset       n     inconsistent%
// ---------------------------------
// BCR=123085    712        2.25
// BI=1         7653       24.06
// BI=3        13010       13.11
// BI=4        27766       17.02

use bio_edit::alignment::pairwise::*;
use bio_edit::alignment::AlignmentMode;

pub fn align_to_vdj_ref(
    seq: &[u8],
    vref: &[u8],
    dref: &[u8],
    jref: &[u8],
) -> bio_edit::alignment::Alignment {
    // Define penalties.

    let mismatch = -2;
    let gap_open = -12;
    let gap_extend = -2;
    let gap_open_at_boundary = -6;
    let gap_extend_at_boundary = -1;

    // Build concatenation.

    let mut concat = Vec::<u8>::new();
    concat.append(&mut vref.to_vec());
    concat.append(&mut dref.to_vec());
    concat.append(&mut jref.to_vec());

    // Make alignment.

    let mut scoring = Scoring::from_scores(gap_open, gap_extend, -mismatch, mismatch);
    scoring.xclip_prefix = MIN_SCORE;
    scoring.xclip_suffix = MIN_SCORE;
    scoring.yclip_prefix = 0;
    scoring.yclip_suffix = 0;
    let mut aligner = Aligner::with_scoring(scoring);
    let mut gap_open_fn = vec![0_i32; concat.len() + 1];
    for j in 1..=concat.len() {
        if j as usize == vref.len() || j as usize == vref.len() + dref.len() {
            gap_open_fn[j] = gap_open_at_boundary;
        } else {
            gap_open_fn[j] = gap_open;
        }
    }
    let mut gap_extend_fn = vec![0_i32; concat.len() + 1];
    for j in 1..=concat.len() {
        if j as usize == vref.len() || j as usize == vref.len() + dref.len() {
            gap_extend_fn[j] = gap_extend_at_boundary;
        } else {
            gap_extend_fn[j] = gap_extend;
        }
    }
    let mut al = aligner.custom_with_gap_fns(&seq, &concat, &gap_open_fn, &gap_extend_fn);
    al.mode = AlignmentMode::Semiglobal;
    al.filter_clip_operations();
    al
}
