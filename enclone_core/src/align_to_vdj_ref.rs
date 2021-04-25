// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Align a sequence to a concatenated V(D)J reference, encouraging gaps at the junction points.
//
// results of various parameter choices
//
// enclone BI=1-4,9 BUILT_IN GVARS=d_inconsistent_%,d_inconsistent_n NOPRINT
//
// d_inconsistent_n = 53373
// d_inconsistent_% = 14.40

use bio_edit::alignment::pairwise::*;
use bio_edit::alignment::AlignmentMode;
use bio_edit::alignment::AlignmentOperation::*;

pub fn align_to_vdj_ref(
    seq: &[u8],
    vref: &[u8],
    dref: &[u8],
    jref: &[u8],
    drefname: &str, // useful for debugging
) -> (bio_edit::alignment::Alignment, f64) {
    // Define penalties.

    let matchp = 2;
    let mismatch = -2;
    let gap_open = -12;
    let gap_extend = -2;
    let gap_open_at_boundary = -4;
    let gap_extend_at_boundary = -1;

    // Build concatenation.

    let mut concat = Vec::<u8>::new();
    concat.append(&mut vref.to_vec());
    concat.append(&mut dref.to_vec());
    concat.append(&mut jref.to_vec());

    // Set clip penalties.  Note that yclip_suffix was set to zero.   This was
    // accompanied by a change to bio_edit in commit ccabb0dd1768738bdeee5b62458048d74f6dcfab,
    // and the entire commit is very flaky.  The alignment of these two sequences illustrates
    // the problem if one does not make the commit:
    // TTACTGTAAAGTCATGCTCTATGATAGTCGTGGTTCTGACTACTACTACGTTATGGACGTCTGGGGC
    // TTACTGTGCGAGACAGTATTACTATGATAGTAGTGGTTATTACTACATTACTACTACTACTACGGTATGGACGTCTGGGGC.

    // Make alignment.

    let mut scoring = Scoring::from_scores(gap_open, gap_extend, matchp, mismatch);
    scoring.xclip_prefix = MIN_SCORE;
    scoring.xclip_suffix = MIN_SCORE;
    scoring.yclip_prefix = MIN_SCORE;
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

    // Create zero-one vectors corresponding to indel-free aligned parts of the D gene; a zero
    // denotes a mismatch.

    let mut zos = Vec::<Vec<u8>>::new();
    {
        let mut rpos = 0;
        let mut zo = Vec::<u8>::new();
        for m in 0..al.operations.len() {
            match al.operations[m] {
                Match => {
                    if rpos >= vref.len() && rpos < vref.len() + dref.len() {
                        zo.push(1);
                    }
                    rpos += 1;
                }
                Subst => {
                    if rpos >= vref.len() && rpos < vref.len() + dref.len() {
                        zo.push(0);
                    }
                    rpos += 1;
                }
                Del => {
                    if zo.len() > 0 {
                        zos.push(zo.clone());
                        zo.clear();
                    }
                    rpos += 1;
                }
                Ins => {
                    if zo.len() > 0 {
                        zos.push(zo.clone());
                        zo.clear();
                    }
                }
                _ => {}
            };
        }
        if zo.len() > 0 {
            zos.push(zo.clone());
        }
    }

    // Compute a match bit score.

    let mut bits = 0.0_f64;
    for i in 0..zos.len() {
        let zo = &zos[i];
        for start in 0..zo.len() {
            for stop in start + 1..=zo.len() {
                let b = &zo[start..stop];
                let n = b.len();
                let mut k = 0; // number of mismatches
                for z in 0..n {
                    if b[z] == 0 {
                        k += 1;
                    }
                }

                // Let p be the probability that a random DNA sequence of length n will match a
                // given DNA sequence with ≤ k mismatches = sum{l=0..=k} (n choose l) * 3^l / 4^n.

                let mut sum = 0.0;
                let mut choose = 1.0;
                for l in 0..=k {
                    sum += choose;
                    choose *= (n - l) as f64;
                    choose /= (l + 1) as f64;
                    choose *= 3.0;
                }
                let p = sum / 4.0_f64.powi(n as i32);

                // Update bits.

                bits = bits.max(-p.log2());
            }
        }
    }

    // Add a constant times bits to the alignment score (null case handled differently).

    if dref.is_empty() {
        bits = 10.0;
    }
    let verbose = false;
    if verbose {
        use string_utils::*;
        println!("{} ==> {:.1}, ref = {}", drefname, bits, strme(&concat));
    }
    let score = al.score as f64 + 1.2 * bits;

    // Return the alignment and score.

    (al, score)
}
