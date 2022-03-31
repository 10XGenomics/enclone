// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Align a sequence to a concatenated V(D)J reference, encouraging insertions exactly at
// junction points and deletions that bridge junction points.
//
// ASSESSMENT
//
// enclone BI=1-4,9 BUILT_IN SUBSET_JSON=subset/outs/all_contig_annotations.json
//         MIN_EXACTS=2 MAX_EXACTS=10 NOPRINT
//
// (slow)
//
// enclone BCR=subset GVARS=d_inconsistent_%,d_inconsistent_n NOPRINT
//
// (fast)
//
// inconsistency rate from this = 13.33%
// sample size = 53,373
//
// If you mess with this, you can test your changes with "cargo t test_enclone_d" and
// "merge_html BUILD" and then manually examine the D gene page.  Note carefully that we do not
// want to worsen the placement of indels.  Also run the above big test.

use bio_edit::alignment::pairwise::{Aligner, Scoring, MIN_SCORE};
use bio_edit::alignment::AlignmentMode;
use bio_edit::alignment::AlignmentOperation::*;
use string_utils::strme;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Create zero-one vectors corresponding to indel-free aligned parts of the D gene; a zero denotes
// a mismatch.

pub fn zero_one(
    ops: &Vec<bio_edit::alignment::AlignmentOperation>,
    start: usize,
    stop: usize,
) -> Vec<Vec<u8>> {
    let mut zos = Vec::<Vec<u8>>::new();
    {
        let mut rpos = 0;
        let mut zo = Vec::<u8>::new();
        for m in 0..ops.len() {
            match ops[m] {
                Match => {
                    if rpos >= start && rpos < stop {
                        zo.push(1);
                    }
                    rpos += 1;
                }
                Subst => {
                    if rpos >= start && rpos < stop {
                        zo.push(0);
                    }
                    rpos += 1;
                }
                Del => {
                    if !zo.is_empty() {
                        zos.push(zo.clone());
                        zo.clear();
                    }
                    rpos += 1;
                }
                Ins => {
                    if !zo.is_empty() {
                        zos.push(zo.clone());
                        zo.clear();
                    }
                }
                _ => {}
            };
        }
        if !zo.is_empty() {
            zos.push(zo.clone());
        }
    }
    zos
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Compute a match bit score.

pub fn match_bit_score(zos: &Vec<Vec<u8>>) -> f64 {
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
    bits
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This is the same as bio::Alignment::cigar(&self, hard_clip: bool), but works on a vector
// of alignment operations, and assumes hard_clip is false.  Code copied from bio source.

pub fn cigar(
    ops: &Vec<bio_edit::alignment::AlignmentOperation>,
    xstart: usize,
    xend: usize,
    xlen: usize,
) -> String {
    use bio_edit::alignment::AlignmentOperation;
    let clip_str = "S";
    let add_op = |op: AlignmentOperation, k, cigar: &mut String| match op {
        AlignmentOperation::Match => cigar.push_str(&format!("{}{}", k, "=")),
        AlignmentOperation::Subst => cigar.push_str(&format!("{}{}", k, "X")),
        AlignmentOperation::Del => cigar.push_str(&format!("{}{}", k, "D")),
        AlignmentOperation::Ins => cigar.push_str(&format!("{}{}", k, "I")),
        _ => {}
    };

    let mut cigar = "".to_owned();
    if ops.is_empty() {
        return cigar;
    }

    let mut last = ops[0];
    if xstart > 0 {
        cigar.push_str(&format!("{}{}", xstart, clip_str))
    }
    let mut k = 1;
    for &op in ops[1..].iter() {
        if op == last {
            k += 1;
        } else {
            add_op(last, k, &mut cigar);
            k = 1;
        }
        last = op;
    }
    add_op(last, k, &mut cigar);
    if xlen > xend {
        cigar.push_str(&format!("{}{}", xlen - xend, clip_str))
    }
    cigar
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn align_to_vdj_ref(
    seq: &[u8],
    vref: &[u8],
    dref: &[u8],
    d2ref: &[u8],
    jref: &[u8],
    drefname: &str, // useful for debugging
    left: bool,
    jscore_match: i32,
    jscore_mismatch: i32,
    jscore_gap_open: i32,
    jscore_gap_extend: i32,
    jscore_bits_multiplier: f64,
) -> (Vec<bio_edit::alignment::AlignmentOperation>, f64) {
    // Define penalties.

    let matchp = jscore_match;
    let mismatch = jscore_mismatch;
    let gap_open = jscore_gap_open;
    let gap_extend = jscore_gap_extend;
    let gap_open_at_boundary = -40_i32;
    let gap_extend_at_boundary = -10_i32;
    let del_gap_extend_at_boundary = -20_i32;
    let align_div = 10.0;
    let bits_multiplier = jscore_bits_multiplier;
    const MIN_BITS_FOR_D2: f64 = 14.0;
    const D2_PENALTY: f64 = -15.0;

    // Define scoring function.  It does not appear that the aligner is scoring in exactly the
    // intended fashion in all cases.  This is likely more of an issue for alv and alj.  The
    // problem has to do with the interpretation of being at the boundary.  This is still somewhat
    // of a problem since although we're rescoring to "fix" the problem, the aligner might have
    // chosen a suboptimal placement in the first place.
    //
    // Note handling of deletions that bridge boundaries.

    let rescore = |ops: &Vec<bio_edit::alignment::AlignmentOperation>| -> f64 {
        let mut score = 0_i32;
        let mut i = 0;
        let mut rpos = 0;
        let b1 = vref.len();
        let b2 = vref.len() + dref.len();
        let b3 = vref.len() + dref.len() + d2ref.len();
        while i < ops.len() {
            if ops[i] == Match {
                rpos += 1;
                score += matchp;
                i += 1;
            } else if ops[i] == Subst {
                rpos += 1;
                score += mismatch;
                i += 1;
            } else if ops[i] == Ins {
                let mut j = i + 1;
                while j < ops.len() && ops[j] == Ins {
                    j += 1;
                }
                if (rpos == vref.len() + dref.len() + d2ref.len())
                    || (rpos == vref.len() || rpos == vref.len() + dref.len())
                {
                    score += gap_open_at_boundary + (j - i - 1) as i32 * gap_extend_at_boundary;
                } else {
                    score += gap_open + (j - i - 1) as i32 * gap_extend;
                }
                i = j;
            } else if ops[i] == Del {
                let mut j = i + 1;
                while j < ops.len() && ops[j] == Del {
                    j += 1;
                }
                if (rpos <= b1 && rpos + j - i >= b1)
                    || (rpos <= b2 && rpos + j - i >= b2)
                    || (rpos <= b3 && rpos + j - i >= b3)
                {
                    score += gap_open_at_boundary + (j - i - 1) as i32 * del_gap_extend_at_boundary;
                } else {
                    score += gap_open + (j - i - 1) as i32 * gap_extend;
                }
                rpos += j - i;
                i = j;
            }
        }
        score as f64 / align_div
    };

    // Build concatenation.

    let mut concat = Vec::<u8>::new();
    concat.append(&mut vref.to_vec());
    concat.append(&mut dref.to_vec());
    concat.append(&mut d2ref.to_vec());
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
        if j as usize == vref.len()
            || j as usize == vref.len() + dref.len()
            || j as usize == vref.len() + dref.len() + d2ref.len()
        {
            gap_open_fn[j] = gap_open_at_boundary;
        } else {
            gap_open_fn[j] = gap_open;
        }
    }
    let mut gap_extend_fn = vec![0_i32; concat.len() + 1];
    for j in 1..=concat.len() {
        if j as usize == vref.len()
            || j as usize == vref.len() + dref.len()
            || j as usize == vref.len() + dref.len() + d2ref.len()
        {
            gap_extend_fn[j] = gap_extend_at_boundary;
        } else {
            gap_extend_fn[j] = gap_extend;
        }
    }
    let mut al = aligner.custom_with_gap_fns(seq, &concat, &gap_open_fn, &gap_extend_fn);
    al.mode = AlignmentMode::Semiglobal;
    let mut ops = al.operations;

    // Fix alignments.

    let mut edits = Vec::<(usize, bio_edit::alignment::AlignmentOperation)>::new();
    let mut i = 0;
    let mut pos = 0;
    let mut rpos = 0;
    let b1 = vref.len();
    let b2 = vref.len() + dref.len();
    let b3 = vref.len() + dref.len() + d2ref.len();
    let mut edited = vec![false; ops.len()];
    while i < ops.len() {
        if ops[i] == Match || ops[i] == Subst {
            pos += 1;
            rpos += 1;
            i += 1;
        } else if ops[i] == Ins {
            let mut j = i + 1;
            while j < ops.len() && ops[j] == Ins {
                j += 1;
            }
            pos += j - i;
            i = j;
        } else if ops[i] == Del {
            let mut j = i + 1;
            while j < ops.len() && ops[j] == Del {
                j += 1;
            }
            let k = j - i;
            for bi in [b1, b2, b3].iter() {
                let bi = *bi;

                // Maybe can shift right one.

                if rpos < bi
                    && rpos + k >= bi
                    && j < ops.len()
                    && pos < seq.len()
                    && rpos < concat.len()
                    && ops[j] == Subst
                    && seq[pos] == concat[rpos]
                {
                    edits.push((i, Match));
                    edits.push((j, Del));
                    edited[i] = true;
                    edited[j] = true;
                    break;

                // Maybe can shift left one.
                } else if rpos + k > bi
                    && i > 0
                    && ops[i - 1] == Subst
                    && seq[pos - 1] == concat[rpos + k - 1]
                    && !edited[i - 1]
                {
                    edits.push((i - 1, Del));
                    edits.push((j - 1, Match));
                    edited[i - 1] = true;
                    edited[j - 1] = true;
                    break;
                }
            }
            rpos += j - i;
            i = j;
        }
    }
    for x in edits.iter() {
        ops[x.0] = x.1;
    }

    // Create zero-one vectors corresponding to indel-free aligned parts of the D gene; a zero
    // denotes a mismatch.  Then compute a match bit score.

    let zos1 = zero_one(&ops, vref.len(), vref.len() + dref.len());
    let zos2 = zero_one(
        &ops,
        vref.len() + dref.len(),
        vref.len() + dref.len() + d2ref.len(),
    );
    let bits1 = match_bit_score(&zos1);
    let bits2 = match_bit_score(&zos2);
    let mut bits = bits1.max(bits2);
    if !d2ref.is_empty() && bits1.min(bits2) < MIN_BITS_FOR_D2 {
        bits = 0.0;
    }

    // Possibly emit verbose logging.

    let verbose = false;
    if verbose {
        let full_score = rescore(&ops) + bits_multiplier * bits;
        println!(
            "\n{} ==> score = {:.1}, bits = {:.1}, full_score = {:.1}",
            drefname,
            rescore(&ops),
            bits,
            full_score,
        );
        println!("seq = {}", strme(seq));
        println!("ref = {}", strme(&concat));
        use itertools::Itertools;
        for zo in zos1.iter() {
            print!("{}", zo.iter().format(""));
        }
        println!();
        println!("ops = {:?}", ops.iter().format(","));
    }

    // Add a constant times bits to the alignment score (null case handled differently).
    //
    // Note that we do not allow the null case if there is an insertion in the alignment.  In an
    // an earlier version, we allowed many more null cases, and we believe that this was distorting
    // our inconsistency scoring.  This is because calling null makes it easier to be consistent.

    if left && dref.is_empty() {
        if !ops.contains(&Ins) {
            bits = 10.0;
        } else {
            bits = -1000.0;
        }
    }
    let mut full_score = rescore(&ops) + bits_multiplier * bits;
    if drefname.contains(':') {
        full_score += D2_PENALTY;
    }

    // Return the alignment and score.

    (ops, full_score)
}
