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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Create zero-one vectors corresponding to indel-free aligned parts of the D gene; a zero denotes a mismatch.

fn zero_one(al: &bio_edit::alignment::Alignment, vref: &[u8], dref: &[u8]) -> Vec<Vec<u8>> {
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
    zos
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Compute a match bit score.

fn match_bit_score(zos: &Vec<Vec<u8>>) -> f64 {
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

pub fn align_to_vdj_ref(
    seq: &[u8],
    vref: &[u8],
    dref: &[u8],
    jref: &[u8],
    drefname: &str, // useful for debugging
) -> (Vec<bio::alignment::AlignmentOperation>, f64) {
    // Define penalties.

    let matchp = 2 as i32;
    let mismatch = -2 as i32;
    let gap_open = -12 as i32;
    let gap_extend = -2 as i32;
    let gap_open_at_boundary = -4 as i32;
    let gap_extend_at_boundary = -1 as i32;

    // Define scoring function.  It does not appear that the aligner is scoring in exactly the intended 
    // fashion in all cases.  This is likely more of an issue for alv and alj.  The problem has to do with the 
    // interpretation of being at the boundary.  This is still somewhat of a problem since although we're 
    // rescoring to "fix" the problem, the aligner might have chosen a suboptimal placement in the first place.

    let rescore = |ops: &Vec<bio_edit::alignment::AlignmentOperation>| -> i32 {
        let mut score = 0 as i32;
        let mut i = 0;
        let mut rpos = 0;
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
                if rpos == vref.len() || rpos == vref.len() + dref.len() {
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
                if rpos + j - i == vref.len() || rpos == vref.len() + dref.len() {
                    score += gap_open_at_boundary + (j - i - 1) as i32 * gap_extend_at_boundary;
                } else {
                    score += gap_open + (j - i - 1) as i32 * gap_extend;
                }
                rpos += j - i;
                i = j;
            }
        }
        score
    };

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
    // denotes a mismatch.  Then compute a match bit score.

    let zos = zero_one(&al, &vref, &dref);
    let mut bits = match_bit_score(&zos);

    // Possibly emit verbose logging.

    let verbose = false;
    if verbose {
        use string_utils::*;
        println!("\n{} ==> {:.1}", drefname, bits);
        println!("seq = {}", strme(&seq));
        println!("ref = {}", strme(&concat));
        use itertools::Itertools;
        for zo in zos.iter() {
            print!("{}", zo.iter().format(""));
        }
        println!("");
        println!("ops = {:?}", al.operations.iter().format(","));
    }

    // In the non-null case, where the D segment is placed without indels, see if the placement 
    // seems like it might be wrong.

    if zos.len() == 1 {

        // Set pos to the start of D on seq.

        let (mut pos, mut rpos) = (0, 0);
        for m in 0..al.operations.len() {
            match al.operations[m] {
                Match => {
                    if rpos == vref.len() {
                        break;
                    }
                    pos += 1;
                    rpos += 1;
                }
                Subst => {
                    if rpos == vref.len() {
                        break;
                    }
                    pos += 1;
                    rpos += 1;
                }
                Del => {
                    rpos += 1;
                }
                Ins => {
                    pos += 1;
                }
                _ => {}
            };
        }

        // Look for alternate starts near the given start.
        // This is in the range pos-DWOBBLE..=pos+DWOBBLE, but we are a bit careful in case the
        //  reference is badly busted.

        const D_WOBBLE: usize = 4;
        let mut dstart = 0;
        if pos >= D_WOBBLE {
            dstart = pos - D_WOBBLE;
        }
        let mut dstop = pos + D_WOBBLE;
        if dstop + dref.len() > seq.len() {
            if seq.len() >= dref.len() {
                dstop = seq.len() - dref.len();
            } else {
                dstop = dstart - 1;
            }
        }
        for dpos in dstart..=dstop {
            if dpos != pos {

                // Compute bits for the alternate start.

                let mut zos = vec![vec![1; dref.len()]];
                for i in 0..dref.len() {
                    if dref[i] != seq[dpos + i] {
                        zos[0][i] = 0;
                    }
                }
                let bits_alt = match_bit_score(&zos);

                // Does the alt placement appear to be better?

                if bits_alt > bits {

                    // Realign, forcing the D segment to be in the given place, without indels.

                    let mut gap_open_fn = vec![gap_open; vref.len() + 1];
                    gap_open_fn[vref.len()] = gap_open_at_boundary;
                    let mut gap_extend_fn = vec![gap_extend; vref.len() + 1];
                    gap_extend_fn[vref.len()] = gap_extend_at_boundary;
                    let mut alv = aligner.custom_with_gap_fns(
                        &seq[0..dpos], 
                        &vref, 
                        &gap_open_fn, 
                        &gap_extend_fn
                    );
                    alv.mode = AlignmentMode::Semiglobal;
                    let mut gap_open_fn = vec![gap_open; jref.len() + 1];
                    gap_open_fn[1] = gap_open_at_boundary;
                    let mut gap_extend_fn = vec![gap_extend; jref.len() + 1];
                    gap_extend_fn[1] = gap_extend_at_boundary;
                    let mut alj = aligner.custom_with_gap_fns(
                        &seq[dpos + dref.len()..], 
                        &jref, 
                        &gap_open_fn, 
                        &gap_extend_fn
                    );
                    alj.mode = AlignmentMode::Semiglobal;

                    // Compute the new score and see if its better.

                    let mut ops_new = Vec::<bio_edit::alignment::AlignmentOperation>::new();
                    ops_new.append(&mut alv.operations.clone());
                    for i in 0..zos[0].len() {
                        if zos[0][i] == 0 {
                            ops_new.push(Subst);
                        } else {
                            ops_new.push(Match);
                        }
                    }
                    ops_new.append(&mut alj.operations.clone());
                    let full_score = al.score as f64 + 1.2 * bits;
                    let score_alt = rescore(&ops_new);
                    let full_score_alt = score_alt as f64 + 1.2 * bits_alt;
                    if full_score_alt > full_score {
                        return (ops_new, full_score_alt);
                    }
                }
            }
        }
    }

    // Add a constant times bits to the alignment score (null case handled differently).

    if dref.is_empty() {
        bits = 5.5;
    }
    let full_score = rescore(&al.operations) as f64 + 1.2 * bits;

    // Return the alignment and score.

    (al.operations, full_score)
}
