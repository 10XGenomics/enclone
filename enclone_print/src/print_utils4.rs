// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use amino::*;
use debruijn::dna_string::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use equiv::EquivRel;
// use itertools::Itertools;
use std::cmp::max;
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Note confusing notation.  The object cdr3 contains pairs (String,usize) consisting of
// the cdr3_aa and the length of seq_del.

pub fn define_mat(
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    cdr3s: &Vec<Vec<(String, usize)>>,
    js: &Vec<usize>,
    od: &Vec<(Vec<usize>, usize, i32)>,
    info: &Vec<CloneInfo>,
) -> Vec<Vec<Option<usize>>> {
    // Form the flattened list of all CDR3_AAs.

    let nexacts = cdr3s.len();
    let mut all_cdr3s = Vec::<(Vec<u8>, usize)>::new();
    for j in 0..nexacts {
        for k in 0..cdr3s[j].len() {
            all_cdr3s.push((cdr3s[j][k].0.as_bytes().to_vec(), cdr3s[j][k].1));
        }
    }

    // Sort the CDR3s.

    unique_sort(&mut all_cdr3s);

    // Form an equivalence relation on the CDR3_AAs, requiring that they are "close enough":
    // 1. They have the same length and differ at no more than 4 positions.
    // 2. Each has a V..J sequence such that the two differ by no more than 50 positions.

    let mut ec: EquivRel = EquivRel::new(all_cdr3s.len() as i32);
    for m1 in 0..all_cdr3s.len() {
        for m2 in m1 + 1..all_cdr3s.len() {
            let (x1, x2) = (&all_cdr3s[m1].0, &all_cdr3s[m2].0);
            let (y1, y2) = (all_cdr3s[m1].1, all_cdr3s[m2].1);
            if x1.len() == x2.len() && y1 == y2 {
                let mut diffs = 0;
                for u in 0..x1.len() {
                    if x1[u] != x2[u] {
                        diffs += 1;
                    }
                }
                if diffs <= 4 {
                    'outer: for l1 in 0..od.len() {
                        let y1: &CloneInfo = &info[od[l1].2 as usize];
                        for u1 in 0..y1.cdr3_aa.len() {
                            if y1.cdr3_aa[u1] == strme(&x1).after(":") {
                                for l2 in 0..od.len() {
                                    let y2: &CloneInfo = &info[od[l2].2 as usize];
                                    for u2 in 0..y2.cdr3_aa.len() {
                                        if y2.cdr3_aa[u2] == strme(&x2).after(":") {
                                            if y1.tigs[u1].len() == y2.tigs[u2].len() {
                                                // Could be we're spending a lot of time
                                                // finding diffs.

                                                let mut vj_diffs = 0;
                                                if !y1.has_del[u1] && !y2.has_del[u2] {
                                                    vj_diffs = ndiffs(&y1.tigsp[u1], &y2.tigsp[u2]);
                                                } else {
                                                    for j in 0..y1.tigs[u1].len() {
                                                        if y1.tigs[u1][j] != y2.tigs[u2][j] {
                                                            vj_diffs += 1;
                                                        }
                                                    }
                                                }
                                                if vj_diffs <= ctl.heur.max_diffs {
                                                    ec.join(m1 as i32, m2 as i32);
                                                    break 'outer;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Create rmap, that sends
    // (index into exact subclonotypes for this clonotype,
    //  index into chains for one of these exact subclonotypes)
    // to an index into the orbit reps for the CDR3s.

    let mut r = Vec::<i32>::new();
    ec.orbit_reps(&mut r);
    let mut rpos = HashMap::<(usize, usize), usize>::new();
    for u in 0..nexacts {
        let x = &cdr3s[u];
        for (iy, y) in x.iter().enumerate() {
            let p = bin_position(&all_cdr3s, &(y.0.as_bytes().to_vec(), y.1));
            let c = ec.class_id(p);
            let q = bin_position(&r, &c) as usize;
            rpos.insert((u, iy), q);
        }
    }

    // Find the maximum multiplicity of each orbit, and the number of columns.

    let mut mm = vec![0; r.len()];
    for (u, x) in cdr3s.iter().enumerate() {
        let mut mm0 = vec![0; r.len()];
        for iy in 0..x.len() {
            let q = rpos[&(u, iy)];
            mm0[q] += 1;
        }
        for i in 0..r.len() {
            mm[i] = max(mm[i], mm0[i]);
        }
    }
    let cols = mm.iter().sum();

    // Define a matrix mat[col][ex] which is the column of the exact subclonotype ex
    // corresponding to the given column col of the clonotype, which may or may not be
    // defined.
    // ◼ This should be propagated so we don't compute something equivalent
    //   over and over.

    let mut mat = vec![vec![None; nexacts]; cols];
    for cx in 0..cols {
        // for every column
        'exact: for u in 0..nexacts {
            // for every exact subclonotype
            let clonotype_id = od[js[u]].1;
            let ex = &exact_clonotypes[clonotype_id];
            let x = &cdr3s[u];
            let mut mm0 = vec![0; r.len()];
            // for every chain in the exact subclonotype:
            for (iy, y) in x.iter().enumerate() {
                let q = rpos[&(u, iy)];
                let mut col = mm0[q];
                for j in 0..q {
                    col += mm[j];
                }
                mm0[q] += 1;
                if col != cx {
                    continue;
                }
                for m in 0..ex.share.len() {
                    if ex.share[m].cdr3_aa == y.0.after(":") && ex.share[m].seq_del.len() == y.1 {
                        mat[cx][u] = Some(m);
                        continue 'exact;
                    }
                }
            }
        }
    }
    mat
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Define amino acid positions to show.

pub fn build_show_aa(
    ctl: &EncloneControl,
    rsi: &ColInfo,
    vars_amino: &Vec<Vec<usize>>,
    shares_amino: &Vec<Vec<usize>>,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    _exacts: &Vec<usize>,
    _exact_clonotypes: &Vec<ExactClonotype>,
) -> Vec<Vec<usize>> {
    let cols = rsi.vids.len();
    let mut show_aa = vec![Vec::<usize>::new(); cols];
    for cx in 0..cols {
        for x in ctl.clono_print_opt.amino.iter() {
            if x.contains('-') {
                let (start, stop) = (x.before("-").force_usize(), x.after("-").force_usize());
                for p in start..=stop {
                    show_aa[cx].push(p);
                }
            }
        }
        if ctl.clono_print_opt.amino.contains(&"cdr1".to_string()) {
            if rsi.cdr1_starts[cx].is_some()
                && rsi.fr2_starts[cx].is_some()
                && rsi.cdr1_starts[cx].unwrap() <= rsi.fr2_starts[cx].unwrap()
            {
                for j in (rsi.cdr1_starts[cx].unwrap()..rsi.fr2_starts[cx].unwrap()).step_by(3) {
                    show_aa[cx].push(j / 3);
                }
            }
        }
        if ctl.clono_print_opt.amino.contains(&"cdr2".to_string()) {
            if rsi.cdr2_starts[cx].is_some()
                && rsi.fr3_starts[cx].is_some()
                && rsi.cdr2_starts[cx].unwrap() <= rsi.fr3_starts[cx].unwrap()
            {
                for j in (rsi.cdr2_starts[cx].unwrap()..rsi.fr3_starts[cx].unwrap()).step_by(3) {
                    show_aa[cx].push(j / 3);
                }
            }
        }
        if ctl.clono_print_opt.amino.contains(&"cdr3".to_string()) {
            for j in 0..rsi.cdr3_lens[cx] {
                let p = rsi.cdr3_starts[cx] / 3 + j;
                show_aa[cx].push(p);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"fwr1".to_string()) {
            if rsi.cdr1_starts[cx].is_some() && rsi.fr1_starts[cx] <= rsi.cdr1_starts[cx].unwrap() {
                for j in (rsi.fr1_starts[cx]..rsi.cdr1_starts[cx].unwrap()).step_by(3) {
                    show_aa[cx].push(j / 3);
                }
            }
        }
        if ctl.clono_print_opt.amino.contains(&"fwr2".to_string()) {
            if rsi.fr2_starts[cx].is_some()
                && rsi.cdr2_starts[cx].is_some()
                && rsi.fr2_starts[cx].unwrap() <= rsi.cdr2_starts[cx].unwrap()
            {
                for j in (rsi.fr2_starts[cx].unwrap()..rsi.cdr2_starts[cx].unwrap()).step_by(3) {
                    show_aa[cx].push(j / 3);
                }
            }
        }
        if ctl.clono_print_opt.amino.contains(&"fwr3".to_string()) {
            if rsi.fr3_starts[cx].is_some() && rsi.fr3_starts[cx].unwrap() <= rsi.cdr3_starts[cx] {
                for j in (rsi.fr3_starts[cx].unwrap()..rsi.cdr3_starts[cx]).step_by(3) {
                    show_aa[cx].push(j / 3);
                }
            }
        }
        if ctl.clono_print_opt.amino.contains(&"fwr4".to_string()) {
            for j in (rsi.cdr3_starts[cx] + 3 * rsi.cdr3_lens[cx]..rsi.seq_lens[cx]).step_by(3) {
                show_aa[cx].push(j / 3);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"var".to_string()) {
            for j in 0..vars_amino[cx].len() {
                let p = vars_amino[cx][j];
                show_aa[cx].push(p / 3);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"share".to_string()) {
            for j in 0..shares_amino[cx].len() {
                let p = shares_amino[cx][j];
                show_aa[cx].push(p / 3);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"donor".to_string()) {
            let vseq1 = refdata.refs[rsi.vids[cx]].to_ascii_vec();
            let jseq1 = refdata.refs[rsi.jids[cx]].to_ascii_vec();
            let vseq2: Vec<u8>;
            if rsi.vpids[cx].is_some() {
                vseq2 = dref[rsi.vpids[cx].unwrap()].nt_sequence.clone();
            } else {
                vseq2 = vseq1.clone();
            }
            let jseq2 = &jseq1;
            let vlen = vseq2.len() - ctl.heur.ref_v_trim;
            let jlen = jseq2.len() - ctl.heur.ref_j_trim;
            // This test must be here for a reason, but we encountered examples where it was
            // triggered, and there was nothing obviously wrong.  In the event that an internal
            // error is encountered elsewhere in the code, we might wish to turn this back on.
            /*
            let gap = rsi.seq_lens[cx] as isize - vlen as isize - jlen as isize;
            if gap < 0 {
                let mut bcs = Vec::<String>::new();
                for u in 0..exacts.len() {
                    let ex = &exact_clonotypes[exacts[u]];
                    for i in 0..ex.clones.len() {
                        bcs.push(ex.clones[i][0].barcode.clone());
                    }
                }
                bcs.sort();
                panic!(
                    "Something is wrong because gap is {}, which is negative.\n\
                    This is happening for chain {} of {} of the clonotype with \
                    these barcodes:\n{}\nand with first V..J sequence\n{}.",
                    gap,
                    cx + 1,
                    cols,
                    bcs.iter().format(","),
                    strme(&exact_clonotypes[exacts[0]].share[cx].seq)
                );
            }
            */
            for j in 0..vlen {
                if j < vseq1.len() && vseq1[j] != vseq2[j] {
                    show_aa[cx].push(j / 3);
                }
            }
            for j in ctl.heur.ref_j_trim..jlen {
                if jseq1[j] != jseq2[j] {
                    let pos = rsi.seq_lens[cx] - jlen + j;
                    show_aa[cx].push(pos / 3);
                }
            }
        }
        if ctl.clono_print_opt.amino.contains(&"donorn".to_string()) {
            let vseq1 = refdata.refs[rsi.vids[cx]].to_ascii_vec();
            let jseq1 = refdata.refs[rsi.jids[cx]].to_ascii_vec();
            let vseq2: Vec<u8>;
            if rsi.vpids[cx].is_some() {
                vseq2 = dref[rsi.vpids[cx].unwrap()].nt_sequence.clone();
            } else {
                vseq2 = vseq1.clone();
            }
            let jseq2 = &jseq1;
            let vlen = vseq2.len() - ctl.heur.ref_v_trim;
            let jlen = jseq2.len() - ctl.heur.ref_j_trim;
            /*
            let gap = rsi.seq_lens[cx] as isize - vlen as isize - jlen as isize;
            assert!(gap >= 0);
            */
            for j in 0..vlen {
                if j < vseq1.len() && vseq1[j] != vseq2[j] {
                    let n = 3 * (j / 3);
                    if n + 3 <= vseq1.len()
                        && codon_to_aa(&vseq1[n..n + 3]) != codon_to_aa(&vseq2[n..n + 3])
                    {
                        show_aa[cx].push(j / 3);
                    }
                }
            }
            for j in ctl.heur.ref_j_trim..jlen - 1 {
                if jseq1[j] != jseq2[j] {
                    let mut p = j;
                    let d = jlen - 1 - j;
                    if d % 3 == 0 {
                    } else if d % 3 == 1 {
                        p -= 2;
                    } else if d % 3 == 2 {
                        p -= 1;
                    }
                    if p + 3 <= jlen
                        && codon_to_aa(&jseq1[p..p + 3]) != codon_to_aa(&jseq2[p..p + 3])
                    {
                        let pos = rsi.seq_lens[cx] - jlen + j;
                        show_aa[cx].push(pos / 3);
                    }
                }
            }
        }
        unique_sort(&mut show_aa[cx]);

        // Remove an amino acid position that is too high.  The way we would expect
        // to get this is that we started from the last base position on the J segment,
        // which is part of a codon whose other two bases lie in a C segment.

        if !show_aa[cx].is_empty() {
            let p = show_aa[cx][show_aa[cx].len() - 1];
            if 3 * p + 3 > rsi.seq_del_lens[cx] {
                show_aa[cx].pop();
            }
        }
    }
    show_aa
}
