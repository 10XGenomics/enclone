// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use amino::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use equiv::EquivRel;
use std::cmp::max;
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn define_mat(
    _ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    exacts: &Vec<usize>,
    _cdr3s: &Vec<Vec<(String, usize)>>,
    _js: &Vec<usize>,
    _ks: &Vec<usize>,
    od: &Vec<(Vec<usize>, usize, i32)>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
) -> Vec<Vec<Option<usize>>> {

    /*
    // XXX:
    let mut tracking = false;
    for u in 0..exacts.len() {
        let ex = &exact_clonotypes[exacts[u]];
        for m in 0..ex.share.len() {
            if ex.share[m].cdr3_aa == "CARHVGGKYSSGWYEFDYW" {
                tracking = true;
            }
        }
    }
    if tracking {
        eprintln!("\ntracking");
        eprintln!("there are {} exact subclonotypes", exacts.len());
        let ex = &exact_clonotypes[exacts[0]];
        eprintln!("the first has {} chains", ex.share.len());
    }
    */

    // Define map to indices into exacts.

    let nexacts = exacts.len();
    let mut to_exacts = HashMap::<usize, usize>::new();
    for u in 0..nexacts {
        to_exacts.insert(exacts[u], u);
    }

    // Get the info indices corresponding to this clonotype.

    let mut infos = Vec::<usize>::new();
    for i in 0..od.len() {
        let x = od[i].2 as usize;
        if to_exacts.contains_key(&info[x].clonotype_index) {
            infos.push(x);
        }
    }
    infos.sort();
    // if tracking { eprintln!("there are {} info entries", infos.len()); } // XXXXXXXXXXXXXXXXXXXXXXX

    // Form the set of all chains that appear in an exact subclonotypes of this clonotype, and
    // also track the V..J sequences for the chains.

    let mut chains = Vec::<(usize, usize)>::new();
    let mut seq_chains = Vec::<(Vec<u8>, usize, usize)>::new();
    for u in 0..nexacts {
        let ex = &exact_clonotypes[exacts[u]];
        for m in 0..ex.share.len() {
            chains.push((u, m));
            seq_chains.push((ex.share[m].seq.clone(), u, m));
        }
    }
    seq_chains.sort();

    // Define an equivalence relation on the chains, introducing connections defined by the
    // raw joins.  Also join where there are identical V..J sequences.

    let mut e = EquivRel::new(chains.len() as i32);
    for i1 in 0..infos.len() {
        let j1 = infos[i1];
        let u1 = info[j1].clonotype_index;
        let v1 = to_exacts[&u1];
        let m1s = &info[j1].exact_cols;
        for x in raw_joins[j1].iter() {
            let i2 = bin_position(&infos, &*x);
            if i2 >= 0 {
                let i2 = i2 as usize;
                let j2 = infos[i2];
                let u2 = info[j2].clonotype_index;
                let v2 = to_exacts[&u2];
                let m2s = &info[j2].exact_cols;
                for j in 0..2 {
                    let z1 = bin_position(&chains, &(v1, m1s[j]));
                    let z2 = bin_position(&chains, &(v2, m2s[j]));
                    e.join(z1, z2);
                }
            }
        }
    }
    let mut i = 0;
    while i < seq_chains.len() {
        let j = next_diff1_3(&seq_chains, i as i32) as usize;
        for k in i+1..j {
            let (x1, x2) = (&seq_chains[i], &seq_chains[k]);
            let z1 = bin_position(&chains, &(x1.1, x1.2));
            let z2 = bin_position(&chains, &(x2.1, x2.2));
            e.join(z1, z2);
        }
        i = j;
    }

    // Get representatives for the chain orbits.

    let mut r = Vec::<i32>::new();
    e.orbit_reps(&mut r);
    // if tracking { eprintln!("there are {} orbits", r.len()); } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Reorder the chains.  This is done to get the heavy chains before the light chains and also
    // to mimic the behavior of the previous version of this algorithm, to minimiize churn.  Then
    // update the representatives.

    let mut chainsp = Vec::<(String, usize, usize, usize)>::new();
    for u in 0..nexacts {
        let ex = &exact_clonotypes[exacts[u]];
        for m in 0..ex.share.len() {
            let mut c = ex.share[m].chain_type.clone();
            if c.starts_with("TRB") {
                c = c.replacen("TRB", "TRX", 1);
            } else if c.starts_with("TRA") {
                c = c.replacen("TRA", "TRY", 1);
            }
            chainsp.push((format!("{}:{}", c, ex.share[m].cdr3_aa), ex.share[m].seq.len(), u, m));
        }
    }
    chainsp.sort();
    let mut chainso = Vec::<(usize, usize)>::new();
    let mut chainsox = Vec::<(usize, usize, usize)>::new();
    for i in 0..chainsp.len() {
        chainso.push((chainsp[i].2, chainsp[i].3));
        chainsox.push((chainsp[i].2, chainsp[i].3, i));
    }
    chainsox.sort();
    for i in 0..r.len() {
        r[i] = chainsox[r[i] as usize].2 as i32;
    }
    r.sort();

    // Create rmap, that sends
    // (index into exact subclonotypes for this clonotype,
    //  index into chains for one of these exact subclonotypes)
    // to an index into the orbit reps for the chains.

    let mut rpos = HashMap::<(usize, usize), usize>::new();
    for i in 0..chains.len() {
        let c = e.class_id(i as i32);
        let f = chainsox[c as usize].2 as i32;
        let q = bin_position(&r, &f) as usize;
        rpos.insert((chains[i].0, chains[i].1), q);
    }

    // Find the maximum multiplicity of each orbit, and the number of columns.

    let mut mm = vec![0; r.len()];
    for u in 0..nexacts {
        let ex = &exact_clonotypes[exacts[u]];
        let mut mm0 = vec![0; r.len()];
        for m in 0..ex.share.len() {
            mm0[rpos[&(u, m)]] += 1;
        }
        for i in 0..r.len() {
            mm[i] = max(mm[i], mm0[i]);
        }
    }
    let cols = mm.iter().sum();

    // Define a matrix mat[col][ex] which is the column of the exact subclonotype ex corresponding
    // to the given column col of the clonotype, which may or may not be defined.

    let mut mat = vec![vec![None; nexacts]; cols];
    for cx in 0..cols {
        // for every column
        'exact: for u in 0..nexacts {
            // for every exact subclonotype
            let ex = &exact_clonotypes[exacts[u]];
            let mut mm0 = vec![0; r.len()];
            for m in 0..ex.share.len() {
                // for every chain in the exact subclonotype:
                let q = rpos[&(u, m)];
                let mut col = mm0[q];
                for j in 0..q {
                    col += mm[j];
                }
                mm0[q] += 1;
                if col == cx {
                    mat[cx][u] = Some(m);
                    continue 'exact;
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
