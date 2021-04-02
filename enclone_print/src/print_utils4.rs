// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::print_utils1::*;
use amino::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use equiv::EquivRel;
use itertools::Itertools;
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn compute_some_stats(
    ctl: &EncloneControl,
    lvars: &Vec<String>,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    n_vdj_gex: &Vec<usize>,
    cred: &mut Vec<Vec<String>>,
    pe: &mut Vec<Vec<String>>,
    ppe: &mut Vec<Vec<String>>,
    npe: &mut Vec<Vec<String>>,
) {
    let nexacts = exacts.len();

    // Compute "cred" stats (credibility/# of neighboring cells that are also B cells).

    *cred = vec![Vec::<String>::new(); lvars.len()];
    for k in 0..lvars.len() {
        if lvars[k] == "cred".to_string() {
            for u in 0..nexacts {
                let clonotype_id = exacts[u];
                let ex = &exact_clonotypes[clonotype_id];
                for l in 0..ex.clones.len() {
                    let bc = &ex.clones[l][0].barcode;
                    let li = ex.clones[l][0].dataset_index;
                    if gex_info.pca[li].contains_key(&bc.clone()) {
                        let mut creds = 0;
                        let mut z = Vec::<(f64, String)>::new();
                        let x = &gex_info.pca[li][&bc.clone()];
                        for y in gex_info.pca[li].iter() {
                            let mut dist2 = 0.0;
                            for m in 0..x.len() {
                                dist2 += (y.1[m] - x[m]) * (y.1[m] - x[m]);
                            }
                            z.push((dist2, y.0.clone()));
                        }
                        z.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        let top = n_vdj_gex[li];
                        for i in 0..top {
                            if bin_member(&vdj_cells[li], &z[i].1) {
                                creds += 1;
                            }
                        }
                        let pc = 100.0 * creds as f64 / top as f64;
                        cred[k].push(format!("{:.1}", pc));
                    } else {
                        cred[k].push("".to_string());
                    }
                }
            }
        }
    }

    // Compute pe (PCA distance).

    *pe = vec![Vec::<String>::new(); lvars.len()];
    for k in 0..lvars.len() {
        if lvars[k].starts_with("pe") {
            let n = lvars[k].after("pe").force_usize();
            let mut bcs = Vec::<String>::new();
            let mut lis = Vec::<usize>::new();
            let mut count = 0;
            let mut to_index = Vec::<usize>::new();
            for u in 0..nexacts {
                let clonotype_id = exacts[u];
                let ex = &exact_clonotypes[clonotype_id];
                for l in 0..ex.clones.len() {
                    let bc = &ex.clones[l][0].barcode;
                    let li = ex.clones[l][0].dataset_index;
                    if gex_info.pca[li].contains_key(&bc.clone()) {
                        bcs.push(bc.to_string());
                        lis.push(li);
                        to_index.push(count);
                    }
                    count += 1;
                }
            }
            let mut e: EquivRel = EquivRel::new(bcs.len() as i32);
            let mut mat = vec![Vec::<f64>::new(); bcs.len()];
            for i in 0..bcs.len() {
                mat[i] = gex_info.pca[lis[i]][&bcs[i].clone()].clone();
            }
            for i1 in 0..bcs.len() {
                for i2 in i1 + 1..bcs.len() {
                    if e.class_id(i1 as i32) != e.class_id(i2 as i32) {
                        let mut d = 0.0;
                        for j in 0..mat[i1].len() {
                            d += (mat[i1][j] - mat[i2][j]) * (mat[i1][j] - mat[i2][j]);
                        }
                        d = d.sqrt();
                        if d <= n as f64 {
                            e.join(i1 as i32, i2 as i32);
                        }
                    }
                }
            }
            pe[k] = vec![String::new(); count];
            let mut ids = Vec::<i32>::new();
            for i in 0..bcs.len() {
                ids.push(e.class_id(i as i32));
            }
            unique_sort(&mut ids);
            let mut reps = Vec::<i32>::new();
            e.orbit_reps(&mut reps);
            reps.sort();
            for i in 0..bcs.len() {
                pe[k][to_index[i]] = format!("{}", bin_position(&ids, &e.class_id(i as i32)));
            }
        }
    }

    // Compute ppe (PCA distance).

    *ppe = vec![Vec::<String>::new(); lvars.len()];
    for k in 0..lvars.len() {
        if lvars[k].starts_with("ppe") {
            let n = lvars[k].after("ppe").force_usize();
            let mut bcs = Vec::<String>::new();
            let mut lis = Vec::<usize>::new();
            let mut count = 0;
            let mut to_index = Vec::<usize>::new();
            for u in 0..nexacts {
                let clonotype_id = exacts[u];
                let ex = &exact_clonotypes[clonotype_id];
                for l in 0..ex.clones.len() {
                    let bc = &ex.clones[l][0].barcode;
                    let li = ex.clones[l][0].dataset_index;
                    if gex_info.pca[li].contains_key(&bc.clone()) {
                        bcs.push(bc.to_string());
                        lis.push(li);
                        to_index.push(count);
                    }
                    count += 1;
                }
            }
            let mut mat = vec![Vec::<f64>::new(); bcs.len()];
            for i in 0..bcs.len() {
                mat[i] = gex_info.pca[lis[i]][&bcs[i].clone()].clone();
            }
            let mut matg = Vec::<Vec<f64>>::new();
            for li in 0..ctl.origin_info.n() {
                for i in gex_info.pca[li].iter() {
                    matg.push(i.1.to_vec());
                }
            }
            let mut x = vec![0; bcs.len()];
            for i1 in 0..mat.len() {
                for i2 in 0..matg.len() {
                    let m1 = &mat[i1];
                    let m2 = &matg[i2];
                    let mut d = 0.0;
                    for j in 0..m1.len() {
                        d += (m1[j] - m2[j]) * (m1[j] - m2[j]);
                    }
                    d = d.sqrt();
                    if d <= n as f64 {
                        x[i1] += 1;
                    }
                }
            }
            let mut y = vec![0; bcs.len()];
            for i1 in 0..mat.len() {
                for i2 in 0..mat.len() {
                    let m1 = &mat[i1];
                    let m2 = &mat[i2];
                    let mut d = 0.0;
                    for j in 0..m1.len() {
                        d += (m1[j] - m2[j]) * (m1[j] - m2[j]);
                    }
                    d = d.sqrt();
                    if d <= n as f64 {
                        y[i1] += 1;
                    }
                }
            }
            ppe[k] = vec![String::new(); count];
            for i in 0..bcs.len() {
                ppe[k][to_index[i]] = format!("{:.1}", 100.0 * y[i] as f64 / x[i] as f64);
            }
        }
    }

    // Compute npe (PCA distance).

    *npe = vec![Vec::<String>::new(); lvars.len()];
    for k in 0..lvars.len() {
        if lvars[k].starts_with("npe") {
            let n = lvars[k].after("npe").force_usize();
            let mut bcs = Vec::<String>::new();
            let mut lis = Vec::<usize>::new();
            let mut count = 0;
            let mut to_index = Vec::<usize>::new();
            for u in 0..nexacts {
                let clonotype_id = exacts[u];
                let ex = &exact_clonotypes[clonotype_id];
                for l in 0..ex.clones.len() {
                    let bc = &ex.clones[l][0].barcode;
                    let li = ex.clones[l][0].dataset_index;
                    if gex_info.pca[li].contains_key(&bc.clone()) {
                        bcs.push(bc.to_string());
                        lis.push(li);
                        to_index.push(count);
                    }
                    count += 1;
                }
            }
            let mut mat = vec![Vec::<f64>::new(); bcs.len()];
            for i in 0..bcs.len() {
                mat[i] = gex_info.pca[lis[i]][&bcs[i].clone()].clone();
            }
            let mut y = vec![0; bcs.len()];
            for i1 in 0..mat.len() {
                for i2 in 0..mat.len() {
                    let m1 = &mat[i1];
                    let m2 = &mat[i2];
                    let mut d = 0.0;
                    for j in 0..m1.len() {
                        d += (m1[j] - m2[j]) * (m1[j] - m2[j]);
                    }
                    d = d.sqrt();
                    if d <= n as f64 {
                        y[i1] += 1;
                    }
                }
            }
            npe[k] = vec![String::new(); count];
            for i in 0..bcs.len() {
                npe[k][to_index[i]] = format!("{}", y[i]);
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn compute_bu(
    u: usize,
    cell_count: usize,
    exacts: &Vec<usize>,
    lvars: &Vec<String>,
    ctl: &EncloneControl,
    bli: &Vec<(String, usize, usize)>,
    ex: &ExactClonotype,
    exact_clonotypes: &Vec<ExactClonotype>,
    row: &mut Vec<String>,
    subrows: &mut Vec<Vec<String>>,
    varmat: &Vec<Vec<Vec<u8>>>,
    have_gex: bool,
    gex_info: &GexInfo,
    rsi: &ColInfo,
    sr: &mut Vec<(Vec<String>, Vec<Vec<String>>, Vec<Vec<u8>>, usize)>,
    fate: &Vec<HashMap<String, String>>,
    nd_fields: &Vec<String>,
    alt_bcs: &Vec<String>,
    cred: &Vec<Vec<String>>,
    pe: &Vec<Vec<String>>,
    ppe: &Vec<Vec<String>>,
    npe: &Vec<Vec<String>>,
    d_all: &Vec<Vec<u32>>,
    ind_all: &Vec<Vec<u32>>,
    mat: &Vec<Vec<Option<usize>>>,
) {
    // Very bad computation because of embedded binary search.

    let cols = mat.len();
    *subrows = Vec::<Vec<String>>::new();
    if ctl.clono_print_opt.bu {
        for bcl in bli.iter() {
            let mut row = Vec::<String>::new();
            let bc = &bcl.0;
            let li = bcl.1;
            let di = ex.clones[bcl.2][0].dataset_index;
            row.push(format!("$  {}", bc.clone()));
            let ex = &exact_clonotypes[exacts[u]];
            for k in 0..lvars.len() {
                let nr = row.len();
                let mut filled = false;
                for l in 0..ctl.origin_info.n() {
                    if lvars[k] == format!("n_{}", ctl.origin_info.dataset_id[l]) {
                        let mut n = 0;
                        if di == l {
                            n = 1;
                        }
                        row.push(format!("{}", n));
                        filled = true;
                    }
                }
                if filled {
                } else if lvars[k] == "n_b".to_string() {
                    let mut n = 0;
                    let li = ex.clones[bcl.2][0].dataset_index;
                    if gex_info.cell_type[li].contains_key(&bc.clone()) {
                        if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                            n = 1;
                        }
                    }
                    row.push(format!("{}", n));
                } else if lvars[k] == "filter".to_string() {
                    let mut f = String::new();
                    if fate[li].contains_key(&bc.clone()) {
                        f = fate[li][&bc.clone()].clone();
                        f = f.between(" ", " ").to_string();
                    }
                    row.push(f);
                } else if lvars[k] == "n_other".to_string() {
                    let mut n = 0;
                    let di = ex.clones[bcl.2][0].dataset_index;
                    let f = format!("n_{}", ctl.origin_info.dataset_id[di]);
                    let mut found = false;
                    for i in 0..nd_fields.len() {
                        if f == nd_fields[i] {
                            found = true;
                        }
                    }
                    if !found {
                        n = 1;
                    }
                    row.push(format!("{}", n));
                } else if lvars[k] == "sec".to_string() {
                    let mut n = 0;
                    if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                        n = ctl.origin_info.secmem[li][&bc.clone()].0;
                    }
                    row.push(format!("{}", n));
                } else if lvars[k] == "mem".to_string() {
                    let mut n = 0;
                    if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                        n = ctl.origin_info.secmem[li][&bc.clone()].1;
                    }
                    row.push(format!("{}", n));
                } else if bin_member(&alt_bcs, &lvars[k]) {
                    let mut val = String::new();
                    let alt = &ctl.origin_info.alt_bc_fields[li];
                    for j in 0..alt.len() {
                        if alt[j].0 == lvars[k] {
                            if alt[j].1.contains_key(&bc.clone()) {
                                val = alt[j].1[&bc.clone()].clone();
                            }
                        }
                    }
                    row.push(val);
                } else if lvars[k] == "datasets".to_string() {
                    row.push(format!("{}", ctl.origin_info.dataset_id[li].clone()));
                } else if lvars[k] == "origins".to_string() {
                    row.push(format!("{}", ctl.origin_info.origin_id[li].clone()));
                } else if lvars[k] == "donors".to_string() {
                    row.push(format!("{}", ctl.origin_info.donor_id[li].clone()));
                } else if lvars[k] == "clust".to_string() && have_gex {
                    let mut cid = 0;
                    if gex_info.cluster[li].contains_key(&bc.clone()) {
                        cid = gex_info.cluster[li][&bc.clone()];
                    }
                    row.push(format!("{}", cid));
                } else if lvars[k].starts_with("pe") && have_gex {
                    row.push(format!("{}", pe[k][cell_count + bcl.2]));
                } else if lvars[k].starts_with("npe") && have_gex {
                    row.push(format!("{}", npe[k][cell_count + bcl.2]));
                } else if lvars[k].starts_with("ppe") && have_gex {
                    row.push(format!("{}", ppe[k][cell_count + bcl.2]));
                } else if lvars[k] == "cred".to_string() && have_gex {
                    row.push(format!("{}", cred[k][cell_count + bcl.2]));
                } else if lvars[k] == "type".to_string() && have_gex {
                    let mut cell_type = "".to_string();
                    if gex_info.cell_type[li].contains_key(&bc.clone()) {
                        cell_type = gex_info.cell_type[li][&bc.clone()].clone();
                    }
                    row.push(cell_type);
                } else if lvars[k] == "n_gex".to_string() && have_gex {
                    let mut n_gex = 0;
                    if bin_member(&gex_info.gex_cell_barcodes[li], &bc) {
                        n_gex = 1;
                    }
                    row.push(format!("{}", n_gex));
                } else if lvars[k] == "mark".to_string() {
                    let mut mark = String::new();
                    if ex.clones[bcl.2][0].marked {
                        mark = "x".to_string();
                    }
                    row.push(mark);
                } else if lvars[k] == "entropy".to_string() && have_gex {
                    // NOTE DUPLICATION WITH CODE BELOW.
                    let mut gex_count = 0;
                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                    if p >= 0 {
                        let mut raw_count = 0;
                        if gex_info.gex_matrices[li].initialized() {
                            let row = gex_info.gex_matrices[li].row(p as usize);
                            for j in 0..row.len() {
                                let f = row[j].0;
                                let n = row[j].1;
                                if gex_info.is_gex[li][f] {
                                    raw_count += n;
                                }
                            }
                        } else {
                            let l = bcl.2;
                            for j in 0..d_all[l].len() {
                                if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                    raw_count += d_all[l][j] as usize;
                                }
                            }
                        }
                        gex_count = raw_count;
                    }
                    let mut entropy = 0.0;
                    if p >= 0 {
                        if gex_info.gex_matrices[li].initialized() {
                            let row = gex_info.gex_matrices[li].row(p as usize);
                            for j in 0..row.len() {
                                let f = row[j].0;
                                let n = row[j].1;
                                if gex_info.is_gex[li][f] {
                                    let q = n as f64 / gex_count as f64;
                                    entropy -= q * q.log2();
                                }
                            }
                        } else {
                            let l = bcl.2;
                            for j in 0..d_all[l].len() {
                                if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                    let n = d_all[l][j] as usize;
                                    let q = n as f64 / gex_count as f64;
                                    entropy -= q * q.log2();
                                }
                            }
                        }
                    }
                    row.push(format!("{:.2}", entropy));
                } else if have_gex {
                    // this calc isn't needed except in _% case below
                    // TODO: ELIMINATE UNNEEDED CALC
                    let mut gex_count = 0.0;
                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                    if p >= 0 {
                        let mut raw_count = 0 as f64;
                        if gex_info.gex_matrices[li].initialized() {
                            let row = gex_info.gex_matrices[li].row(p as usize);
                            for j in 0..row.len() {
                                let f = row[j].0;
                                let n = row[j].1;
                                if gex_info.is_gex[li][f] {
                                    raw_count += n as f64;
                                }
                            }
                        } else {
                            let l = bcl.2;
                            for j in 0..d_all[l].len() {
                                if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                    raw_count += d_all[l][j] as f64;
                                }
                            }
                        }
                        if !ctl.gen_opt.full_counts {
                            gex_count = raw_count * gex_info.gex_mults[li];
                        } else {
                            gex_count = raw_count;
                        }
                    }
                    if lvars[k] == "gex".to_string() {
                        row.push(format!("{}", gex_count.round()));
                    } else {
                        let mut y = lvars[k].clone();
                        if y.contains(':') {
                            y = y.after(":").to_string();
                        }
                        let y0 = y.clone();
                        let suffixes = ["_min", "_max", "_μ", "_Σ", "_cell", "_%"];
                        for s in suffixes.iter() {
                            if y.ends_with(s) {
                                y = y.rev_before(&s).to_string();
                                break;
                            }
                        }
                        let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                        let mut computed = false;
                        let mut count = 0.0;
                        let l = bcl.2;
                        if p >= 0 {
                            let mut ux = Vec::<usize>::new();
                            if ctl.clono_print_opt.regex_match[li].contains_key(&y) {
                                ux = ctl.clono_print_opt.regex_match[li][&y].clone();
                            }
                            if ux.len() > 0 {
                                computed = true;
                                for fid in ux.iter() {
                                    let counti = get_gex_matrix_entry(
                                        &ctl, &gex_info, *fid, &d_all, &ind_all, li, l, p as usize,
                                        &y,
                                    );
                                    count += counti;
                                }
                            } else if gex_info.feature_id[li].contains_key(&y) {
                                computed = true;
                                let fid = gex_info.feature_id[li][&y];
                                count = get_gex_matrix_entry(
                                    &ctl, &gex_info, fid, &d_all, &ind_all, li, l, p as usize, &y,
                                );
                            }
                        }
                        if computed {
                            // note unneeded calculation above in certain cases
                            // TODO: ELIMINATE!
                            if y0.ends_with("_min") {
                            } else if y0.ends_with("_max") {
                            } else if y0.ends_with("_μ") {
                            } else if y0.ends_with("_Σ") {
                            } else if y0.ends_with("_%") {
                                row.push(format!("{:.2}", (100.0 * count) / gex_count));
                            } else {
                                row.push(format!("{}", count.round()));
                            }
                        }
                    }
                }
                if row.len() == nr {
                    row.push("".to_string());
                }
            }
            let mut ncall = 0;
            for k in 0..cols {
                ncall += rsi.cvars[k].len();
            }
            let mut cx = vec!["".to_string(); ncall];
            let mut cp = 0;
            for col in 0..cols {
                let m = mat[col][u];
                if m.is_some() {
                    let m = m.unwrap();
                    for p in 0..rsi.cvars[col].len() {
                        if rsi.cvars[col][p] == "u".to_string() {
                            let numi = ex.clones[bcl.2][m].umi_count;
                            cx[cp + p] = format!("{}", numi);
                        } else if rsi.cvars[col][p] == "r".to_string() {
                            let r = ex.clones[bcl.2][m].read_count;
                            cx[cp + p] = format!("{}", r);
                        } else if rsi.cvars[col][p] == "nval".to_string() {
                            let mut n = 0;
                            if ex.clones[bcl.2][m].validated_umis.is_some() {
                                n = ex.clones[bcl.2][m].validated_umis.as_ref().unwrap().len();
                            }
                            cx[cp + p] = format!("{}", n);
                        } else if rsi.cvars[col][p] == "nnval".to_string() {
                            let mut n = 0;
                            if ex.clones[bcl.2][m].non_validated_umis.is_some() {
                                n = ex.clones[bcl.2][m]
                                    .non_validated_umis
                                    .as_ref()
                                    .unwrap()
                                    .len();
                            }
                            cx[cp + p] = format!("{}", n);
                        } else if rsi.cvars[col][p] == "nival".to_string() {
                            let mut n = 0;
                            if ex.clones[bcl.2][m].invalidated_umis.is_some() {
                                n = ex.clones[bcl.2][m].invalidated_umis.as_ref().unwrap().len();
                            }
                            cx[cp + p] = format!("{}", n);
                        } else if rsi.cvars[col][p] == "valumis".to_string() {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].validated_umis.is_some() {
                                n = ex.clones[bcl.2][m].non_validated_umis.clone().unwrap();
                            }
                            cx[cp + p] = format!("{}", n.iter().format(","));
                        } else if rsi.cvars[col][p] == "valbcumis".to_string() {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].validated_umis.is_some() {
                                n = ex.clones[bcl.2][m].validated_umis.clone().unwrap();
                                for i in 0..n.len() {
                                    n[i] = format!(
                                        "{}{}",
                                        ex.clones[bcl.2][m].barcode.before("-"),
                                        n[i]
                                    );
                                }
                            }
                            cx[cp + p] = format!("{}", n.iter().format(","));
                        } else if rsi.cvars[col][p] == "nvalumis".to_string() {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].non_validated_umis.is_some() {
                                n = ex.clones[bcl.2][m].non_validated_umis.clone().unwrap();
                            }
                            cx[cp + p] = format!("{}", n.iter().format(","));
                        } else if rsi.cvars[col][p] == "nvalbcumis".to_string() {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].non_validated_umis.is_some() {
                                n = ex.clones[bcl.2][m].non_validated_umis.clone().unwrap();
                                for i in 0..n.len() {
                                    n[i] = format!(
                                        "{}{}",
                                        ex.clones[bcl.2][m].barcode.before("-"),
                                        n[i]
                                    );
                                }
                            }
                            cx[cp + p] = format!("{}", n.iter().format(","));
                        } else if rsi.cvars[col][p] == "ivalumis".to_string() {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].invalidated_umis.is_some() {
                                n = ex.clones[bcl.2][m].invalidated_umis.clone().unwrap();
                            }
                            cx[cp + p] = format!("{}", n.iter().format(","));
                        } else if rsi.cvars[col][p] == "ivalbcumis".to_string() {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].invalidated_umis.is_some() {
                                n = ex.clones[bcl.2][m].invalidated_umis.clone().unwrap();
                                for i in 0..n.len() {
                                    n[i] = format!(
                                        "{}{}",
                                        ex.clones[bcl.2][m].barcode.before("-"),
                                        n[i]
                                    );
                                }
                            }
                            cx[cp + p] = format!("{}", n.iter().format(","));
                        }
                    }
                }
                cp += rsi.cvars[col].len();
            }
            row.append(&mut cx);
            subrows.push(row);
        }
    }
    sr.push((row.to_vec(), subrows.to_vec(), varmat[u].clone(), u));
}
