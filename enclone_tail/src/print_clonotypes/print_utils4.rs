// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use amino::codon_to_aa;
use enclone_core::{
    defs::{ColInfo, EncloneControl, ExactClonotype, GexInfo},
    enclone_structs::BarcodeFates,
};
use enclone_proto::types::DonorReferenceItem;
use equiv::EquivRel;
use itertools::Itertools;

use super::finish_table::Sr;
use enclone_process::get_gex_matrix_entry::get_gex_matrix_entry;
use string_utils::TextUtils;
use vdj_ann::refx::RefData;
use vector_utils::{bin_member, bin_position, bin_position1_2, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Define amino acid positions to show.

pub fn build_show_aa(
    ctl: &EncloneControl,
    rsi: &ColInfo,
    vars_amino: &[Vec<usize>],
    shares_amino: &[Vec<usize>],
    refdata: &RefData,
    dref: &[DonorReferenceItem],
    _exacts: &[usize],
    _exact_clonotypes: &[ExactClonotype],
) -> Vec<Vec<usize>> {
    let cols = rsi.vids.len();
    let mut show_aa = vec![Vec::<usize>::new(); cols];
    for cx in 0..cols {
        for x in &ctl.clono_print_opt.amino {
            if x.contains('-') {
                let (start, stop) = (x.before("-").force_usize(), x.after("-").force_usize());
                for p in start..=stop {
                    show_aa[cx].push(p);
                }
            }
        }
        if ctl.clono_print_opt.amino.contains(&"cdr1".to_string())
            && rsi.cdr1_starts[cx].is_some()
            && rsi.fr2_starts[cx].is_some()
            && rsi.cdr1_starts[cx].unwrap() <= rsi.fr2_starts[cx].unwrap()
        {
            for j in (rsi.cdr1_starts[cx].unwrap()..rsi.fr2_starts[cx].unwrap()).step_by(3) {
                show_aa[cx].push(j / 3);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"cdr2".to_string())
            && rsi.cdr2_starts[cx].is_some()
            && rsi.fr3_starts[cx].is_some()
            && rsi.cdr2_starts[cx].unwrap() <= rsi.fr3_starts[cx].unwrap()
        {
            for j in (rsi.cdr2_starts[cx].unwrap()..rsi.fr3_starts[cx].unwrap()).step_by(3) {
                show_aa[cx].push(j / 3);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"cdr3".to_string()) {
            for j in 0..rsi.cdr3_lens[cx] {
                let p = rsi.cdr3_starts[cx] / 3 + j;
                show_aa[cx].push(p);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"fwr1".to_string())
            && rsi.cdr1_starts[cx].is_some()
            && rsi.fr1_starts[cx] <= rsi.cdr1_starts[cx].unwrap()
        {
            for j in (rsi.fr1_starts[cx]..rsi.cdr1_starts[cx].unwrap()).step_by(3) {
                show_aa[cx].push(j / 3);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"fwr2".to_string())
            && rsi.fr2_starts[cx].is_some()
            && rsi.cdr2_starts[cx].is_some()
            && rsi.fr2_starts[cx].unwrap() <= rsi.cdr2_starts[cx].unwrap()
        {
            for j in (rsi.fr2_starts[cx].unwrap()..rsi.cdr2_starts[cx].unwrap()).step_by(3) {
                show_aa[cx].push(j / 3);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"fwr3".to_string())
            && rsi.fr3_starts[cx].is_some()
            && rsi.fr3_starts[cx].unwrap() <= rsi.cdr3_starts[cx]
        {
            for j in (rsi.fr3_starts[cx].unwrap()..rsi.cdr3_starts[cx]).step_by(3) {
                show_aa[cx].push(j / 3);
            }
        }
        if ctl.clono_print_opt.amino.contains(&"fwr4".to_string()) {
            for j in
                (rsi.cdr3_starts[cx] + 3 * rsi.cdr3_lens[cx]..rsi.seq_del_lens[cx] - 1).step_by(3)
            {
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
            let vseq2 = if let Some(vpid) = rsi.vpids[cx] {
                &dref[vpid].nt_sequence
            } else {
                &vseq1
            };
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
            let vseq2 = if let Some(vpid) = rsi.vpids[cx] {
                &dref[vpid].nt_sequence
            } else {
                &vseq1
            };
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

pub struct SomeStats {
    pub cred: Vec<Vec<String>>,
    pub pe: Vec<Vec<String>>,
    pub ppe: Vec<Vec<String>>,
    pub npe: Vec<Vec<String>>,
}

pub fn compute_some_stats(
    ctl: &EncloneControl,
    lvars: &[String],
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    gex_info: &GexInfo,
    vdj_cells: &[Vec<String>],
    n_vdj_gex: &[usize],
) -> SomeStats {
    let nexacts = exacts.len();

    // Compute "cred" stats (credibility/# of neighboring cells that are also B cells).

    let mut cred = vec![Vec::<String>::new(); lvars.len()];
    for (lvar, cred) in lvars.iter().zip(cred.iter_mut()) {
        if lvar == "cred" {
            for &clonotype_id in exacts.iter().take(nexacts) {
                let ex = &exact_clonotypes[clonotype_id];
                for clone in &ex.clones {
                    let bc = &clone[0].barcode;
                    let li = clone[0].dataset_index;
                    if gex_info.pca[li].contains_key(&bc.clone()) {
                        let mut creds = 0;
                        let mut z = Vec::<(f64, String)>::new();
                        let x = &gex_info.pca[li][&bc.clone()];
                        for y in &gex_info.pca[li] {
                            let dist2 = x
                                .iter()
                                .zip(y.1.iter())
                                .map(|(&xm, &ym)| (ym - xm) * (ym - xm))
                                .sum::<f64>();
                            z.push((dist2, y.0.clone()));
                        }
                        z.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        let top = n_vdj_gex[li];
                        for zi in z.iter().take(top) {
                            if bin_member(&vdj_cells[li], &zi.1) {
                                creds += 1;
                            }
                        }
                        let pc = 100.0 * creds as f64 / top as f64;
                        cred.push(format!("{pc:.1}"));
                    } else {
                        cred.push(String::new());
                    }
                }
            }
        }
    }

    // Compute pe (PCA distance).

    let mut pe = vec![Vec::<String>::new(); lvars.len()];
    for k in 0..lvars.len() {
        if lvars[k].starts_with("pe") {
            let n = lvars[k].after("pe").force_usize();
            let mut bcs = Vec::<String>::new();
            let mut lis = Vec::<usize>::new();
            let mut count = 0;
            let mut to_index = Vec::<usize>::new();
            for &clonotype_id in exacts.iter().take(nexacts) {
                let ex = &exact_clonotypes[clonotype_id];
                for clone in &ex.clones {
                    let bc = &clone[0].barcode;
                    let li = clone[0].dataset_index;
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
            reps.sort_unstable();
            for (i, idx) in to_index.into_iter().take(bcs.len()).enumerate() {
                pe[k][idx] = format!("{}", bin_position(&ids, &e.class_id(i as i32)));
            }
        }
    }

    // Compute ppe (PCA distance).

    let mut ppe = vec![Vec::<String>::new(); lvars.len()];
    for k in 0..lvars.len() {
        if lvars[k].starts_with("ppe") {
            let n = lvars[k].after("ppe").force_usize();
            let mut bcs = Vec::<&str>::new();
            let mut lis = Vec::<usize>::new();
            let mut count = 0;
            let mut to_index = Vec::<usize>::new();
            for &clonotype_id in exacts.iter().take(nexacts) {
                let ex = &exact_clonotypes[clonotype_id];
                for clone in &ex.clones {
                    let bc = &clone[0].barcode;
                    let li = clone[0].dataset_index;
                    if gex_info.pca[li].contains_key(&bc.clone()) {
                        bcs.push(bc.as_str());
                        lis.push(li);
                        to_index.push(count);
                    }
                    count += 1;
                }
            }
            let mat = bcs
                .iter()
                .zip(lis.iter())
                .map(|(bc, &li)| gex_info.pca[li][&(*bc).to_string()].clone())
                .collect::<Vec<_>>();
            let mut matg = Vec::<Vec<f64>>::new();
            for li in 0..ctl.origin_info.n() {
                for i in &gex_info.pca[li] {
                    matg.push(i.1.clone());
                }
            }
            let mut x = vec![0; bcs.len()];
            for (m1, x) in mat.iter().zip(x.iter_mut()) {
                for m2 in &matg {
                    let d = m1
                        .iter()
                        .zip(m2.iter())
                        .map(|(&m1, &m2)| (m1 - m2) * (m1 - m2))
                        .sum::<f64>()
                        .sqrt();
                    if d <= n as f64 {
                        *x += 1;
                    }
                }
            }
            let mut y = vec![0; bcs.len()];
            for i1 in 0..mat.len() {
                for i2 in 0..mat.len() {
                    let m1 = &mat[i1];
                    let m2 = &mat[i2];
                    let d = m1
                        .iter()
                        .zip(m2.iter())
                        .map(|(&m1, &m2)| (m1 - m2) * (m1 - m2))
                        .sum::<f64>()
                        .sqrt();
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

    let mut npe = vec![Vec::<String>::new(); lvars.len()];
    for k in 0..lvars.len() {
        if lvars[k].starts_with("npe") {
            let n = lvars[k].after("npe").force_usize();
            let mut bcs = Vec::<&str>::new();
            let mut lis = Vec::<usize>::new();
            let mut count = 0;
            let mut to_index = Vec::<usize>::new();
            for &clonotype_id in exacts.iter().take(nexacts) {
                let ex = &exact_clonotypes[clonotype_id];
                for clone in &ex.clones {
                    let bc = &clone[0].barcode;
                    let li = clone[0].dataset_index;
                    if gex_info.pca[li].contains_key(&bc.clone()) {
                        bcs.push(bc.as_str());
                        lis.push(li);
                        to_index.push(count);
                    }
                    count += 1;
                }
            }
            let mat = bcs
                .iter()
                .zip(lis.iter())
                .map(|(&bc, &li)| gex_info.pca[li][&bc.to_string()].clone())
                .collect::<Vec<_>>();
            let y = mat
                .iter()
                .map(|m1| {
                    mat.iter()
                        .filter(|m2| {
                            m1.iter()
                                .zip(m2.iter())
                                .map(|(&m1, &m2)| (m1 - m2) * (m1 - m2))
                                .sum::<f64>()
                                .sqrt()
                                <= n as f64
                        })
                        .count()
                })
                .collect::<Vec<_>>();
            npe[k] = vec![String::new(); count];
            for i in 0..bcs.len() {
                npe[k][to_index[i]] = format!("{}", y[i]);
            }
        }
    }

    SomeStats { cred, pe, ppe, npe }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn compute_bu(
    u: usize,
    cell_count: usize,
    exacts: &[usize],
    lvars: &[String],
    ctl: &EncloneControl,
    bli: &[(String, usize, usize)],
    ex: &ExactClonotype,
    exact_clonotypes: &[ExactClonotype],
    row: &mut [String],
    subrows: &mut Vec<Vec<String>>,
    have_gex: bool,
    gex_info: &GexInfo,
    rsi: &ColInfo,
    sr: &mut Vec<Sr>,
    fate: &[BarcodeFates],
    nd_fields: &[String],
    alt_bcs: &[String],
    cred: &[Vec<String>],
    pe: &[Vec<String>],
    ppe: &[Vec<String>],
    npe: &[Vec<String>],
    d_all: &[Vec<u32>],
    ind_all: &[Vec<u32>],
    mat: &[Vec<Option<usize>>],
    these_stats: &[(String, Vec<String>)],
    refdata: &RefData,
) {
    // Very bad computation because of embedded binary search.

    let cols = mat.len();
    *subrows = Vec::<Vec<String>>::new();
    if ctl.clono_print_opt.bu {
        for bcl in bli {
            let mut row = Vec::<String>::new();
            let bc = &bcl.0;
            let li = bcl.1;
            let di = ex.clones[bcl.2][0].dataset_index;
            row.push(format!("$  {bc}"));
            let ex = &exact_clonotypes[exacts[u]];
            for (k, var) in lvars.iter().enumerate() {
                let p = bin_position1_2(these_stats, var);
                let nr = row.len();
                let mut filled = false;
                for l in 0..ctl.origin_info.n() {
                    if var.starts_with("n_") && var[2..] == ctl.origin_info.dataset_id[l] {
                        let n = if di == l { "1" } else { "0" };
                        row.push(n.to_string());
                        filled = true;
                    }
                }

                // Not sure it makes sense to do the below for just var_def variables.
                // Why not all?

                let mut in_var_def = false;
                for i in 0..ctl.gen_opt.var_def.len() {
                    if ctl.gen_opt.var_def[i].0 == *var {
                        in_var_def = true;
                    }
                }
                if var == "nbc" {
                    in_var_def = true;
                }
                if in_var_def && p >= 0 {
                    let stats_me = &these_stats[p as usize].1;
                    row.push(stats_me[bcl.2].clone());
                    filled = true;
                }

                // Proceed.

                if filled {

                    // Many of the variables (presumably a lot) should be computed by just
                    // looking up their values in these_stats.  This should allow elimination of a
                    // lot of the code below.  We give one example here.
                } else if p >= 0
                    && var.starts_with("fb")
                    && var.after("fb").ends_with("_n")
                    && var.after("fb").rev_before("_n").parse::<usize>().is_ok()
                    && var.after("fb").rev_before("_n").force_usize() >= 1
                {
                    let stats_me = &these_stats[p as usize].1;
                    row.push(stats_me[bcl.2].clone());
                } else if var == "n_b" {
                    let mut n = 0;
                    let li = ex.clones[bcl.2][0].dataset_index;
                    if gex_info.cell_type[li].contains_key(bc)
                        && gex_info.cell_type[li][&bc.clone()].starts_with('B')
                    {
                        n = 1;
                    }
                    row.push(format!("{n}"));
                } else if var == "filter" {
                    let mut f = String::new();
                    if fate[li].contains_key(bc) {
                        f = fate[li][bc].label().to_string();
                    }
                    row.push(f);
                } else if var == "n_other" {
                    let mut n = 0;
                    let di = ex.clones[bcl.2][0].dataset_index;
                    let f = format!("n_{}", ctl.origin_info.dataset_id[di]);
                    let found = nd_fields.iter().any(|ff| *ff == f);
                    if !found {
                        n = 1;
                    }
                    row.push(format!("{n}"));
                } else if var == "sec" {
                    let mut n = 0;
                    if ctl.origin_info.secmem[li].contains_key(bc) {
                        n = ctl.origin_info.secmem[li][bc].0;
                    }
                    row.push(format!("{n}"));
                } else if var == "mem" {
                    let mut n = 0;
                    if ctl.origin_info.secmem[li].contains_key(bc) {
                        n = ctl.origin_info.secmem[li][bc].1;
                    }
                    row.push(format!("{n}"));
                } else if bin_member(alt_bcs, var) {
                    let mut val = String::new();
                    let alt = &ctl.origin_info.alt_bc_fields[li];
                    for aj in alt {
                        if aj.0 == *var && aj.1.contains_key(bc) {
                            val = aj.1[bc].clone();
                        }
                    }
                    row.push(val);
                } else if var == "datasets" {
                    row.push(ctl.origin_info.dataset_id[li].clone());
                } else if var == "origins" {
                    row.push(ctl.origin_info.origin_id[li].clone());
                } else if var == "donors" {
                    row.push(ctl.origin_info.donor_id[li].clone());
                } else if var == "clust" && have_gex {
                    let mut cid = 0;
                    if gex_info.cluster[li].contains_key(bc) {
                        cid = gex_info.cluster[li][bc];
                    }
                    row.push(format!("{cid}"));
                } else if var.starts_with("pe") && have_gex {
                    row.push(pe[k][cell_count + bcl.2].to_string());
                } else if var.starts_with("npe") && have_gex {
                    row.push(npe[k][cell_count + bcl.2].to_string());
                } else if var.starts_with("ppe") && have_gex {
                    row.push(ppe[k][cell_count + bcl.2].to_string());
                } else if var == "cred" && have_gex {
                    row.push(cred[k][cell_count + bcl.2].to_string());
                } else if var == "type" && have_gex {
                    let cell_type = if gex_info.cell_type[li].contains_key(bc) {
                        gex_info.cell_type[li][bc].clone()
                    } else {
                        String::default()
                    };
                    row.push(cell_type);
                } else if var == "n_gex" && have_gex {
                    let mut n_gex = 0;
                    if bin_member(&gex_info.gex_cell_barcodes[li], bc) {
                        n_gex = 1;
                    }
                    row.push(format!("{n_gex}"));
                } else if var == "mark" {
                    let mut mark = String::new();
                    if ex.clones[bcl.2][0].marked {
                        mark = "x".to_string();
                    }
                    row.push(mark);
                } else if var == "entropy" && have_gex {
                    // NOTE DUPLICATION WITH CODE BELOW.
                    let mut gex_count = 0;
                    let p = bin_position(&gex_info.gex_barcodes[li], bc);
                    if p >= 0 {
                        let mut raw_count = 0;

                        let l = bcl.2;
                        for j in 0..d_all[l].len() {
                            if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                raw_count += d_all[l][j] as usize;
                            }
                        }

                        gex_count = raw_count;
                    }
                    let mut entropy = 0.0;
                    if p >= 0 {
                        let l = bcl.2;
                        for j in 0..d_all[l].len() {
                            if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                let n = d_all[l][j] as usize;
                                let q = n as f64 / gex_count as f64;
                                entropy -= q * q.log2();
                            }
                        }
                    }
                    row.push(format!("{entropy:.2}"));
                } else if have_gex {
                    // this calc isn't needed except in _% case below
                    // TODO: ELIMINATE UNNEEDED CALC
                    let mut gex_count = 0.0;
                    let p = bin_position(&gex_info.gex_barcodes[li], bc);
                    if p >= 0 {
                        let mut raw_count = 0 as f64;
                        let l = bcl.2;
                        for j in 0..d_all[l].len() {
                            if gex_info.is_gex[li][ind_all[l][j] as usize] {
                                raw_count += d_all[l][j] as f64;
                            }
                        }

                        if !ctl.gen_opt.full_counts {
                            gex_count = raw_count * gex_info.gex_mults[li];
                        } else {
                            gex_count = raw_count;
                        }
                    }
                    if var == "gex" {
                        row.push(format!("{}", gex_count.round()));
                    } else {
                        let mut y = var.as_str();
                        if y.contains(':') {
                            y = y.after(":");
                        }
                        let y0 = y;
                        let suffixes = ["_min", "_max", "_μ", "_Σ", "_cell", "_%"];
                        for &s in &suffixes {
                            if y.ends_with(s) {
                                y = y.rev_before(s);
                                break;
                            }
                        }
                        let p = bin_position(&gex_info.gex_barcodes[li], bc);
                        let mut computed = false;
                        let mut count = 0.0;
                        let l = bcl.2;
                        if p >= 0 {
                            let ux = ctl.clono_print_opt.regex_match[li]
                                .get(&y.to_string())
                                .cloned()
                                .unwrap_or_default();
                            if !ux.is_empty() {
                                computed = true;
                                for fid in &ux {
                                    let counti = get_gex_matrix_entry(
                                        ctl, gex_info, *fid, d_all, ind_all, li, l, y,
                                    );
                                    count += counti;
                                }
                            } else if let Some(&fid) = gex_info.feature_id[li].get(&y.to_string()) {
                                computed = true;
                                count = get_gex_matrix_entry(
                                    ctl, gex_info, fid, d_all, ind_all, li, l, y,
                                );
                            }
                        }
                        if computed {
                            // note unneeded calculation above in certain cases
                            // TODO: ELIMINATE!
                            if y0.ends_with("_min")
                                || y0.ends_with("_max")
                                || y0.ends_with("_μ")
                                || y0.ends_with("_Σ")
                            {
                            } else if y0.ends_with("_%") {
                                row.push(format!("{:.2}", (100.0 * count) / gex_count));
                            } else {
                                row.push(format!("{}", count.round()));
                            }
                        }
                    }
                }
                if row.len() == nr {
                    row.push(String::new());
                }
            }
            let mut ncall = 0;
            for k in 0..cols {
                ncall += rsi.cvars[k].len();
            }
            let mut cx = vec![String::new(); ncall];
            let mut cp = 0;
            for (m, cvars) in mat.iter().take(cols).zip(rsi.cvars.iter()) {
                if let Some(m) = m[u] {
                    for (cvar, cxp) in cvars.iter().zip(cx.iter_mut().skip(cp)) {
                        if cvar == "v_name_orig" {
                            let v = &refdata.name[ex.clones[bcl.2][m].v_ref_id];
                            *cxp = v.to_string();
                        } else if cvar == "u" {
                            let numi = ex.clones[bcl.2][m].umi_count;
                            *cxp = format!("{numi}");
                        } else if cvar == "r" {
                            let r = ex.clones[bcl.2][m].read_count;
                            *cxp = format!("{r}");
                        } else if cvar == "nval" {
                            let mut n = 0;
                            if ex.clones[bcl.2][m].validated_umis.is_some() {
                                n = ex.clones[bcl.2][m].validated_umis.as_ref().unwrap().len();
                            }
                            *cxp = format!("{n}");
                        } else if cvar == "nnval" {
                            let mut n = 0;
                            if ex.clones[bcl.2][m].non_validated_umis.is_some() {
                                n = ex.clones[bcl.2][m]
                                    .non_validated_umis
                                    .as_ref()
                                    .unwrap()
                                    .len();
                            }
                            *cxp = format!("{n}");
                        } else if cvar == "nival" {
                            let mut n = 0;
                            if ex.clones[bcl.2][m].invalidated_umis.is_some() {
                                n = ex.clones[bcl.2][m].invalidated_umis.as_ref().unwrap().len();
                            }
                            *cxp = format!("{n}");
                        } else if cvar == "valumis" {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].validated_umis.is_some() {
                                n = ex.clones[bcl.2][m].non_validated_umis.clone().unwrap();
                            }
                            *cxp = format!("{}", n.iter().format(","));
                        } else if cvar == "valbcumis" {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].validated_umis.is_some() {
                                n = ex.clones[bcl.2][m]
                                    .validated_umis
                                    .as_ref()
                                    .unwrap()
                                    .iter()
                                    .map(|ni| {
                                        format!("{}{ni}", ex.clones[bcl.2][m].barcode.before("-"))
                                    })
                                    .collect();
                            }
                            *cxp = format!("{}", n.iter().format(","));
                        } else if cvar == "nvalumis" {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].non_validated_umis.is_some() {
                                n = ex.clones[bcl.2][m].non_validated_umis.clone().unwrap();
                            }
                            *cxp = format!("{}", n.iter().format(","));
                        } else if cvar == "nvalbcumis" {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].non_validated_umis.is_some() {
                                n = ex.clones[bcl.2][m]
                                    .non_validated_umis
                                    .as_ref()
                                    .unwrap()
                                    .iter()
                                    .map(|ni| {
                                        format!("{}{ni}", ex.clones[bcl.2][m].barcode.before("-"))
                                    })
                                    .collect();
                            }
                            *cxp = format!("{}", n.iter().format(","));
                        } else if cvar == "ivalumis" {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].invalidated_umis.is_some() {
                                n = ex.clones[bcl.2][m].invalidated_umis.clone().unwrap();
                            }
                            *cxp = format!("{}", n.iter().format(","));
                        } else if cvar == "ivalbcumis" {
                            let mut n = Vec::<String>::new();
                            if ex.clones[bcl.2][m].invalidated_umis.is_some() {
                                n = ex.clones[bcl.2][m]
                                    .invalidated_umis
                                    .as_ref()
                                    .unwrap()
                                    .iter()
                                    .map(|ni| {
                                        format!("{}{ni}", ex.clones[bcl.2][m].barcode.before("-"))
                                    })
                                    .collect::<Vec<_>>();
                            }
                            *cxp = format!("{}", n.iter().format(","));
                        }
                    }
                }
                cp += cvars.len();
            }
            row.append(&mut cx);
            subrows.push(row);
        }
    }
    sr.push(Sr {
        row: row.to_vec(),
        subrows: subrows.clone(),
    });
}
