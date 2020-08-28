// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::*;
use enclone_proto::types::*;
use itertools::*;
use std::cmp::max;
use std::collections::HashMap;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Find variant positions.  And some other things.

pub fn vars_and_shares(
    pass: usize,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    rsi: &ColInfo,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    vars: &mut Vec<Vec<usize>>,
    vars_amino: &mut Vec<Vec<usize>>,
    shares_amino: &mut Vec<Vec<usize>>,
    out_data: &mut Vec<HashMap<String, String>>,
) {
    // Copied stuff.

    let pcols_sort = &ctl.parseable_opt.pcols_sort;
    let nexacts = exacts.len();
    let cols = rsi.vids.len();

    // Go through the columns.

    for cx in 0..cols {
        // go through each column
        let (mut vref, mut jref) = (Vec::<u8>::new(), Vec::<u8>::new());
        let mut vseq2 = Vec::<u8>::new();
        for u in 0..nexacts {
            let m = rsi.mat[cx][u];
            if m.is_some() {
                let m = m.unwrap();
                // Reference assigned multiple times here, is wrong.
                // Also where using allelized reference, need to explain.
                // vref is supposed to be the donor reference, but seems like it isn't
                vref = exact_clonotypes[exacts[u]].share[m].vs.to_ascii_vec();
                jref = exact_clonotypes[exacts[u]].share[m].js.to_ascii_vec();
            }
            let vseq1 = refdata.refs[rsi.vids[cx]].to_ascii_vec();
            if rsi.vpids[cx].is_some() {
                vseq2 = dref[rsi.vpids[cx].unwrap()].nt_sequence.clone();
            } else {
                vseq2 = vseq1.clone();
            }
        }
        let mut n = 0;
        for z in 0..rsi.seqss[cx].len() {
            n = max(n, rsi.seqss[cx][z].len());
        }
        let (mut v, mut s) = (Vec::<usize>::new(), Vec::<usize>::new());
        let (mut v_amino, mut s_amino) = (Vec::<usize>::new(), Vec::<usize>::new());
        for p in 0..n {
            let mut bases = Vec::<u8>::new();
            let mut bases_amino = Vec::<u8>::new();
            for s in 0..rsi.seqss[cx].len() {
                // ◼ Hideous workaround for the problem that a productive pair
                // ◼ could have two contigs with identical CDR3_AA sequences.
                // (but also because we now have some null seq entries?)
                if p >= rsi.seqss[cx][s].len() {
                    // if pass == 2 { fwriteln!( &mut mlog, "DIFFERENT LENGTHS" ); }
                    continue;
                }
                bases.push(rsi.seqss[cx][s][p]);
                bases_amino.push(rsi.seqss_amino[cx][s][p]);
            }
            unique_sort(&mut bases);
            unique_sort(&mut bases_amino);
            if bases.len() > 1 {
                v.push(p);
            } else if pass == 2 {
                let b = bases[0];
                if p < vref.len() - ctl.heur.ref_v_trim && b != vref[p] {
                    s.push(p);
                }
                if p >= n - (jref.len() - ctl.heur.ref_j_trim) && b != jref[jref.len() - (n - p)] {
                    s.push(p);
                }
            }
            if bases_amino.len() > 1 {
                v_amino.push(p);
            } else if pass == 2 {
                let b = bases_amino[0];
                if p < vseq2.len() - ctl.heur.ref_v_trim && b != vseq2[p] {
                    // if p < vref.len() - ctl.heur.ref_v_trim && b != vref[p] {
                    s_amino.push(p);
                }
                if p >= n - (jref.len() - ctl.heur.ref_j_trim) && b != jref[jref.len() - (n - p)] {
                    s_amino.push(p);
                }
            }
        }
        let mut va = Vec::<usize>::new();
        for x in v_amino.iter() {
            va.push(*x / 3);
        }
        unique_sort(&mut va);
        let mut sa = Vec::<usize>::new();
        for x in s_amino.iter() {
            sa.push(*x / 3);
        }
        unique_sort(&mut sa);
        for u in 0..nexacts {
            macro_rules! speakc {
                ($u:expr, $col:expr, $var:expr, $val:expr) => {
                    if ctl.parseable_opt.pout.len() > 0 && $col + 1 <= ctl.parseable_opt.pchains {
                        let varc = format!("{}{}", $var, $col + 1);
                        if pass == 2 && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc)) {
                            out_data[$u].insert(varc, $val);
                        }
                    }
                };
            }
            speakc![
                u,
                cx,
                "var_indices_dna",
                format!("{}", v.iter().format(","))
            ];
            speakc![
                u,
                cx,
                "share_indices_dna",
                format!("{}", s.iter().format(","))
            ];
            speakc![
                u,
                cx,
                "var_indices_aa",
                format!("{}", va.iter().format(","))
            ];
            speakc![
                u,
                cx,
                "share_indices_aa",
                format!("{}", sa.iter().format(","))
            ];
        }
        vars.push(v);
        vars_amino.push(v_amino);
        shares_amino.push(s_amino);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn delete_weaks(
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    mults: &Vec<usize>, // should eliminate
    exact_clonotypes: &Vec<ExactClonotype>,
    _total_cells: usize,
    mat: &Vec<Vec<Option<usize>>>,
    refdata: &RefData,
    vars: &Vec<Vec<usize>>,
    bads: &mut Vec<bool>,
) {
    // Mark for deletion exact subclonotypes that fail the MIN_CELLS_EXACT or MIN_CHAINS_EXACT
    // or CHAINS_EXACT tests.

    let nexacts = exacts.len();
    for u in 0..nexacts {
        if exact_clonotypes[exacts[u]].ncells() < ctl.gen_opt.min_cells_exact {
            bads[u] = true;
        }
        if exact_clonotypes[exacts[u]].share.len() < ctl.gen_opt.min_chains_exact {
            bads[u] = true;
        }
        if ctl.gen_opt.chains_exact > 0
            && exact_clonotypes[exacts[u]].share.len() != ctl.gen_opt.chains_exact
        {
            bads[u] = true;
        }
    }

    // Mark for deletion exact subclonotypes, based on CONST_IGH and CONST_IGKL
    // (see enclone help special).

    if ctl.gen_opt.const_igh.is_some() {
        for u in 0..nexacts {
            let mut ok = false;
            let ex = &exact_clonotypes[exacts[u]];
            for m in 0..ex.share.len() {
                if ex.share[m].left && ex.share[m].c_ref_id.is_some() {
                    let id = ex.share[m].c_ref_id.unwrap();
                    if ctl
                        .gen_opt
                        .const_igh
                        .as_ref()
                        .unwrap()
                        .is_match(&refdata.name[id])
                    {
                        ok = true;
                    }
                }
            }
            if !ok {
                bads[u] = true;
            }
        }
    }
    if ctl.gen_opt.const_igkl.is_some() {
        for u in 0..nexacts {
            let mut ok = false;
            let ex = &exact_clonotypes[exacts[u]];
            for m in 0..ex.share.len() {
                if !ex.share[m].left && ex.share[m].c_ref_id.is_some() {
                    let id = ex.share[m].c_ref_id.unwrap();
                    if ctl
                        .gen_opt
                        .const_igkl
                        .as_ref()
                        .unwrap()
                        .is_match(&refdata.name[id])
                    {
                        ok = true;
                    }
                }
            }
            if !ok {
                bads[u] = true;
            }
        }
    }

    // Find and mark for deletion exact subclonotypes having a variant base in V..J that,
    // accounting for all the cells in all the exact subclonotypes, never occurs as Q60
    // doesn't occur as Q40 twice, and disagrees with the reference.

    let cols = mat.len();
    // (column, pos, base, qual, row)
    if ctl.clono_filt_opt.qual_filter {
        let mut vquals = Vec::<(usize, usize, u8, u8, usize)>::new();
        for u in 0..nexacts {
            let clonotype_id = exacts[u];
            let ex = &exact_clonotypes[clonotype_id];
            for col in 0..cols {
                let m = mat[col][u];
                if m.is_some() {
                    let m = m.unwrap();
                    if ex.share[m].annv.len() > 1 {
                        continue;
                    }
                    let n = ex.share[m].seq_del.len();
                    let vref = &exact_clonotypes[exacts[u]].share[m].vs.to_ascii_vec();
                    let jref = &exact_clonotypes[exacts[u]].share[m].js.to_ascii_vec();
                    for z in 0..vars[col].len() {
                        let p = vars[col][z];
                        let b = ex.share[m].seq_del[p];
                        let mut refdiff = false;
                        if p < vref.len() - ctl.heur.ref_v_trim && b != vref[p] {
                            refdiff = true;
                        }
                        if p >= n - (jref.len() - ctl.heur.ref_j_trim)
                            && b != jref[jref.len() - (n - p)]
                        {
                            refdiff = true;
                        }
                        if refdiff {
                            for j in 0..ex.clones.len() {
                                let qual = ex.clones[j][m].quals[p];
                                vquals.push((col, p, b, qual, u));
                            }
                        }
                    }
                }
            }
        }
        vquals.sort();
        let mut j = 0;
        while j < vquals.len() {
            let mut k = j + 1;
            while k < vquals.len() {
                if vquals[k].0 != vquals[j].0
                    || vquals[k].1 != vquals[j].1
                    || vquals[k].2 != vquals[j].2
                {
                    break;
                }
                k += 1;
            }
            let mut q60 = false;
            let mut q40 = 0;
            for m in j..k {
                if vquals[m].3 >= 60 {
                    q60 = true;
                } else if vquals[m].3 >= 40 {
                    q40 += 1;
                }
            }
            if !q60 && q40 < 2 {
                let u = vquals[j].4;
                bads[u] = true;
            }
            j = k;
        }
    }

    // Based on the number of cells in each column, decide which exact subclonotypes
    // look like junk.  Preliminary heuristic.

    if cols > 2 && ctl.clono_filt_opt.weak_chains {
        let mut ncells = vec![0; cols];
        let mut col_entries = vec![Vec::<usize>::new(); cols];
        for u in 0..nexacts {
            for col in 0..cols {
                let mid = mat[col][u];
                if mid.is_some() {
                    col_entries[col].push(u);
                    let clonotype_id = exacts[u];
                    ncells[col] += exact_clonotypes[clonotype_id].clones.len();
                }
            }
        }
        let total_cells: usize = mults.iter().sum();
        for j in 0..cols {
            if ncells[j] <= 5 && 8 * ncells[j] < total_cells {
                for d in col_entries[j].iter() {
                    bads[*d] = true;
                }
            }
        }
    }
}
