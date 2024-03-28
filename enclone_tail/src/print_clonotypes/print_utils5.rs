// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use amino::codon_to_aa;
use ansi_escape::emit_end_escape;
use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use enclone_core::print_tools::emit_codon_color_escape;
use enclone_proto::types::DonorReferenceItem;
use itertools::Itertools;
use std::collections::HashMap;
use string_utils::{strme, TextUtils};
use vdj_ann::refx::RefData;
use vector_utils::{bin_member, meet_size, unique_sort, VecUtils};

use super::print_utils1::aa_classes;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Find variant positions.  And some other things.

pub fn vars_and_shares(
    pass: usize,
    ctl: &EncloneControl,
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    rsi: &ColInfo,
    refdata: &RefData,
    dref: &[DonorReferenceItem],
    vars: &mut Vec<Vec<usize>>,
    vars_amino: &mut Vec<Vec<usize>>,
    shares_amino: &mut Vec<Vec<usize>>,
    ref_diff_pos: &mut Vec<Vec<Vec<usize>>>,
    out_data: &mut [HashMap<String, String>],
) {
    // Copied stuff.

    let pcols_sort = &ctl.parseable_opt.pcols_sort;
    let nexacts = exacts.len();
    let cols = rsi.vids.len();
    *ref_diff_pos = vec![vec![Vec::<usize>::new(); nexacts]; cols];

    // Go through the columns.

    for (cx, ((mat, (&vid, &vpid)), (ref_diff, (seqss, seqss_amino)))) in rsi
        .mat
        .iter()
        .zip(rsi.vids.iter().zip(rsi.vpids.iter()))
        .zip(
            ref_diff_pos
                .iter_mut()
                .zip(rsi.seqss.iter().zip(rsi.seqss_amino.iter())),
        )
        .take(cols)
        .enumerate()
    {
        // go through each column
        let (mut vref, mut jref) = (Vec::<u8>::new(), Vec::<u8>::new());
        let mut vseq2 = Vec::<u8>::new();
        for (&exact, &m) in exacts.iter().zip(mat.iter()) {
            if let Some(m) = m {
                // Reference assigned multiple times here, is wrong.
                // Also where using allelized reference, need to explain.
                // vref is supposed to be the donor reference, but seems like it isn't
                vref = exact_clonotypes[exact].share[m].vs.to_ascii_vec();
                jref = exact_clonotypes[exact].share[m].js.to_ascii_vec();
            }
            let vseq1 = refdata.refs[vid].to_ascii_vec();
            if let Some(vpid) = vpid {
                vseq2 = dref[vpid].nt_sequence.clone();
            } else {
                vseq2 = vseq1;
            }
        }

        for (&exact, (&m, ref_diff)) in exacts.iter().zip(mat.iter().zip(ref_diff.iter_mut())) {
            if let Some(m) = m {
                let seq = &exact_clonotypes[exact].share[m].seq_del_amino;
                let n = seq.len();
                for (p, &b) in seq.iter().enumerate() {
                    if (p < vref.len() - ctl.heur.ref_v_trim && b != vref[p])
                        || (p >= n - (jref.len() - ctl.heur.ref_j_trim)
                            && b != jref[jref.len() - (n - p)])
                    {
                        ref_diff.push(p);
                    }
                }
            }
        }

        let n = seqss.iter().map(std::vec::Vec::len).max().unwrap_or(0);
        let (mut v, mut s) = (Vec::<usize>::new(), Vec::<usize>::new());
        let (mut v_amino, mut s_amino) = (Vec::<usize>::new(), Vec::<usize>::new());
        for p in 0..n {
            let mut bases = Vec::<u8>::new();
            let mut bases_amino = Vec::<u8>::new();
            for (seqss, seqss_amino) in seqss.iter().zip(seqss_amino.iter()) {
                // ◼ Hideous workaround for the problem that a productive pair
                // ◼ could have two contigs with identical CDR3_AA sequences.
                // (but also because we now have some null seq entries?)
                if p >= seqss.len() {
                    // if pass == 2 { fwriteln!( &mut mlog, "DIFFERENT LENGTHS" ); }
                    continue;
                }
                bases.push(seqss[p]);
                bases_amino.push(seqss_amino[p]);
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
        let mut va = v_amino.iter().map(|&x| x / 3).collect::<Vec<_>>();
        unique_sort(&mut va);
        let mut sa = s_amino.iter().map(|&x| x / 3).collect::<Vec<_>>();
        unique_sort(&mut sa);
        for u in out_data.iter_mut().take(nexacts) {
            macro_rules! speakc {
                ($u:expr, $col:expr, $var:expr, $val:expr) => {
                    if ctl.parseable_opt.pout.len() > 0
                        && (ctl.parseable_opt.pchains == "max"
                            || cx < ctl.parseable_opt.pchains.force_usize())
                    {
                        let varc = format!("{}{}", $var, $col + 1);
                        if pass == 2 && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc)) {
                            $u.insert(varc, $val);
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
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    mat: &[Vec<Option<usize>>],
    refdata: &RefData,
    bads: &mut [bool],
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

    if ctl.clono_filt_opt.const_igh.is_some() {
        for u in 0..nexacts {
            let mut ok = false;
            let ex = &exact_clonotypes[exacts[u]];
            for m in 0..ex.share.len() {
                if ex.share[m].left && ex.share[m].c_ref_id.is_some() {
                    let id = ex.share[m].c_ref_id.unwrap();
                    if ctl
                        .clono_filt_opt
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
    if ctl.clono_filt_opt.const_igkl.is_some() {
        for u in 0..nexacts {
            let mut ok = false;
            let ex = &exact_clonotypes[exacts[u]];
            for m in 0..ex.share.len() {
                if !ex.share[m].left && ex.share[m].c_ref_id.is_some() {
                    let id = ex.share[m].c_ref_id.unwrap();
                    if ctl
                        .clono_filt_opt
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

    // Remove onesies that do not have an exact match.

    let cols = mat.len();
    if cols > 1 {
        for u1 in 0..nexacts {
            let ex1 = &exact_clonotypes[exacts[u1]];
            if ex1.share.len() == 1 && !bads[u1] {
                let mut perf = false;
                'u2: for u2 in 0..nexacts {
                    let ex2 = &exact_clonotypes[exacts[u2]];
                    if ex2.share.len() > 1 && !bads[u2] {
                        for i in 0..ex2.share.len() {
                            if ex1.share[0].seq == ex2.share[i].seq {
                                perf = true;
                                break 'u2;
                            }
                        }
                    }
                }
                if !perf {
                    bads[u1] = true;
                }
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Build the diff row.  This also has the variable names!

pub fn build_diff_row(
    ctl: &EncloneControl,
    rsi: &ColInfo,
    rows: &mut Vec<Vec<String>>,
    drows: &mut Vec<Vec<String>>,
    row1: &[String],
    nexacts: usize,
    field_types: &[Vec<u8>],
    show_aa: &[Vec<usize>],
) {
    let mat = &rsi.mat;
    let cols = mat.len();
    let diff_pos = rows.len();
    if !drows.is_empty() {
        let mut row = row1.to_owned();
        for col in 0..cols {
            for m in 0..rsi.cvars[col].len() {
                if rsi.cvars[col][m] == *"amino" {
                    let mut xdots = String::new();
                    for k in 0..show_aa[col].len() {
                        if k > 0
                            && field_types[col][k] != field_types[col][k - 1]
                            && !ctl.gen_opt.nospaces
                        {
                            xdots.push(' ');
                        }
                        let p = show_aa[col][k];
                        let q = 3 * p;
                        let leader = q < rsi.fr1_starts[col];
                        let mut cdr = false;
                        if rsi.cdr1_starts[col].is_some()
                            && rsi.cdr2_starts[col].is_some()
                            && rsi.fr2_starts[col].is_some()
                            && rsi.fr3_starts[col].is_some()
                            && q >= rsi.cdr1_starts[col].unwrap()
                            && q < rsi.fr2_starts[col].unwrap()
                        {
                            cdr = true;
                        }

                        if rsi.cdr2_starts[col].is_some()
                            && rsi.fr3_starts[col].is_some()
                            && q >= rsi.cdr2_starts[col].unwrap()
                            && q < rsi.fr3_starts[col].unwrap()
                        {
                            cdr = true;
                        }
                        if q >= rsi.cdr3_starts[col]
                            && q < rsi.cdr3_starts[col] + 3 * rsi.cdr3_lens[col]
                        {
                            cdr = true;
                        }
                        let mut codons = Vec::<Vec<u8>>::new();
                        for u in 0..nexacts {
                            if mat[col][u].is_some() {
                                let seq_amino = rsi.seqss_amino[col][u].clone();
                                if 3 * p + 3 <= seq_amino.len() {
                                    codons.push(seq_amino[3 * p..3 * p + 3].to_vec());
                                }
                            }
                        }
                        unique_sort(&mut codons);
                        if codons.len() > 1 {
                            if cdr {
                                if ctl.gen_opt.diff_style == *"C1" {
                                    xdots.push('C');
                                } else if ctl.gen_opt.diff_style == *"C2" {
                                    xdots.push('');
                                    xdots.push('[');
                                    xdots.push('0');
                                    xdots.push('1');
                                    xdots.push('m');
                                    xdots.push('');
                                    xdots.push('[');
                                    xdots.push('3');
                                    xdots.push('1');
                                    xdots.push('m');
                                    xdots.push('◼');
                                    xdots.push('');
                                    xdots.push('[');
                                    xdots.push('0');
                                    xdots.push('1');
                                    xdots.push('m');
                                    xdots.push('');
                                    xdots.push('[');
                                    xdots.push('3');
                                    xdots.push('0');
                                    xdots.push('m');
                                } else {
                                    xdots.push('x');
                                }
                            } else if !leader {
                                if ctl.gen_opt.diff_style == *"C1" {
                                    xdots.push('F');
                                } else if ctl.gen_opt.diff_style == *"C2" {
                                    xdots.push('▮');
                                } else {
                                    xdots.push('x');
                                }
                            } else if ctl.gen_opt.diff_style == *"C1" {
                                xdots.push('L');
                            } else if ctl.gen_opt.diff_style == *"C2" {
                                xdots.push('▮');
                            } else {
                                xdots.push('x');
                            }
                        } else {
                            xdots.push('.');
                        }
                    }
                    row.push(xdots);
                } else {
                    let mut v = rsi.cvars[col][m].clone();
                    if v.contains(':') {
                        v = v.before(":").to_string();
                    }
                    row.push(v);
                }
            }
            for r in &mut row {
                *r = format!("[01m{}[0m", *r);
            }
        }
        rows.push(row);
    } else {
        rows[diff_pos - 1][..row1.len()].clone_from_slice(row1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn insert_consensus_row(
    ctl: &EncloneControl,
    rsi: &ColInfo,
    nexacts: usize,
    field_types: &[Vec<u8>],
    show_aa: &[Vec<usize>],
    row1: &[String],
    rows: &mut Vec<Vec<String>>,
) {
    let mat = &rsi.mat;
    if ctl.clono_print_opt.conx || ctl.clono_print_opt.conp {
        let style = if ctl.clono_print_opt.conx { "x" } else { "p" };
        let mut row = vec!["consensus".to_string()];
        for _ in 1..row1.len() {
            row.push("\\ext".to_string());
        }
        let classes = aa_classes();
        for col in 0..rsi.mat.len() {
            for m in 0..rsi.cvars[col].len() {
                if rsi.cvars[col][m] == *"amino" {
                    let mut xdots = String::new();
                    for k in 0..show_aa[col].len() {
                        if k > 0
                            && field_types[col][k] != field_types[col][k - 1]
                            && !ctl.gen_opt.nospaces
                        {
                            xdots.push(' ');
                        }
                        let p = show_aa[col][k];
                        let mut codons = Vec::<Vec<u8>>::new();
                        for u in 0..nexacts {
                            if mat[col][u].is_some() {
                                let seq_amino = rsi.seqss_amino[col][u].clone();
                                if 3 * p + 3 <= seq_amino.len() {
                                    codons.push(seq_amino[3 * p..3 * p + 3].to_vec());
                                }
                            }
                        }
                        unique_sort(&mut codons);
                        let mut gap = false;
                        for x in &codons {
                            if x.contains(&b'-') {
                                gap = true;
                            }
                        }
                        if codons.solo() && gap {
                            xdots += "g";
                        } else if codons.solo() {
                            let codon = &codons[0];
                            let aa = codon_to_aa(codon);
                            let mut log = Vec::<u8>::new();
                            emit_codon_color_escape(codon, &mut log);
                            log.push(aa);
                            emit_end_escape(&mut log);
                            xdots += strme(&log);
                        } else if gap {
                            xdots += "X";
                        } else {
                            let mut aas = Vec::<u8>::new();
                            for x in &codons {
                                aas.push(codon_to_aa(x));
                            }
                            unique_sort(&mut aas);
                            if aas.solo() {
                                xdots.push(aas[0] as char);
                            } else if style == "x" {
                                xdots += "X";
                            } else {
                                for m in &classes {
                                    if meet_size(&aas, m.1) == aas.len() {
                                        xdots.push(m.0);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    row.push(xdots);
                } else {
                    row.push(String::new());
                }
            }
        }
        rows.push(row);
    }
}
