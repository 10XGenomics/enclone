// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
// This file is auto-generated by the crate enclone_vars, please do not edit.

use crate::print_utils1::*;
use crate::print_utils3::*;
use amino::*;
use enclone_core::align_to_vdj_ref::*;
use enclone_core::defs::*;
use enclone_core::median::*;
use enclone_core::opt_d::*;
use enclone_proto::types::*;
use itertools::Itertools;
use stats_utils::*;
use std::cmp::min;
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::RefData;
use vector_utils::*;

pub fn proc_cvar_auto(
    j: usize,
    pass: usize,
    var: &String,
    ex: &ExactClonotype,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    mid: usize,
    col: usize,
    u: usize,
    rsi: &ColInfo,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    ctl: &EncloneControl,
    extra_args: &Vec<String>,
    pcols_sort: &Vec<String>,
    cx: &mut Vec<Vec<String>>,
    varmat: &Vec<Vec<Vec<u8>>>,
    out_data: &mut Vec<HashMap<String, String>>,
    stats: &mut Vec<(String, Vec<String>)>,
) -> Result<bool, String> {
    let cvars = &ctl.clono_print_opt.cvars;

    // Set up macros.

    macro_rules! speakc {
        ($u:expr, $col:expr, $var:expr, $val:expr) => {
            if pass == 2
                && ((ctl.parseable_opt.pout.len() > 0
                    && (ctl.parseable_opt.pchains == "max"
                        || col < ctl.parseable_opt.pchains.force_usize()))
                    || extra_args.len() > 0)
            {
                let mut v = $var.clone();
                v = v.replace("_Σ", "_sum");
                v = v.replace("_μ", "_mean");

                // Strip escape character sequences from val.  Can happen in notes, maybe
                // other places.

                let mut val_clean = String::new();
                let mut chars = Vec::<char>::new();
                let valx = format!("{}", $val);
                for c in valx.chars() {
                    chars.push(c);
                }
                let mut escaped = false;
                for l in 0..chars.len() {
                    if chars[l] == '' {
                        escaped = true;
                    }
                    if escaped {
                        if chars[l] == 'm' {
                            escaped = false;
                        }
                        continue;
                    }
                    val_clean.push(chars[l]);
                }

                // Proceed.

                let varc = format!("{}{}", v, $col + 1);
                if pcols_sort.is_empty()
                    || bin_member(&pcols_sort, &varc)
                    || bin_member(&extra_args, &varc)
                {
                    out_data[$u].insert(varc, val_clean);
                }
            }
        };
    }

    // Test variable.

    let val = if false {
        (String::new(), Vec::<String>::new())
    } else if var == "aa%" {
        let xm = &ex.share[mid];
        let mut diffs = 0;
        let mut denom = 0;
        let aa_seq = &xm.aa_mod_indel;
        let mut vref = refdata.refs[xm.v_ref_id].to_ascii_vec();
        if xm.v_ref_id_donor_alt_id.is_some() {
            vref = dref[xm.v_ref_id_donor.unwrap()].nt_sequence.clone();
        }
        let jref = refdata.refs[xm.j_ref_id].to_ascii_vec();
        let z = 3 * aa_seq.len() + 1;
        for p in 0..aa_seq.len() {
            if aa_seq[p] == b'-' {
                diffs += 1;
                denom += 1;
                continue;
            }
            if 3 * p + 3 <= vref.len() - ctl.heur.ref_v_trim {
                denom += 1;
                if aa_seq[p] != codon_to_aa(&vref[3 * p..3 * p + 3]) {
                    diffs += 1;
                }
            }
            if 3 * p > z - (jref.len() - ctl.heur.ref_j_trim) + 3 {
                denom += 1;
                if aa_seq[p]
                    != codon_to_aa(&jref[jref.len() - (z - 3 * p)..jref.len() - (z - 3 * p) + 3])
                {
                    diffs += 1;
                }
            }
        }

        (
            format!("{:.1}", percent_ratio(denom - diffs, denom)),
            Vec::new(),
        )
    } else if var == "cdiff" {
        let cstart = ex.share[mid].j_stop;
        let clen = ex.share[mid].full_seq.len() - cstart;
        let cid = ex.share[mid].c_ref_id;
        let mut cdiff = String::new();
        let mut ndiffs = 0;
        if cid.is_some() {
            let r = &refdata.refs[cid.unwrap()];
            let mut extra = 0;
            if clen > r.len() {
                extra = clen - r.len();
            }
            for i in 0..min(clen, r.len()) {
                let tb = ex.share[mid].full_seq[cstart + i];
                let rb = r.to_ascii_vec()[i];
                if tb != rb {
                    ndiffs += 1;
                    if ndiffs <= 5 {
                        cdiff += &format!("{}{}", i, tb as char);
                    }
                }
            }
            if ndiffs > 5 {
                cdiff += "...";
            }
            if extra > 0 {
                cdiff += &format!("+{}", extra);
            }
        } else if clen > 0 {
            cdiff = format!("+{}", clen);
        }

        (cdiff, Vec::new())
    } else if var == "cdr3_aa_conp" {
        (
            cdr3_aa_con("p", col, exacts, exact_clonotypes, rsi),
            Vec::new(),
        )
    } else if var == "cdr3_aa_conx" {
        (
            cdr3_aa_con("x", col, exacts, exact_clonotypes, rsi),
            Vec::new(),
        )
    } else if var == "cdr3_start" {
        (ex.share[mid].cdr3_start.to_string(), Vec::new())
    } else if var.starts_with("cdr")
        && var.ends_with("_aa_ref")
        && var.between2("cdr", "_aa_ref").parse::<usize>().is_ok()
        && var.between2("cdr", "_aa_ref").force_i64() >= 1
        && var.between2("cdr", "_aa_ref").force_usize() <= 2
    {
        let arg1 = var.between2("cdr", "_aa_ref").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if arg1 == 1 {
            if x.cdr1_start.is_some()
                && x.fr2_start.is_some()
                && x.cdr1_start.unwrap() <= x.fr2_start.unwrap()
            {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                    [x.cdr1_start.unwrap()..x.fr2_start.unwrap()]
                    .to_vec();
                y = stringme(&aa_seq(&dna, 0));
            }
        } else {
            if x.cdr2_start.is_some()
                && x.fr3_start.is_some()
                && x.cdr2_start.unwrap() <= x.fr3_start.unwrap()
            {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                    [x.cdr2_start.unwrap()..x.fr3_start.unwrap()]
                    .to_vec();
                y = stringme(&aa_seq(&dna, 0));
            }
        }

        (y, Vec::new())
    } else if var.starts_with("cdr")
        && var.ends_with("_dna_ref")
        && var.between2("cdr", "_dna_ref").parse::<usize>().is_ok()
        && var.between2("cdr", "_dna_ref").force_i64() >= 1
        && var.between2("cdr", "_dna_ref").force_usize() <= 2
    {
        let arg1 = var.between2("cdr", "_dna_ref").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if arg1 == 1 {
            if x.cdr1_start.is_some()
                && x.fr2_start.is_some()
                && x.cdr1_start.unwrap() <= x.fr2_start.unwrap()
            {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                    [x.cdr1_start.unwrap()..x.fr2_start.unwrap()]
                    .to_vec();
                y = stringme(&dna);
            }
        } else {
            if x.cdr2_start.is_some()
                && x.fr3_start.is_some()
                && x.cdr2_start.unwrap() <= x.fr3_start.unwrap()
            {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                    [x.cdr2_start.unwrap()..x.fr3_start.unwrap()]
                    .to_vec();
                y = stringme(&dna);
            }
        }

        (y, Vec::new())
    } else if var.starts_with("cdr")
        && var.ends_with("_aa")
        && var.between2("cdr", "_aa").parse::<usize>().is_ok()
        && var.between2("cdr", "_aa").force_i64() >= 1
        && var.between2("cdr", "_aa").force_usize() <= 3
    {
        let arg1 = var.between2("cdr", "_aa").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        let c;
        if arg1 == 1 {
            c = get_cdr1(&x, 0, 0);
            if c.is_some() {
                y = stringme(&aa_seq(c.unwrap().as_bytes(), 0));
            }
        } else if arg1 == 2 {
            c = get_cdr2(&x, 0, 0);
            if c.is_some() {
                y = stringme(&aa_seq(c.unwrap().as_bytes(), 0));
            }
        } else {
            y = x.cdr3_aa.clone();
        }

        (y, Vec::new())
    } else if var.starts_with("cdr")
        && var.ends_with("_aa_north")
        && var.between2("cdr", "_aa_north").parse::<usize>().is_ok()
        && var.between2("cdr", "_aa_north").force_i64() >= 1
        && var.between2("cdr", "_aa_north").force_usize() <= 3
    {
        let arg1 = var.between2("cdr", "_aa_north").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        let c;
        if arg1 == 1 {
            let (mut left, mut right) = (0, 0);
            if x.left {
                left = 3;
                right = 3;
            }
            c = get_cdr1(&x, left, right);
            if c.is_some() {
                y = stringme(&aa_seq(c.unwrap().as_bytes(), 0));
            }
        } else if arg1 == 2 {
            let (left, right);
            if ex.share[mid].left {
                left = 2;
                right = 3;
            } else {
                left = 1;
                right = 0;
            }
            c = get_cdr2(&x, left, right);
            if c.is_some() {
                y = stringme(&aa_seq(c.unwrap().as_bytes(), 0));
            }
        } else {
            c = get_cdr3(&x, -1, -1);
            if c.is_some() {
                y = stringme(&aa_seq(c.unwrap().as_bytes(), 0));
            }
        }

        (y, Vec::new())
    } else if var.starts_with("cdr")
        && var.ends_with("_dna")
        && var.between2("cdr", "_dna").parse::<usize>().is_ok()
        && var.between2("cdr", "_dna").force_i64() >= 1
        && var.between2("cdr", "_dna").force_usize() <= 3
    {
        let arg1 = var.between2("cdr", "_dna").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        let c;
        if arg1 == 1 {
            c = get_cdr1(&x, 0, 0);
        } else if arg1 == 2 {
            c = get_cdr2(&x, 0, 0);
        } else {
            c = Some(x.cdr3_dna.clone());
        }
        if c.is_some() {
            y = c.unwrap();
        }

        (y, Vec::new())
    } else if var.starts_with("cdr")
        && var.ends_with("_len")
        && var.between2("cdr", "_len").parse::<usize>().is_ok()
        && var.between2("cdr", "_len").force_i64() >= 1
        && var.between2("cdr", "_len").force_usize() <= 3
    {
        let arg1 = var.between2("cdr", "_len").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        let c;
        if arg1 == 1 {
            c = get_cdr1(&x, 0, 0);
        } else if arg1 == 2 {
            c = get_cdr2(&x, 0, 0);
        } else {
            c = Some(x.cdr3_dna.clone());
        }
        if c.is_some() {
            y = format!("{}", c.unwrap().len() / 3);
        }

        (y, Vec::new())
    } else if var == "cigar" {
        let vref = refdata.refs[rsi.vids[col]].to_ascii_vec();
        let mut dref = Vec::<u8>::new();
        if rsi.dids[col].is_some() {
            dref = refdata.refs[rsi.dids[col].unwrap()].to_ascii_vec();
        }
        let d2ref = Vec::<u8>::new();
        let jref = refdata.refs[rsi.jids[col]].to_ascii_vec();
        let td = &ex.share[mid];
        let tig = &td.seq;
        let ops = align_to_vdj_ref(
            tig,
            &vref,
            &dref,
            &d2ref,
            &jref,
            "", // drefname
            ex.share[mid].left,
            ctl,
        )
        .0;

        (cigar(&ops, 0, tig.len(), tig.len()), Vec::new())
    } else if var == "clen" {
        (
            format!("{}", ex.share[mid].full_seq.len() - ex.share[mid].j_stop),
            Vec::new(),
        )
    } else if var == "comp" {
        let (comp, _edit) = comp_edit(&ex, mid, col, &refdata, &dref, &rsi);

        (format!("{}", comp), Vec::new())
    } else if var == "const" {
        let mut constx = Vec::<String>::new();
        let cid = ex.share[mid].c_ref_id;
        if cid.is_some() {
            constx.push(refdata.name[cid.unwrap()].clone());
        } else {
            constx.push("?".to_string());
        }
        unique_sort(&mut constx);
        // This is overcomplicated because there is now at most one
        // const entry per exact subclonotype.

        (format!("{}", constx.iter().format(",")), Vec::new())
    } else if var == "const_id" {
        let mut const_id = String::new();
        if ex.share[mid].c_ref_id.is_some() {
            const_id = format!("{}", refdata.id[ex.share[mid].c_ref_id.unwrap()]);
        }

        (const_id, Vec::new())
    } else if var == "d1_name" {
        let mut opt_name = String::new();
        if ex.share[mid].left {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut opt = Vec::new();
            if !ds.is_empty() {
                opt = ds[0].clone();
            }
            if opt.is_empty() {
                opt_name = "none".to_string();
            } else {
                for i in 0..opt.len() {
                    if i > 0 {
                        opt_name += ":";
                    }
                    opt_name += &refdata.name[opt[i]];
                }
            }
        }

        (opt_name, Vec::new())
    } else if var == "d1_score" {
        let mut score = String::new();
        if ex.share[mid].left {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut delta = 0.0;
            if scores.len() > 1 {
                delta = scores[0] - scores[1];
            }
            score = format!("{:.1}", delta)
        }

        (score, Vec::new())
    } else if var == "d2_name" {
        let mut opt2_name = String::new();
        if ex.share[mid].left {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut opt2 = Vec::new();
            if ds.len() > 1 {
                opt2 = ds[1].clone();
            }
            if opt2.is_empty() {
                opt2_name = "none".to_string();
            } else {
                for i in 0..opt2.len() {
                    if i > 0 {
                        opt2_name += ":";
                    }
                    opt2_name += &refdata.name[opt2[i]];
                }
            }
        }

        (opt2_name, Vec::new())
    } else if var == "d2_score" {
        let mut scorex = String::new();
        if ex.share[mid].left {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut score = 0.0;
            if scores.len() > 1 {
                score = scores[1];
            }
            scorex = format!("{:.1}", score)
        }

        (scorex, Vec::new())
    } else if var == "d_delta" {
        let mut del = String::new();
        if ex.share[mid].left {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut delta = 0.0;
            if scores.len() > 1 {
                delta = scores[0] - scores[1];
            }
            del = format!("{:.1}", delta)
        }

        (del, Vec::new())
    } else if var == "d_donor" {
        let vid = ex.share[mid].v_ref_id;
        let mut vref = refdata.refs[vid].to_ascii_vec();
        if rsi.vpids[col].is_some() {
            vref = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
        }
        let jid = ex.share[mid].j_ref_id;
        let jref = &refdata.refs[jid].to_ascii_vec();
        let tig = &ex.share[mid].seq_del;
        let n = tig.len();
        let mut diffs = 0;
        for p in 0..n {
            if tig[p] == b'-' {
                continue;
            }
            if p < vref.len() - ctl.heur.ref_v_trim && tig[p] != vref[p] {
                diffs += 1;
            } else if p >= n - (jref.len() - ctl.heur.ref_j_trim)
                && tig[p] != jref[jref.len() - (n - p)]
            {
                diffs += 1;
            }
        }

        (format!("{}", diffs), Vec::new())
    } else if var == "d_frame" {
        let mut d_frame = String::new();
        if ex.share[mid].d_start.is_some() {
            d_frame = format!(
                "{}",
                (ex.share[mid].d_start.unwrap() - ex.share[mid].v_start) % 3
            );
        }

        (d_frame, Vec::new())
    } else if var == "d_id" {
        let did = if rsi.dids[col].is_some() {
            format!("{}", refdata.id[rsi.dids[col].unwrap()])
        } else {
            String::new()
        };

        (did, Vec::new())
    } else if var == "d_name" {
        let dname = if rsi.dids[col].is_some() {
            refdata.name[rsi.dids[col].unwrap()].clone()
        } else {
            String::new()
        };

        (dname, Vec::new())
    } else if var == "d_start" {
        let mut d_start = String::new();
        if ex.share[mid].d_start.is_some() {
            d_start = format!("{}", ex.share[mid].d_start.unwrap());
        }

        (d_start, Vec::new())
    } else if var == "d_univ" {
        let vid = ex.share[mid].v_ref_id;
        let vref = &refdata.refs[vid].to_ascii_vec();
        let jid = ex.share[mid].j_ref_id;
        let jref = &refdata.refs[jid].to_ascii_vec();
        let tig = &ex.share[mid].seq_del;
        let n = tig.len();
        let mut diffs = 0;
        for p in 0..n {
            if tig[p] == b'-' {
                continue;
            }
            if p < vref.len() - ctl.heur.ref_v_trim && tig[p] != vref[p] {
                diffs += 1;
            } else if p >= n - (jref.len() - ctl.heur.ref_j_trim)
                && tig[p] != jref[jref.len() - (n - p)]
            {
                diffs += 1;
            }
        }

        (format!("{}", diffs), Vec::new())
    } else if var == "d_Δ" {
        let mut del = String::new();
        if ex.share[mid].left {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut delta = 0.0;
            if scores.len() > 1 {
                delta = scores[0] - scores[1];
            }
            del = format!("{:.1}", delta)
        }

        (del, Vec::new())
    } else if var == "dna%" {
        let xm = &ex.share[mid];
        let mut diffs = 0;
        let mut denom = 0;
        let seq = &xm.seq_del_amino;
        let mut vref = refdata.refs[xm.v_ref_id].to_ascii_vec();
        if xm.v_ref_id_donor_alt_id.is_some() {
            vref = dref[xm.v_ref_id_donor.unwrap()].nt_sequence.clone();
        }
        let jref = refdata.refs[xm.j_ref_id].to_ascii_vec();
        let z = seq.len();
        for p in 0..z {
            let b = seq[p];
            if b == b'-' {
                diffs += 1;
                denom += 1;
                continue;
            }
            if p < vref.len() - ctl.heur.ref_v_trim {
                denom += 1;
                if b != vref[p] {
                    diffs += 1;
                }
            }
            if p >= z - (jref.len() - ctl.heur.ref_j_trim) {
                denom += 1;
                if b != jref[jref.len() - (z - p)] {
                    diffs += 1;
                }
            }
        }

        (
            format!("{:.1}", percent_ratio(denom - diffs, denom)),
            Vec::new(),
        )
    } else if var == "edit" {
        let (_comp, edit) = comp_edit(&ex, mid, col, &refdata, &dref, &rsi);

        (edit, Vec::new())
    } else if var.starts_with("fwr")
        && var.ends_with("_aa")
        && var.between2("fwr", "_aa").parse::<usize>().is_ok()
        && var.between2("fwr", "_aa").force_i64() >= 1
        && var.between2("fwr", "_aa").force_usize() <= 4
    {
        let arg1 = var.between2("fwr", "_aa").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        let c;
        if arg1 == 1 {
            c = get_fwr1(&x);
        } else if arg1 == 2 {
            c = get_fwr2(&x);
        } else if arg1 == 3 {
            c = get_fwr3(&x);
        } else {
            let x = &ex.share[mid];
            let start = rsi.cdr3_starts[col] + 3 * rsi.cdr3_lens[col];
            let stop = rsi.seq_del_lens[col];
            let dna = &x.seq_del_amino[start..stop];
            c = Some(stringme(&dna));
        }
        if c.is_some() {
            y = stringme(&aa_seq(c.unwrap().as_bytes(), 0));
        }

        (y, Vec::new())
    } else if var.starts_with("fwr")
        && var.ends_with("_aa_ref")
        && var.between2("fwr", "_aa_ref").parse::<usize>().is_ok()
        && var.between2("fwr", "_aa_ref").force_i64() >= 1
        && var.between2("fwr", "_aa_ref").force_usize() <= 4
    {
        let arg1 = var.between2("fwr", "_aa_ref").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if arg1 == 1 {
            if x.cdr1_start.is_some() && x.fr1_start <= x.cdr1_start.unwrap() {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                    [x.fr1_start..x.cdr1_start.unwrap()]
                    .to_vec();
                y = stringme(&aa_seq(&dna, 0));
            }
        } else if arg1 == 2 {
            if x.fr2_start.unwrap() <= x.cdr2_start.unwrap() {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                    [x.fr2_start.unwrap()..x.cdr2_start.unwrap()]
                    .to_vec();
                y = stringme(&aa_seq(&dna, 0));
            }
        } else if arg1 == 3 {
            if x.fr3_start.is_some() && x.fr3_start.unwrap() <= x.cdr3_start - x.ins_len() {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec();
                if x.cdr3_start <= dna.len() {
                    let dna = dna[x.fr3_start.unwrap()..x.cdr3_start - x.ins_len()].to_vec();
                    y = stringme(&aa_seq(&dna, 0));
                }
            }
        } else {
            let heavy = refdata.rtype[x.j_ref_id] == 0;
            let aa_len;
            if heavy {
                aa_len = 10;
            } else {
                aa_len = 9;
            }
            let dna = refdata.refs[x.j_ref_id].to_ascii_vec();
            let dna = dna[dna.len() - 1 - 3 * aa_len..dna.len() - 1].to_vec();
            y = stringme(&aa_seq(&dna, 0));
        }

        (y, Vec::new())
    } else if var.starts_with("fwr")
        && var.ends_with("_dna")
        && var.between2("fwr", "_dna").parse::<usize>().is_ok()
        && var.between2("fwr", "_dna").force_i64() >= 1
        && var.between2("fwr", "_dna").force_usize() <= 4
    {
        let arg1 = var.between2("fwr", "_dna").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        let c;
        if arg1 == 1 {
            c = get_fwr1(&x);
        } else if arg1 == 2 {
            c = get_fwr2(&x);
        } else if arg1 == 3 {
            c = get_fwr3(&x);
        } else {
            let x = &ex.share[mid];
            let start = rsi.cdr3_starts[col] + 3 * rsi.cdr3_lens[col];
            let stop = rsi.seq_del_lens[col];
            let dna = &x.seq_del_amino[start..stop];
            c = Some(stringme(&dna));
        }
        if c.is_some() {
            y = c.unwrap();
        }

        (y, Vec::new())
    } else if var.starts_with("fwr")
        && var.ends_with("_dna_ref")
        && var.between2("fwr", "_dna_ref").parse::<usize>().is_ok()
        && var.between2("fwr", "_dna_ref").force_i64() >= 1
        && var.between2("fwr", "_dna_ref").force_usize() <= 4
    {
        let arg1 = var.between2("fwr", "_dna_ref").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if arg1 == 1 {
            if x.cdr1_start.is_some() && x.fr1_start <= x.cdr1_start.unwrap() {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                    [x.fr1_start..x.cdr1_start.unwrap()]
                    .to_vec();
                y = stringme(&dna);
            }
        } else if arg1 == 2 {
            if x.fr2_start.unwrap() <= x.cdr2_start.unwrap() {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                    [x.fr2_start.unwrap()..x.cdr2_start.unwrap()]
                    .to_vec();
                y = stringme(&dna);
            }
        } else if arg1 == 3 {
            if x.fr3_start.is_some() && x.fr3_start.unwrap() <= x.cdr3_start - x.ins_len() {
                let dna = refdata.refs[x.v_ref_id].to_ascii_vec();
                if x.cdr3_start <= dna.len() {
                    let dna = dna[x.fr3_start.unwrap()..x.cdr3_start - x.ins_len()].to_vec();
                    y = stringme(&dna);
                }
            }
        } else {
            let heavy = refdata.rtype[x.j_ref_id] == 0;
            let aa_len;
            if heavy {
                aa_len = 10;
            } else {
                aa_len = 9;
            }
            let dna = refdata.refs[x.j_ref_id].to_ascii_vec();
            let dna = dna[dna.len() - 1 - 3 * aa_len..dna.len() - 1].to_vec();
            y = stringme(&dna);
        }

        (y, Vec::new())
    } else if var.starts_with("fwr")
        && var.ends_with("_len")
        && var.between2("fwr", "_len").parse::<usize>().is_ok()
        && var.between2("fwr", "_len").force_i64() >= 1
        && var.between2("fwr", "_len").force_usize() <= 4
    {
        let arg1 = var.between2("fwr", "_len").force_usize();
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        let c;
        if arg1 == 1 {
            c = get_fwr1(&x);
        } else if arg1 == 2 {
            c = get_fwr2(&x);
        } else if arg1 == 3 {
            c = get_fwr3(&x);
        } else {
            let x = &ex.share[mid];
            let start = rsi.cdr3_starts[col] + 3 * rsi.cdr3_lens[col];
            let stop = rsi.seq_del_lens[col];
            let dna = &x.seq_del_amino[start..stop];
            c = Some(stringme(&dna));
        }
        if c.is_some() {
            y = format!("{}", c.unwrap().len() / 3);
        }

        (y, Vec::new())
    } else if var == "j_id" {
        (format!("{}", refdata.id[rsi.jids[col]]), Vec::new())
    } else if var == "j_name" {
        (refdata.name[rsi.jids[col]].clone(), Vec::new())
    } else if var.starts_with("ndiff")
        && var.ends_with("vj")
        && var.between2("ndiff", "vj").parse::<usize>().is_ok()
        && var.between2("ndiff", "vj").force_i64() >= 1
    {
        let arg1 = var.between2("ndiff", "vj").force_usize();
        let nd;
        let mat = &rsi.mat;
        let u0 = arg1 - 1;
        if u0 < exacts.len() && mat[col][u0].is_some() && mat[col][u].is_some() {
            let m0 = mat[col][u0].unwrap();
            let m = mat[col][u].unwrap();
            let mut ndiff = 0;
            let ex0 = &exact_clonotypes[exacts[u0]];
            let ex = &exact_clonotypes[exacts[u]];
            for p in 0..ex0.share[m0].seq_del.len() {
                if ex0.share[m0].seq_del[p] != ex.share[m].seq_del[p] {
                    ndiff += 1;
                }
            }
            nd = format!("{}", ndiff)
        } else {
            nd = "_".to_string()
        }

        (nd, Vec::new())
    } else if var == "notes" {
        (ex.share[mid].vs_notesx.clone(), Vec::new())
    } else if var.starts_with("q")
        && var.ends_with("_")
        && var.between2("q", "_").parse::<usize>().is_ok()
        && var.between2("q", "_").force_i64() >= 0
    {
        let arg1 = var.between2("q", "_").force_usize();
        let mut val = String::new();
        if arg1 < ex.share[mid].seq.len() {
            let mut quals = Vec::<u8>::new();
            for j in 0..ex.clones.len() {
                quals.push(ex.clones[j][mid].quals[arg1]);
            }
            val = format!("{}", quals.iter().format(","));
        }

        (val, Vec::new())
    } else if var == "r" {
        let mut nreads = Vec::<String>::new();
        let mut nreads_sorted = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            nreads.push(format!("{}", ex.clones[j][mid].read_count));
            nreads_sorted.push(ex.clones[j][mid].read_count);
        }
        nreads_sorted.sort_unstable();

        (format!("{}", rounded_median(&nreads_sorted)), nreads)
    } else if var == "r_cell" {
        let mut nreads = Vec::<String>::new();
        let mut nreads_sorted = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            nreads.push(format!("{}", ex.clones[j][mid].read_count));
            nreads_sorted.push(ex.clones[j][mid].read_count);
        }
        nreads_sorted.sort_unstable();

        let _exact = format!("{}", rounded_median(&nreads_sorted));
        (String::new(), nreads)
    } else if var == "r_max" {
        let mut nreads = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            nreads.push(ex.clones[j][mid].read_count);
        }

        (format!("{}", *nreads.iter().max().unwrap()), Vec::new())
    } else if var == "r_mean" {
        let mut nreads = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            nreads.push(ex.clones[j][mid].read_count);
        }
        let rtot: usize = nreads.iter().sum();
        let r_mean = (rtot as f64 / nreads.len() as f64).round() as usize;

        (format!("{}", r_mean), Vec::new())
    } else if var == "r_min" {
        let mut nreads = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            nreads.push(ex.clones[j][mid].read_count);
        }

        (format!("{}", *nreads.iter().min().unwrap()), Vec::new())
    } else if var == "r_sum" {
        let mut nreads = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            nreads.push(ex.clones[j][mid].read_count);
        }
        let rtot: usize = nreads.iter().sum();

        (format!("{}", rtot), Vec::new())
    } else if var == "r_Σ" {
        let mut nreads = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            nreads.push(ex.clones[j][mid].read_count);
        }
        let rtot: usize = nreads.iter().sum();

        (format!("{}", rtot), Vec::new())
    } else if var == "r_μ" {
        let mut nreads = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            nreads.push(ex.clones[j][mid].read_count);
        }
        let rtot: usize = nreads.iter().sum();
        let r_mean = (rtot as f64 / nreads.len() as f64).round() as usize;

        (format!("{}", r_mean), Vec::new())
    } else if var == "u" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        numis.sort_unstable();
        let median_numis = rounded_median(&numis);
        let mut vals = Vec::<String>::new();
        for k in 0..ex.ncells() {
            vals.push(format!("{}", ex.clones[k][mid].umi_count));
        }

        (format!("{}", median_numis), vals)
    } else if var == "u_cell" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        numis.sort_unstable();
        let median_numis = rounded_median(&numis);
        let mut vals = Vec::<String>::new();
        for k in 0..ex.ncells() {
            vals.push(format!("{}", ex.clones[k][mid].umi_count));
        }

        let _exact = format!("{}", median_numis);
        (String::new(), vals)
    } else if var == "u_max" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        numis.sort();

        (format!("{}", numis.iter().max().unwrap()), Vec::new())
    } else if var == "u_mean" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        numis.sort();
        let utot: usize = numis.iter().sum();
        let u_mean = (utot as f64 / numis.len() as f64).round() as usize;

        (format!("{}", u_mean), Vec::new())
    } else if var == "u_min" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        numis.sort();

        (format!("{}", numis.iter().min().unwrap()), Vec::new())
    } else if var == "u_sum" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        let utot: usize = numis.iter().sum();

        (format!("{}", utot), Vec::new())
    } else if var == "u_Σ" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        let utot: usize = numis.iter().sum();

        (format!("{}", utot), Vec::new())
    } else if var == "u_μ" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        numis.sort();
        let utot: usize = numis.iter().sum();
        let u_mean = (utot as f64 / numis.len() as f64).round() as usize;

        (format!("{}", u_mean), Vec::new())
    } else if var == "udiff" {
        let ulen = ex.share[mid].v_start;
        let uid = ex.share[mid].u_ref_id;
        let mut udiff = String::new();
        let mut ndiffs = 0;
        if uid.is_some() {
            let r = &refdata.refs[uid.unwrap()];
            let mut extra = 0;
            if ulen > r.len() {
                extra = ulen - r.len();
            }
            for i in 0..ulen {
                let mut rpos = i;
                if ulen < r.len() {
                    rpos += r.len() - ulen;
                } else {
                    if i + r.len() < ulen {
                        continue;
                    }
                    rpos -= ulen - r.len();
                }
                let tb = ex.share[mid].full_seq[i];
                let rb = r.to_ascii_vec()[rpos];
                if tb != rb {
                    ndiffs += 1;
                    if ndiffs <= 5 {
                        udiff += &format!("{}{}", rpos, tb as char);
                    }
                }
            }
            if ndiffs > 5 {
                udiff += "...";
            }
            if extra > 0 {
                udiff += &format!("+{}", extra);
            }
        } else if ulen > 0 {
            udiff = format!("+{}", ulen);
        }

        (udiff, Vec::new())
    } else if var == "ulen" {
        (format!("{}", ex.share[mid].v_start), Vec::new())
    } else if var == "utr_id" {
        let mut u = String::new();
        let uid = ex.share[mid].u_ref_id;
        if uid.is_some() {
            u = format!("{}", refdata.id[uid.unwrap()]);
        }

        (u, Vec::new())
    } else if var == "utr_name" {
        let mut u = String::new();
        let uid = ex.share[mid].u_ref_id;
        if uid.is_some() {
            u = refdata.name[uid.unwrap()].clone();
        }

        (u, Vec::new())
    } else if var == "v_id" {
        (format!("{}", refdata.id[rsi.vids[col]]), Vec::new())
    } else if var == "v_name" {
        (refdata.name[rsi.vids[col]].clone(), Vec::new())
    } else if var == "v_start" {
        (format!("{}", ex.share[mid].v_start), Vec::new())
    } else if var == "var" {
        (stringme(&varmat[u][col]), Vec::new())
    } else if var == "vjlen" {
        (
            format!("{}", ex.share[mid].j_stop - ex.share[mid].v_start),
            Vec::new(),
        )
    } else {
        ("$UNDEFINED".to_string(), Vec::<String>::new())
    };
    if val.0 == "$UNDEFINED" {
        return Ok(false);
    } else {
        let (exact, cell) = &val;
        let varc = format!("{}{}", var, col + 1);
        if exact.len() > 0 {
            if j < rsi.cvars[col].len() && cvars.contains(&var) {
                cx[col][j] = exact.clone();
            }
            speakc!(u, col, var, exact);
            if val.1.is_empty() {
                stats.push((varc, vec![exact.to_string(); ex.ncells()]));
            } else {
                stats.push((varc, cell.to_vec()));
            }
        } else if cell.len() > 0 {
            if pass == 2
                && ((ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize())
                    || !extra_args.is_empty())
            {
                if pcols_sort.is_empty() || bin_member(pcols_sort, &varc) {
                    let vals = format!("{}", cell.iter().format(&POUT_SEP));
                    out_data[u].insert(varc, vals);
                }
            }
        }
        return Ok(true);
    }
}
