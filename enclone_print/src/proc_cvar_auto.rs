// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
// This file is auto-generated by the crate enclone_vars, please do not edit.

use crate::print_utils1::*;
use amino::*;
use enclone_core::defs::*;
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

    macro_rules! cvar_stats1 {
        ($i: expr, $var:expr, $val:expr) => {
            if $i < rsi.cvars[col].len() && cvars.contains(&$var) {
                cx[col][$i] = $val.clone();
            }
            speakc!(u, col, $var, $val);
            let varc = format!("{}{}", $var, col + 1);
            stats.push((varc, vec![$val.to_string(); ex.ncells()]));
        };
    }

    // Test variable.

    let val = if false {
        String::new()
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
        format!("{:.1}", percent_ratio(denom - diffs, denom))
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
        cdiff
    } else if var == "cdr3_aa_conp" {
        cdr3_aa_con("p", col, exacts, exact_clonotypes, rsi)
    } else if var == "cdr3_aa_conx" {
        cdr3_aa_con("x", col, exacts, exact_clonotypes, rsi)
    } else if var == "cdr3_start" {
        ex.share[mid].cdr3_start.to_string()
    } else if var == "clen" {
        format!("{}", ex.share[mid].full_seq.len() - ex.share[mid].j_stop)
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
        format!("{}", constx.iter().format(","))
    } else if var == "const_id" {
        let mut const_id = String::new();
        if ex.share[mid].c_ref_id.is_some() {
            const_id = format!("{}", refdata.id[ex.share[mid].c_ref_id.unwrap()]);
        }
        const_id
    } else if var == "d1_name" {
        if !ex.share[mid].left {
            String::new()
        } else {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut opt = Vec::new();
            if !ds.is_empty() {
                opt = ds[0].clone();
            }
            let mut opt_name = String::new();
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
            opt_name
        }
    } else if var == "d1_score" {
        if !ex.share[mid].left {
            String::new()
        } else {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut delta = 0.0;
            if scores.len() > 1 {
                delta = scores[0] - scores[1];
            }
            format!("{:.1}", delta)
        }
    } else if var == "d2_name" {
        if !ex.share[mid].left {
            String::new()
        } else {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut opt2 = Vec::new();
            if ds.len() > 1 {
                opt2 = ds[1].clone();
            }
            let mut opt2_name = String::new();
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
            opt2_name
        }
    } else if var == "d2_score" {
        if !ex.share[mid].left {
            String::new()
        } else {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut score = 0.0;
            if scores.len() > 1 {
                score = scores[1];
            }
            format!("{:.1}", score)
        }
    } else if var == "d_delta" {
        if !ex.share[mid].left {
            String::new()
        } else {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut delta = 0.0;
            if scores.len() > 1 {
                delta = scores[0] - scores[1];
            }
            format!("{:.1}", delta)
        }
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
        format!("{}", diffs)
    } else if var == "d_frame" {
        let mut d_frame = String::new();
        if ex.share[mid].d_start.is_some() {
            d_frame = format!(
                "{}",
                (ex.share[mid].d_start.unwrap() - ex.share[mid].v_start) % 3
            );
        }
        d_frame
    } else if var == "d_id" {
        let did = if rsi.dids[col].is_some() {
            format!("{}", refdata.id[rsi.dids[col].unwrap()])
        } else {
            String::new()
        };
        did
    } else if var == "d_start" {
        let mut d_start = String::new();
        if ex.share[mid].d_start.is_some() {
            d_start = format!("{}", ex.share[mid].d_start.unwrap());
        }
        d_start
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
        format!("{}", diffs)
    } else if var == "d_Δ" {
        if !ex.share[mid].left {
            String::new()
        } else {
            let mut scores = Vec::<f64>::new();
            let mut ds = Vec::<Vec<usize>>::new();
            opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
            let mut delta = 0.0;
            if scores.len() > 1 {
                delta = scores[0] - scores[1];
            }
            format!("{:.1}", delta)
        }
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
        format!("{:.1}", percent_ratio(denom - diffs, denom))
    } else if var == "j_id" {
        format!("{}", refdata.id[rsi.jids[col]])
    } else if var == "j_name" {
        refdata.name[rsi.jids[col]].clone()
    } else if var == "notes" {
        ex.share[mid].vs_notesx.clone()
    } else if var == "u_max" {
        let mut numis = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
        }
        numis.sort();
        format!("{}", numis.iter().max().unwrap())
    } else if var == "v_id" {
        format!("{}", refdata.id[rsi.vids[col]])
    } else if var == "v_name" {
        refdata.name[rsi.vids[col]].clone()
    } else if var == "v_start" {
        format!("{}", ex.share[mid].v_start)
    } else if var == "vjlen" {
        format!("{}", ex.share[mid].j_stop - ex.share[mid].v_start)
    } else {
        "$UNDEFINED".to_string()
    };
    if val == "$UNDEFINED" {
        return Ok(false);
    } else {
        cvar_stats1![j, var, val];
        return Ok(true);
    }
}
