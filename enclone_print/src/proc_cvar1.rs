// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_cvar.

use crate::print_utils1::{color_codon, test_internal_error_seq};
use amino::aa_seq;
use bio_edit::alignment::pairwise::Aligner;
use bio_edit::alignment::AlignmentOperation::*;
use enclone_core::align_to_vdj_ref::{align_to_vdj_ref, cigar};
use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use enclone_proto::types::DonorReferenceItem;
use itertools::Itertools;
use std::collections::HashMap;
use string_utils::{stringme, strme, TextUtils};
use vdj_ann::refx::RefData;
use vector_utils::{bin_member, next_diff, sort_sync2};

pub fn proc_cvar1(
    var: &String,
    j: usize,
    col: usize,
    mid: usize,
    pass: usize,
    u: usize,
    ex: &ExactClonotype,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    refdata: &RefData,
    _varmat: &Vec<Vec<Vec<u8>>>,
    out_data: &mut Vec<HashMap<String, String>>,
    rsi: &ColInfo,
    dref: &Vec<DonorReferenceItem>,
    peer_groups: &Vec<Vec<(usize, u8, u32)>>,
    show_aa: &Vec<Vec<usize>>,
    ref_diff_pos: &Vec<Vec<Vec<usize>>>,
    field_types: &Vec<Vec<u8>>,
    col_var: bool,
    pcols_sort: &Vec<String>,
    _bads: &mut Vec<bool>,
    cx: &mut Vec<Vec<String>>,
    _median_numis: usize,
    _utot: usize,
    _median_nreads: usize,
    _r_min: usize,
    _r_max: usize,
    _r_mean: usize,
    _rtot: usize,
    extra_args: &Vec<String>,
    stats: &mut Vec<(String, Vec<String>)>,
    cdr3_con: &Vec<Vec<u8>>,
) -> Result<bool, String> {
    let seq_amino = &rsi.seqss_amino[col][u];
    let cvars = &ctl.clono_print_opt.cvars;
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

    // Set up chain variable macro.  This is the mechanism for generating
    // both human-readable and parseable output for chain variables.

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

    if var.starts_with('q')
        && var.ends_with('_')
        && var.after("q").rev_before("_").parse::<usize>().is_ok()
    {
        let n = var.between("q", "_").force_usize();
        let mut val = String::new();
        if n < ex.share[mid].seq.len() {
            let mut quals = Vec::<u8>::new();
            for j in 0..ex.clones.len() {
                quals.push(ex.clones[j][mid].quals[n]);
            }
            val = format!("{}", quals.iter().format(","));
        }
        cvar_stats1![j, var, val];
    } else if *var == "amino" && col_var {
        let mut last_color = "black".to_string();
        for k in 0..show_aa[col].len() {
            let p = show_aa[col][k];
            if k > 0 && field_types[col][k] != field_types[col][k - 1] {
                cx[col][j] += " ";
            }
            if 3 * p + 3 <= seq_amino.len()
                && seq_amino[3 * p..3 * p + 3].to_vec() == b"---".to_vec()
            {
                cx[col][j] += "-";
            } else if 3 * p + 3 > seq_amino.len() || seq_amino[3 * p..3 * p + 3].contains(&b'-') {
                cx[col][j] += "*";
            } else {
                let x = &peer_groups[rsi.vids[col]];
                let last = k == show_aa[col].len() - 1;
                let log = color_codon(
                    ctl,
                    seq_amino,
                    ref_diff_pos,
                    x,
                    col,
                    mid,
                    p,
                    u,
                    &mut last_color,
                    last,
                    cdr3_con,
                    exacts,
                    exact_clonotypes,
                );
                cx[col][j] += strme(&log);
            }
        }
    } else if *var == "cigar" {
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
        let c = cigar(&ops, 0, tig.len(), tig.len());
        cvar_stats1![j, var, c];
    } else if *var == "comp" || *var == "edit" {
        let mut comp = 1000000;
        let mut edit = String::new();
        let td = &ex.share[mid];
        let tig = &td.seq;
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = Aligner::new(-6, -1, &score);

        // Go through passes.  If IGH/TRB, we go through every D segment.  Otherwise
        // there is just one pass.

        let mut z = 1;
        if ex.share[mid].left {
            z = refdata.ds.len();
        }
        let mut ds = Vec::<usize>::new();
        let mut counts = Vec::<usize>::new();
        for di in 0..z {
            let mut d = 0;
            if ex.share[mid].left {
                d = refdata.ds[di];
            }

            // Start to build reference concatenation.  First append the V segment.

            let mut concat = Vec::<u8>::new();
            let mut vref = refdata.refs[rsi.vids[col]].to_ascii_vec();
            if rsi.vpids[col].is_none() {
            } else {
                vref = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
            }
            concat.append(&mut vref.clone());

            // Append the D segment if IGH/TRB.

            if ex.share[mid].left {
                let mut x = refdata.refs[d].to_ascii_vec();
                concat.append(&mut x);
            }

            // Append the J segment.

            let mut x = refdata.refs[rsi.jids[col]].to_ascii_vec();
            concat.append(&mut x);

            // Align the V..J sequence on the contig to the reference concatenation.

            let al = aligner.semiglobal(tig, &concat);
            let mut m = 0;
            let mut pos = al.xstart;
            let mut rpos = (al.ystart as isize) - (vref.len() as isize);
            let mut count = 0;
            let start = td.cdr3_start - td.ins_len();
            let stop = td.j_stop - td.v_start;
            let mut edits = Vec::<String>::new();
            while m < al.operations.len() {
                let n = next_diff(&al.operations, m);
                match al.operations[m] {
                    Match => {
                        pos += 1;
                        rpos += 1;
                    }
                    Subst => {
                        if pos >= start && pos < stop {
                            count += 1;
                            edits.push(format!("S{}", rpos));
                        }
                        pos += 1;
                        rpos += 1;
                    }
                    Del => {
                        if pos >= start && pos < stop {
                            if *var == "comp" || *var == "edit" {
                                count += 1;
                            } else {
                                count += n - m;
                            }
                            edits.push(format!("D{}:{}", rpos, n - m));
                        }
                        pos += n - m;
                        m = n - 1;
                    }
                    Ins => {
                        if pos >= start && pos < stop {
                            if *var == "comp" || *var == "edit" {
                                count += 1;
                            } else {
                                count += n - m;
                            }
                            edits.push(format!("I{}:{}", rpos, n - m));
                        }
                        rpos += (n - m) as isize;
                        m = n - 1;
                    }
                    _ => {}
                };
                m += 1;
            }
            counts.push(count);
            ds.push(d);
            if count < comp {
                comp = count;
                edit = format!("{}", edits.iter().format("•"));
            }
        }
        sort_sync2(&mut counts, &mut ds);
        let mut comp = 0;
        if !counts.is_empty() {
            comp = counts[0];
        }
        if *var == "comp".to_string() {
            cvar_stats1![j, var, format!("{}", comp)];
        } else {
            cvar_stats1![j, var, edit];
        }
    } else if var.starts_with("cdr1_aa_") && var.ends_with("_ext") {
        let (mut left, mut right) = (0, 0);
        if var.ends_with("_ext") {
            left = var.between("aa_", "_").force_i64() * 3;
            right = var.after("aa_").between("_", "_").force_i64() * 3;
        }
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.cdr1_start.is_some()
            && x.fr2_start.is_some()
            && x.cdr1_start.unwrap() <= x.fr2_start.unwrap()
        {
            let mut dna = Vec::<u8>::new();
            if x.cdr1_start.unwrap() as i64 - left >= 0
                && x.cdr1_start.unwrap() as i64 - left < x.seq_del_amino.len() as i64
                && x.fr2_start.unwrap() as i64 + right > 0
                && x.fr2_start.unwrap() as i64 + right <= x.seq_del_amino.len() as i64
            {
                for p in x.cdr1_start.unwrap() as i64 - left..x.fr2_start.unwrap() as i64 + right {
                    let p = p as usize;
                    for j in 0..x.ins.len() {
                        if x.ins[j].0 == p {
                            let mut z = x.ins[j].1.clone();
                            dna.append(&mut z);
                        }
                    }
                    if x.seq_del_amino[p] != b'-' {
                        dna.push(x.seq_del_amino[p]);
                    }
                }
                test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa)?;
                y = stringme(&aa_seq(&dna, 0));
            }
        }
        cvar_stats1![j, var, y];
    } else if var.starts_with("cdr2_aa_") && var.ends_with("_ext") {
        let (mut left, mut right) = (0, 0);
        if var.ends_with("_ext") {
            left = var.between("aa_", "_").force_i64() * 3;
            right = var.after("aa_").between("_", "_").force_i64() * 3;
        }
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.cdr2_start.is_some()
            && x.fr3_start.is_some()
            && x.cdr2_start.unwrap() <= x.fr3_start.unwrap()
        {
            let mut dna = Vec::<u8>::new();
            if x.cdr2_start.unwrap() as i64 - left >= 0
                && x.cdr2_start.unwrap() as i64 - left < x.seq_del_amino.len() as i64
                && x.fr3_start.unwrap() as i64 + right > 0
                && x.fr3_start.unwrap() as i64 + right <= x.seq_del_amino.len() as i64
            {
                for p in x.cdr2_start.unwrap() as i64 - left..x.fr3_start.unwrap() as i64 + right {
                    let p = p as usize;
                    for j in 0..x.ins.len() {
                        if x.ins[j].0 == p {
                            let mut z = x.ins[j].1.clone();
                            dna.append(&mut z);
                        }
                    }
                    if x.seq_del_amino[p] != b'-' {
                        dna.push(x.seq_del_amino[p]);
                    }
                }
                test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa)?;
                y = stringme(&aa_seq(&dna, 0));
            }
        }
        cvar_stats1![j, var, y];
    } else if var.starts_with("cdr3_aa_") && var.ends_with("_ext") {
        let mut left = -1 * 3;
        let mut right = -1 * 3;
        if var.ends_with("_ext") {
            left = var.between("aa_", "_").force_i64() * 3;
            right = var.after("aa_").between("_", "_").force_i64() * 3;
        }
        let x = &ex.share[mid];
        let mut dna = Vec::<u8>::new();
        let mut y = "unknown".to_string();
        if x.cdr3_start as i64 - left >= 0
            && x.cdr3_start as i64 - left < x.seq_del_amino.len() as i64
            && x.cdr3_start as i64 + 3 * x.cdr3_aa.len() as i64 + right > 0
            && x.cdr3_start as i64 + 3 * x.cdr3_aa.len() as i64 + right
                <= x.seq_del_amino.len() as i64
        {
            for p in
                x.cdr3_start as i64 - left..x.cdr3_start as i64 + 3 * x.cdr3_aa.len() as i64 + right
            {
                let p = p as usize;
                for j in 0..x.ins.len() {
                    if x.ins[j].0 == p {
                        let mut z = x.ins[j].1.clone();
                        dna.append(&mut z);
                    }
                }
                if x.seq_del_amino[p] != b'-' {
                    dna.push(x.seq_del_amino[p]);
                }
            }
            test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa)?;
            y = stringme(&aa_seq(&dna, 0));
        }
        cvar_stats1![j, var, y];
    } else {
        return Ok(false);
    }
    Ok(true)
}
