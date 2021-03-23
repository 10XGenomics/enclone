// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_cvar.

use crate::print_utils1::*;
use amino::*;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use itertools::*;
use stats_utils::*;
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

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
    field_types: &Vec<Vec<u8>>,
    col_var: bool,
    pcols_sort: &Vec<String>,
    _bads: &mut Vec<bool>,
    cx: &mut Vec<Vec<String>>,
    _u_min: usize,
    _u_max: usize,
    _u_mean: usize,
    _median_numis: usize,
    _utot: usize,
    _median_nreads: usize,
    _r_min: usize,
    _r_max: usize,
    _r_mean: usize,
    _rtot: usize,
    extra_args: &Vec<String>,
) -> bool {
    let seq_amino = &rsi.seqss_amino[col][u];
    let mat = &rsi.mat;
    let cvars = &ctl.clono_print_opt.cvars;
    macro_rules! speakc {
        ($u:expr, $col:expr, $var:expr, $val:expr) => {
            if pass == 2
                && ((ctl.parseable_opt.pout.len() > 0 && $col + 1 <= ctl.parseable_opt.pchains)
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

    macro_rules! cvar {
        ($i: expr, $var:expr, $val:expr) => {
            if $i < rsi.cvars[col].len() && cvars.contains(&$var) {
                cx[col][$i] = $val.clone();
            }
            speakc!(u, col, $var, $val);
        };
    }

    if *var == "amino".to_string() && col_var {
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
                let log = color_codon(&ctl, &seq_amino, &x, p, &mut last_color, last);
                cx[col][j] += strme(&log);
            }
        }
    } else if *var == "comp".to_string() || *var == "edit".to_string() {
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

            let al = aligner.semiglobal(&tig, &concat);
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
                            count += 1;
                            edits.push(format!("D{}:{}", rpos, n - m));
                        }
                        pos += n - m;
                        m = n - 1;
                    }
                    Ins => {
                        if pos >= start && pos < stop {
                            count += 1;
                            edits.push(format!("I{}:{}", rpos, n - m));
                        }
                        rpos += (n - m) as isize;
                        m = n - 1;
                    }
                    _ => {}
                };
                m += 1;
            }
            if count < comp {
                comp = count;
                edit = format!("{}", edits.iter().format("•"));
            }
        }
        if *var == "comp".to_string() {
            cvar![j, var, format!("{}", comp)];
        } else {
            cvar![j, var, format!("{}", edit)];
        }
    } else if *var == "cdr1_dna".to_string()
        || *var == "cdr1_aa".to_string()
        || *var == "cdr1_aa_north".to_string()
        || *var == "cdr1_len".to_string()
        || (var.starts_with("cdr1_aa_") && var.ends_with("_ext"))
    {
        let (mut left, mut right) = (0, 0);
        if var.ends_with("_ext") {
            left = var.between("aa_", "_").force_i64() * 3;
            right = var.after("aa_").between("_", "_").force_i64() * 3;
        } else if var.ends_with("_north") {
            if ex.share[mid].left {
                left = 3 * 3;
                right = 3 * 3;
            }
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

                // Test for internal error.

                let mut found = false;
                for i in 0..x.seq.len() {
                    if x.seq[i..].starts_with(&dna) {
                        found = true;
                    }
                }
                if !found {
                    eprintln!(
                        "\nInternal error, failed to find {}, CDR3 = {}.\n",
                        strme(&dna),
                        x.cdr3_aa
                    );
                    std::process::exit(1);
                }
                if *var == "cdr1_dna".to_string() {
                    y = stringme(&dna);
                } else if var.starts_with("cdr1_aa") {
                    y = stringme(&aa_seq(&dna, 0));
                } else {
                    y = format!("{}", dna.len() / 3);
                }
            }
        }
        cvar![j, var, y];
    } else if *var == "cdr1_dna_ref".to_string() || *var == "cdr1_aa_ref".to_string() {
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.cdr1_start.is_some()
            && x.fr2_start.is_some()
            && x.cdr1_start.unwrap() <= x.fr2_start.unwrap()
        {
            let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                [x.cdr1_start.unwrap()..x.fr2_start.unwrap()]
                .to_vec();
            if *var == "cdr1_dna_ref".to_string() {
                y = stringme(&dna);
            } else {
                y = stringme(&aa_seq(&dna, 0));
            }
        }
        cvar![j, var, y];
    } else if *var == "cdr2_dna".to_string()
        || *var == "cdr2_aa".to_string()
        || *var == "cdr2_aa_north".to_string()
        || *var == "cdr2_len".to_string()
        || (var.starts_with("cdr2_aa_") && var.ends_with("_ext"))
    {
        let (mut left, mut right) = (0, 0);
        if var.ends_with("_ext") {
            left = var.between("aa_", "_").force_i64() * 3;
            right = var.after("aa_").between("_", "_").force_i64() * 3;
        } else if var.ends_with("_north") {
            if ex.share[mid].left {
                left = 2 * 3;
                right = 3 * 3;
            } else {
                left = 1 * 3;
            }
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

                // Test for internal error.

                let mut found = false;
                for i in 0..x.seq.len() {
                    if x.seq[i..].starts_with(&dna) {
                        found = true;
                    }
                }
                if !found {
                    eprintln!(
                        "\nInternal error, failed to find {}, CDR3 = {}.\n",
                        strme(&dna),
                        x.cdr3_aa
                    );
                    std::process::exit(1);
                }
                if *var == "cdr2_dna".to_string() {
                    y = stringme(&dna);
                } else if var.starts_with("cdr2_aa") {
                    y = stringme(&aa_seq(&dna, 0));
                } else {
                    y = format!("{}", dna.len() / 3);
                }
            }
        }
        cvar![j, var, y];
    } else if *var == "cdr2_dna_ref".to_string() || *var == "cdr2_aa_ref".to_string() {
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.cdr2_start.is_some()
            && x.fr3_start.is_some()
            && x.cdr2_start.unwrap() <= x.fr3_start.unwrap()
        {
            let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                [x.cdr2_start.unwrap()..x.fr3_start.unwrap()]
                .to_vec();
            if *var == "cdr2_dna_ref".to_string() {
                y = stringme(&dna);
            } else {
                y = stringme(&aa_seq(&dna, 0));
            }
        }
        cvar![j, var, y];
    } else if *var == "cdr3_aa".to_string() {
        cvar![j, var, ex.share[mid].cdr3_aa.clone()];
    } else if (var.starts_with("cdr3_aa_") && var.ends_with("_ext")) || var == "cdr3_aa_north" {
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

            // Test for internal error.

            let mut found = false;
            for i in 0..x.seq.len() {
                if x.seq[i..].starts_with(&dna) {
                    found = true;
                }
            }
            if !found {
                eprintln!(
                    "\nInternal error, failed to find {}, CDR3 = {}.\n",
                    strme(&dna),
                    x.cdr3_aa
                );
                std::process::exit(1);
            }
            y = stringme(&aa_seq(&dna, 0));
        }
        cvar![j, var, y];
    } else if *var == "cdr3_dna".to_string() {
        cvar![j, var, ex.share[mid].cdr3_dna.clone()];
    } else if *var == "cdr3_len".to_string() {
        cvar![j, var, ex.share[mid].cdr3_aa.len().to_string()];
    } else if *var == "cdr3_aa_conx".to_string() || *var == "cdr3_aa_conp".to_string() {
        let c;
        if *var == "cdr3_aa_conx" {
            c = cdr3_aa_con("x", col, &exacts, &exact_clonotypes, &rsi);
        } else {
            c = cdr3_aa_con("p", col, &exacts, &exact_clonotypes, &rsi);
        }
        cvar![j, var, c];
    } else if *var == "fwr1_dna".to_string()
        || *var == "fwr1_aa".to_string()
        || *var == "fwr1_len".to_string()
    {
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.cdr1_start.is_some() && x.fr1_start <= x.cdr1_start.unwrap() {
            let mut dna = Vec::<u8>::new();
            for p in x.fr1_start..x.cdr1_start.unwrap() {
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

            // Test for internal error.

            let mut found = false;
            for i in 0..x.seq.len() {
                if x.seq[i..].starts_with(&dna) {
                    found = true;
                }
            }
            if !found {
                eprintln!(
                    "\nInternal error, failed to find {}, CDR3 = {}.\n",
                    strme(&dna),
                    x.cdr3_aa
                );
                std::process::exit(1);
            }
            if *var == "fwr1_dna".to_string() {
                y = stringme(&dna);
            } else if *var == "fwr1_aa".to_string() {
                y = stringme(&aa_seq(&dna, 0));
            } else {
                y = format!("{}", dna.len() / 3);
            }
        }
        cvar![j, var, y];
    } else if *var == "fwr1_dna_ref".to_string() || *var == "fwr1_aa_ref".to_string() {
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.cdr1_start.is_some() && x.fr1_start <= x.cdr1_start.unwrap() {
            let dna = refdata.refs[x.v_ref_id].to_ascii_vec()[x.fr1_start..x.cdr1_start.unwrap()]
                .to_vec();
            if *var == "fwr1_dna_ref".to_string() {
                y = stringme(&dna);
            } else {
                y = stringme(&aa_seq(&dna, 0));
            }
        }
        cvar![j, var, y];
    } else if *var == "fwr2_dna".to_string()
        || *var == "fwr2_aa".to_string()
        || *var == "fwr2_len".to_string()
    {
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.fr2_start.unwrap() <= x.cdr2_start.unwrap() {
            let mut dna = Vec::<u8>::new();
            for p in x.fr2_start.unwrap()..x.cdr2_start.unwrap() {
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

            // Test for internal error.

            let mut found = false;
            for i in 0..x.seq.len() {
                if x.seq[i..].starts_with(&dna) {
                    found = true;
                }
            }
            if !found {
                eprintln!(
                    "\nInternal error, failed to find {}, CDR3 = {}.\n",
                    strme(&dna),
                    x.cdr3_aa
                );
                std::process::exit(1);
            }
            if *var == "fwr2_dna".to_string() {
                y = stringme(&dna);
            } else if *var == "fwr2_aa".to_string() {
                y = stringme(&aa_seq(&dna, 0));
            } else {
                y = format!("{}", dna.len() / 3);
            }
        }
        cvar![j, var, y];
    } else if *var == "fwr2_dna_ref".to_string() || *var == "fwr2_aa_ref".to_string() {
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.fr2_start.unwrap() <= x.cdr2_start.unwrap() {
            let dna = refdata.refs[x.v_ref_id].to_ascii_vec()
                [x.fr2_start.unwrap()..x.cdr2_start.unwrap()]
                .to_vec();
            if *var == "fwr2_dna_ref".to_string() {
                y = stringme(&dna);
            } else {
                y = stringme(&aa_seq(&dna, 0));
            }
        }
        cvar![j, var, y];
    } else if *var == "fwr3_dna".to_string()
        || *var == "fwr3_aa".to_string()
        || *var == "fwr3_len".to_string()
    {
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.fr3_start.is_some() && x.fr3_start.unwrap() <= x.cdr3_start - x.ins_len() {
            let mut dna = Vec::<u8>::new();
            for p in x.fr3_start.unwrap()..x.cdr3_start - x.ins_len() {
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

            // Test for internal error.

            let mut found = false;
            for i in 0..x.seq.len() {
                if x.seq[i..].starts_with(&dna) {
                    found = true;
                }
            }
            if !found {
                eprintln!(
                    "\nInternal error, failed to find {}, CDR3 = {}.\n",
                    strme(&dna),
                    x.cdr3_aa
                );
                std::process::exit(1);
            }
            if *var == "fwr3_dna".to_string() {
                y = stringme(&dna);
            } else if *var == "fwr3_aa".to_string() {
                y = stringme(&aa_seq(&dna, 0));
            } else {
                y = format!("{}", dna.len() / 3);
            }
        }
        cvar![j, var, y];
    } else if *var == "fwr3_dna_ref".to_string() || *var == "fwr3_aa_ref".to_string() {
        let x = &ex.share[mid];
        let mut y = "unknown".to_string();
        if x.fr3_start.is_some() && x.fr3_start.unwrap() <= x.cdr3_start - x.ins_len() {
            let dna = refdata.refs[x.v_ref_id].to_ascii_vec();
            if x.cdr3_start <= dna.len() {
                let dna = dna[x.fr3_start.unwrap()..x.cdr3_start - x.ins_len()].to_vec();
                if *var == "fwr3_dna_ref".to_string() {
                    y = stringme(&dna);
                } else {
                    y = stringme(&aa_seq(&dna, 0));
                }
            }
        }
        cvar![j, var, y];
    } else if *var == "fwr4_dna".to_string()
        || *var == "fwr4_aa".to_string()
        || *var == "fwr4_len".to_string()
    {
        let x = &ex.share[mid];
        let start = x.cdr3_start - x.ins_len() + 3 * x.cdr3_aa.len();
        let stop = x.seq_del_amino.len();
        let dna = &x.seq_del_amino[start..stop];
        let y;
        if *var == "fwr4_dna".to_string() {
            y = stringme(&dna);
        } else if *var == "fwr4_aa".to_string() {
            y = stringme(&aa_seq(&dna.to_vec(), 0));
        } else {
            y = format!("{}", dna.len() / 3);
        }
        cvar![j, var, y];
    } else if *var == "fwr4_dna_ref".to_string() || *var == "fwr4_aa_ref".to_string() {
        let x = &ex.share[mid];
        let heavy = refdata.rtype[x.j_ref_id] == 0;
        let aa_len;
        if heavy {
            aa_len = 10;
        } else {
            aa_len = 9;
        }
        let dna = refdata.refs[x.j_ref_id].to_ascii_vec();
        let dna = dna[dna.len() - 1 - 3 * aa_len..dna.len() - 1].to_vec();
        let y;
        if *var == "fwr4_dna_ref".to_string() {
            y = stringme(&dna);
        } else {
            y = stringme(&aa_seq(&dna.to_vec(), 0));
        }
        cvar![j, var, y];
    } else if *var == "ulen".to_string() {
        cvar![j, *var, format!("{}", ex.share[mid].v_start)];
    } else if *var == "clen".to_string() {
        cvar![
            j,
            var,
            format!("{}", ex.share[mid].full_seq.len() - ex.share[mid].j_stop)
        ];
    } else if *var == "aa%".to_string() {
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
        cvar![
            j,
            *var,
            format!("{:.1}", percent_ratio(denom - diffs, denom))
        ];
    } else if *var == "dna%".to_string() {
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
        cvar![
            j,
            *var,
            format!("{:.1}", percent_ratio(denom - diffs, denom))
        ];
    } else if *var == "vjlen".to_string() {
        cvar![
            j,
            var,
            format!("{}", ex.share[mid].j_stop - ex.share[mid].v_start)
        ];
    } else if var.starts_with("ndiff") {
        let u0 = var.between("ndiff", "vj").force_usize() - 1;
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
            cvar![j, *var, format!("{}", ndiff)];
        } else {
            cvar![j, *var, "_".to_string()];
        }
    } else {
        return false;
    }
    return true;
}
