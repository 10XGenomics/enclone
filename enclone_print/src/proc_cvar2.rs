// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_cvar.

use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype, POUT_SEP};
use itertools::Itertools;
use stats_utils::percent_ratio;
use std::collections::HashMap;
use string_utils::TextUtils;
use vector_utils::{bin_member, next_diff12_4};

pub fn proc_cvar2(
    var: &String,
    j: usize,
    col: usize,
    mid: usize,
    pass: usize,
    u: usize,
    ex: &ExactClonotype,
    ctl: &EncloneControl,
    _exacts: &Vec<usize>,
    _exact_clonotypes: &Vec<ExactClonotype>,
    _varmat: &Vec<Vec<Vec<u8>>>,
    out_data: &mut Vec<HashMap<String, String>>,
    rsi: &ColInfo,
    _peer_groups: &Vec<Vec<(usize, u8, u32)>>,
    _show_aa: &Vec<Vec<usize>>,
    _field_types: &Vec<Vec<u8>>,
    col_var: bool,
    pcols_sort: &Vec<String>,
    bads: &mut Vec<bool>,
    cx: &mut Vec<Vec<String>>,
    _median_numis: usize,
    median_nreads: usize,
    r_min: usize,
    r_max: usize,
    r_mean: usize,
    rtot: usize,
    extra_args: &Vec<String>,
    stats: &mut Vec<(String, Vec<String>)>,
) -> bool {
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
            stats.push((varc, vec![$val; ex.ncells()]));
        };
    }
    macro_rules! cvar_stats {
        ($i: expr, $var:expr, $val:expr, $stats:expr) => {
            if $i < rsi.cvars[col].len() && cvars.contains(&$var) {
                cx[col][$i] = $val.clone();
            }
            speakc!(u, col, $var, $val);
            let varc = format!("{}{}", $var, col + 1);
            stats.push((varc, $stats.clone()));
        };
    }

    // Proceed.

    if var == "nval" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut vals = String::new();
                let mut valsx = Vec::<String>::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        vals += POUT_SEP;
                    }
                    let mut n = 0;
                    if ex.clones[k][mid].validated_umis.is_some() {
                        n = ex.clones[k][mid].validated_umis.as_ref().unwrap().len();
                    }
                    vals += &format!("{}", n);
                    valsx.push(format!("{}", n));
                }
                out_data[u].insert(varc.clone(), vals.to_string());
                stats.push((varc, valsx));
            }
        }
    } else if var == "nnval" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut nvals = String::new();
                let mut nvalsx = Vec::<String>::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        nvals += POUT_SEP;
                    }
                    let mut n = 0;
                    if ex.clones[k][mid].non_validated_umis.is_some() {
                        n = ex.clones[k][mid].non_validated_umis.as_ref().unwrap().len();
                    }
                    nvals += &format!("{}", n);
                    nvalsx.push(format!("{}", n));
                }
                out_data[u].insert(varc.clone(), nvals.to_string());
                stats.push((varc, nvalsx));
            }
        }
    } else if var == "nival" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut nvals = String::new();
                let mut nvalsx = Vec::<String>::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        nvals += POUT_SEP;
                    }
                    let mut n = 0;
                    if ex.clones[k][mid].invalidated_umis.is_some() {
                        n = ex.clones[k][mid].invalidated_umis.as_ref().unwrap().len();
                    }
                    nvals += &format!("{}", n);
                    nvalsx.push(format!("{}", n));
                }
                out_data[u].insert(varc.clone(), nvals.to_string());
                stats.push((varc, nvalsx));
            }
        }
    } else if var == "valumis" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut vals = String::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        vals += POUT_SEP;
                    }
                    let mut n = String::new();
                    if ex.clones[k][mid].validated_umis.is_some() {
                        n = format!(
                            "{}",
                            ex.clones[k][mid]
                                .validated_umis
                                .as_ref()
                                .unwrap()
                                .iter()
                                .format(",")
                        );
                    }
                    vals += &n.to_string();
                }
                out_data[u].insert(varc, vals.to_string());
            }
        }
    } else if var == "valbcumis" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut vals = String::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        vals += POUT_SEP;
                    }
                    let mut n = String::new();
                    if ex.clones[k][mid].validated_umis.is_some() {
                        let mut bc_umis = ex.clones[k][mid].validated_umis.clone().unwrap();
                        for i in 0..bc_umis.len() {
                            bc_umis[i] =
                                format!("{}{}", ex.clones[k][mid].barcode.before("-"), bc_umis[i]);
                        }
                        n = format!("{}", bc_umis.iter().format(","));
                    }
                    vals += &n.to_string();
                }
                out_data[u].insert(varc, vals.to_string());
            }
        }
    } else if var == "nvalumis" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut nvals = String::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        nvals += POUT_SEP;
                    }
                    let mut n = String::new();
                    if ex.clones[k][mid].non_validated_umis.is_some() {
                        n = format!(
                            "{}",
                            ex.clones[k][mid]
                                .non_validated_umis
                                .as_ref()
                                .unwrap()
                                .iter()
                                .format(",")
                        );
                    }
                    nvals += &n.to_string();
                }
                out_data[u].insert(varc, nvals.to_string());
            }
        }
    } else if var == "nvalbcumis" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut vals = String::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        vals += POUT_SEP;
                    }
                    let mut n = String::new();
                    if ex.clones[k][mid].non_validated_umis.is_some() {
                        let mut bc_umis = ex.clones[k][mid].non_validated_umis.clone().unwrap();
                        for i in 0..bc_umis.len() {
                            bc_umis[i] =
                                format!("{}{}", ex.clones[k][mid].barcode.before("-"), bc_umis[i]);
                        }
                        n = format!("{}", bc_umis.iter().format(","));
                    }
                    vals += &n.to_string();
                }
                out_data[u].insert(varc, vals.to_string());
            }
        }
    } else if var == "ivalumis" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut nvals = String::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        nvals += POUT_SEP;
                    }
                    let mut n = String::new();
                    if ex.clones[k][mid].invalidated_umis.is_some() {
                        n = format!(
                            "{}",
                            ex.clones[k][mid]
                                .invalidated_umis
                                .as_ref()
                                .unwrap()
                                .iter()
                                .format(",")
                        );
                    }
                    nvals += &n.to_string();
                }
                out_data[u].insert(varc, nvals.to_string());
            }
        }
    } else if var == "ivalbcumis" {
        cvar_stats1![j, *var, "".to_string()];
        if pass == 2
            && ((!ctl.parseable_opt.pout.is_empty()
                && (ctl.parseable_opt.pchains == "max"
                    || col < ctl.parseable_opt.pchains.force_usize()))
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut vals = String::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        vals += POUT_SEP;
                    }
                    let mut n = String::new();
                    if ex.clones[k][mid].invalidated_umis.is_some() {
                        let mut bc_umis = ex.clones[k][mid].invalidated_umis.clone().unwrap();
                        for i in 0..bc_umis.len() {
                            bc_umis[i] =
                                format!("{}{}", ex.clones[k][mid].barcode.before("-"), bc_umis[i]);
                        }
                        n = format!("{}", bc_umis.iter().format(","));
                    }
                    vals += &n.to_string();
                }
                out_data[u].insert(varc, vals.to_string());
            }
        }
    } else if *var == "u_cell" {
        let var = var.clone();
        if pass == 2
            && ((ctl.parseable_opt.pchains == "max"
                || col < ctl.parseable_opt.pchains.force_usize())
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty() || bin_member(pcols_sort, &varc) {
                let mut vals = String::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        vals += POUT_SEP;
                    }
                    vals += &format!("{}", ex.clones[k][mid].umi_count);
                }
                out_data[u].insert(varc, vals.to_string());
            }
        }
    } else if *var == "r" {
        let mut nreads = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            nreads.push(format!("{}", ex.clones[j][mid].read_count));
        }
        cvar_stats![j, var, format!("{}", median_nreads), nreads];
    } else if *var == "r_min" {
        cvar_stats1![j, var, format!("{}", r_min)];
    } else if *var == "r_max" {
        cvar_stats1![j, var, format!("{}", r_max)];
    } else if *var == "r_μ" {
        cvar_stats1![j, var, format!("{}", r_mean)];
    } else if *var == "r_Σ" {
        cvar_stats1![j, var, format!("{}", rtot)];
    } else if *var == "r_cell" {
        let var = var.clone();
        if pass == 2
            && ((ctl.parseable_opt.pchains == "max"
                || col < ctl.parseable_opt.pchains.force_usize())
                || !extra_args.is_empty())
        {
            let varc = format!("{}{}", var, col + 1);
            if pcols_sort.is_empty()
                || bin_member(pcols_sort, &varc)
                || bin_member(extra_args, &varc)
            {
                let mut vals = String::new();
                for k in 0..ex.ncells() {
                    if k > 0 {
                        vals += POUT_SEP;
                    }
                    vals += &format!("{}", ex.clones[k][mid].read_count);
                }
                out_data[u].insert(varc, vals.to_string());
            }
        }

    // Compute potential whitelist contamination percent and filter.
    // This is an undocumented option.
    } else if *var == "white" || ctl.clono_filt_opt_def.whitef {
        let mut bch = vec![Vec::<(usize, String, usize, usize)>::new(); 2];
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = &ex.clones[l][0].barcode;
            let mut numi = 0;
            for j in 0..ex.clones[l].len() {
                numi += ex.clones[l][j].umi_count;
            }
            bch[0].push((li, bc[0..8].to_string(), numi, l));
            bch[1].push((li, bc[8..16].to_string(), numi, l));
        }
        let mut junk = 0;
        let mut bad = vec![false; ex.clones.len()];
        for l in 0..2 {
            bch[l].sort();
            let mut m = 0;
            while m < bch[l].len() {
                let n = next_diff12_4(&bch[l], m as i32) as usize;
                for u1 in m..n {
                    for u2 in m..n {
                        if bch[l][u1].2 >= 10 * bch[l][u2].2 {
                            bad[bch[l][u2].3] = true;
                        }
                    }
                }
                m = n;
            }
        }
        for u in 0..bad.len() {
            if bad[u] {
                junk += 1;
            }
        }
        // Don't look at very large clones because of course they
        // show overlap.
        /* // BROKEN AND WAS UGLY ANYWAY
        const MAX_WHITELIST_CLONE: usize = 100;
        if ex.clones.len() <= MAX_WHITELIST_CLONE {
            res.3 += junk;
            res.4 += ex.clones.len();
        }
        */
        let junk_rate = percent_ratio(junk, ex.clones.len());
        if *var == "white".to_string() && col_var {
            cx[col][j] = format!("{:.1}", junk_rate);
        }
        // WRONG!  THIS IS SUPPOSED TO BE EXECUTED ON PASS 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if ctl.clono_filt_opt_def.whitef && junk_rate == 0.0
        /* && pass == 1 */
        {
            bads[u] = true;
        }
    } else {
        return false;
    }
    true
}
