// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_lvar,
// plus a small helper function get_gex_matrix_entry.

use crate::print_utils4::get_gex_matrix_entry;
use enclone_core::defs::{EncloneControl, ExactClonotype, GexInfo, POUT_SEP};
use enclone_core::median::rounded_median;
use itertools::Itertools;
use std::collections::HashMap;
use string_utils::TextUtils;
use vector_utils::{bin_member, bin_position};

pub fn proc_lvar2(
    i: usize,
    x: &String,
    pass: usize,
    u: usize,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    gex_info: &GexInfo,
    row: &mut Vec<String>,
    out_data: &mut Vec<HashMap<String, String>>,
    d_all: &mut Vec<Vec<u32>>,
    ind_all: &mut Vec<Vec<u32>>,
    stats: &mut Vec<(String, Vec<String>)>,
    lvars: &Vec<String>,
    _alt_bcs: &Vec<String>,
    gex_mean: f64,
    gex_sum: f64,
    gex_fcounts_unsorted: &Vec<f64>,
    extra_args: &Vec<String>,
) -> bool {
    let clonotype_id = exacts[u];
    let ex = &exact_clonotypes[clonotype_id];
    let verbose = ctl.gen_opt.row_fill_verbose;

    // Set up speak macro.

    macro_rules! speak {
        ($u:expr, $var:expr, $val:expr) => {
            if pass == 2 && (ctl.parseable_opt.pout.len() > 0 || extra_args.len() > 0) {
                let mut v = $var.to_string();
                if ctl.parseable_opt.pcols.is_empty()
                    || bin_member(&ctl.parseable_opt.pcols_sortx, &v)
                    || bin_member(&extra_args, &v)
                {
                    v = v.replace("_Σ", "_sum");
                    v = v.replace("_μ", "_mean");
                    out_data[$u].insert(v, $val);
                }
            }
        };
    }

    // Set up lead variable macros.  This is the mechanism for generating
    // both human-readable and parseable output for lead variables.

    macro_rules! lvar {
        ($i: expr, $var:expr, $val:expr) => {
            if verbose {
                eprint!("lvar {} ==> {}; ", $var, $val);
                eprintln!("$i = {}, lvars.len() = {}", $i, lvars.len());
            }
            if $i < lvars.len() {
                row.push($val)
            }
            if pass == 2 {
                speak!(u, $var.to_string(), $val);
            }
        };
    }
    macro_rules! lvar_stats1 {
        ($i: expr, $var:expr, $val:expr) => {
            if verbose {
                eprint!("lvar {} ==> {}; ", $var, $val);
                eprintln!("$i = {}, lvars.len() = {}", $i, lvars.len());
            }
            if $i < lvars.len() {
                row.push($val)
            }
            if pass == 2 {
                speak!(u, $var.to_string(), $val);
            }
            stats.push(($var.to_string(), vec![$val; ex.ncells()]));
        };
    }

    // Proceed.

    let (mut counts_sub, mut fcounts_sub) = (Vec::<usize>::new(), Vec::<f64>::new());
    let xorig = x.clone();
    let (mut x, mut y) = (x.to_string(), x.to_string());
    if x.contains(':') {
        x = x.before(":").to_string();
    }
    if y.contains(':') {
        y = y.after(":").to_string();
    }
    let y0 = y.clone();
    for _ in 1..=2 {
        let suffixes = ["_min", "_max", "_μ", "_Σ", "_cell", "_%"];
        for s in suffixes.iter() {
            if y.ends_with(s) {
                y = y.rev_before(s).to_string();
                break;
            }
        }
    }
    let mut computed = false;
    for l in 0..ex.clones.len() {
        let li = ex.clones[l][0].dataset_index;
        let bc = ex.clones[l][0].barcode.clone();
        let mut ux = Vec::<usize>::new();
        if ctl.clono_print_opt.regex_match[li].contains_key(&y) {
            ux = ctl.clono_print_opt.regex_match[li][&y].clone();
        }
        if !ux.is_empty() {
            let p = bin_position(&gex_info.gex_barcodes[li], &bc);
            if p >= 0 {
                computed = true;
                let mut raw_count = 0.0;
                for fid in ux.iter() {
                    let raw_counti = get_gex_matrix_entry(
                        ctl, gex_info, *fid, d_all, ind_all, li, l, p as usize, &y,
                    );
                    raw_count += raw_counti;
                }
                counts_sub.push(raw_count.round() as usize);
                fcounts_sub.push(raw_count);
            }
        } else if gex_info.feature_id[li].contains_key(&y) {
            computed = true;
            let p = bin_position(&gex_info.gex_barcodes[li], &bc);
            if p >= 0 {
                let fid = gex_info.feature_id[li][&y];
                let raw_count =
                    get_gex_matrix_entry(ctl, gex_info, fid, d_all, ind_all, li, l, p as usize, &y);
                counts_sub.push(raw_count.round() as usize);
                fcounts_sub.push(raw_count);
            }
        } else {
            counts_sub.push(0);
            fcounts_sub.push(0.0);
        }
    }
    if computed {
        let mut f = Vec::<String>::new();
        for x in fcounts_sub.iter() {
            f.push(format!("{}", x));
        }
        if !y0.ends_with("_%") {
            stats.push((x.clone(), f));
        } else {
            let mut f = Vec::<String>::new();
            for i in 0..fcounts_sub.len() {
                let mut x = 0.0;
                if gex_mean > 0.0 {
                    x = 100.0 * fcounts_sub[i] / gex_mean;
                }
                f.push(format!("{}", x));
            }
            stats.push((x.clone(), f));
        }
        let mut counts_sub_sorted = counts_sub.clone();
        counts_sub_sorted.sort_unstable();
        let sum = fcounts_sub.iter().sum::<f64>();
        let mean = sum / counts_sub.len() as f64;

        if xorig.ends_with("_%_cell") {
            if pass == 2 {
                let mut c = Vec::<String>::new();
                for j in 0..counts_sub.len() {
                    c.push(format!(
                        "{:.2}",
                        100.0 * counts_sub[j] as f64 / gex_fcounts_unsorted[j]
                    ));
                }
                let val = format!("{}", c.iter().format(POUT_SEP));
                speak!(u, x, val);
            }
        } else if xorig.ends_with("_cell") {
            if pass == 2 {
                let val = format!("{}", counts_sub.iter().format(POUT_SEP));
                speak!(u, x, val);
            }
        } else if y0.ends_with("_min") {
            lvar![i, x, format!("{}", counts_sub_sorted[0])];
        } else if y0.ends_with("_max") {
            lvar![i, x, format!("{}", counts_sub_sorted[counts_sub.len() - 1])];
        } else if y0.ends_with("_μ") {
            lvar![i, x, format!("{}", mean.round())];
        } else if y0.ends_with("_Σ") {
            lvar![i, x, format!("{}", sum.round())];
        } else if y0.ends_with("_%") {
            lvar![i, x, format!("{:.2}", (100.0 * sum) / gex_sum)];
        } else {
            let mut median = 0;
            if !counts_sub_sorted.is_empty() {
                median = rounded_median(&counts_sub_sorted);
            }
            lvar![i, x, format!("{}", median)];
        }
    } else if i < lvars.len() {
        lvar_stats1![i, x, "".to_string()];
    }
    true
}
