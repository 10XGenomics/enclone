// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_lvar,
// plus a small helper function get_gex_matrix_entry.

use amino::aa_seq;
use enclone_core::defs::{EncloneControl, ExactClonotype, GexInfo, POUT_SEP};
use enclone_core::median::rounded_median;
use itertools::Itertools;
use regex::Regex;
use std::collections::HashMap;
use string_utils::{strme, TextUtils};
use vector_utils::{bin_member, bin_position};

pub fn get_gex_matrix_entry(
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    fid: usize,
    d_all: &Vec<Vec<u32>>,
    ind_all: &Vec<Vec<u32>>,
    li: usize,
    l: usize,
    p: usize,
    y: &str,
) -> f64 {
    let mut raw_count = 0 as f64;
    if gex_info.gex_matrices[li].initialized() {
        raw_count = gex_info.gex_matrices[li].value(p as usize, fid) as f64;
    } else {
        for j in 0..d_all[l].len() {
            if ind_all[l][j] == fid as u32 {
                raw_count = d_all[l][j] as f64;
                break;
            }
        }
    }
    let mult: f64;
    if y.ends_with("_g") {
        mult = gex_info.gex_mults[li];
    } else {
        mult = gex_info.fb_mults[li];
    }
    if !ctl.gen_opt.full_counts {
        raw_count *= mult;
    }
    raw_count
}

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
    macro_rules! lvar_stats {
        ($i: expr, $var:expr, $val:expr, $stats: expr) => {
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
            stats.push(($var.to_string(), $stats.clone()));
        };
    }

    // Proceed.

    if x.starts_with("count_cdr1_") || x.contains(":count_cdr1_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_cdr1_") {
            x = x.before(":count_cdr1_").to_string();
        }
        y = y.after("count_cdr1_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr1_start.is_some() && ex.share[j].fr2_start.is_some() {
                let cdr1 = ex.share[j].cdr1_start.unwrap();
                let fwr2 = ex.share[j].fr2_start.unwrap();
                if cdr1 < fwr2 {
                    let aa = aa_seq(&ex.share[j].seq[cdr1..fwr2], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_cdr2_") || x.contains(":count_cdr2_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_cdr2_") {
            x = x.before(":count_cdr2_").to_string();
        }
        y = y.after("count_cdr2_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr2_start.is_some() && ex.share[j].fr3_start.is_some() {
                let cdr2 = ex.share[j].cdr2_start.unwrap();
                let fwr3 = ex.share[j].fr3_start.unwrap();
                if cdr2 < fwr3 {
                    let aa = aa_seq(&ex.share[j].seq[cdr2..fwr3], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_cdr3_") || x.contains(":count_cdr3_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_cdr3_") {
            x = x.before(":count_cdr3_").to_string();
        }
        y = y.after("count_cdr3_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            let cdr3 = ex.share[j].cdr3_start;
            let fwr4 = cdr3 + 3 * ex.share[j].cdr3_aa.len();
            let aa = aa_seq(&ex.share[j].seq[cdr3..fwr4], 0);
            n += reg.find_iter(strme(&aa)).count();
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_fwr1_") || x.contains(":count_fwr1_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr1_") {
            x = x.before(":count_fwr1_").to_string();
        }
        y = y.after("count_fwr1_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr1_start.is_some() {
                let fwr1 = ex.share[j].fr1_start;
                let cdr1 = ex.share[j].cdr1_start.unwrap();
                if fwr1 < cdr1 {
                    let aa = aa_seq(&ex.share[j].seq[fwr1..cdr1], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_fwr2_") || x.contains(":count_fwr2_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr2_") {
            x = x.before(":count_fwr2_").to_string();
        }
        y = y.after("count_fwr2_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].fr2_start.is_some() && ex.share[j].cdr2_start.is_some() {
                let fwr2 = ex.share[j].fr2_start.unwrap();
                let cdr2 = ex.share[j].cdr2_start.unwrap();
                if fwr2 < cdr2 {
                    let aa = aa_seq(&ex.share[j].seq[fwr2..cdr2], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_fwr3_") || x.contains(":count_fwr3_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr3_") {
            x = x.before(":count_fwr3_").to_string();
        }
        y = y.after("count_fwr3_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].fr3_start.is_some() {
                let fwr3 = ex.share[j].fr3_start.unwrap();
                let cdr3 = ex.share[j].cdr3_start;
                if fwr3 < cdr3 {
                    let aa = aa_seq(&ex.share[j].seq[fwr3..cdr3], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_fwr4_") || x.contains(":count_fwr4_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr4_") {
            x = x.before(":count_fwr4_").to_string();
        }
        y = y.after("count_fwr4_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            let fwr4 = ex.share[j].cdr3_start + 3 * ex.share[j].cdr3_aa.len();
            let aa = aa_seq(&ex.share[j].seq[fwr4..], 0);
            n += reg.find_iter(strme(&aa)).count();
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_cdr_") || x.contains(":count_cdr_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_cdr_") {
            x = x.before(":count_cdr_").to_string();
        }
        y = y.after("count_cdr_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr1_start.is_some() && ex.share[j].fr2_start.is_some() {
                let cdr1 = ex.share[j].cdr1_start.unwrap();
                let fwr2 = ex.share[j].fr2_start.unwrap();
                if cdr1 < fwr2 {
                    let aa = aa_seq(&ex.share[j].seq[cdr1..fwr2], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
            if ex.share[j].cdr2_start.is_some() && ex.share[j].fr3_start.is_some() {
                let cdr2 = ex.share[j].cdr2_start.unwrap();
                let fwr3 = ex.share[j].fr3_start.unwrap();
                if cdr2 < fwr3 {
                    let aa = aa_seq(&ex.share[j].seq[cdr2..fwr3], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
            let cdr3 = ex.share[j].cdr3_start;
            let fwr4 = cdr3 + 3 * ex.share[j].cdr3_aa.len();
            let aa = aa_seq(&ex.share[j].seq[cdr3..fwr4], 0);
            n += reg.find_iter(strme(&aa)).count();
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_fwr_") || x.contains(":count_fwr_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr_") {
            x = x.before(":count_fwr_").to_string();
        }
        y = y.after("count_fwr_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr1_start.is_some() {
                let fwr1 = ex.share[j].fr1_start;
                let cdr1 = ex.share[j].cdr1_start.unwrap();
                if fwr1 < cdr1 {
                    let aa = aa_seq(&ex.share[j].seq[fwr1..cdr1], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
            if ex.share[j].fr2_start.is_some() && ex.share[j].cdr2_start.is_some() {
                let fwr2 = ex.share[j].fr2_start.unwrap();
                let cdr2 = ex.share[j].cdr2_start.unwrap();
                if fwr2 < cdr2 {
                    let aa = aa_seq(&ex.share[j].seq[fwr2..cdr2], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
            if ex.share[j].fr3_start.is_some() {
                let fwr3 = ex.share[j].fr3_start.unwrap();
                let cdr3 = ex.share[j].cdr3_start;
                if fwr3 < cdr3 {
                    let aa = aa_seq(&ex.share[j].seq[fwr3..cdr3], 0);
                    n += reg.find_iter(strme(&aa)).count();
                }
            }
            let fwr4 = ex.share[j].cdr3_start + 3 * ex.share[j].cdr3_aa.len();
            let aa = aa_seq(&ex.share[j].seq[fwr4..], 0);
            n += reg.find_iter(strme(&aa)).count();
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else if x.starts_with("count_") || x.contains(":count_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_") {
            x = x.before(":count_").to_string();
        }
        y = y.after("count_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            let aa = aa_seq(&ex.share[j].seq, 0); // seems inefficient
            n += reg.find_iter(strme(&aa)).count();
        }
        lvar_stats![i, x, format!("{}", n), vec![format!("{}", n); ex.ncells()]];
    } else {
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
                    let raw_count = get_gex_matrix_entry(
                        ctl, gex_info, fid, d_all, ind_all, li, l, p as usize, &y,
                    );
                    counts_sub.push(raw_count.round() as usize);
                    fcounts_sub.push(raw_count);
                }
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
    }
    true
}
