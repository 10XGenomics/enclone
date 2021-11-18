// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_lvar,
// plus a small helper function get_gex_matrix_entry.

use enclone_core::defs::{EncloneControl, ExactClonotype, GexInfo, POUT_SEP};
use itertools::Itertools;
use std::collections::HashMap;
use string_utils::{abbrev_list, strme, TextUtils};
use vector_utils::bin_member;

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

pub fn proc_lvar1(
    i: usize,
    x: &String,
    pass: usize,
    u: usize,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    mults: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    gex_info: &GexInfo,
    row: &mut Vec<String>,
    out_data: &mut Vec<HashMap<String, String>>,
    _d_all: &mut Vec<Vec<u32>>,
    _ind_all: &mut Vec<Vec<u32>>,
    groups: &HashMap<usize, Vec<usize>>,
    stats: &mut Vec<(String, Vec<String>)>,
    _vdj_cells: &Vec<Vec<String>>,
    _n_vdj_gex: &Vec<usize>,
    nd_fields: &Vec<String>,
    lvars: &Vec<String>,
    alt_bcs: &Vec<String>,
    _n_gex: usize,
    _n_gexs: &Vec<usize>,
    _gex_min: usize,
    _gex_max: usize,
    _gex_mean: f64,
    _gex_sum: f64,
    _gex_median: usize,
    _count_unsorted: &Vec<usize>,
    _entropy: f64,
    _entropies_unsorted: &Vec<f64>,
    _fcounts: &Vec<f64>,
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
                v = v.replace("_Σ", "_sum");
                v = v.replace("_μ", "_mean");
                if ctl.parseable_opt.pcols.is_empty()
                    || bin_member(&ctl.parseable_opt.pcols_sortx, &v)
                    || bin_member(&extra_args, &v)
                {
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

    // Check INFO.

    if ctl.gen_opt.info.is_some() {
        for q in 0..ctl.gen_opt.info_fields.len() {
            if *x == ctl.gen_opt.info_fields[q] {
                let mut found = false;
                let mut lvarred = false;
                if ex.share.len() == 2 && ex.share[0].left != ex.share[1].left {
                    let mut tag = String::new();
                    for j in 0..ex.share.len() {
                        if ex.share[j].left {
                            tag += strme(&ex.share[j].seq);
                        }
                    }
                    tag += "_";
                    for j in 0..ex.share.len() {
                        if !ex.share[j].left {
                            tag += strme(&ex.share[j].seq);
                        }
                    }
                    if ctl.gen_opt.info_data.contains_key(&tag) {
                        let val = &ctl.gen_opt.info_data[&tag][q];
                        lvar_stats1![i, x, val.clone()];
                        lvarred = true;
                        found = true;
                    }
                }
                if !lvarred {
                    lvar_stats1![i, x, String::new()];
                }
                if !found {
                    stats.push((x.to_string(), vec![]));
                }
                return true;
            }
        }
    }

    // Proceed.

    if x.starts_with('g') && x.after("g").parse::<usize>().is_ok() {
        let d = x.after("g").force_usize();
        if groups.contains_key(&d) {
            lvar_stats1![i, x, format!("{}", groups[&d][u] + 1)];
            return true;
        }
    }

    if x == "n" {
        let counts = vec!["1.0".to_string(); mults[u]];
        lvar_stats![i, x, format!("{}", mults[u]), counts];
    } else if x == "clust" {
        let mut clust = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            let mut cid = 0;
            let bc = &ex.clones[j][0].barcode;
            let li = ex.clones[j][0].dataset_index;
            if gex_info.cluster[li].contains_key(&bc.clone()) {
                cid = gex_info.cluster[li][&bc.clone()];
            }
            clust.push(cid);
        }
        let mut clustf = Vec::<String>::new();
        for x in clust.iter() {
            clustf.push(format!("{}", x));
        }
        clust.sort_unstable();
        lvar_stats![i, x, abbrev_list(&clust), clustf];
    } else if x == "n_other" {
        let mut n = 0;
        let mut ns = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            let mut found = false;
            let di = ex.clones[j][0].dataset_index;
            let f = format!("n_{}", ctl.origin_info.dataset_id[di]);
            for i in 0..nd_fields.len() {
                if f == nd_fields[i] {
                    found = true;
                }
            }
            if !found {
                n += 1;
                ns.push("1.0".to_string());
            } else {
                ns.push("0.0".to_string());
            }
        }
        lvar_stats![i, x, format!("{}", n), ns];
    } else if x == "n_b" {
        let mut n_b = 0;
        let mut ns = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            let bc = &ex.clones[j][0].barcode;
            let li = ex.clones[j][0].dataset_index;
            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                    n_b += 1;
                    ns.push("1.0".to_string());
                } else {
                    ns.push("0.0".to_string());
                }
            }
        }
        lvar_stats![i, x, format!("{}", n_b), ns];
    } else if x == "type" {
        let mut cell_types = Vec::<String>::new();
        /*
        for j in 0..ex.clones.len() {
            let mut cell_type = "".to_string();
            let bc = &ex.clones[j][0].barcode;
            let li = ex.clones[j][0].dataset_index;
            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                cell_type = gex_info.cell_type[li][&bc.clone()].clone();
            }
            cell_types.push(cell_type);
        }
        */
        cell_types.sort();
        lvar![i, x, abbrev_list(&cell_types)];
    } else if x == "mark" {
        let mut n = 0;
        for j in 0..ex.clones.len() {
            if ex.clones[j][0].marked {
                n += 1;
            }
        }
        lvar![i, x, format!("{}", n)];
    } else if x.starts_with("pe") {
        lvar_stats1![i, x, format!("")];
    } else if x.starts_with("npe") {
        lvar_stats1![i, x, format!("")];
    } else if x.starts_with("ppe") {
        lvar_stats1![i, x, format!("")];
    } else if bin_member(alt_bcs, x) {
        let mut r = Vec::<String>::new();
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = ex.clones[l][0].barcode.clone();
            let mut val = String::new();
            let alt = &ctl.origin_info.alt_bc_fields[li];
            for j in 0..alt.len() {
                if alt[j].0 == *x && alt[j].1.contains_key(&bc.clone()) {
                    val = alt[j].1[&bc.clone()].clone();
                }
            }
            r.push(val);
        }
        lvar_stats![i, x, format!(""), r];
        if pass == 2 {
            speak!(u, x, format!("{}", r.iter().format(POUT_SEP)));
        }
    } else if x.starts_with("n_") && !x.starts_with("n_gex") {
        let name = x.after("n_");
        let mut count = 0;
        let mut counts = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            let x = &ex.clones[j][0];
            if ctl.origin_info.dataset_id[x.dataset_index] == name {
                count += 1;
                counts.push("1.0".to_string());
            } else if x.origin_index.is_some()
                && ctl.origin_info.origin_list[x.origin_index.unwrap()] == name
            {
                count += 1;
                counts.push("1.0".to_string());
            } else if x.donor_index.is_some()
                && ctl.origin_info.donor_list[x.donor_index.unwrap()] == name
            {
                count += 1;
                counts.push("1.0".to_string());
            } else if x.tag_index.is_some()
                && ctl.origin_info.tag_list[x.tag_index.unwrap()] == name
            {
                count += 1;
                counts.push("1.0".to_string());
            }
        }
        lvar_stats![i, x, format!("{}", count), counts];
    } else {
        return false;
    }
    true
}
