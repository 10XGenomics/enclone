// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_lvar,
// plus a small helper function get_gex_matrix_entry.

use enclone_core::defs::{EncloneControl, ExactClonotype, GexInfo};
use std::collections::HashMap;

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
    _i: usize,
    _x: &String,
    _pass: usize,
    _u: usize,
    _ctl: &EncloneControl,
    _exacts: &Vec<usize>,
    _exact_clonotypes: &Vec<ExactClonotype>,
    _gex_info: &GexInfo,
    _row: &mut Vec<String>,
    _out_data: &mut Vec<HashMap<String, String>>,
    _d_all: &mut Vec<Vec<u32>>,
    _ind_all: &mut Vec<Vec<u32>>,
    _stats: &mut Vec<(String, Vec<String>)>,
    _lvars: &Vec<String>,
    _alt_bcs: &Vec<String>,
    _gex_mean: f64,
    _gex_sum: f64,
    _gex_fcounts_unsorted: &Vec<f64>,
    _extra_args: &Vec<String>,
) -> bool {
    false
}
