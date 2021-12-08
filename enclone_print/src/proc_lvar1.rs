// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_lvar,
// plus a small helper function get_gex_matrix_entry.

use enclone_core::defs::{EncloneControl, ExactClonotype, GexInfo};
use std::collections::HashMap;

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
