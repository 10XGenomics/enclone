// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_cvar.

use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use std::collections::HashMap;

pub fn proc_cvar2(
    _var: &String,
    _j: usize,
    _col: usize,
    _mid: usize,
    _pass: usize,
    _u: usize,
    _ex: &ExactClonotype,
    _ctl: &EncloneControl,
    _exacts: &Vec<usize>,
    _exact_clonotypes: &Vec<ExactClonotype>,
    _out_data: &mut Vec<HashMap<String, String>>,
    _rsi: &ColInfo,
    _peer_groups: &Vec<Vec<(usize, u8, u32)>>,
    _show_aa: &Vec<Vec<usize>>,
    _field_types: &Vec<Vec<u8>>,
    _col_var: bool,
    _pcols_sort: &Vec<String>,
    _bads: &mut Vec<bool>,
    _cx: &mut Vec<Vec<String>>,
    _extra_args: &Vec<String>,
    _stats: &mut Vec<(String, Vec<String>)>,
) -> bool {
    false
}
