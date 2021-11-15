// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_cvar.

use crate::print_utils1::color_codon;
use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use std::collections::HashMap;
use string_utils::strme;

pub fn proc_cvar1(
    var: &String,
    j: usize,
    col: usize,
    mid: usize,
    _pass: usize,
    u: usize,
    _ex: &ExactClonotype,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    _out_data: &mut Vec<HashMap<String, String>>,
    rsi: &ColInfo,
    peer_groups: &Vec<Vec<(usize, u8, u32)>>,
    show_aa: &Vec<Vec<usize>>,
    ref_diff_pos: &Vec<Vec<Vec<usize>>>,
    field_types: &Vec<Vec<u8>>,
    col_var: bool,
    _pcols_sort: &Vec<String>,
    cx: &mut Vec<Vec<String>>,
    _extra_args: &Vec<String>,
    _stats: &mut Vec<(String, Vec<String>)>,
    cdr3_con: &Vec<Vec<u8>>,
) -> Result<bool, String> {
    let seq_amino = &rsi.seqss_amino[col][u];
    if *var == "amino" && col_var {
        let mut last_color = "black".to_string();
        for k in 0..show_aa[col].len() {
            let p = show_aa[col][k];
            if k > 0 && field_types[col][k] != field_types[col][k - 1] && !ctl.gen_opt.nospaces {
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
    } else {
        return Ok(false);
    }
    Ok(true)
}
