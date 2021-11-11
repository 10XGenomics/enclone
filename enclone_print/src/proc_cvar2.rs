// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_cvar.

use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use stats_utils::percent_ratio;
use std::collections::HashMap;
use vector_utils::next_diff12_4;

pub fn proc_cvar2(
    var: &String,
    j: usize,
    col: usize,
    _mid: usize,
    _pass: usize,
    _u: usize,
    ex: &ExactClonotype,
    _ctl: &EncloneControl,
    _exacts: &Vec<usize>,
    _exact_clonotypes: &Vec<ExactClonotype>,
    _out_data: &mut Vec<HashMap<String, String>>,
    _rsi: &ColInfo,
    _peer_groups: &Vec<Vec<(usize, u8, u32)>>,
    _show_aa: &Vec<Vec<usize>>,
    _field_types: &Vec<Vec<u8>>,
    col_var: bool,
    _pcols_sort: &Vec<String>,
    _bads: &mut Vec<bool>,
    cx: &mut Vec<Vec<String>>,
    _extra_args: &Vec<String>,
    _stats: &mut Vec<(String, Vec<String>)>,
) -> bool {
    // Compute potential whitelist contamination percent and filter.
    // This is an undocumented option.

    if *var == "white" {
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
    } else {
        return false;
    }
    true
}
