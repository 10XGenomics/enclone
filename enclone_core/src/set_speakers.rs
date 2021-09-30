// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::allowed_vars::{CVARS_ALLOWED, CVARS_ALLOWED_PCELL, LVARS_ALLOWED};
use crate::defs::EncloneControl;
use string_utils::TextUtils;
use vector_utils::bin_member;

// Define the set "parseable_fields" of fields that could occur in parseable output.
//
// The overlap with code in proc_args_check.rs is not nice.

pub fn set_speakers(ctl: &EncloneControl, parseable_fields: &mut Vec<String>, max_chains: usize) {
    // Make some abbreviations.

    let lvars = &ctl.clono_print_opt.lvars;

    // Define parseable output columns.  The entire machinery for parseable output is controlled
    // by macros that begin with "speak".

    let pcols_sort = &ctl.parseable_opt.pcols_sort;
    macro_rules! speaker {
        ($var:expr) => {
            if ctl.parseable_opt.pcols.is_empty() || bin_member(&pcols_sort, &$var.to_string()) {
                parseable_fields.push($var.to_string());
            }
        };
    }
    let mut have_gex = false;
    for i in 0..ctl.origin_info.gex_path.len() {
        if !ctl.origin_info.gex_path[i].is_empty() {
            have_gex = true;
        }
    }
    let mut all_lvars = lvars.clone();
    for i in 0..LVARS_ALLOWED.len() {
        let x = &LVARS_ALLOWED[i];
        if !have_gex
            && (*x == "gex".to_string()
                || x.starts_with("gex_")
                || x.ends_with("_g")
                || x.ends_with("_g_μ")
                || *x == "n_gex_cell".to_string()
                || *x == "n_gex".to_string()
                || *x == "n_b".to_string()
                || *x == "clust".to_string()
                || *x == "type".to_string()
                || *x == "entropy".to_string()
                || *x == "cred".to_string()
                || *x == "cred_cell".to_string())
        {
            continue;
        }
        if !lvars.contains(&x.to_string()) {
            all_lvars.push(x.to_string());
        }
    }
    for x in all_lvars.iter() {
        if (*x == "sec" || *x == "mem") && !ctl.gen_opt.using_secmem {
            continue;
        }
        speaker!(x);
    }

    // Define chain variables for parseable output.

    macro_rules! speakerc {
        ($col:expr, $var:expr) => {
            let varc = format!("{}{}", $var, $col + 1);
            if ctl.parseable_opt.pcols.is_empty() || bin_member(&pcols_sort, &varc) {
                parseable_fields.push(format!("{}{}", $var, $col + 1));
            }
        };
    }
    let pchains;
    if ctl.parseable_opt.pchains == "max" {
        pchains = max_chains;
    } else {
        pchains = ctl.parseable_opt.pchains.force_usize();
    }
    for col in 0..pchains {
        for x in CVARS_ALLOWED.iter() {
            speakerc!(col, x);
        }
        if ctl.parseable_opt.pbarcode {
            for x in CVARS_ALLOWED_PCELL.iter() {
                speakerc!(col, x);
            }
        }
        for x in &[
            "var_indices_dna",
            "var_indices_aa",
            "share_indices_dna",
            "share_indices_aa",
            "seq",
            "vj_seq",
            "vj_seq_nl",
            "vj_aa",
            "vj_aa_nl",
            "var_aa",
        ] {
            speakerc!(col, x);
        }
        for i in 0..pcols_sort.len() {
            if pcols_sort[i].starts_with('q') && pcols_sort[i].ends_with(&format!("_{}", col + 1)) {
                let x = pcols_sort[i].after("q").rev_before("_");
                if x.parse::<usize>().is_ok() {
                    parseable_fields.push(pcols_sort[i].clone());
                }
            }
        }
    }

    // Define more lead variables for parseable output.

    speaker!("group_id");
    speaker!("group_ncells");
    speaker!("clonotype_id");
    speaker!("exact_subclonotype_id");
    speaker!("barcodes");
    for x in ctl.origin_info.dataset_list.iter() {
        if !x.is_empty() {
            speaker!(&format!("{}_barcodes", x));
        }
    }
    if ctl.parseable_opt.pbarcode {
        speaker!("barcode");
        for x in ctl.origin_info.dataset_list.iter() {
            speaker!(&format!("{}_barcode", x));
        }
    }
}
