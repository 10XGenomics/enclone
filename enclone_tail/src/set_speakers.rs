// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::allowed_vars::{CVARS_ALLOWED, CVARS_ALLOWED_PCELL, LVARS_ALLOWED};
use enclone_core::defs::EncloneControl;
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
    for x in LVARS_ALLOWED {
        if !have_gex
            && (x == "gex"
                || x.starts_with("gex_")
                || x.ends_with("_g")
                || x.ends_with("_g_μ")
                || x == "n_gex_cell"
                || x == "n_gex"
                || x == "n_b"
                || x == "clust"
                || x == "type"
                || x == "entropy"
                || x == "cred"
                || x == "cred_cell")
        {
            continue;
        }
        if !lvars.contains(&x.to_string()) {
            all_lvars.push(x.to_string());
        }
    }
    for x in &all_lvars {
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
    let pchains = if ctl.parseable_opt.pchains == "max" {
        max_chains
    } else {
        ctl.parseable_opt.pchains.force_usize()
    };
    for col in 0..pchains {
        for x in CVARS_ALLOWED {
            speakerc!(col, x);
        }
        if ctl.parseable_opt.pbarcode {
            for x in CVARS_ALLOWED_PCELL {
                speakerc!(col, x);
            }
        }
        for x in [
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
        for pcol in pcols_sort {
            if pcol.starts_with('q') && pcol.ends_with(&format!("_{}", col + 1)) {
                let x = pcol.after("q").rev_before("_");
                if x.parse::<usize>().is_ok() {
                    parseable_fields.push(pcol.clone());
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
    for x in &ctl.origin_info.dataset_list {
        if !x.is_empty() {
            speaker!(&format!("{x}_barcodes"));
        }
    }
    if ctl.parseable_opt.pbarcode {
        speaker!("barcode");
        for x in &ctl.origin_info.dataset_list {
            speaker!(&format!("{x}_barcode"));
        }
    }
}
