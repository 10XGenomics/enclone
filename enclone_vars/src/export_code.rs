// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Read the vars file and export code.  This is a partial implementation.

use crate::var::parse_variables;
use io_utils::*;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn export_code() {

    // Define code start/stop for cvar_vdj.

    let cvar_vdj_start =
        r###"

        // Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
        // This file is auto-generated, please do not edit.

        macro_rules! speakc {
            ($u:expr, $col:expr, $var:expr, $val:expr) => {
                if pass == 2
                    && ((ctl.parseable_opt.pout.len() > 0
                        && (ctl.parseable_opt.pchains == "max"
                            || col < ctl.parseable_opt.pchains.force_usize()))
                        || extra_args.len() > 0)
                {
                    let mut v = $var.clone();
                    v = v.replace("_Σ", "_sum");
                    v = v.replace("_μ", "_mean");
    
                    // Strip escape character sequences from val.  Can happen in notes, maybe
                    // other places.
    
                    let mut val_clean = String::new();
                    let mut chars = Vec::<char>::new();
                    let valx = format!("{}", $val);
                    for c in valx.chars() {
                        chars.push(c);
                    }
                    let mut escaped = false;
                    for l in 0..chars.len() {
                        if chars[l] == '^[' {
                            escaped = true;
                        }
                        if escaped {
                            if chars[l] == 'm' {
                                escaped = false;
                            }
                            continue;
                        }
                        val_clean.push(chars[l]);
                    }
    
                    // Proceed.
    
                    let varc = format!("{}{}", v, $col + 1);
                    if pcols_sort.is_empty()
                        || bin_member(&pcols_sort, &varc)
                        || bin_member(&extra_args, &varc)
                    {
                        out_data[$u].insert(varc, val_clean);
                    }
                }
            };
        }

        macro_rules! cvar_stats1 {
            ($i: expr, $var:expr, $val:expr) => {
                if $i < rsi.cvars[col].len() && cvars.contains(&$var) {
                    cx[col][$i] = $val.clone();
                }
                speakc!(u, col, $var, $val);
                let varc = format!("{}{}", $var, col + 1);
                stats.push((varc, vec![$val.to_string(); ex.ncells()]));
            };
        }

        use enclone_core::*;
        use enclone_proto::types::*;
        use vdj_ann::refx::RefData;
        pub fn proc_cvar_auto(
            var: &str,
            ex: &ExactClonotype,
            mid: usize,
            col: usize,
            u: usize,
            rsi: &ColInfo,
            refdata: &RefData,
            dref: &Vec<DonorReferenceItem>,
            ctl: &EncloneControl,
            out_data: &mut Vec<HashMap<String, String>>,
            stats: &mut Vec<(String, Vec<String>)>,
        ) -> Result<bool, String> {
            let val =
            if false {

        "###;

    let cvar_vdj_stop =

        r###"

            else {
                "$UNDEFINED".to_string()
            };
            if val == "$UNDEFINED" {
                return Ok(false);
            } else {
                cvar_stats1![j, var, val];
                return Ok(true);
            }
        }

        "###;

    // Fetch variables and traverse.

    let mut f = open_for_write_new!["enclone_vars/src/proc_cvar_auto.rs"];
    fwrite!(f, "{}", cvar_vdj_start);
    let vars = std::fs::read_to_string("enclone_vars/src/vars").unwrap();
    let vars = parse_variables(&vars);
    for v in vars.iter() {
        if v.inputs == "cvar_vdj" { // RESTRICTION 1
            let mut upper = false;
            let var = &v.name;
            for c in var.chars() {
                if c.is_ascii_uppercase() {
                    upper = true;
                }
            }
            if !upper { // RESTRICTION 2
                if !var.contains('{') { // RESTRICTION 3
                    fwriteln!(f, "}} else if var == \"{}\"", var);
                    fwriteln!(f, "{}", v.code);
                }
            }
        }
    }

    // Add tail.

    fwrite!(f, "{}", cvar_vdj_stop);
}

// candidates:
// aa%
// cdiff
// d1_name
