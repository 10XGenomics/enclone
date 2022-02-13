// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use enclone_core::defs::EncloneControl;
use io_utils::*;
use std::collections::HashMap;
use std::io::Write;
use tables::*;
use vector_utils::*;

pub fn print_fate(ctl: &EncloneControl, fate: &Vec<HashMap<String, String>>, logx: &mut Vec<u8>) {
    // Print barcode fate.

    fwriteln!(logx, "2. barcode fate");
    let mut fates = Vec::<String>::new();
    for i in 0..fate.len() {
        for f in fate[i].iter() {
            if f.1.contains(" GEX ") && ctl.clono_filt_opt_def.ngex {
                continue;
            }
            if f.1.contains(" CROSS ") && ctl.clono_filt_opt_def.ncross {
                continue;
            }
            if f.1.contains(" UMI ") && !ctl.clono_filt_opt_def.umi_filt {
                continue;
            }
            if f.1.contains(" UMI_RATIO ") && !ctl.clono_filt_opt_def.umi_ratio_filt {
                continue;
            }
            if f.1.contains(" GRAPH_FILTER ") && ctl.gen_opt.ngraph_filter {
                continue;
            }
            if f.1.contains(" QUAL") && !ctl.clono_filt_opt.qual_filter {
                continue;
            }
            if f.1.contains(" WEAK_CHAINS ") && !ctl.clono_filt_opt_def.weak_chains {
                continue;
            }
            if f.1.contains(" FOURSIE_KILL ") && !ctl.clono_filt_opt_def.weak_foursies {
                continue;
            }
            if f.1.contains(" WHITEF ") && ctl.gen_opt.nwhitef {
                continue;
            }
            if f.1.contains(" BC_DUP ") && !ctl.clono_filt_opt_def.bc_dup {
                continue;
            }
            if f.1.contains(" IMPROPER ") && ctl.merge_all_impropers {
                continue;
            }
            fates.push(f.1.clone());
        }
    }
    fates.sort();
    let mut freq = Vec::<(u32, String)>::new();
    make_freq(&fates, &mut freq);
    let mut rows = Vec::<Vec<String>>::new();
    rows.push(vec!["barcodes".to_string(), "why deleted".to_string()]);
    rows.push(vec!["\\hline".to_string(); 2]);
    for i in 0..freq.len() {
        rows.push(vec![format!("{}", freq[i].0), freq[i].1.clone()]);
    }
    rows.push(vec![format!("{}", fates.len()), "total".to_string()]);
    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows, 2, &b"r|l".to_vec(), false, false);
    log.truncate(log.len() - 1);
    log = log.replace("\n", "\n   ");
    fwrite!(logx, "   {}\n", log);
}
