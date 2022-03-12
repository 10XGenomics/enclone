// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use enclone_core::defs::{EncloneControl, GexInfo};
use io_utils::{fwrite, fwriteln};
use std::io::Write;
use string_utils::TextUtils;
use tables::print_tabular_vbox;

// Print dataset-level variable values.

pub fn print_dataset_vars(ctl: &EncloneControl, gex_info: &GexInfo, logx: &mut Vec<u8>) {
    if !ctl.gen_opt.dvars.is_empty() {
        fwriteln!(logx, "\nDATASET-LEVEL METRICS");
        let mut row = vec!["dataset".to_string()];
        for j in 0..ctl.gen_opt.dvars.len() {
            let var = ctl.gen_opt.dvars[j].clone();
            let mut display_var = var.clone();
            if var.contains(':') {
                display_var = var.before(":").to_string();
            }
            row.push(display_var);
        }
        let mut rows = vec![row];
        for i in 0..ctl.origin_info.n() {
            let mut row = Vec::<String>::new();
            let dataset_name = &ctl.origin_info.dataset_id[i];
            row.push(dataset_name.clone());
            for j in 0..ctl.gen_opt.dvars.len() {
                let mut var = ctl.gen_opt.dvars[j].clone();
                if var.contains(':') {
                    var = var.after(":").to_string();
                }
                let mut value = String::new();
                if gex_info.json_metrics[i].contains_key(&var.to_string()) {
                    value = format!("{:.2}", gex_info.json_metrics[i][&var.to_string()]);
                }
                if value.is_empty() {
                    let mut feature = String::new();
                    let mut typex = String::new();
                    let mut fail = false;
                    if var.ends_with("_cellular_r") {
                        feature = var.before("_cellular_r").to_string();
                        typex = "r".to_string();
                    } else if var.ends_with("_cellular_u") {
                        feature = var.before("_cellular_u").to_string();
                        typex = "u".to_string();
                    } else {
                        fail = true;
                    }
                    if fail {
                        value = "undefined".to_string();
                    } else if typex == "r" {
                        if !gex_info.feature_metrics[i]
                            .contains_key(&(feature.clone(), "num_reads".to_string()))
                            || !gex_info.feature_metrics[i]
                                .contains_key(&(feature.clone(), "num_reads_cells".to_string()))
                        {
                            value = "undefined".to_string();
                        } else {
                            let num = gex_info.feature_metrics[i]
                                [&(feature.clone(), "num_reads_cells".to_string())]
                                .force_usize();
                            let den = gex_info.feature_metrics[i]
                                [&(feature.clone(), "num_reads".to_string())]
                                .force_usize();
                            if den == 0 {
                                value = "0/0".to_string();
                            } else {
                                value = format!("{:.1}", 100.0 * num as f64 / den as f64);
                            }
                        }
                    } else if !gex_info.feature_metrics[i]
                        .contains_key(&(feature.clone(), "num_umis".to_string()))
                    {
                        value = "undefined".to_string();
                    } else if !gex_info.feature_metrics[i]
                        .contains_key(&(feature.clone(), "num_umis_cells".to_string()))
                    {
                        value = "undefined".to_string();
                    } else {
                        let num = gex_info.feature_metrics[i]
                            [&(feature.clone(), "num_umis_cells".to_string())]
                            .force_usize();
                        let den = gex_info.feature_metrics[i]
                            [&(feature.clone(), "num_umis".to_string())]
                            .force_usize();
                        if den == 0 {
                            value = "0/0".to_string();
                        } else {
                            value = format!("{:.1}", 100.0 * num as f64 / den as f64);
                        }
                    }
                }
                row.push(value);
            }
            rows.push(vec!["\\hline".to_string(); row.len()]);
            rows.push(row);
        }
        let mut just = vec![b'l'];
        for _ in 0..ctl.gen_opt.dvars.len() {
            just.push(b'|');
            just.push(b'r');
        }
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &just, false, false);
        fwrite!(logx, "{}", log);
    }
}
