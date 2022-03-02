// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.

use enclone_core::cell_color::CellColor;
use enclone_core::defs::{ExactClonotype, PlotOpt, POUT_SEP};
use itertools::Itertools;
use std::collections::HashMap;
use vector_utils::make_freq;

pub fn setup_cat_var(
    plot_opt: &PlotOpt, // overrides ctl
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    out_datas: &Vec<Vec<HashMap<String, String>>>,
    by_cat_var: &mut bool,
    barcode_to_cat_var_color: &mut HashMap<(usize, String), String>,
    cat_var_labels: &mut Vec<String>,
) {
    let mut vars = Vec::<String>::new();
    let mut maxcat = 0;
    match plot_opt.cell_color {
        CellColor::ByCategoricalVariableValue(ref x) => {
            *by_cat_var = true;
            vars = x.vars.clone();
            maxcat = x.maxcat;
        }
        _ => {}
    };
    let mut barcode_to_cat_var_value = HashMap::<(usize, String), String>::new();
    if *by_cat_var {
        let mut bvv = Vec::<(usize, String, usize, String)>::new();
        for i in 0..out_datas.len() {
            for j in 0..out_datas[i].len() {
                for (z, var) in vars.iter().enumerate() {
                    if out_datas[i][j].contains_key(&*var) {
                        let val_list = &out_datas[i][j][&*var];
                        let vals = val_list.split(POUT_SEP).collect::<Vec<&str>>();
                        let ex = &exact_clonotypes[exacts[i][j]];
                        for k in 0..ex.ncells() {
                            let val = if vals.len() > 1 { &vals[k] } else { &vals[0] };
                            let li = ex.clones[k][0].dataset_index;
                            let bc = &ex.clones[k][0].barcode;
                            bvv.push((li, bc.clone(), z, val.to_string()));
                        }
                    } else {
                        let ex = &exact_clonotypes[exacts[i][j]];
                        for k in 0..ex.ncells() {
                            let li = ex.clones[k][0].dataset_index;
                            let bc = &ex.clones[k][0].barcode;
                            bvv.push((li, bc.clone(), z, String::new()));
                        }
                    }
                }
            }
        }
        bvv.sort();
        let mut values = Vec::<String>::new();
        for i in (0..bvv.len()).step_by(vars.len()) {
            let mut vals = Vec::<String>::new();
            for j in i..i + vars.len() {
                vals.push(bvv[j].3.clone());
            }
            let v = format!("{}", vals.iter().format(", "));
            barcode_to_cat_var_value.insert((bvv[i].0, bvv[i].1.clone()), v.clone());
            values.push(v);
        }
        values.sort();
        let mut freq = Vec::<(u32, String)>::new();
        make_freq(&values, &mut freq);
        let mut trunc = false;
        if freq.len() > maxcat {
            trunc = true;
            freq.truncate(maxcat - 1);
        }
        for i in 0..freq.len() {
            cat_var_labels.push(freq[i].1.clone());
        }
        if trunc {
            cat_var_labels.push("other".to_string());
        }
        for i in (0..bvv.len()).step_by(vars.len()) {
            let mut vals = Vec::<String>::new();
            for j in i..i + vars.len() {
                vals.push(bvv[j].3.clone());
            }
            let v = format!("{}", vals.iter().format(", "));
            let mut color = String::new();
            let mut found = false;
            for j in 0..freq.len() {
                if v == freq[j].1 {
                    let c = j % 256;
                    color = format!("default-pre-{}", c);
                    found = true;
                    break;
                }
            }
            if !found {
                let c = freq.len();
                color = format!("default-pre-{}", c);
            }
            barcode_to_cat_var_color.insert((bvv[i].0, bvv[i].1.clone()), color);
        }
    }
}
