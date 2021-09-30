// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Assign the color to a cell in a honeycomb plot.

use crate::TextUtils;
use ansi_escape::print_color13;
use enclone_core::cell_color::CellColor;
use enclone_core::defs::{EncloneControl, ExactClonotype, PlotOpt, POUT_SEP};
use lazy_static::lazy_static;
use std::collections::HashMap;
use std::sync::Mutex;
use vdj_ann::refx::RefData;
use vector_utils::{bin_position, unique_sort, VecUtils};

// The use of these global variables is horrible.  It already lead to one bug that was very
// difficult to track down.

lazy_static! {
    pub static ref VAR_LOW: Mutex<Vec<(String, f64)>> = Mutex::new(Vec::<(String, f64)>::new());
    pub static ref VAR_HIGH: Mutex<Vec<(String, f64)>> = Mutex::new(Vec::<(String, f64)>::new());
}

pub fn assign_cell_color(
    ctl: &EncloneControl,
    plot_opt: &PlotOpt, // overrides ctl
    refdata: &RefData,
    const_names: &Vec<String>,
    dsx: usize,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    out_datas: &Vec<Vec<HashMap<String, String>>>,
    // three variables that specify the cell:
    i: usize, // index into exacts    = clonotype index
    j: usize, // index into exacts[i] = exact subclonotype index
    k: usize, // index into clones    = cell index
) -> String {
    let ex = &exact_clonotypes[exacts[i][j]];
    let mut color = "black".to_string();

    // Determine color for coloring by variable.

    let mut by_var = false;
    match plot_opt.cell_color {
        CellColor::ByVariableValue(_) => {
            by_var = true;
        }
        _ => {}
    };
    if by_var {
        match plot_opt.cell_color {
            CellColor::ByVariableValue(ref x) => {
                color = "undefined".to_string();
                if out_datas[i][j].contains_key(&x.var) {
                    let n = VAR_LOW.lock().unwrap().len();
                    let mut computed = false;
                    for z in 0..n {
                        if VAR_LOW.lock().unwrap()[z].0 == x.var {
                            computed = true;
                        }
                    }
                    if !computed {
                        let (mut low, mut high) = (f64::MAX, f64::MIN);
                        for i in 0..out_datas.len() {
                            for j in 0..out_datas[i].len() {
                                if out_datas[i][j].contains_key(&x.var) {
                                    let val_list = &out_datas[i][j][&x.var];
                                    let vals = val_list.split(POUT_SEP).collect::<Vec<&str>>();
                                    for k in 0..vals.len() {
                                        let val = &vals[k];
                                        if val.parse::<f64>().is_ok() {
                                            let v = val.force_f64();
                                            if v.is_finite() {
                                                low = low.min(v);
                                                high = high.max(v);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        VAR_LOW.lock().unwrap().push((x.var.clone(), low));
                        VAR_HIGH.lock().unwrap().push((x.var.clone(), high));
                    }
                    let (mut low, mut high) = (0.0, 0.0);
                    let n = VAR_LOW.lock().unwrap().len();
                    for z in 0..n {
                        if VAR_LOW.lock().unwrap()[z].0 == x.var {
                            low = VAR_LOW.lock().unwrap()[z].1;
                        }
                    }
                    let n = VAR_HIGH.lock().unwrap().len();
                    for z in 0..n {
                        if VAR_HIGH.lock().unwrap()[z].0 == x.var {
                            high = VAR_HIGH.lock().unwrap()[z].1;
                        }
                    }
                    let val_list = &out_datas[i][j][&x.var];
                    let vals = val_list.split(POUT_SEP).collect::<Vec<&str>>();
                    let val;
                    if vals.solo() {
                        val = &vals[0];
                    } else {
                        val = &vals[k];
                    }
                    if val.parse::<f64>().is_ok() {
                        let mut v = val.force_f64();
                        if v.is_finite() {
                            let xmin;
                            if x.min.is_some() {
                                xmin = x.min.unwrap();
                            } else {
                                xmin = low;
                            }
                            let xmax;
                            if x.max.is_some() {
                                xmax = x.max.unwrap();
                            } else {
                                xmax = high;
                            }
                            v = v.min(xmax);
                            v = v.max(xmin);
                            let vnorm = (v - xmin) / (xmax - xmin);
                            // let c = &TURBO_SRGB_BYTES[(vnorm * 255.0).round() as usize];
                            // color = format!("rgb({},{},{})", c[0], c[1], c[2]);
                            color = format!("turbo-pre-{}", (vnorm * 255.0).round() as usize);
                        }
                    }
                }
            }
            _ => {}
        };

    // Determine color for PLOT_BY_ISOTYPE.
    } else if plot_opt.plot_by_isotype {
        let mut crefs = Vec::<Option<usize>>::new();
        for l in 0..ex.share.len() {
            if ex.share[l].left {
                crefs.push(ex.share[l].c_ref_id);
            }
        }
        unique_sort(&mut crefs);
        let mut color_id = 0;
        if crefs.solo() && crefs[0].is_some() {
            let c = &refdata.name[crefs[0].unwrap()];
            // Note that is possible for p to be -1 in the following.  This is known
            // to happen if a heavy chain V gene is on the same contig as a light
            // chain C gene (which may be an artifact).  There is an example in
            // enclone_main/testx/inputs/flaky.
            let p = bin_position(const_names, c);
            color_id = (1 + p) as usize;
        }
        if plot_opt.plot_by_isotype_color.is_empty() {
            let x = print_color13(color_id);
            color = format!("rgb({},{},{})", x.0, x.1, x.2);
        } else {
            color = plot_opt.plot_by_isotype_color[color_id].clone();
        }

    // Determine color for PLOT_BY_MARK.
    } else if plot_opt.plot_by_mark {
        let dom = ex.clones[k][0].dataset_index == dsx;
        let marked = ex.clones[k][0].marked;
        if dom {
            if !marked {
                color = "red".to_string();
            } else {
                color = "rgb(255,200,200)".to_string();
            }
        } else if !marked {
            color = "blue".to_string();
        } else {
            color = "rgb(200,200,255)".to_string();
        }

    // Determine color in other cases.
    } else {
        if ex.clones[k][0].origin_index.is_some() {
            let s = &ctl.origin_info.origin_list[ex.clones[k][0].origin_index.unwrap()];
            if ctl.gen_opt.origin_color_map.contains_key(&s.clone()) {
                color = ctl.gen_opt.origin_color_map[s].clone();
            }
        }
        if ctl.gen_opt.origin_color_map.is_empty() {
            let mut dataset_colors = false;
            for c in ctl.origin_info.color.iter() {
                if !c.is_empty() {
                    dataset_colors = true;
                }
            }
            let di = ex.clones[k][0].dataset_index;
            if dataset_colors {
                color = ctl.origin_info.color[di].clone();
            } else {
                let bc = &ex.clones[k][0].barcode;
                if ctl.origin_info.barcode_color[di].contains_key(bc) {
                    color = ctl.origin_info.barcode_color[di][bc].clone();
                }
            }
        }
    }
    color
}
