// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Assign the color to a cell in a honeycomb plot.

use crate::*;
use ansi_escape::*;
use enclone_core::defs::*;
use vdj_ann::refx::*;
use vector_utils::*;

pub fn assign_cell_color(
    ctl: &EncloneControl,
    plot_opt: &PlotOpt, // overrides ctl
    refdata: &RefData,
    const_names: &Vec<String>,
    dsx: usize,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    // three variables that specify the cell:
    i: usize, // index into exacts
    j: usize, // index into exacts[i]
    k: usize, // index into clones
) -> String {
    let ex = &exact_clonotypes[exacts[i][j]];
    let mut color = "black".to_string();

    // Determine color for PLOT_BY_ISOTYPE.

    if plot_opt.plot_by_isotype {
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
            let p = bin_position(&const_names, &c);
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
        } else {
            if !marked {
                color = "blue".to_string();
            } else {
                color = "rgb(200,200,255)".to_string();
            }
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
