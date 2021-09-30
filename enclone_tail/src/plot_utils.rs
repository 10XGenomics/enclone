// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::assign_cell_color::*;
use crate::colors::*;
use crate::hex::*;
use crate::*;
use enclone_core::defs::*;
use std::collections::HashMap;
use vdj_ann::refx::*;
use vector_utils::*;

#[derive(Clone)]
pub struct PlotCluster {
    pub clonotype_index: usize,         // index of the clonotype
    pub colors: Vec<String>,            // color of a cell
    pub coords: Vec<(f64, f64)>,        // coordinates of a cell
    pub barcodes: Vec<(usize, String)>, // (dataset index, barcode) of a cell
    pub radius: f64,                    // radius of the group of circles defining the clonotype
}

// Traverse the clonotypes, building one cluster for each.

pub fn build_clusters(
    ctl: &EncloneControl,
    plot_opt: &PlotOpt, // overrides ctl
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    out_datas: &Vec<Vec<HashMap<String, String>>>,
    const_names: &Vec<String>,
) -> Vec<PlotCluster> {
    let mut clusters = Vec::<PlotCluster>::new();
    const SEP: f64 = 1.0; // separation between clusters
    let mut passes = 1;
    if plot_opt.split_plot_by_origin {
        passes = ctl.origin_info.origin_list.len();
    }
    let tcn = turbo_color_names();
    for i in 0..exacts.len() {
        for pass in 0..passes {
            let mut colors = Vec::<String>::new();
            let mut coords = Vec::<(f64, f64)>::new();
            let mut barcodes = Vec::<(usize, String)>::new();
            let mut n = 0;

            // For PLOT_BY_MARK, find the dataset having the largest number of cells.
            // Ignoring SPLIT_PLOT_BY_ORIGIN.

            let mut dsx = 0;
            if plot_opt.plot_by_mark {
                let mut ds = Vec::<usize>::new();
                for j in 0..exacts[i].len() {
                    let ex = &exact_clonotypes[exacts[i][j]];
                    for j in 0..ex.clones.len() {
                        ds.push(ex.clones[j][0].dataset_index);
                    }
                }
                ds.sort_unstable();
                let mut freq = Vec::<(u32, usize)>::new();
                make_freq(&ds, &mut freq);
                dsx = freq[0].1;
            }

            // Go through the exact subclonotypes in a clonotype.

            for j in 0..exacts[i].len() {
                let ex = &exact_clonotypes[exacts[i][j]];

                // Traverse the cells in the exact subclonotype.

                for k in 0..ex.clones.len() {
                    if passes > 1 {
                        let li = ex.clones[k][0].dataset_index;
                        let p = bin_position(
                            &ctl.origin_info.origin_list,
                            &ctl.origin_info.origin_id[li],
                        );
                        if pass != p as usize {
                            continue;
                        }
                    }
                    barcodes.push((
                        ex.clones[k][0].dataset_index,
                        ex.clones[k][0].barcode.clone(),
                    ));
                    let mut color = assign_cell_color(
                        ctl,
                        plot_opt,
                        refdata,
                        const_names,
                        dsx,
                        exacts,
                        exact_clonotypes,
                        out_datas,
                        i,
                        j,
                        k,
                    );

                    // Partially translate turbo colors.

                    if color.starts_with("turbo-pre-") {
                        let n = color.after("turbo-pre-").force_usize();
                        color = tcn[n].clone();
                    }

                    // Save.

                    colors.push(color);
                    coords.push(hex_coord(n, 1.0));
                    n += 1;
                }
            }

            // Move colors around to get vertical separation within a given clonotype,
            // e.g. blues on left, reds on right.

            coords.sort_by(|a, b| a.partial_cmp(b).unwrap());
            sort_sync2(&mut colors, &mut barcodes);

            // Substitute enclone colors.

            for j in 0..colors.len() {
                substitute_enclone_color(&mut colors[j]);
            }

            // Save.

            if !barcodes.is_empty() {
                let mut radius = 0.0f64;
                for j in 0..coords.len() {
                    radius = radius
                        .max(1.0 + (coords[j].0 * coords[j].0 + coords[j].1 * coords[j].1).sqrt());
                }
                radius += SEP;
                clusters.push(PlotCluster {
                    colors,
                    coords,
                    clonotype_index: i,
                    barcodes,
                    radius,
                });
            }
        }
    }
    clusters
}
