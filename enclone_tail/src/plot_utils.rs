// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::assign_cell_color::*;
use crate::hex::*;
use crate::*;
use enclone_core::defs::*;
use vdj_ann::refx::*;
use vector_utils::*;

// Traverse the clonotypes, building one cluster for each.

pub fn build_clusters(
    ctl: &EncloneControl,
    plot_opt: &PlotOpt, // overrides ctl
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    const_names: &Vec<String>,
) -> Vec<(
    Vec<String>,
    Vec<(f64, f64)>,
    usize,
    Vec<(usize, String)>,
    f64,
)> {
    let mut clusters = Vec::<(
        Vec<String>,
        Vec<(f64, f64)>,
        usize,
        Vec<(usize, String)>,
        f64,
    )>::new();
    const SEP: f64 = 1.0; // separation between clusters
    for i in 0..exacts.len() {
        let mut colors = Vec::<String>::new();
        let mut coords = Vec::<(f64, f64)>::new();
        let mut barcodes = Vec::<(usize, String)>::new();
        let mut n = 0;

        // For PLOT_BY_MARK, find the dataset having the largest number of cells.

        let mut dsx = 0;
        if plot_opt.plot_by_mark {
            let mut ds = Vec::<usize>::new();
            for j in 0..exacts[i].len() {
                let ex = &exact_clonotypes[exacts[i][j]];
                for j in 0..ex.clones.len() {
                    ds.push(ex.clones[j][0].dataset_index);
                }
            }
            ds.sort();
            let mut freq = Vec::<(u32, usize)>::new();
            make_freq(&ds, &mut freq);
            dsx = freq[0].1;
        }

        // Go through the exact subclonotypes in a clonotype.

        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];

            // Traverse the cells in the exact subclonotype.

            for k in 0..ex.clones.len() {
                barcodes.push((
                    ex.clones[k][0].dataset_index,
                    ex.clones[k][0].barcode.clone(),
                ));
                let color = assign_cell_color(
                    &ctl,
                    &plot_opt,
                    &refdata,
                    &const_names,
                    dsx,
                    &exacts,
                    &exact_clonotypes,
                    i,
                    j,
                    k,
                );
                colors.push(color);
                coords.push(hex_coord(n, 1.0));
                n += 1;
            }
        }

        // Move colors around to get vertical separation, e.g. blues on left, reds on right.

        coords.sort_by(|a, b| a.partial_cmp(b).unwrap());
        sort_sync2(&mut colors, &mut barcodes);

        // Substitute enclone colors.

        for j in 0..colors.len() {
            substitute_enclone_color(&mut colors[j]);
        }

        // Save.

        let mut radius = 0.0f64;
        for j in 0..coords.len() {
            radius =
                radius.max(1.0 + (coords[j].0 * coords[j].0 + coords[j].1 * coords[j].1).sqrt());
        }
        radius += SEP;
        clusters.push((colors, coords, i, barcodes, radius));
    }
    clusters
}
