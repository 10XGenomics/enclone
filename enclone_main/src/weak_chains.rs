// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Based on the number of cells in each column, decide which exact subclonotypes
// look like junk.  Preliminary heuristic.

use enclone_core::defs::*;
use enclone_print::define_mat::*;
use rayon::prelude::*;
use std::collections::HashMap;
use vector_utils::*;

pub fn weak_chains(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<f64>>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
    fate: &mut Vec<HashMap<String, String>>,
) {
    // Note mat calculation duplicated with print_clonotypes and also doublet detection.

    let mut results = Vec::<(usize, Vec<(usize, String, String)>, Vec<usize>)>::new();
    for i in 0..orbits.len() {
        results.push((i, Vec::new(), Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let o = orbits[i].clone();
        let mut od = Vec::<(Vec<usize>, usize, i32)>::new();
        for id in o.iter() {
            let x: &CloneInfo = &info[*id as usize];
            od.push((x.origin.clone(), x.clonotype_id, *id));
        }
        od.sort();
        let mut exacts = Vec::<usize>::new();
        let mut j = 0;
        while j < od.len() {
            let k = next_diff12_3(&od, j as i32) as usize;
            exacts.push(od[j].1);
            j = k;
        }
        let mat = define_mat(
            is_bcr,
            &to_bc,
            &sr,
            &ctl,
            &exact_clonotypes,
            &exacts,
            &od,
            &info,
            &raw_joins,
        );
        let cols = mat.len();
        if cols > 2 {
            let nexacts = exacts.len();
            let mut ncells = vec![0; cols];
            let mut col_entries = vec![Vec::<usize>::new(); cols];
            for u in 0..nexacts {
                for col in 0..cols {
                    let mid = mat[col][u];
                    if mid.is_some() {
                        col_entries[col].push(u);
                        let clonotype_id = exacts[u];
                        ncells[col] += exact_clonotypes[clonotype_id].clones.len();
                    }
                }
            }
            let mut total_cells = 0;
            for j in 0..exacts.len() {
                total_cells += exact_clonotypes[exacts[j]].ncells();
            }
            for j in 0..cols {
                if ncells[j] <= 20 && 8 * ncells[j] < total_cells {
                    for d in col_entries[j].iter() {
                        if ctl.clono_filt_opt_def.weak_chains {
                            res.2.push(exacts[*d]);
                        }
                        let ex = &exact_clonotypes[exacts[*d]];
                        for i in 0..ex.ncells() {
                            res.1.push((
                                ex.clones[i][0].dataset_index,
                                ex.clones[i][0].barcode.clone(),
                                "failed WEAK_CHAINS filter".to_string(),
                            ));
                        }
                    }
                }
            }
        }
    });
    let mut to_delete = vec![false; exact_clonotypes.len()];
    let mut dels = Vec::<i32>::new();
    for i in 0..results.len() {
        for j in 0..results[i].1.len() {
            fate[results[i].1[j].0].insert(results[i].1[j].1.clone(), results[i].1[j].2.clone());
        }
        for x in results[i].2.iter() {
            to_delete[*x] = true;
        }
    }
    dels.sort();
    let mut orbits2 = Vec::<Vec<i32>>::new();
    for i in 0..orbits.len() {
        let mut o = orbits[i].clone();
        let mut del = vec![false; o.len()];
        for j in 0..o.len() {
            let id = info[o[j] as usize].clonotype_index;
            if to_delete[id] {
                del[j] = true;
            }
        }
        erase_if(&mut o, &del);
        orbits2.push(o);
    }
    *orbits = orbits2;
}
