// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Filter B cells based on UMI counts.

use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use equiv::EquivRel;
use stats_utils::binomial_sum;
use std::cmp::max;
use std::collections::HashMap;
use vector_utils::{erase_if, next_diff1_5, reverse_sort, VecUtils};

pub fn filter_umi(
    eq: &EquivRel,
    orbits: &mut Vec<Vec<i32>>,
    ctl: &EncloneControl,
    exact_clonotypes: &mut Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    fate: &mut Vec<HashMap<String, String>>,
) {
    let (mut is_tcr, mut is_bcr) = (true, true);
    if ctl.gen_opt.tcr {
        is_bcr = false;
    }
    if ctl.gen_opt.bcr {
        is_tcr = false;
    }

    // For B cells, filter based on UMI counts.  More details in heuristics.html.
    // Find all clonotypes having one cell which has two chains,
    // one heavy and one light.  Get the sum of the chain UMI counts for this cell.
    //
    // For each cell, let umish be the umi count for its heavy chain having the most umis, and
    // similarly define umisl.  Let umitot = umish + umisl.
    //
    // If every cell in a clonotype would have been deleted, first find the exact subclonotype for
    // which the sum of its umitot values is greatest, and then in it, find the cell having
    // highest umitot value.  Protect this cell, so long as it has at least two chains.

    *orbits = Vec::<Vec<i32>>::new();
    let mut reps = Vec::<i32>::new();
    eq.orbit_reps(&mut reps);
    if is_tcr {
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            eq.orbit(reps[i], &mut o);
            orbits.push(o);
        }
    } else {
        let mut umis = vec![Vec::<usize>::new(); ctl.origin_info.n()];
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            eq.orbit(reps[i], &mut o);
            if o.solo() {
                let x: &CloneInfo = &info[o[0] as usize];
                let ex = &exact_clonotypes[x.clonotype_index];
                if ex.ncells() == 1 && ex.share.duo() && ex.share[0].left != ex.share[1].left {
                    umis[ex.clones[0][0].dataset_index]
                        .push(ex.clones[0][0].umi_count + ex.clones[0][1].umi_count);
                }
            }
        }
        let mut nu = vec![0; ctl.origin_info.n()];
        let mut umin = vec![0.0; ctl.origin_info.n()];
        for l in 0..ctl.origin_info.n() {
            umis[l].sort_unstable();
            nu[l] = umis[l].len();
            if ctl.gen_opt.baseline {
                println!(
                    "\n{} umi counts for dataset {} = {}",
                    nu[l],
                    l + 1,
                    ctl.origin_info.dataset_id[l]
                );
            }
            if nu[l] > 0 {
                let n10 = umis[l][nu[l] / 10] as f64;
                let n50 = umis[l][nu[l] / 2] as f64;
                umin[l] = n10.min(n50 - (4.0 * n50.sqrt()));
            }
            if nu[l] > 0 && ctl.gen_opt.baseline {
                println!("1% ==> {}", umis[l][umis[l].len() / 100]);
                println!("2% ==> {}", umis[l][umis[l].len() / 50]);
                println!("5% ==> {}", umis[l][umis[l].len() / 20]);
                println!("10% ==> {}", umis[l][umis[l].len() / 10]);
                println!("20% ==> {}", umis[l][umis[l].len() / 5]);
                println!("50% ==> {}", umis[l][umis[l].len() / 2]);
                println!("umin = {:.2}", umin[l]);
            }
        }
        // if ctl.clono_filt_opt_def.umi_filt || ctl.clono_filt_opt_def.umi_filt_mark {
        const MIN_BASELINE_CELLS: usize = 20;
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            eq.orbit(reps[i], &mut o);
            let mut ncells = 0;
            for j in 0..o.len() {
                let x: &CloneInfo = &info[o[j] as usize];
                let ex = &exact_clonotypes[x.clonotype_index];
                ncells += ex.ncells();
            }
            let mut nbads = 0;
            if ncells >= 2 {
                let mut to_deletex = vec![false; o.len()];
                let (mut best_ex, mut best_ex_sum) = (0, 0);
                let (mut best_cell, mut best_cell_count) = (0, 0);
                let mut baselined = true;
                let mut protected = false;
                for pass in 1..=3 {
                    if pass == 2 {
                        if nbads == 0 {
                            protected = true;
                        } else {
                            let p = 0.1;
                            let bound = 0.01;

                            // Find probability of observing nbads or more events of probability
                            // p in a sample of size ncells, and if that is at least bound,
                            // don't delete any cells (except onesies).

                            if binomial_sum(ncells, ncells - nbads, 1.0 - p) >= bound {
                                protected = true;
                            }
                        }
                    }
                    for j in 0..o.len() {
                        let x: &CloneInfo = &info[o[j] as usize];
                        let ex = &mut exact_clonotypes[x.clonotype_index];
                        let mut to_delete = vec![false; ex.ncells()];
                        let mut ex_sum = 0;
                        for k in 0..ex.ncells() {
                            let li = ex.clones[k][0].dataset_index;
                            if nu[li] >= MIN_BASELINE_CELLS {
                                let (mut umish, mut umisl) = (0, 0);
                                for l in 0..ex.share.len() {
                                    if ex.share[l].left {
                                        umish = max(umish, ex.clones[k][l].umi_count);
                                    } else {
                                        umisl = max(umish, ex.clones[k][l].umi_count);
                                    }
                                }
                                let umitot = umish + umisl;
                                if pass == 1 {
                                    ex_sum += umitot;
                                }
                                if pass == 2
                                    && j == best_ex
                                    && umitot > best_cell_count
                                    && ex.share.len() > 1
                                {
                                    best_cell = k;
                                    best_cell_count = umitot;
                                }
                                if (umitot as f64) < umin[li] {
                                    if pass == 1 {
                                        nbads += 1;
                                    } else if pass == 3 && protected {
                                        if ex.share.len() == 1 {
                                            to_delete[k] = true;
                                            if ctl.clono_filt_opt_def.umi_filt_mark {
                                                ex.clones[k][0].marked = true;
                                            }
                                        }
                                    } else if pass == 3
                                        && (!baselined
                                            || (best_ex, best_cell) != (j, k)
                                            || ex.share.len() == 1)
                                    {
                                        to_delete[k] = true;
                                        if ctl.clono_filt_opt_def.umi_filt_mark {
                                            ex.clones[k][0].marked = true;
                                        }
                                    }
                                }
                            } else {
                                baselined = false;
                            }
                        }
                        if pass == 1 && ex_sum > best_ex_sum {
                            best_ex = j;
                            best_ex_sum = ex_sum;
                        }
                        if pass == 3 {
                            for i in 0..ex.clones.len() {
                                if to_delete[i] {
                                    fate[ex.clones[i][0].dataset_index].insert(
                                        ex.clones[i][0].barcode.clone(),
                                        "failed UMI filter".to_string(),
                                    );
                                }
                            }
                            if ctl.clono_filt_opt_def.umi_filt {
                                erase_if(&mut ex.clones, &to_delete);
                            }
                        }
                    }
                }
                for j in 0..o.len() {
                    let x: &CloneInfo = &info[o[j] as usize];
                    let ex = &mut exact_clonotypes[x.clonotype_index];
                    if ex.ncells() == 0 {
                        to_deletex[j] = true;
                    }
                }
                erase_if(&mut o, &to_deletex);
            }
            if !o.is_empty() {
                orbits.push(o.clone());
            }
        }
        // }
    }

    // Filter B cells based on UMI count ratios.  This assumes V..J identity to filter.

    if is_bcr {
        const MIN_UMI_RATIO: usize = 500;
        let mut orbits2 = Vec::<Vec<i32>>::new();
        'orbit: for i in 0..orbits.len() {
            let mut ncells = 0;
            let mut o = orbits[i].clone();
            for j in 0..o.len() {
                let x: &CloneInfo = &info[o[j] as usize];
                let ex = &exact_clonotypes[x.clonotype_index];
                ncells += ex.ncells();
            }
            let mut nbads = 0;
            for pass in 1..=2 {
                if pass == 2 {
                    if nbads == 0 {
                        orbits2.push(o.clone());
                        continue 'orbit;
                    } else {
                        let p = 0.1;
                        let bound = 0.01;

                        // Find probability of observing nbads or more events of probability
                        // p in a sample of size ncells, and if that is at least bound,
                        // don't delete any cells.

                        if binomial_sum(ncells, ncells - nbads, 1.0 - p) >= bound {
                            orbits2.push(o.clone());
                            continue 'orbit;
                        }
                    }
                }
                let mut to_deletex = vec![false; o.len()];
                let mut z = Vec::<(Vec<u8>, usize, usize, usize, usize)>::new();
                let mut to_delete = Vec::<Vec<bool>>::new();
                for j in 0..o.len() {
                    let x: &CloneInfo = &info[o[j] as usize];
                    let ex = &mut exact_clonotypes[x.clonotype_index];
                    to_delete.push(vec![false; ex.ncells()]);
                    for k in 0..ex.ncells() {
                        let mut tot = 0;
                        for m in 0..ex.clones[k].len() {
                            tot += ex.clones[k][m].umi_count;
                        }
                        for m in 0..ex.clones[k].len() {
                            z.push((
                                ex.share[m].seq.clone(),
                                ex.clones[k][m].umi_count,
                                j,
                                k,
                                tot,
                            ));
                        }
                    }
                }
                reverse_sort(&mut z);
                let mut j = 0;
                while j < z.len() {
                    let k = next_diff1_5(&z, j as i32) as usize;
                    for l in j..k {
                        if z[j].1 >= MIN_UMI_RATIO * z[l].4 {
                            to_delete[z[l].2][z[l].3] = true;
                        }
                    }
                    j = k;
                }
                for j in 0..o.len() {
                    let x: &CloneInfo = &info[o[j] as usize];
                    let ex = &mut exact_clonotypes[x.clonotype_index];
                    for l in 0..ex.ncells() {
                        if to_delete[j][l] {
                            if ctl.clono_filt_opt_def.umi_ratio_filt_mark {
                                ex.clones[l][0].marked = true;
                            }
                            nbads += 1;
                        }
                    }

                    if pass == 2 {
                        for i in 0..ex.clones.len() {
                            if to_delete[j][i] {
                                fate[ex.clones[i][0].dataset_index].insert(
                                    ex.clones[i][0].barcode.clone(),
                                    "failed UMI_RATIO filter".to_string(),
                                );
                            }
                        }
                        if ctl.clono_filt_opt_def.umi_ratio_filt {
                            erase_if(&mut ex.clones, &to_delete[j]);
                            if ex.ncells() == 0 {
                                to_deletex[j] = true;
                            }
                        }
                    }
                }
                if pass == 2 {
                    if ctl.clono_filt_opt_def.umi_ratio_filt {
                        erase_if(&mut o, &to_deletex);
                    }
                    if !o.is_empty() {
                        orbits2.push(o.clone());
                    }
                }
            }
        }
        *orbits = orbits2;
    }
}
