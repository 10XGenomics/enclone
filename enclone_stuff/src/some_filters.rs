// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::doublets::delete_doublets;
use crate::merge_onesies::merge_onesies;
use crate::split_orbits::split_orbits;
use crate::weak_chains::weak_chains;
use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use enclone_print::define_mat::define_mat;
use enclone_print::print_utils3::define_column_info;
use enclone_proto::types::DonorReferenceItem;
use equiv::EquivRel;
use qd::Double;
use rayon::prelude::*;
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use vdj_ann::refx::RefData;
use vector_utils::{erase_if, next_diff12_3, next_diff1_2, unique_sort};

pub fn some_filters(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<Double>>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
    eq: &EquivRel,
    disintegrated: &Vec<bool>,
    fate: &mut Vec<HashMap<String, String>>,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
) {
    // Delete exact subclonotypes that appear to represent doublets.

    delete_doublets(
        orbits,
        is_bcr,
        to_bc,
        sr,
        ctl,
        exact_clonotypes,
        info,
        raw_joins,
        &refdata,
        dref,
        fate,
    );

    // Given a signature s having at least two chains, if the total cells in the two-chain
    // signatures that are different from it but share a chain with it is at least 20 times
    // greater, delete s.
    //
    // Note duplication of calls to define_mat with other code.  This is expensive.

    let tsig = Instant::now();
    const SIG_MULT: usize = 20;
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
            to_bc,
            sr,
            ctl,
            exact_clonotypes,
            &exacts,
            &od,
            info,
            raw_joins,
            &refdata,
            dref,
        );

        // Find all the signatures and cell counts associated to each.

        let mut freq = Vec::<(usize, Vec<usize>)>::new();
        {
            let mut types = Vec::<(Vec<usize>, usize)>::new();
            for u in 0..exacts.len() {
                let mut t = Vec::<usize>::new();
                for col in 0..mat.len() {
                    if mat[col][u].is_some() {
                        t.push(col);
                    }
                }
                if t.len() >= 2 {
                    types.push((t, exact_clonotypes[exacts[u]].ncells()));
                }
            }
            types.sort();
            let mut i = 0;
            while i < types.len() {
                let j = next_diff1_2(&types, i as i32) as usize;
                let mut mult = 0;
                for k in i..j {
                    mult += types[k].1;
                }
                freq.push((mult, types[i].0.clone()));
                i = j;
            }
        }
        /*
        let mut msg = "\nfrequencies:\n".to_string(); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        use itertools::Itertools; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        for i in 0..freq.len() { // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            msg += &mut format!("{} ==> {}\n", freq[i].0, freq[i].1.iter().format(",")); // XXX
        } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        */

        // Decide which signatures to delete.

        let mut dels = HashSet::<Vec<usize>>::new();
        for i in 0..freq.len() {
            let mut n2 = 0;
            for j in 0..freq.len() {
                if j != i && freq[j].1.len() == 2 {
                    let mut share = false;
                    for x in freq[j].1.iter() {
                        if freq[i].1.contains(x) {
                            share = true;
                        }
                    }
                    if share {
                        n2 += freq[j].0;
                    }
                }
            }
            if n2 > SIG_MULT * freq[i].0 {
                dels.insert(freq[i].1.clone());
                /*
                msg += &mut format!("delete {}\n", freq[i].1.iter().format(",")); // XXXXXXXXXX
                */
            }
        }
        /*
        if dels.len() > 0 { println!("{}", msg); } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        */
        for u in 0..exacts.len() {
            let mut t = Vec::<usize>::new();
            for col in 0..mat.len() {
                if mat[col][u].is_some() {
                    t.push(col);
                }
            }
            if dels.contains(&t) && ctl.clono_filt_opt_def.signature {
                res.2.push(exacts[u]);
                let ex = &exact_clonotypes[exacts[u]];
                for i in 0..ex.ncells() {
                    res.1.push((
                        ex.clones[i][0].dataset_index,
                        ex.clones[i][0].barcode.clone(),
                        "failed SIGNATURE filter".to_string(),
                    ));
                }
            }
        }
    });
    let mut to_delete = vec![false; exact_clonotypes.len()];
    for i in 0..results.len() {
        for j in 0..results[i].1.len() {
            fate[results[i].1[j].0].insert(results[i].1[j].1.clone(), results[i].1[j].2.clone());
        }
        for j in 0..results[i].2.len() {
            to_delete[results[i].2[j]] = true;
        }
    }
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
    ctl.perf_stats(&tsig, "signature filtering");

    // Merge onesies where totally unambiguous.

    let tmerge = Instant::now();
    merge_onesies(orbits, ctl, exact_clonotypes, info, eq, disintegrated);
    ctl.perf_stats(&tmerge, "merging onesies");

    // Check for disjoint orbits.

    let tsplit = Instant::now();
    split_orbits(
        orbits,
        is_bcr,
        to_bc,
        sr,
        ctl,
        exact_clonotypes,
        info,
        raw_joins,
        &refdata,
        dref,
    );
    ctl.perf_stats(&tsplit, "splitting orbits 1");

    // Test for weak chains.

    let tweak = Instant::now();
    weak_chains(
        orbits,
        is_bcr,
        to_bc,
        sr,
        ctl,
        exact_clonotypes,
        info,
        raw_joins,
        fate,
        &refdata,
        dref,
    );
    ctl.perf_stats(&tweak, "weak chain filtering");

    // Check for disjoint orbits (again).

    let tsplit = Instant::now();
    split_orbits(
        orbits,
        is_bcr,
        to_bc,
        sr,
        ctl,
        exact_clonotypes,
        info,
        raw_joins,
        &refdata,
        dref,
    );
    ctl.perf_stats(&tsplit, "splitting orbits 2");

    // Find and mark for deletion exact subclonotypes having a variant base in V..J that,
    // accounting for all the cells in all the exact subclonotypes, never occurs as Q60
    // doesn't occur as Q40 twice, and disagrees with the reference.

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
            to_bc,
            sr,
            ctl,
            exact_clonotypes,
            &exacts,
            &od,
            info,
            raw_joins,
            &refdata,
            dref,
        );
        let cols = mat.len();
        let rsi = define_column_info(ctl, &exacts, exact_clonotypes, &mat, refdata);

        // Create vars, copied from vars_and_shares.

        let mut vars = Vec::<Vec<usize>>::new();
        for cx in 0..cols {
            let mut n = 0;
            for z in 0..rsi.seqss[cx].len() {
                n = max(n, rsi.seqss[cx][z].len());
            }
            let mut v = Vec::<usize>::new();
            for p in 0..n {
                let mut bases = Vec::<u8>::new();
                for s in 0..rsi.seqss[cx].len() {
                    if p >= rsi.seqss[cx][s].len() {
                        continue;
                    }
                    bases.push(rsi.seqss[cx][s][p]);
                }
                unique_sort(&mut bases);
                if bases.len() > 1 {
                    v.push(p);
                }
            }
            vars.push(v);
        }

        // Pretest if using JOIN_BASIC_H.  The code would crash without this.

        let mut neuter = false;
        if ctl.join_alg_opt.basic_h.is_some() {
            let mut ns = vec![Vec::<usize>::new(); cols];
            for u in 0..exacts.len() {
                let clonotype_id = exacts[u];
                let ex = &exact_clonotypes[clonotype_id];
                for col in 0..cols {
                    let m = mat[col][u];
                    if m.is_some() {
                        let m = m.unwrap();
                        if ex.share[m].annv.len() > 1 {
                            continue;
                        }
                        let n = ex.share[m].seq_del.len();
                        ns[col].push(n);
                    }
                }
            }
            for col in 0..cols {
                unique_sort(&mut ns[col]);
                if ns[col].len() > 1 {
                    neuter = true;
                }
            }
        }

        // Proceed.

        // (column, pos, base, qual, row)
        let mut vquals = Vec::<(usize, usize, u8, u8, usize)>::new();
        for u in 0..exacts.len() {
            if neuter {
                continue;
            }
            let clonotype_id = exacts[u];
            let ex = &exact_clonotypes[clonotype_id];
            for col in 0..cols {
                let m = mat[col][u];
                if m.is_some() {
                    let m = m.unwrap();
                    if ex.share[m].annv.len() > 1 {
                        continue;
                    }
                    let n = ex.share[m].seq_del.len();
                    let vref = &exact_clonotypes[exacts[u]].share[m].vs.to_ascii_vec();
                    let jref = &exact_clonotypes[exacts[u]].share[m].js.to_ascii_vec();
                    for z in 0..vars[col].len() {
                        let p = vars[col][z];
                        // not sure how this can happen
                        if ctl.join_alg_opt.basic_h.is_some() && p >= ex.share[m].seq_del.len() {
                            neuter = true;
                            continue;
                        }
                        let b = ex.share[m].seq_del[p];
                        let mut refdiff = false;
                        if p < vref.len() - ctl.heur.ref_v_trim && b != vref[p] {
                            refdiff = true;
                        }
                        if p >= n - (jref.len() - ctl.heur.ref_j_trim)
                            && b != jref[jref.len() - (n - p)]
                        {
                            refdiff = true;
                        }
                        if refdiff {
                            for j in 0..ex.clones.len() {
                                let qual = ex.clones[j][m].quals[p];
                                vquals.push((col, p, b, qual, u));
                            }
                        }
                    }
                }
            }
        }
        if neuter {
            vquals.clear();
        }
        vquals.sort_unstable();
        let mut j = 0;
        while j < vquals.len() {
            let mut k = j + 1;
            while k < vquals.len() {
                if vquals[k].0 != vquals[j].0
                    || vquals[k].1 != vquals[j].1
                    || vquals[k].2 != vquals[j].2
                {
                    break;
                }
                k += 1;
            }
            let mut q60 = false;
            let mut q40 = 0;
            for m in j..k {
                if vquals[m].3 >= 60 {
                    q60 = true;
                } else if vquals[m].3 >= 40 {
                    q40 += 1;
                }
            }
            if !q60 && q40 < 2 {
                let u = vquals[j].4;
                if ctl.clono_filt_opt.qual_filter {
                    res.2.push(exacts[u]);
                }
                let ex = &exact_clonotypes[exacts[u]];
                for i in 0..ex.ncells() {
                    res.1.push((
                        ex.clones[i][0].dataset_index,
                        ex.clones[i][0].barcode.clone(),
                        "failed QUAL filter".to_string(),
                    ));
                }
            }
            j = k;
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
    dels.sort_unstable();
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

    // Check for disjoint orbits (again again).

    let tsplit = Instant::now();
    split_orbits(
        orbits,
        is_bcr,
        to_bc,
        sr,
        ctl,
        exact_clonotypes,
        info,
        raw_joins,
        &refdata,
        dref,
    );
    ctl.perf_stats(&tsplit, "splitting orbits 3");
}
