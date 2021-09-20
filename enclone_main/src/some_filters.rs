// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::doublets::*;
use crate::merge_onesies::*;
use crate::split_orbits::*;
use crate::weak_chains::*;
use enclone_core::defs::*;
use enclone_print::define_mat::*;
use equiv::EquivRel;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use vector_utils::*;

pub fn some_filters(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<f64>>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
    eq: &EquivRel,
    disintegrated: &Vec<bool>,
    fate: &mut Vec<HashMap<String, String>>,
) {
    // Delete exact subclonotypes that appear to represent doublets.

    let tdoublet = Instant::now();
    delete_doublets(
        orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
    );
    ctl.perf_stats(&tdoublet, "doublet filtering");

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
            &to_bc,
            &sr,
            &ctl,
            &exact_clonotypes,
            &exacts,
            &od,
            &info,
            &raw_joins,
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
            if dels.contains(&t) {
                if ctl.clono_filt_opt_def.signature {
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
    merge_onesies(orbits, &ctl, &exact_clonotypes, &info, &eq, &disintegrated);
    ctl.perf_stats(&tmerge, "merging onesies");

    // Check for disjoint orbits.

    let tsplit = Instant::now();
    split_orbits(
        orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
    );
    ctl.perf_stats(&tsplit, "splitting orbits 1");

    // Test for weak chains.

    let tweak = Instant::now();
    weak_chains(
        orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
        fate,
    );
    ctl.perf_stats(&tweak, "weak chain filtering");

    // Check for disjoint orbits (again).

    let tsplit = Instant::now();
    split_orbits(
        orbits,
        is_bcr,
        &to_bc,
        &sr,
        &ctl,
        &exact_clonotypes,
        &info,
        &raw_joins,
    );
    ctl.perf_stats(&tsplit, "splitting orbits 2");
}
