// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Delete exact subclonotypes that appear to represent doublets.

use enclone_core::defs::*;
use enclone_print::define_mat::*;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::HashMap;
use vector_utils::*;

pub fn delete_doublets(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<f64>>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
) {
    if ctl.clono_filt_opt.doublet {
        // Define pure subclonotypes.  To do this we break each clonotype up by chain signature.
        // Note duplication of code with print_clonotypes.rs.  And this is doing some
        // superfluous compute.

        let mut results = Vec::<(usize, Vec<Vec<usize>>)>::new();
        for i in 0..orbits.len() {
            results.push((i, Vec::new()));
        }
        let mut pures = Vec::<Vec<usize>>::new();
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
            let nexacts = mat[0].len();
            let mut priority = Vec::<Vec<bool>>::new();
            for u in 0..nexacts {
                let mut typex = vec![false; mat.len()];
                for col in 0..mat.len() {
                    if mat[col][u].is_some() {
                        typex[col] = true;
                    }
                }
                priority.push(typex.clone());
            }
            sort_sync2(&mut priority, &mut exacts);
            let mut j = 0;
            while j < priority.len() {
                let k = next_diff(&priority, j);
                let mut p = Vec::<usize>::new();
                for l in j..k {
                    p.push(exacts[l] as usize);
                }
                res.1.push(p);
                j = k;
            }
        });
        for i in 0..results.len() {
            pures.append(&mut results[i].1);
        }

        // Define the number of cells in each pure subclonotype.

        let mut npure = vec![0; pures.len()];
        for j in 0..pures.len() {
            for id in pures[j].iter() {
                npure[j] += exact_clonotypes[*id].ncells();
            }
        }

        // Find the pairs of pure subclonotypes that share identical CDR3 sequences.

        let mut shares = Vec::<(usize, usize)>::new();
        {
            let mut content = Vec::<(String, usize)>::new();
            for j in 0..pures.len() {
                for id in pures[j].iter() {
                    let ex = &exact_clonotypes[*id];
                    for k in 0..ex.share.len() {
                        content.push((ex.share[k].cdr3_dna.clone(), j));
                    }
                }
            }
            content.par_sort();
            content.dedup();
            let mut j = 0;
            while j < content.len() {
                let k = next_diff1_2(&content, j as i32) as usize;
                for l1 in j..k {
                    for l2 in j + 1..k {
                        shares.push((content[l1].1, content[l2].1));
                        shares.push((content[l2].1, content[l1].1));
                    }
                }
                j = k;
            }
            shares.par_sort();
            shares.dedup();
        }

        // Find triples of pure subclonotypes in which the first two have no share, but both
        // of the first two share with the third.

        const MIN_MULT_DOUBLET: usize = 5;
        let mut trips = Vec::<(usize, usize, usize)>::new();
        {
            let mut us = Vec::<usize>::new();
            let mut vs = Vec::<Vec<usize>>::new();
            let mut j = 0;
            while j < shares.len() {
                let k = next_diff1_2(&shares, j as i32) as usize;
                let u = shares[j].0;
                us.push(u);
                let mut x = Vec::<usize>::new();
                for l in j..k {
                    let v = shares[l].1;
                    if MIN_MULT_DOUBLET * npure[u] <= npure[v] {
                        x.push(v);
                    }
                }
                vs.push(x);
                j = k;
            }
            let mut results = Vec::<(usize, Vec<(usize, usize, usize)>)>::new();
            for i in 0..us.len() {
                results.push((i, Vec::new()));
            }
            results.par_iter_mut().for_each(|res| {
                let i = res.0;
                let u = us[i];
                let vs = &vs[i];
                for l1 in 0..vs.len() {
                    for l2 in l1 + 1..vs.len() {
                        let v1 = vs[l1];
                        let v2 = vs[l2];
                        if !bin_member(&shares, &(v1, v2)) {
                            res.1.push((v1, v2, u));
                        }
                    }
                }
            });
            for j in 0..results.len() {
                trips.append(&mut results[j].1.clone());
            }
        }

        // Delete some of the third members of the triples.

        let mut to_delete = vec![false; exact_clonotypes.len()];
        for j in 0..trips.len() {
            let (v0, v1, v2) = (trips[j].2, trips[j].0, trips[j].1);
            {
                let verbose = false;
                if verbose {
                    println!("\n{}, {}, {}", v0, v1, v2);
                    println!("DELETING");
                    for (u, m) in pures[v0].iter().enumerate() {
                        let ex = &exact_clonotypes[*m];
                        let mut cdrs = Vec::<String>::new();
                        for k in 0..ex.share.len() {
                            cdrs.push(ex.share[k].cdr3_aa.clone());
                        }
                        println!("[{}] {}", u + 1, cdrs.iter().format(","));
                    }
                    println!("USING");
                    for (u, m) in pures[v1].iter().enumerate() {
                        let ex = &exact_clonotypes[*m];
                        let mut cdrs = Vec::<String>::new();
                        for k in 0..ex.share.len() {
                            cdrs.push(ex.share[k].cdr3_aa.clone());
                        }
                        println!("[{}] {}", u + 1, cdrs.iter().format(","));
                    }
                    println!("AND");
                    for (u, m) in pures[v2].iter().enumerate() {
                        let ex = &exact_clonotypes[*m];
                        let mut cdrs = Vec::<String>::new();
                        for k in 0..ex.share.len() {
                            cdrs.push(ex.share[k].cdr3_aa.clone());
                        }
                        println!("[{}] {}", u + 1, cdrs.iter().format(","));
                    }
                }
                for m in pures[v0].iter() {
                    to_delete[*m] = true;
                }
            }
        }
        let mut orbits2 = Vec::<Vec<i32>>::new();
        for i in 0..orbits.len() {
            let mut o = orbits[i].clone();

            let mut del2 = vec![false; o.len()];
            for j in 0..o.len() {
                let id = info[o[j] as usize].clonotype_index;
                if to_delete[id] {
                    del2[j] = true;
                }
            }
            erase_if(&mut o, &del2);
            orbits2.push(o);
        }
        *orbits = orbits2;
    }
}
