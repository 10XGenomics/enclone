// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use enclone_print::define_mat::define_mat;
use enclone_proto::types::DonorReferenceItem;
use equiv::EquivRel;
use qd::Double;
use std::collections::HashMap;
use vdj_ann::refx::RefData;
use vector_utils::{bin_position, next_diff12_3, unique_sort, VecUtils};

// Check for disjoint orbits.

pub fn split_orbits(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<Double>>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
) {
    let mut orbits2 = Vec::<Vec<i32>>::new();
    for i in 0..orbits.len() {
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

        // Define map of indices into exacts.

        let nexacts = exacts.len();
        let mut to_exacts = HashMap::<usize, usize>::new();
        for u in 0..nexacts {
            to_exacts.insert(exacts[u], u);
        }

        // Get the info indices corresponding to this clonotype.

        let mut infos = Vec::<usize>::new();
        for i in 0..o.len() {
            infos.push(o[i] as usize);
        }

        // Define map of exacts to infos.

        let mut to_infos = vec![Vec::<usize>::new(); nexacts];
        for i in 0..infos.len() {
            let u = to_exacts[&info[infos[i]].clonotype_index];
            to_infos[u].push(i);
        }

        // Determine which columns are "left", meaning IGH or TRB.

        let mut left = vec![false; cols];
        for m in 0..cols {
            for u in 0..mat[0].len() {
                if mat[m][u].is_some() {
                    let c = mat[m][u].unwrap();
                    let ex = &exact_clonotypes[exacts[u]];
                    if ex.share[c].left {
                        left[m] = true;
                    }
                    break;
                }
            }
        }

        // Determine which pairs of configurations share both chain types, and if so, call
        // them joined.

        let mut matu = Vec::<Vec<Option<usize>>>::new();
        for u in 0..nexacts {
            let mut m = Vec::<Option<usize>>::new();
            for c in 0..cols {
                m.push(mat[c][u]);
            }
            matu.push(m);
        }
        unique_sort(&mut matu);
        let mut eqm = vec![vec![false; matu.len()]; matu.len()];
        for j1 in 0..matu.len() {
            for j2 in 0..matu.len() {
                let (mut l, mut r) = (false, false);
                for m in 0..cols {
                    if matu[j1][m].is_some() && matu[j2][m].is_some() {
                        if left[m] {
                            l = true;
                        } else {
                            r = true;
                        }
                    }
                }
                if l && r {
                    eqm[j1][j2] = true;
                }
            }
        }

        // Propagate this to an equivalence relation on the orbit elements.

        let mut eqx = EquivRel::new(o.len() as i32);
        let mut lists = vec![Vec::<usize>::new(); matu.len()];
        for u in 0..nexacts {
            let mut m = Vec::<Option<usize>>::new();
            for c in 0..cols {
                m.push(mat[c][u]);
            }
            lists[bin_position(&matu, &m) as usize].push(u);
        }
        for j1 in 0..matu.len() {
            for j2 in j1..matu.len() {
                if eqm[j1][j2] {
                    let u1 = lists[j1][0];
                    for u2 in lists[j2].iter() {
                        for i1 in to_infos[u1].iter() {
                            for i2 in to_infos[*u2].iter() {
                                eqx.join(*i1 as i32, *i2 as i32);
                            }
                        }
                    }
                }
            }
        }
        let mut reps = Vec::<i32>::new();
        eqx.orbit_reps(&mut reps);

        // Join onesies where possible.  This should probably be more efficient.

        for u1 in 0..nexacts {
            let ex1 = &exact_clonotypes[exacts[u1]];
            if ex1.share.solo() {
                let mut is = Vec::<usize>::new();
                for u2 in 0..nexacts {
                    let ex2 = &exact_clonotypes[exacts[u2]];
                    if ex2.share.solo() {
                        if ex1.share[0].seq == ex2.share[0].seq {
                            eqx.join(to_infos[u1][0] as i32, to_infos[u2][0] as i32);
                        }
                    } else {
                        for j in 0..ex2.share.len() {
                            if ex2.share[j].seq == ex1.share[0].seq {
                                is.push(to_infos[u2][0]);
                            }
                        }
                    }
                }
                let mut rs = Vec::<usize>::new();
                for j in 0..is.len() {
                    rs.push(eqx.class_id(is[j] as i32) as usize);
                }
                unique_sort(&mut rs);
                if rs.solo() {
                    eqx.join(to_infos[u1][0] as i32, is[0] as i32);
                }
            }
        }

        // Divide the orbit if needed.

        if eqx.norbits() == 1 {
            orbits2.push(o.clone());
        } else {
            let mut repsx = Vec::<i32>::new();
            eqx.orbit_reps(&mut repsx);
            for j in 0..repsx.len() {
                let mut ox = Vec::<i32>::new();
                eqx.orbit(repsx[j], &mut ox);
                let mut o2 = Vec::<i32>::new();
                for k in 0..ox.len() {
                    o2.push(o[ox[k] as usize]);
                }
                orbits2.push(o2);
            }
        }
    }
    *orbits = orbits2;
}
