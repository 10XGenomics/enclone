// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// If NWEAK_ONESIES is not specified, disintegrate certain onesie clonotypes into single cell
// clonotypes.  This requires editing of exact_clonotypes, info, eq, join_info and raw_joins.

use enclone_core::defs::*;
use equiv::EquivRel;
use std::collections::HashMap;
use std::time::Instant;
use vector_utils::*;

pub fn disintegrate_onesies(
    ctl: &EncloneControl,
    disintegrated: &mut Vec<bool>,
    eq: &mut EquivRel,
    exact_clonotypes: &mut Vec<ExactClonotype>,
    info: &mut Vec<CloneInfo>,
    join_info: &mut Vec<(usize, usize, bool, Vec<u8>)>,
    raw_joins: &mut Vec<(i32, i32)>,
) {
    if ctl.clono_filt_opt.weak_onesies {
        let t = Instant::now();
        let ncells_total = exact_clonotypes.iter().map(|x| x.ncells()).sum();
        let mut to_info = HashMap::<usize, usize>::new();
        let mut exacts2 = Vec::<ExactClonotype>::new();
        for i in 0..info.len() {
            to_info.insert(info[i].clonotype_index, i);
        }
        let mut to_exact_new = Vec::<Vec<usize>>::new();
        for i in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[i];
            let mut enew = Vec::<usize>::new();
            if ex.share.len() == 1
                && ex.ncells() > 1
                && ex.ncells() * 1000 < ncells_total
                && to_info.contains_key(&i)
                && eq.orbit_size(to_info[&i] as i32) == 1
            {
                for j in 0..ex.clones.len() {
                    enew.push(exacts2.len());
                    exacts2.push(ExactClonotype {
                        share: ex.share.clone(),
                        clones: vec![ex.clones[j].clone()],
                    });
                    disintegrated.push(true);
                }
            } else {
                enew.push(exacts2.len());
                exacts2.push(exact_clonotypes[i].clone());
                disintegrated.push(false);
            }
            to_exact_new.push(enew);
        }
        let mut join_info2 = Vec::new();
        for i in 0..join_info.len() {
            let (u1, u2) = (join_info[i].0, join_info[i].1);
            for v1 in to_exact_new[u1].iter() {
                for v2 in to_exact_new[u2].iter() {
                    let mut x = join_info[i].clone();
                    x.0 = *v1;
                    x.1 = *v2;
                    join_info2.push(x);
                }
            }
        }
        ctl.perf_stats(&t, "disintegrating onesies 1");
        let t = Instant::now();
        *join_info = join_info2;
        *exact_clonotypes = exacts2;
        let mut info2 = Vec::<CloneInfo>::new();
        let mut to_info2 = Vec::<Vec<usize>>::new();
        for i in 0..info.len() {
            let j = info[i].clonotype_index;
            let mut x = Vec::<usize>::new();
            for k in 0..to_exact_new[j].len() {
                info[i].clonotype_index = to_exact_new[j][k];
                info[i].clonotype_id = to_exact_new[j][k];
                let mut origins = Vec::<usize>::new();
                let ex = &exact_clonotypes[info[i].clonotype_index];
                for i in 0..ex.clones.len() {
                    origins.push(ex.clones[i][0].dataset_index);
                }
                unique_sort(&mut origins);
                info[i].origin = origins;
                x.push(info2.len());
                info2.push(info[i].clone());
            }
            to_info2.push(x);
        }
        ctl.perf_stats(&t, "disintegrating onesies 2");
        let t = Instant::now();
        *info = info2;
        let mut raw_joins2 = Vec::<(i32, i32)>::new();
        for i in 0..raw_joins.len() {
            let (j1, j2) = (
                &to_info2[raw_joins[i].0 as usize],
                &to_info2[raw_joins[i].1 as usize],
            );
            raw_joins2.push((j1[0] as i32, j2[0] as i32));
        }
        *raw_joins = raw_joins2;
        let mut reps = Vec::<i32>::new();
        eq.orbit_reps(&mut reps);
        let mut eq2 = EquivRel::new(info.len() as i32);
        for i in 0..reps.len() {
            let mut o = Vec::<i32>::new();
            eq.orbit(reps[i], &mut o);
            if o.len() > 1 {
                for j in 0..o.len() - 1 {
                    eq2.join(
                        to_info2[o[j] as usize][0] as i32,
                        to_info2[o[j + 1] as usize][0] as i32,
                    );
                }
            }
        }
        *eq = eq2;
        ctl.perf_stats(&t, "disintegrating onesies 3");
    }
}
