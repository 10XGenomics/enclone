// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Merge onesies where totally unambiguous.  Possibly inefficient and should optimize.

use enclone_core::defs::*;
use equiv::EquivRel;
use vector_utils::*;

pub fn merge_onesies(
    orbits: &mut Vec<Vec<i32>>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    eq: &EquivRel,
    disintegrated: &Vec<bool>,
) {
    if ctl.join_alg_opt.merge_onesies {
        // ctl.join_alg_opt.merge_onesies is always true
        let mut eqo = EquivRel::new(orbits.len() as i32);
        let mut to_orbit = vec![None; info.len()];
        for i in 0..orbits.len() {
            for j in 0..orbits[i].len() {
                to_orbit[orbits[i][j] as usize] = Some(i);
            }
        }
        let mut ncells_total = 0;
        for i in 0..exact_clonotypes.len() {
            ncells_total += exact_clonotypes[i].ncells();
        }
        let mut onesies = Vec::<usize>::new();
        for i in 0..info.len() {
            if to_orbit[i].is_some() && info[i].tigs.len() == 1 {
                if !ctl.clono_filt_opt.weak_onesies || !disintegrated[info[i].clonotype_index] {
                    onesies.push(i);
                }
            }
        }
        let mut alltigs2 = Vec::<(Vec<u8>, usize)>::new();
        for i in 0..info.len() {
            if to_orbit[i].is_some() && info[i].tigs.len() >= 2 {
                for j in 0..info[i].tigs.len() {
                    alltigs2.push((info[i].tigs[j].clone(), i));
                }
            }
        }
        alltigs2.sort();
        for x in onesies.iter() {
            let low = lower_bound1_2(&alltigs2, &info[*x].tigs[0]);
            let high = upper_bound1_2(&alltigs2, &info[*x].tigs[0]);
            let mut ms = Vec::<usize>::new();
            for m in low..high {
                if alltigs2[m as usize].0 == info[*x].tigs[0] {
                    ms.push(m as usize);
                }
            }
            let mut ok = ms.len() > 0;
            let mut exacts = Vec::<usize>::new();
            for j in 0..ms.len() {
                if eq.class_id(alltigs2[ms[j]].1 as i32) != eq.class_id(alltigs2[ms[0]].1 as i32) {
                    ok = false;
                }
                let mut o = Vec::<i32>::new();
                eq.orbit(alltigs2[ms[j]].1 as i32, &mut o);
                for z in o.iter() {
                    exacts.push(info[*z as usize].clonotype_index);
                }
            }
            unique_sort(&mut exacts);
            if ctl.join_alg_opt.merge_onesies_ctl {
                let ncells0 = exact_clonotypes[info[*x].clonotype_index].ncells();
                if ncells0 * 10000 < ncells_total {
                    ok = false;
                }
            }
            if ok {
                let orb1 = to_orbit[*x].unwrap();
                let orb2 = to_orbit[alltigs2[ms[0]].1].unwrap();

                // Test for donor mixing.

                if !ctl.clono_filt_opt.donor {
                    let mut donors = vec![Vec::<Option<usize>>::new(); 2];
                    let orbs = [&orb1, &orb2];
                    for (pass, orb) in orbs.iter().enumerate() {
                        for id in orbits[**orb as usize].iter() {
                            let ex = &exact_clonotypes[info[*id as usize].clonotype_id];
                            for i in 0..ex.clones.len() {
                                donors[pass].push(ex.clones[i][0].donor_index);
                            }
                        }
                        unique_sort(&mut donors[pass]);
                    }
                    if donors[0] != donors[1] {
                        continue;
                    }
                }

                // Make join.

                eqo.join(orb1 as i32, orb2 as i32);
            }
        }
        let mut orbits2 = Vec::<Vec<i32>>::new();
        let mut repsx = Vec::<i32>::new();
        eqo.orbit_reps(&mut repsx);
        for i in 0..repsx.len() {
            let mut ox = Vec::<i32>::new();
            eqo.orbit(repsx[i], &mut ox);
            let mut o = Vec::<i32>::new();
            for j in 0..ox.len() {
                o.append(&mut orbits[ox[j] as usize].clone());
            }
            orbits2.push(o);
        }
        *orbits = orbits2;
    }
}
