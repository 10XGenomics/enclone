// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This is the actual grouping code.

use enclone_core::defs::*;
use equiv::EquivRel;
use vdj_ann::refx::*;
use vector_utils::*;

pub fn grouper(
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
) -> Vec<Vec<i32>> {
    // Group clonotypes.

    let mut e: EquivRel = EquivRel::new(exacts.len() as i32);
    if ctl.clono_group_opt.heavy_cdr3_aa {
        let mut all = Vec::<(String, usize)>::new();
        for i in 0..exacts.len() {
            for x in exacts[i].iter() {
                for m in 0..exact_clonotypes[*x].share.len() {
                    let y = &exact_clonotypes[*x].share[m];
                    if y.left {
                        all.push((y.cdr3_aa.clone(), i));
                    }
                }
            }
        }
        all.sort();
        let mut i = 0;
        while i < all.len() {
            let j = next_diff1_2(&all, i as i32) as usize;
            for k in i + 1..j {
                e.join(all[i].1 as i32, all[k].1 as i32);
            }
            i = j;
        }
    }
    if ctl.clono_group_opt.vj_refname || ctl.clono_group_opt.vj_refname_strong {
        let mut all = Vec::<(Vec<String>, usize)>::new();
        for i in 0..exacts.len() {
            let ex = &exact_clonotypes[exacts[i][0]];
            let mut s = Vec::<String>::new();
            for j in 0..ex.share.len() {
                s.push(refdata.name[ex.share[j].v_ref_id].clone());
                s.push(refdata.name[ex.share[j].j_ref_id].clone());
            }
            s.sort();
            all.push((s, i));
        }
        // Note duplication with above code.
        all.sort();
        let mut i = 0;
        while i < all.len() {
            let j = next_diff1_2(&all, i as i32) as usize;
            for k in i + 1..j {
                let m1 = all[i].1;
                let m2 = all[k].1;
                if ctl.clono_group_opt.vj_refname_strong {
                    let ex1 = &exact_clonotypes[exacts[m1][0]];
                    let ex2 = &exact_clonotypes[exacts[m2][0]];
                    let mut lens1 = Vec::<(usize, usize)>::new();
                    let mut lens2 = Vec::<(usize, usize)>::new();
                    for j in 0..ex1.share.len() {
                        lens1.push((ex1.share[j].seq_del.len(), ex1.share[j].cdr3_aa.len()));
                    }
                    for j in 0..ex2.share.len() {
                        lens2.push((ex2.share[j].seq_del.len(), ex2.share[j].cdr3_aa.len()));
                    }
                    lens1.sort();
                    lens2.sort();
                    if lens1 != lens2 {
                        continue;
                    }
                }
                e.join(m1 as i32, m2 as i32);
            }
            i = j;
        }
    }
    let mut greps = Vec::<i32>::new();
    e.orbit_reps(&mut greps);

    // Gather groups and sort so that larger groups (as measured by cells) come first.

    let mut groups = Vec::<Vec<i32>>::new();
    let mut grepsn = Vec::<usize>::new();
    for i in 0..greps.len() {
        let mut o = Vec::<i32>::new();
        e.orbit(greps[i], &mut o);
        if o.len() < ctl.clono_group_opt.min_group {
            continue;
        }
        groups.push(o.clone());
        let mut n = 0;
        for j in 0..o.len() {
            let x = o[j] as usize;
            let s = &exacts[x];
            for k in 0..s.len() {
                n += exact_clonotypes[s[k]].clones.len();
            }
        }
        grepsn.push(n);
    }
    sort_sync2(&mut grepsn, &mut groups);
    groups.reverse();
    groups
}
