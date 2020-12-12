// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This is the actual grouping code.

use edit_distance::edit_distance;
use enclone_core::defs::*;
use equiv::EquivRel;
use rayon::prelude::*;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

pub fn grouper(
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    in_center: &Vec<bool>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
) -> Vec<Vec<(i32, String)>> {
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Case 1: symmetric grouping.

    if !ctl.clono_group_opt.asymmetric {
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

        let mut groups = Vec::<Vec<(i32, String)>>::new();
        let mut grepsn = Vec::<usize>::new();
        for i in 0..greps.len() {
            let mut o = Vec::<i32>::new();
            e.orbit(greps[i], &mut o);
            if o.len() < ctl.clono_group_opt.min_group {
                continue;
            }
            let mut z = Vec::<(i32, String)>::new();
            for j in 0..o.len() {
                z.push((o[j], String::new()));
            }
            groups.push(z);
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
        return groups;

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Case 2: asymmetric grouping.
    } else {
        // Define the center.

        let mut center = Vec::<usize>::new();
        for i in 0..exacts.len() {
            if in_center[i] {
                let mut d1 = false;
                let s = &exacts[i];
                for k in 0..s.len() {
                    let ex = &exact_clonotypes[s[k]];
                    for u in 0..ex.clones.len() {
                        if ctl.origin_info.donor_id[ex.clones[u][0].dataset_index] == "d1" {
                            d1 = true;
                        }
                    }
                }
                if d1 {
                    center.push(i);
                }
            }
        }

        // Define the distance from center clonotypes to other clonotypes.

        let infinity = 1_000_000_000.0_f64;
        let mut dist = Vec::<(usize, Vec<f64>)>::new();
        for i in 0..center.len() {
            dist.push((i, vec![infinity; exacts.len()]));
        }
        dist.par_iter_mut().for_each(|res| {
            let i = res.0;
            for j in 0..exacts.len() {
                let (c1, c2) = (&exacts[center[i]], &exacts[j]);
                for k1 in 0..c1.len() {
                    for k2 in 0..c2.len() {
                        let (ex1, ex2) = (&exact_clonotypes[c1[k1]], &exact_clonotypes[c2[k2]]);
                        let (mut heavy, mut light) = (infinity, infinity);
                        for m1 in 0..ex1.share.len() {
                            let cdr3_aa1 = &ex1.share[m1].cdr3_aa;
                            for m2 in 0..ex2.share.len() {
                                let cdr3_aa2 = &ex2.share[m2].cdr3_aa;
                                if ex1.share[m1].left && ex2.share[m2].left {
                                    let x = edit_distance(&cdr3_aa1, &cdr3_aa2) as f64;
                                    heavy = heavy.min(x);
                                }
                                if !ex1.share[m1].left && !ex2.share[m2].left {
                                    let x = edit_distance(&cdr3_aa1, &cdr3_aa2) as f64;
                                    light = light.min(x);
                                }
                            }
                        }
                        res.1[j] = res.1[j].min(heavy + light);
                    }
                }
            }
        });

        // Compute the groups.

        let mut groups = Vec::<Vec<(i32, String)>>::new();
        let bound = &ctl.clono_group_opt.asymmetric_dist_bound;
        let mut top_dist = None;
        if bound.starts_with("top=") {
            top_dist = Some(bound.after("top=").force_usize());
        }
        let mut max_dist = None;
        if bound.starts_with("max=") {
            max_dist = Some(bound.after("max=").force_f64());
        }
        for i in 0..center.len() {
            let mut g = Vec::<(i32, String)>::new();
            g.push((center[i] as i32, "distance = 0".to_string()));
            let mut id = Vec::<(f64, usize)>::new();
            for j in 0..exacts.len() {
                id.push((dist[i].1[j], j));
            }
            id.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for j in 0..id.len() {
                if top_dist.is_some() && j > top_dist.unwrap() {
                    break;
                }
                let d = id[j].0;
                if max_dist.is_some() && d > max_dist.unwrap() {
                    break;
                }
                if id[j].1 != center[i] {
                    g.push((id[j].1 as i32, format!("distance = {}", d)));
                }
            }
            if g.len() >= ctl.clono_group_opt.min_group {
                groups.push(g);
            }
        }
        groups
    }
}
