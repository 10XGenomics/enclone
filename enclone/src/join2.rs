// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file provides the tail end code for join.rs, plus a small function used there.

use enclone_core::defs::*;
use equiv::EquivRel;
use stats_utils::*;
use std::time::Instant;
use vector_utils::*;

// partial_bernoulli_sum( n, k ): return sum( choose(n,i), i = 0..=k ).
//
// Beware of overflow.

pub fn partial_bernoulli_sum(n: usize, k: usize) -> f64 {
    assert!(n >= 1);
    assert!(k <= n);
    let mut sum = 0.0;
    let mut choose = 1.0;
    for i in 0..=k {
        sum += choose;
        choose *= (n - i) as f64;
        choose /= (i + 1) as f64;
    }
    sum
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn finish_join(
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    results: &Vec<(
        usize,
        usize,
        usize,
        usize,
        Vec<(usize, usize, bool, Vec<u8>)>,
        Vec<(usize, usize)>,
    )>,
    join_info: &mut Vec<(usize, usize, bool, Vec<u8>)>,
) -> EquivRel {
    // Tally results.

    let (mut joins, mut errors) = (0, 0);
    let timer3 = Instant::now();
    for l in 0..results.len() {
        joins += results[l].2;
        errors += results[l].3;
        for i in 0..results[l].4.len() {
            let u1 = results[l].4[i].0;
            let u2 = results[l].4[i].1;
            let err = results[l].4[i].2;
            let log = results[l].4[i].3.clone();
            join_info.push((u1, u2, err, log));
        }
    }
    if !ctl.silent {
        println!("{} joins", joins);
        if ctl.origin_info.donors > 1 {
            println!("{} errors", errors);
        }
    }

    // Make equivalence relation.

    let mut eq: EquivRel = EquivRel::new(info.len() as i32);
    for l in 0..results.len() {
        for j in 0..results[l].5.len() {
            eq.join(results[l].5[j].0 as i32, results[l].5[j].1 as i32);
        }
    }

    // Join orbits that cross subclones of a clone.  This arose because we split up multi-chain
    // clonotypes into two-chain clonotypes.

    let mut ox = Vec::<(usize, i32)>::new();
    for i in 0..info.len() {
        let x: &CloneInfo = &info[i];
        ox.push((x.clonotype_id, eq.class_id(i as i32)));
    }
    ox.sort();
    let mut i = 0;
    while i < ox.len() {
        let j = next_diff1_2(&ox, i as i32) as usize;
        for k in i..j - 1 {
            eq.join(ox[k].1, ox[k + 1].1);
        }
        i = j;
    }

    // Merge onesies where totally unambiguous.  Possibly inefficient and should optimize.

    if ctl.join_alg_opt.merge_onesies {
        let mut ncells_total = 0;
        for i in 0..exact_clonotypes.len() {
            ncells_total += exact_clonotypes[i].ncells();
        }
        let mut onesies = Vec::<usize>::new();
        for i in 0..info.len() {
            if info[i].tigs.len() == 1 {
                onesies.push(i);
            }
        }
        let mut alltigs2 = Vec::<(Vec<u8>, usize)>::new();
        for i in 0..info.len() {
            if info[i].tigs.len() >= 2 {
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
                eq.join(*x as i32, alltigs2[ms[0]].1 as i32);
            }
        }
    }

    // Tally whitelist contamination.
    // WARNING: THIS ONLY WORKS IF YOU RUN WITH CLONES=1 AND NO OTHER FILTERS.

    let mut white = ctl.clono_filt_opt.whitef;
    for j in 0..ctl.clono_print_opt.cvars.len() {
        if ctl.clono_print_opt.cvars[j] == "white" {
            white = true;
        }
    }
    if white {
        let mut bads = 0;
        let mut denom = 0;
        for i in 0..results.len() {
            bads += results[i].2;
            denom += results[i].3;
        }
        let bad_rate = percent_ratio(bads, denom);
        println!("whitelist contamination rate = {:.2}%", bad_rate);
    }
    ctl.perf_stats(&timer3, "in tail of join");
    eq
}
