// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file provides the single function join_exacts, which computes the equivalence relation
// on exact subclonotypes.
//
// Note that in principle the specificity of joining might be increased by using nonproductive
// contigs that represent the sequence of the "other" allele.  This does not look easy to
// execute.

use vdj_ann::{annotate, refx};

use self::annotate::print_annotations;
use self::refx::RefData;
use crate::join2::finish_join;
use crate::join_core::join_core;
use debruijn::dna_string::DnaString;
use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype, PotentialJoin};
use equiv::EquivRel;
use io_utils::{fwrite, fwriteln};
use itertools::Itertools;
use rayon::prelude::*;
use std::cmp::min;
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;
use vector_utils::{bin_member, erase_if, next_diff1_2};

pub fn join_exacts(
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    refdata: &RefData,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    mut join_info: &mut Vec<(usize, usize, bool, Vec<u8>)>,
    raw_joins: &mut Vec<(i32, i32)>,
    sr: &Vec<Vec<f64>>,
) -> EquivRel {
    // Run special option for joining by barcode identity.

    let timer1 = Instant::now();
    if ctl.join_alg_opt.bcjoin {
        let mut eq: EquivRel = EquivRel::new(info.len() as i32);
        let mut bcx = Vec::<(String, usize)>::new(); // {(barcode, info_index)}
        for i in 0..info.len() {
            let ex = &exact_clonotypes[info[i].clonotype_id];
            for j in 0..ex.clones.len() {
                bcx.push((ex.clones[j][0].barcode.clone(), i));
            }
        }
        bcx.sort();
        let mut i = 0;
        while i < bcx.len() {
            let j = next_diff1_2(&bcx, i as i32) as usize;
            for k in i + 1..j {
                eq.join(bcx[i].1 as i32, bcx[k].1 as i32);
            }
            i = j;
        }
        return eq;
    }

    // Find potential joins.

    let mut i = 0;
    let mut results = Vec::<(
        usize,                              // i
        usize,                              // j
        usize,                              // joins
        usize,                              // errors
        Vec<(usize, usize, bool, Vec<u8>)>, // log+ (index1, index2, err?, log)
        Vec<(usize, usize)>,                // joinlist
    )>::new();
    while i < info.len() {
        let mut j = i + 1;
        while j < info.len() {
            if info[j].lens != info[i].lens {
                break;
            }
            j += 1;
        }
        results.push((
            i,
            j,
            0,
            0,
            Vec::<(usize, usize, bool, Vec<u8>)>::new(),
            Vec::<(usize, usize)>::new(),
        ));
        i = j;
    }
    if !ctl.silent {
        println!("comparing {} simple clonotypes", info.len());
    }
    ctl.perf_stats(&timer1, "join setup");
    let timer2 = Instant::now();

    let joinf = |r: &mut (
        usize,
        usize,
        usize,
        usize,
        Vec<(usize, usize, bool, Vec<u8>)>,
        Vec<(usize, usize)>,
    )| {
        let (i, j) = (r.0, r.1);
        let joins = &mut r.2;
        let errors = &mut r.3;
        let logplus = &mut r.4;
        let mut pot = Vec::<PotentialJoin>::new();

        // Main join logic.  If you change par_iter_mut to iter_mut above, and run happening,
        // a lot of time shows up on the following line.  If further you manually inline join_core
        // here, then instead the time shows up on the results.iter_mut line.  It's not clear
        // what this means.

        join_core(
            is_bcr,
            i,
            j,
            ctl,
            exact_clonotypes,
            info,
            to_bc,
            sr,
            &mut pot,
        );

        // Run two passes.

        for _ in 0..2 {
            // Form the equivalence relation implied by the potential joins.

            let mut eq: EquivRel = EquivRel::new((j - i) as i32);
            for pj in 0..pot.len() {
                let (k1, k2) = (pot[pj].k1, pot[pj].k2);
                eq.join((k1 - i) as i32, (k2 - i) as i32);
            }

            // Impose a higher bar on joins that involve only two cells. (not documented)

            let mut to_pot = vec![Vec::<usize>::new(); j - i];
            for pj in 0..pot.len() {
                let k1 = pot[pj].k1;
                to_pot[k1 - i].push(pj);
            }
            let mut to_delete = vec![false; pot.len()];
            let mut reps = Vec::<i32>::new();
            eq.orbit_reps(&mut reps);
            for s in 0..reps.len() {
                // Examine a potential orbit.

                let mut x = Vec::<i32>::new();
                eq.orbit(reps[s], &mut x);

                // Count the number of cells in the orbit.

                let mut ncells = 0;
                for t in 0..x.len() {
                    let k = x[t] as usize + i;
                    let mult = exact_clonotypes[info[k].clonotype_index].ncells();
                    ncells += mult;
                }

                // Impose more stringent conditions if number of cells is two.

                if ncells == 2 && x.len() == 2 {
                    let (k1, k2) = (x[0] as usize + i, x[1] as usize + i);
                    let k = min(k1, k2);
                    for pj in to_pot[k - i].iter() {
                        let cd = pot[*pj].cd;
                        let shares = &pot[*pj].shares;

                        // Impose higher bar.

                        let min_shares = shares.iter().min().unwrap();
                        if cd > *min_shares / 2 && !ctl.join_alg_opt.easy {
                            to_delete[*pj] = true;
                        }
                    }
                }
            }
            erase_if(&mut pot, &to_delete);
        }

        // Analyze potential joins.

        let mut eq: EquivRel = EquivRel::new((j - i) as i32);
        for pj in 0..pot.len() {
            let k1 = pot[pj].k1;
            let k2 = pot[pj].k2;
            let nrefs = pot[pj].nrefs;
            let cd = pot[pj].cd;
            let diffs = pot[pj].diffs;
            let bcs1 = &pot[pj].bcs1;
            let bcs2 = &pot[pj].bcs2;
            let shares = &pot[pj].shares;
            let indeps = &pot[pj].indeps;
            let shares_details = &pot[pj].shares_details;
            let share_pos_v = &pot[pj].share_pos_v;
            let share_pos_j = &pot[pj].share_pos_j;
            let score = pot[pj].score;
            let err = pot[pj].err;
            let p1 = pot[pj].p1;
            let mult = pot[pj].mult;

            // Do nothing if join could have no effect on equivalence relation.

            if !ctl.force && (eq.class_id((k1 - i) as i32) == eq.class_id((k2 - i) as i32)) {
                continue;
            }

            // Save join and tally stats.

            r.5.push((k1, k2));
            *joins += 1;
            if err {
                *errors += 1;
            }
            eq.join((k1 - i) as i32, (k2 - i) as i32);

            // Decide if we're going to print.

            if ctl.join_print_opt.quiet || (!err && *joins % ctl.join_print_opt.pfreq != 0) {
                continue;
            }

            // Print join.

            let mut log = Vec::<u8>::new();
            fwriteln!(
                log,
                "\n====================================================\
                 ==============================================="
            );
            if err {
                fwriteln!(log, "\nJOIN ERROR");
            }
            let (mut lena1, mut lena2) = (Vec::<String>::new(), Vec::<String>::new());
            for l1 in info[k1].origin.iter() {
                lena1.push(ctl.origin_info.dataset_id[*l1].clone());
            }
            for l2 in info[k2].origin.iter() {
                lena2.push(ctl.origin_info.dataset_id[*l2].clone());
            }
            fwriteln!(
                log,
                "\n[{}] = {}.clono{} * {}.clono{}",
                diffs,
                lena1.iter().format(","),
                info[k1].clonotype_id,
                lena2.iter().format(","),
                info[k2].clonotype_id
            );

            fwrite!(log, "ref v segs1 = ");
            for (j, t) in info[k1].vsids.iter().enumerate() {
                if j > 0 {
                    fwrite!(log, ", ");
                }
                fwrite!(log, "{}={}", t, refdata.transcript[*t]);
            }
            fwrite!(log, "\nref v segs2 = ");
            for (j, t) in info[k2].vsids.iter().enumerate() {
                if j > 0 {
                    fwrite!(log, ", ");
                }
                fwrite!(log, "{}={}", t, refdata.transcript[*t]);
            }
            fwriteln!(log, "");
            for l1 in info[k1].origin.iter() {
                fwriteln!(
                    log,
                    "{} = {}",
                    ctl.origin_info.dataset_id[*l1],
                    ctl.origin_info.descrips[*l1]
                );
            }
            for l2 in info[k2].origin.iter() {
                fwriteln!(
                    log,
                    "{} = {}",
                    ctl.origin_info.dataset_id[*l2],
                    ctl.origin_info.descrips[*l2]
                );
            }
            let ci1 = info[k1].clonotype_index;
            let ci2 = info[k2].clonotype_index;
            let mut mega1 = String::new();
            for j in 0..exact_clonotypes[ci1].share.len() {
                let x = &exact_clonotypes[ci1].share[j];
                if j > 0 {
                    mega1 += ";";
                }
                mega1 += format!("{}:{}", x.chain_type, x.cdr3_aa).as_str();
            }
            let mut mega2 = String::new();
            for j in 0..exact_clonotypes[ci2].share.len() {
                let x = &exact_clonotypes[ci2].share[j];
                if j > 0 {
                    mega2 += ";";
                }
                mega2 += format!("{}:{}", x.chain_type, x.cdr3_aa).as_str();
            }
            fwriteln!(
                log,
                "{}, mult = {}",
                mega1,
                exact_clonotypes[info[k1].clonotype_index].ncells()
            );
            if ctl.join_print_opt.show_bc {
                fwriteln!(log, "bcs = {}", bcs1.iter().format(" "));
            }
            fwriteln!(
                log,
                "{}, mult = {}",
                mega2,
                exact_clonotypes[info[k2].clonotype_index].ncells()
            );
            if ctl.join_print_opt.show_bc {
                fwriteln!(log, "bcs = {}", bcs2.iter().format(" "));
            }
            fwriteln!(log, "score = {:.1}", score);
            for u in 0..nrefs {
                fwrite!(
                    log,
                    "versus ref {}, shares = {}, indeps = {}",
                    u + 1,
                    shares[u],
                    indeps[u]
                );
                fwriteln!(
                    log,
                    " (shares breakdown = V{}, J{}; V{}, J{})",
                    shares_details[u][0],
                    shares_details[u][1],
                    shares_details[u][2],
                    shares_details[u][3]
                );
                fwriteln!(
                    log,
                    "              share pos v = {}",
                    share_pos_v[u].iter().format(",")
                );
            }
            if cd < 0 {
                fwriteln!(log, "CDR3 length difference");
            } else {
                fwriteln!(log, "CDR3 nucleotide diffs = {}", cd);
                fwriteln!(log, "in amino acid space:");
                fwriteln!(log, "{}", mega1);
                fwriteln!(log, "{}", mega2);
            }
            fwriteln!(
                log,
                "p1 = prob of getting so many shares by accident = {}",
                p1
            );
            fwriteln!(
                log,
                "mult = CDR3: partial_bernoulli_sum(3 * cn, cd as usize) = {}",
                mult
            );
            fwriteln!(log, "score = p1 * mult = {}", p1 * mult);

            // Show difference patterns.  And x denotes a different base.  A ▓ denotes an
            // equal base that differs from the reference.  Otherwise - is shown.

            let nchains = info[k1].lens.len();
            for m in 0..nchains {
                let (tig1, tig2) = (&info[k1].tigs[m], &info[k2].tigs[m]);
                fwriteln!(log, "difference pattern for chain {}", m + 1);
                for i in 0..tig1.len() {
                    if i > 0 && i % 80 == 0 {
                        fwrite!(log, "\n");
                    }
                    if tig1[i] == tig2[i] {
                        if bin_member(&share_pos_v[m], &i)
                            || bin_member(&share_pos_j[m], &(tig1.len() - i))
                        {
                            fwrite!(log, "▓");
                        } else {
                            fwrite!(log, "-");
                        }
                    } else {
                        fwrite!(log, "x");
                    }
                }
                fwriteln!(log, "");
            }

            // Show sequences and annotations.

            if ctl.join_print_opt.seq || ctl.join_print_opt.ann || ctl.join_print_opt.ann0 {
                let nchains = info[k1].lens.len();
                for m in 0..nchains {
                    let (tig1, tig2) = (&info[k1].tigs[m], &info[k2].tigs[m]);
                    let ex1 = &exact_clonotypes[info[k1].clonotype_index];
                    let otig1 =
                        DnaString::from_acgt_bytes(&ex1.share[info[k1].exact_cols[m]].full_seq);
                    let ex2 = &exact_clonotypes[info[k2].clonotype_index];
                    let otig2 =
                        DnaString::from_acgt_bytes(&ex2.share[info[k2].exact_cols[m]].full_seq);
                    if ctl.join_print_opt.seq {
                        fwriteln!(log, "\nchain {}, tig 1 = {}", m + 1, otig1.to_string());
                    }
                    if ctl.join_print_opt.ann0 {
                        // somewhat broken for the moment, because tig1 could have - characters
                        if !info[k1].has_del[m] {
                            fwriteln!(log, "chain {}, tig 1", m + 1);
                            let t1 = DnaString::from_acgt_bytes(tig1);
                            print_annotations(&t1, refdata, &mut log, false, true, false);
                        }
                    }
                    if ctl.join_print_opt.ann {
                        fwriteln!(log, "chain {}, tig 1", m + 1);
                        print_annotations(&otig1, refdata, &mut log, false, true, false);
                    }
                    if ctl.join_print_opt.seq {
                        fwriteln!(log, "\nchain {}, tig 2 = {}", m + 1, otig2.to_string());
                    }
                    if ctl.join_print_opt.ann0 {
                        // somewhat broken for the moment, because tig2 could have - characters
                        if !info[k2].has_del[m] {
                            fwriteln!(log, "chain {}, tig 2", m + 1);
                            let t2 = DnaString::from_acgt_bytes(tig2);
                            print_annotations(&t2, refdata, &mut log, false, true, false);
                        }
                    }
                    if ctl.join_print_opt.ann {
                        fwriteln!(log, "chain {}, tig 2", m + 1);
                        print_annotations(&otig2, refdata, &mut log, false, true, false);
                    }
                }
            }
            // not sure why this logging is here, so turned off for now
            /*
            if ctl.join_print_opt.seq {
                for x in 0..info[k1].lens.len() {
                    fwriteln!(log, "{}", strme(&info[k1].tigs[x]));
                }
                fwriteln!(log, "----------------------------------------");
                for x in 0..info[k2].lens.len() {
                    fwriteln!(log, "{:?}", strme(&info[k2].tigs[x]));
                }
            }
            */
            logplus.push((info[k1].clonotype_index, info[k2].clonotype_index, err, log));
        }
    };

    results.par_iter_mut().for_each(joinf);

    ctl.perf_stats(&timer2, "in main part of join");
    for l in 0..results.len() {
        for j in 0..results[l].5.len() {
            raw_joins.push((results[l].5[j].0 as i32, results[l].5[j].1 as i32));
        }
    }
    finish_join(ctl, info, &results, &mut join_info)
}
