// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::defs::{CloneInfo, EncloneControl, ExactClonotype, PotentialJoin};
use debruijn::{
    dna_string::{ndiffs, DnaString},
    Mer,
};
use stats_utils::abs_diff;
use std::collections::HashMap;
use stirling_numbers::p_at_most_m_distinct_in_sample_of_x_from_n;
use vector_utils::{meet, unique_sort};

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

pub fn join_one(
    is_bcr: bool,
    k1: usize,
    k2: usize,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<f64>>,
    pot: &mut Vec<PotentialJoin>,
) -> bool {
    // Do not merge onesies or foursies with anything.  Deferred until later.
    // Note that perhaps some foursies should be declared doublets and deleted.
    // Note onesies merging above is turned off so this appears to be moot.

    let (clono1, clono2) = (info[k1].clonotype_id, info[k2].clonotype_id);
    let chains1 = exact_clonotypes[clono1].share.len();
    let chains2 = exact_clonotypes[clono2].share.len();
    if !(2..=3).contains(&chains1) || chains2 < 2 || chains2 > 3 {
        return false;
    }
    // NEED FOR THIS SEEMS LIKE A BUG:
    if info[k1].vs.len() == 1 || info[k2].vs.len() == 4 {
        return false;
    }
    if info[k1].vs.len() > 2 {
        return false;
    }

    // Require that CDR3s have the same length.  Ugly.
    // First part should be a tautology.

    let (x1, x2) = (&info[k1].cdr3s, &info[k2].cdr3s);
    if x1.len() != x2.len() {
        return false;
    }
    for i in 0..x1.len() {
        if x1[i].len() != x2[i].len() {
            return false;
        }
    }

    // Compute number of differences.

    let mut diffs = 0_usize;
    for x in 0..info[k1].lens.len() {
        if !info[k1].has_del[x] && !info[k2].has_del[x] {
            // A great deal of time is spent in the call to ndiffs.  Notes on this:
            // 1. It is slower than if the computation is done outside
            //    the ndiffs function.  This is mysterious but must have something to
            //    do with the storage of the 256-byte lookup table.
            // 2. Adding #[inline(always)] in front of the ndiffs function definition
            //    doesn't help.
            // 3. Adding a bounds test for diffs > ctl.heur.max_diffs inside the ndiffs
            //    function doesn't help, whether placed in the inner loop or the other
            //    loop.
            diffs += ndiffs(&info[k1].tigsp[x], &info[k2].tigsp[x]);
        } else {
            for j in 0..info[k1].tigs[x].len() {
                if info[k1].tigs[x][j] != info[k2].tigs[x][j] {
                    diffs += 1;
                }
            }
        }
    }

    // Another test for acceptable join.

    if diffs > ctl.heur.max_diffs {
        return false;
    }
    if !is_bcr && diffs > 5 {
        return false;
    }

    // Unless MIX_DONORS specified, do not join across donors.
    // And test for error.
    //
    // WARNING!  There are actually two cases: where an individual exact subclonotype
    // itself crosses donors, and where we cross donors in making a join.  Note that
    // the former case is most improbable, unless there is cross-sample contamination.
    // And if that did happen, the output would be confusing and might have a greatly
    // exaggerated number of fails.

    let (mut donors1, mut donors2) = (Vec::<usize>::new(), Vec::<usize>::new());
    let ex1 = &exact_clonotypes[info[k1].clonotype_index];
    let ex2 = &exact_clonotypes[info[k2].clonotype_index];
    for j in 0..ex1.clones.len() {
        if ex1.clones[j][0].donor_index.is_some() {
            donors1.push(ex1.clones[j][0].donor_index.unwrap());
        }
    }
    for j in 0..ex2.clones.len() {
        if ex2.clones[j][0].donor_index.is_some() {
            donors2.push(ex2.clones[j][0].donor_index.unwrap());
        }
    }
    unique_sort(&mut donors1);
    unique_sort(&mut donors2);
    if !ctl.clono_filt_opt_def.donor
        && !donors1.is_empty()
        && !donors2.is_empty()
        && donors1 != donors2
    {
        return false;
    }
    let err = donors1 != donors2 || donors1.len() != 1 || donors2.len() != 1;

    // Analyze the two clonotypes versus the reference.  First traverse the reference
    // sequences.  Either we use the references for k1 or the references for k2, but
    // these are nearly always the same.

    let mut nrefs = 1;
    for m in 0..2 {
        if info[k1].vs[m] != info[k2].vs[m] || info[k1].js[m] != info[k2].js[m] {
            nrefs = 2;
        }
    }
    let mut shares = vec![0; nrefs]; // shared mutations from reference
    let mut indeps = vec![0; nrefs]; // independent mutations from reference
    let mut total = vec![vec![0; 2]; nrefs]; // total differences from reference
    let mut shares_details = vec![vec![0; 4]; nrefs];
    let mut share_pos_v = vec![Vec::<usize>::new(); 2];
    let mut share_pos_j = vec![Vec::<usize>::new(); 2];
    for u in 0..nrefs {
        let k: usize;
        if u == 0 {
            k = k1;
        } else {
            k = k2;
        }

        // Traverse the chains in the clonotype.

        let nchains = info[k1].lens.len();
        for m in 0..nchains {
            let (tig1, tig2) = (&info[k1].tigs[m], &info[k2].tigs[m]);

            // Traverse the two segments (V and J).

            for si in 0..2 {
                let seg: &DnaString;
                if si == 0 {
                    seg = &info[k].vs[m];
                } else {
                    seg = &info[k].js[m];
                }
                let mut ref_trim = ctl.heur.ref_v_trim;
                if si == 1 {
                    ref_trim = ctl.heur.ref_j_trim;
                }
                for p in 0..seg.len() - ref_trim {
                    let (t1, t2);
                    let r;
                    if si == 0 {
                        // Ugly bailout arising very rarely if the two reference
                        // sequences have different lengths.
                        if p >= tig1.len() || p >= tig2.len() {
                            return false;
                        }
                        t1 = tig1[p];
                        t2 = tig2[p];
                        // r = seg.get(p);
                        let rx = seg.get(p);
                        if rx == 0 {
                            r = b'A';
                        } else if rx == 1 {
                            r = b'C';
                        } else if rx == 2 {
                            r = b'G';
                        } else {
                            r = b'T';
                        }
                    } else {
                        t1 = tig1[tig1.len() - p - 1];
                        t2 = tig2[tig2.len() - p - 1];
                        // r = seg.get( seg.len() - p - 1 );
                        let rx = seg.get(seg.len() - p - 1);
                        if rx == 0 {
                            r = b'A';
                        } else if rx == 1 {
                            r = b'C';
                        } else if rx == 2 {
                            r = b'G';
                        } else {
                            r = b'T';
                        }
                    }
                    if t1 == t2 && t1 != r {
                        shares[u] += 1;
                        shares_details[u][2 * m + si] += 1;
                        if si == 0 {
                            share_pos_v[m].push(p);
                        } else {
                            share_pos_j[m].push(p);
                        }
                    } else if (t1 == r && t2 != r) || (t2 == r && t1 != r) {
                        indeps[u] += 1;
                    } else if t1 != r && t2 != r {
                        indeps[u] += 2;
                    }
                    if t1 != r {
                        total[u][0] += 1;
                    }
                    if t2 != r {
                        total[u][1] += 1;
                    }
                }
            }
        }
    }

    // Don't allow different references if one is strongly favored.
    // (not documented)

    if nrefs == 2 {
        for m in 0..2 {
            if abs_diff(total[0][m], total[1][m]) > ctl.heur.max_degradation {
                return false;
            }
        }
    }

    // Compute junction diffs.  Ugly.

    let mut cd = 0_isize;
    for l in 0..x1.len() {
        for m in 0..x1[l].len() {
            if x1[l].as_bytes()[m] != x2[l].as_bytes()[m] {
                cd += 1;
            }
        }
    }

    // Cap CDR3 diffs.

    if cd > ctl.join_alg_opt.max_cdr3_diffs as isize || (!is_bcr && cd > 0) {
        return false;
    }

    // Another test for acceptable join.  (not fully documented)

    let min_shares = shares.iter().min().unwrap();
    let min_indeps = indeps.iter().min().unwrap();

    // Reject if barcode overlap. (not documented)

    let (mut bcs1, mut bcs2) = (Vec::<String>::new(), Vec::<String>::new());
    for origin in info[k1].origin.iter() {
        bcs1.append(&mut to_bc[&(*origin, info[k1].clonotype_id)].clone());
    }
    for origin in info[k2].origin.iter() {
        bcs2.append(&mut to_bc[&(*origin, info[k2].clonotype_id)].clone());
    }
    unique_sort(&mut bcs1);
    unique_sort(&mut bcs2);
    if meet(&bcs1, &bcs2) {
        return false;
    }

    // Estimate the probability p1 that drawing k = min_indeps + 2 * min_shares
    // objects from n = 3 * (sum of VJ contig lengths) yields d = min_shares or
    // more duplicates.

    let n = 3 * (info[k1].tigs[0].len() + info[k1].tigs[1].len());
    let k = *min_indeps + 2 * *min_shares;
    let d = *min_shares;
    let p1 = p_at_most_m_distinct_in_sample_of_x_from_n((k - d) as usize, k as usize, n, sr);
    assert!(!p1.is_infinite()); // TODO: IS THIS SAFE?

    // Multiply by 80^cd, or if using old version, the number of DNA sequences that differ from
    // the given CDR3 sequences on <= cd bases.  This is sum( choose(3cn, m), m = 0..=cd ).

    let mult;
    if ctl.join_alg_opt.old_mult {
        let mut cn = 0;
        for l in 0..x1.len() {
            cn += x1[l].len();
        }
        mult = partial_bernoulli_sum(3 * cn, cd as usize);
        assert!(!mult.is_infinite()); // TODO: IS THIS SAFE?
    } else {
        mult = ctl.join_alg_opt.mult_pow.powi(cd as i32);
    }

    // Compute score.

    let mut score = p1 * mult;

    // Test for concentration of SHM in the junction regions.

    if cd as f64 >= ctl.join_alg_opt.cdr3_mult * std::cmp::max(1, *min_indeps) as f64 {
        score = ctl.join_alg_opt.max_score + 1.0;
    }

    // Do not merge cells if they were assigned different light chain constant regions.

    if !ctl.join_alg_opt.old_light {
        for i in 0..info[k1].cdr3s.len() {
            let (j1, j2) = (info[k1].exact_cols[i], info[k2].exact_cols[i]);
            if !ex1.share[j1].left {
                if ex1.share[j1].c_ref_id.is_some() && ex2.share[j2].c_ref_id.is_some() {
                    if ex1.share[j1].c_ref_id.unwrap() != ex2.share[j2].c_ref_id.unwrap() {
                        score = ctl.join_alg_opt.max_score + 1.0;
                    }
                }
            }
        }
    }

    // Threshold on score.

    if score > ctl.join_alg_opt.max_score {
        return false;
    }

    // Save potential joins.  Note that this jacks up memory usage significantly,
    // so it would likely be more efficient to duplicate some of the computations
    // during the analysis phase.

    if !ctl.join_print_opt.show_bc {
        bcs1.clear();
        bcs2.clear();
    }
    pot.push(PotentialJoin {
        k1,
        k2,
        nrefs,
        cd,
        diffs,
        bcs1,
        bcs2,
        shares,
        indeps,
        shares_details,
        share_pos_v,
        share_pos_j,
        score,
        err,
        p1,
        mult,
    });
    true
}
