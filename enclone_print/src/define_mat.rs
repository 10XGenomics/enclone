// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype, PotentialJoin};
use enclone_core::join_one::join_one;
use equiv::EquivRel;
use qd::Double;
use std::cmp::max;
use std::collections::HashMap;
use vdj_ann::refx::RefData;
use vector_utils::{bin_position, next_diff12_3, next_diff1_3, unique_sort};

// Define an equivalence relation on the chains, introducing connections defined by the
// raw joins.  Also join where there are identical V..J sequences.

fn joiner(
    infos: &Vec<usize>,
    info: &Vec<CloneInfo>,
    to_exacts: &HashMap<usize, usize>,
    raw_joinsx: &Vec<Vec<usize>>,
    chains: &Vec<(usize, usize)>,
    seq_chains: &Vec<(Vec<u8>, usize, usize)>,
) -> EquivRel {
    let mut e = EquivRel::new(chains.len() as i32);
    for i1 in 0..infos.len() {
        let j1 = infos[i1];
        let u1 = info[j1].clonotype_index;
        let v1 = to_exacts[&u1];
        let m1s = &info[j1].exact_cols;
        for i2 in raw_joinsx[i1].iter() {
            let j2 = infos[*i2];
            let u2 = info[j2].clonotype_index;
            let v2 = to_exacts[&u2];
            let m2s = &info[j2].exact_cols;
            for j in 0..2 {
                let z1 = bin_position(chains, &(v1, m1s[j]));
                let z2 = bin_position(chains, &(v2, m2s[j]));
                e.join(z1, z2);
            }
        }
    }
    let mut i = 0;
    while i < seq_chains.len() {
        let j = next_diff1_3(seq_chains, i as i32) as usize;
        for k in i + 1..j {
            let (x1, x2) = (&seq_chains[i], &seq_chains[k]);
            let z1 = bin_position(chains, &(x1.1, x1.2));
            let z2 = bin_position(chains, &(x2.1, x2.2));
            e.join(z1, z2);
        }
        i = j;
    }
    e
}

// This generates a cols x nexacts matrices for a given clonotype, where cols is defined by the
// algorithm, and is the number of columns (chains) in the clonotype table.

pub fn define_mat(
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &Vec<Vec<Double>>,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    exacts: &Vec<usize>,
    od: &Vec<(Vec<usize>, usize, i32)>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
    refdata: &RefData,
) -> Vec<Vec<Option<usize>>> {
    // Define map of indices into exacts.

    let nexacts = exacts.len();
    let mut to_exacts = HashMap::<usize, usize>::new();
    for u in 0..nexacts {
        to_exacts.insert(exacts[u], u);
    }

    // Get the info indices corresponding to this clonotype.

    let mut infos = Vec::<usize>::new();
    for i in 0..od.len() {
        let x = od[i].2 as usize;
        if to_exacts.contains_key(&info[x].clonotype_index) {
            infos.push(x);
        }
    }
    infos.sort_unstable();

    // Define map of exacts to infos.

    let mut to_infos = vec![Vec::<usize>::new(); nexacts];
    for i in 0..infos.len() {
        let u = to_exacts[&info[infos[i]].clonotype_index];
        to_infos[u].push(i);
    }

    // Form the set of all chains that appear in an exact subclonotypes of this clonotype, and
    // also track the V..J sequences for the chains.

    let mut chains = Vec::<(usize, usize)>::new();
    let mut seq_chains = Vec::<(Vec<u8>, usize, usize)>::new();
    for u in 0..nexacts {
        let ex = &exact_clonotypes[exacts[u]];
        for m in 0..ex.share.len() {
            chains.push((u, m));
            seq_chains.push((ex.share[m].seq.clone(), u, m));
        }
    }
    seq_chains.sort();

    // Gather the raw joins.

    let mut raw_joinsx = vec![Vec::<usize>::new(); infos.len()];
    for i1 in 0..infos.len() {
        let j1 = infos[i1];
        for x in raw_joins[j1].iter() {
            let i2 = bin_position(&infos, &*x);
            if i2 >= 0 {
                raw_joinsx[i1].push(i2 as usize);
            }
        }
    }

    // Look for additional raw joins.  In the join stage, for efficiency reasons, we don't try to
    // make a join if two info elements are already in the same equivalence class.  This causes
    // us to miss raw joins for two reasons:
    // • One reason is that where we have two exact subclonotypes, and at least one one has more
    //   than two chains, then we may miss a raw join that is needed here.
    // • Another reason is that exact subclonotypes may have been deleted since the original
    //   equivalence relation was built.  This we address partially.

    let mut extras = Vec::<(usize, usize)>::new();
    for i1 in 0..raw_joinsx.len() {
        for i2 in raw_joinsx[i1].iter() {
            let i2 = *i2;
            let (j1, j2) = (infos[i1], infos[i2]);
            let (u1, u2) = (
                info[j1].clonotype_index as usize,
                info[j2].clonotype_index as usize,
            );
            let (ex1, ex2) = (&exact_clonotypes[u1], &exact_clonotypes[u2]);
            let (v1, v2) = (to_exacts[&u1], to_exacts[&u2]);
            if ex1.share.len() > 2 || ex2.share.len() > 2 {
                let (s1, s2) = (&to_infos[v1], &to_infos[v2]);
                for k1 in s1.iter() {
                    for k2 in s2.iter() {
                        let (k1, k2) = (*k1, *k2);
                        if (k1 == i1 && k2 == i2) || (k1 == i2 && k2 == i1) {
                            continue;
                        }
                        let (l1, l2) = (infos[k1], infos[k2]);
                        if info[l1].lens == info[l2].lens {
                            let mut pot = Vec::<PotentialJoin>::new();
                            if join_one(
                                is_bcr,
                                l1,
                                l2,
                                ctl,
                                exact_clonotypes,
                                info,
                                to_bc,
                                sr,
                                &mut pot,
                                &refdata,
                            ) {
                                extras.push((k1, k2));
                            }
                        }
                    }
                }
            }
        }
    }
    for x in extras.iter() {
        raw_joinsx[x.0].push(x.1);
    }

    // Define an initial equivalence relation on the chains, and get orbit representatives.

    let mut e = joiner(&infos, info, &to_exacts, &raw_joinsx, &chains, &seq_chains);
    let mut r = Vec::<i32>::new();
    e.orbit_reps(&mut r);

    // First for each pair of chain orbits with one "heavy" and one "light", pick some info
    // entries, if there are any.  This is effectively at random.  A parameter governs how
    // much we pick.

    let mut rxi = Vec::<(usize, usize, usize)>::new(); // (heavy orbit, light orbit, infos index)
    for i in 0..infos.len() {
        let z = &info[infos[i]];
        let u = z.clonotype_index;
        let v = to_exacts[&u];
        if z.exact_cols.len() != 2 {
            continue;
        }
        let (m1, m2) = (z.exact_cols[0], z.exact_cols[1]);
        let ex = &exact_clonotypes[u];
        if !ex.share[m1].left || ex.share[m2].left {
            continue; // maybe never happens
        }
        let p1 = e.class_id(bin_position(&chains, &(v, m1)));
        let p2 = e.class_id(bin_position(&chains, &(v, m2)));
        let q1 = bin_position(&r, &p1) as usize;
        let q2 = bin_position(&r, &p2) as usize;
        rxi.push((q1, q2, i));
    }
    rxi.sort_unstable();
    const MAX_USE: usize = 5; // knob set empirically
    let mut rxir = Vec::<(usize, usize, usize)>::new(); // (heavy orbit, light orbit, info index)
    let mut i = 0;
    while i < rxi.len() {
        let j = next_diff12_3(&rxi, i as i32) as usize;
        for k in i..j {
            if k < i + MAX_USE {
                rxir.push(rxi[k]);
            }
        }
        i = j;
    }

    // Now for each pair of these, if they are not effectively joined, attempt to join them.
    // This partially addresses the "second reason" described above.  It is partial because we
    // picked an info entry above at random, rather than trying them all.

    for f1 in rxir.iter() {
        for f2 in rxir.iter() {
            if f1.0 != f2.0 || f1.1 != f2.1 {
                let (i1, i2) = (infos[f1.2], infos[f2.2]);
                if info[i1].lens != info[i2].lens {
                    continue;
                }
                let mut pot = Vec::<PotentialJoin>::new();
                if join_one(
                    is_bcr,
                    i1,
                    i2,
                    ctl,
                    exact_clonotypes,
                    info,
                    to_bc,
                    sr,
                    &mut pot,
                    &refdata,
                ) {
                    e.join(r[f1.0], r[f2.0]);
                    e.join(r[f1.1], r[f2.1]);
                }
            }
        }
    }

    // Find the exact subclonotypes having three chains and list their orbits, allowing for order
    // to vary.

    let mut r = Vec::<i32>::new();
    e.orbit_reps(&mut r);
    let mut threes = Vec::<Vec<usize>>::new();
    let mut threesp = HashMap::<Vec<usize>, usize>::new();
    for u in 0..nexacts {
        let ex = &exact_clonotypes[exacts[u]];
        if ex.share.len() == 3 {
            let zs = [
                [0, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ];
            for z in zs.iter() {
                if ex.share[z[0]].left && !ex.share[z[2]].left {
                    let p1 = e.class_id(bin_position(&chains, &(u, z[0])));
                    let p2 = e.class_id(bin_position(&chains, &(u, z[1])));
                    let p3 = e.class_id(bin_position(&chains, &(u, z[2])));
                    let q1 = bin_position(&r, &p1) as usize;
                    let q2 = bin_position(&r, &p2) as usize;
                    let q3 = bin_position(&r, &p3) as usize;
                    threes.push(vec![q1, q2, q3]);
                    threesp.insert(vec![q1, q2, q3], u);
                }
            }
        }
    }
    unique_sort(&mut threes);

    // There is one more case to deal with.  This is where we have two exact subclonotypes, each
    // with three chains, and we joined two of their chains, but not the third.  And where the
    // join algorithm would not have joined the third.  In this case, if the third chains are
    // "close enough", we join them anyway.  As before, we only test representatives.

    for t1 in threes.iter() {
        't2_loop: for t2 in threes.iter() {
            if t1 == t2 {
                continue;
            }
            let (mut matches, mut mismatch) = (0, 0);
            for i in 0..3 {
                if t1[i] == t2[i] {
                    matches += 1;
                } else {
                    mismatch = i;
                }
            }
            if matches != 2 {
                continue;
            }
            let (u1, u2) = (threesp[t1], threesp[t2]);
            let (ex1, ex2) = (&exact_clonotypes[exacts[u1]], &exact_clonotypes[exacts[u2]]);
            for m1 in 0..ex1.share.len() {
                let p1 = bin_position(&chains, &(u1, m1));
                let q1 = bin_position(&r, &p1) as usize;
                if q1 == t1[mismatch] {
                    for m2 in 0..ex2.share.len() {
                        let p2 = bin_position(&chains, &(u2, m2));
                        let q2 = bin_position(&r, &p2) as usize;
                        if q2 == t2[mismatch] {
                            let (seq1, seq2) = (&ex1.share[m1].seq, &ex2.share[m2].seq);
                            if seq1.len() == seq2.len() {
                                const MAX_DIFFS: usize = 10;
                                let mut diffs = 0;
                                for j in 0..seq1.len() {
                                    if seq1[j] != seq2[j] {
                                        diffs += 1;
                                    }
                                }
                                if diffs <= MAX_DIFFS {
                                    e.join(p1, p2);
                                    break 't2_loop;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Get representatives for the chain orbits.

    let mut r = Vec::<i32>::new();
    e.orbit_reps(&mut r);

    // Reorder the chains.  This is done to get the heavy chains before the light chains and also
    // to mimic the behavior of the previous version of this algorithm, to minimiize churn.  Then
    // update the representatives.

    let mut chainsp = Vec::<(String, usize, usize, usize)>::new();
    for u in 0..nexacts {
        let ex = &exact_clonotypes[exacts[u]];
        for m in 0..ex.share.len() {
            let mut c = ex.share[m].chain_type.clone();
            if c.starts_with("TRB") {
                c = c.replacen("TRB", "TRX", 1);
            } else if c.starts_with("TRA") {
                c = c.replacen("TRA", "TRY", 1);
            }
            chainsp.push((
                format!("{}:{}", c, ex.share[m].cdr3_aa),
                ex.share[m].seq.len(),
                u,
                m,
            ));
        }
    }
    chainsp.sort();
    let mut chainso = Vec::<(usize, usize)>::new();
    let mut chainsox = Vec::<(usize, usize, usize)>::new();
    for i in 0..chainsp.len() {
        chainso.push((chainsp[i].2, chainsp[i].3));
        chainsox.push((chainsp[i].2, chainsp[i].3, i));
    }
    chainsox.sort_unstable();
    for i in 0..r.len() {
        r[i] = chainsox[r[i] as usize].2 as i32;
    }
    r.sort_unstable();

    // Create rmap, that sends
    // (index into exact subclonotypes for this clonotype,
    //  index into chains for one of these exact subclonotypes)
    // to an index into the orbit reps for the chains.

    let mut rpos = HashMap::<(usize, usize), usize>::new();
    for i in 0..chains.len() {
        let c = e.class_id(i as i32);
        let f = chainsox[c as usize].2 as i32;
        let q = bin_position(&r, &f) as usize;
        rpos.insert((chains[i].0, chains[i].1), q);
    }

    // Find the maximum multiplicity of each orbit, and the number of columns.

    let mut mm = vec![0; r.len()];
    for u in 0..nexacts {
        let ex = &exact_clonotypes[exacts[u]];
        let mut mm0 = vec![0; r.len()];
        for m in 0..ex.share.len() {
            mm0[rpos[&(u, m)]] += 1;
        }
        for i in 0..r.len() {
            mm[i] = max(mm[i], mm0[i]);
        }
    }
    let cols = mm.iter().sum();

    // Define a matrix mat[col][ex] which is the column of the exact subclonotype ex corresponding
    // to the given column col of the clonotype, which may or may not be defined.

    let mut mat = vec![vec![None; nexacts]; cols];
    for cx in 0..cols {
        // for every column
        'exact: for u in 0..nexacts {
            // for every exact subclonotype
            let ex = &exact_clonotypes[exacts[u]];
            let mut mm0 = vec![0; r.len()];
            for m in 0..ex.share.len() {
                // for every chain in the exact subclonotype:
                let q = rpos[&(u, m)];
                let mut col = mm0[q];
                for j in 0..q {
                    col += mm[j];
                }
                mm0[q] += 1;
                if col == cx {
                    mat[cx][u] = Some(m);
                    continue 'exact;
                }
            }
        }
    }
    mat
}
