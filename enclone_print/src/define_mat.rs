// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::*;
use equiv::EquivRel;
use std::cmp::max;
use std::collections::HashMap;
use vector_utils::*;

pub fn define_mat(
    _ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    exacts: &Vec<usize>,
    _cdr3s: &Vec<Vec<(String, usize)>>,
    _js: &Vec<usize>,
    _ks: &Vec<usize>,
    od: &Vec<(Vec<usize>, usize, i32)>,
    info: &Vec<CloneInfo>,
    raw_joins: &Vec<Vec<usize>>,
) -> Vec<Vec<Option<usize>>> {

    // Define map to indices into exacts.

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
    infos.sort();

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

    // Define an equivalence relation on the chains, introducing connections defined by the
    // raw joins.  Also join where there are identical V..J sequences.

    let mut e = EquivRel::new(chains.len() as i32);
    for i1 in 0..infos.len() {
        let j1 = infos[i1];
        let u1 = info[j1].clonotype_index;
        let v1 = to_exacts[&u1];
        let m1s = &info[j1].exact_cols;
        for x in raw_joins[j1].iter() {
            let i2 = bin_position(&infos, &*x);
            if i2 >= 0 {
                let i2 = i2 as usize;
                let j2 = infos[i2];
                let u2 = info[j2].clonotype_index;
                let v2 = to_exacts[&u2];
                let m2s = &info[j2].exact_cols;
                for j in 0..2 {
                    let z1 = bin_position(&chains, &(v1, m1s[j]));
                    let z2 = bin_position(&chains, &(v2, m2s[j]));
                    e.join(z1, z2);
                }
            }
        }
    }
    let mut i = 0;
    while i < seq_chains.len() {
        let j = next_diff1_3(&seq_chains, i as i32) as usize;
        for k in i+1..j {
            let (x1, x2) = (&seq_chains[i], &seq_chains[k]);
            let z1 = bin_position(&chains, &(x1.1, x1.2));
            let z2 = bin_position(&chains, &(x2.1, x2.2));
            e.join(z1, z2);
        }
        i = j;
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
            chainsp.push((format!("{}:{}", c, ex.share[m].cdr3_aa), ex.share[m].seq.len(), u, m));
        }
    }
    chainsp.sort();
    let mut chainso = Vec::<(usize, usize)>::new();
    let mut chainsox = Vec::<(usize, usize, usize)>::new();
    for i in 0..chainsp.len() {
        chainso.push((chainsp[i].2, chainsp[i].3));
        chainsox.push((chainsp[i].2, chainsp[i].3, i));
    }
    chainsox.sort();
    for i in 0..r.len() {
        r[i] = chainsox[r[i] as usize].2 as i32;
    }
    r.sort();

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
