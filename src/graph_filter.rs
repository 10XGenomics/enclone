// Copyright (c) 2019 10X Genomics, Inc. All rights reserved.

// This file provides the single function graph_filter.

use crate::defs::*;
use graph_simple::GraphSimple;
use io_utils::*;
use petgraph::prelude::*;
use std::cmp::*;
use std::io::Write;
use string_utils::*;
use vector_utils::*;

// Create a digraph which has one vertex for each V..J that appears in a productive
// pair, and for a given light chain and a given heavy chain vertex, a weighted edge
// from the light to the heavy, so long as they co-occur in some cell; the weight is
// a pair (numi, ncells), where numi is the sum over all cells for which the co-occur,
// of the minimum of the number of UMIs that support the two chains.
// seqs = { (V..J, is_igh) }
// Also tracking CDR3_AAs for now for exploratory purposes at least.  This makes
// things less efficient.
//
// As part of this, identify weak links and kill them.
//
// This code tiptoes around the fact that whereas these calculations should be at the clonotype
// level, we compute here at the exact subclonotype level.  The use of V segments is part of this.
//
// Hmm, seems like the edges go from heavy to light.

pub fn graph_filter(mut tig_bc: &mut Vec<Vec<TigData>>, graph: bool) {
    let mut seqs = Vec::<(Vec<u8>, bool, String, usize)>::new();
    for i in 0..tig_bc.len() {
        for j in 0..tig_bc[i].len() {
            let x = &tig_bc[i][j];
            seqs.push((x.seq.clone(), x.left, x.cdr3_aa.clone(), x.v_ref_id));
        }
    }
    unique_sort(&mut seqs);
    let mut edges0 = Vec::<(usize, usize, usize)>::new();
    for i in 0..tig_bc.len() {
        for j1 in 0..tig_bc[i].len() {
            if tig_bc[i][j1].left {
                let x1 = &tig_bc[i][j1];
                let p1 = bin_position(
                    &seqs,
                    &(x1.seq.clone(), true, x1.cdr3_aa.clone(), x1.v_ref_id),
                ) as usize;
                for j2 in 0..tig_bc[i].len() {
                    if !tig_bc[i][j2].left {
                        let x2 = &tig_bc[i][j2];
                        let p2 = bin_position(
                            &seqs,
                            &(x2.seq.clone(), false, x2.cdr3_aa.clone(), x2.v_ref_id),
                        ) as usize;
                        edges0.push((p1, p2, min(x1.umi_count, x2.umi_count)));
                    }
                }
            }
        }
    }
    edges0.sort();
    let mut edges1 = Vec::<(usize, usize, (usize, usize))>::new();
    let mut i = 0;
    while i < edges0.len() {
        let j = next_diff12_3(&edges0, i as i32) as usize;
        let mut weight = 0;
        for k in i..j {
            weight += edges0[k].2;
        }
        edges1.push((edges0[i].0, edges0[i].1, (weight, j - i)));
        i = j;
    }
    let mut g = Graph::<u32, (usize, usize), Directed>::new();
    g.reserve_exact_nodes(seqs.len());
    g.reserve_exact_edges(edges1.len());
    for i in 0..seqs.len() {
        g.add_node(i as u32);
    }
    for e in 0..edges1.len() {
        let v = edges1[e].0;
        let w = edges1[e].1;
        let weight = edges1[e].2;
        g.add_edge(NodeIndex::<u32>::new(v), NodeIndex::<u32>::new(w), weight);
    }

    // Kill weak branches from light to heavy chains.  Also kill light chain onesies that
    // have too many heavy chain partners.

    let mut log = Vec::<u8>::new();
    fwriteln!(log, "\nBRANCHING FROM LIGHT CHAINS");
    const MIN_RATIO_KILL: usize = 8;
    const MAX_KILL: usize = 5;
    const MAX_KILL_CELLS: usize = 2;
    const MAX_PARTNERS: usize = 50;
    let mut kills = Vec::<(usize, usize)>::new();
    let mut badones = Vec::<usize>::new();
    for v in 0..g.node_count() {
        if g.n_to(v) > 1 {
            let mut stats = Vec::<(usize, usize, usize)>::new();
            fwriteln!(log, "\nlight chain {} = {}", v, seqs[v].2);
            for i in 0..g.n_to(v) {
                let (w, e) = (g.v_to(v, i), g.e_to(v, i));
                let numi = g.edge_obj(e as u32).0;
                let ncells = g.edge_obj(e as u32).1;
                fwriteln!(
                    log,
                    "• heavy chain {} = {}, weight = {}/{}",
                    w,
                    seqs[w].2,
                    numi,
                    ncells
                );
                stats.push((numi, ncells, w));
            }
            reverse_sort(&mut stats);
            for i in 1..stats.len() {
                // Below, the part seqs[stats[i].2].3 != seqs[stats[0].2].3
                // is requiring that the partner V segments are different.  See discussion at
                // the top around the roles of clonotypes versus exact subclonotypes.

                if seqs[stats[i].2].3 != seqs[stats[0].2].3 {
                    let numi = stats[i].0;
                    let ncells = stats[i].1;
                    let numi_best = stats[0].0;
                    let ncells_best = stats[0].1;
                    if numi_best >= MIN_RATIO_KILL * max(1, numi) && numi <= MAX_KILL {
                        kills.push((v, stats[i].2));
                    } else if numi_best >= numi && ncells_best >= MIN_RATIO_KILL * max(1, ncells) {
                        if ncells <= MAX_KILL_CELLS {
                            let w = stats[i].2;
                            if graph {
                                println!(
                                    "\nkill type 2, from {} to {}, ncells = {}, numi = {}",
                                    seqs[v].2, seqs[w].2, ncells, numi
                                );
                                println!(
                                    "killed by {} to {}, ncells = {}, numi = {}",
                                    seqs[v].2, seqs[stats[0].2].2, ncells_best, numi_best
                                );
                            }
                            kills.push((v, w));
                        } else {
                            let w = stats[i].2;
                            for j in 0..g.n_from(w) {
                                let (vx, ex) = (g.v_from(w, j), g.e_from(w, j));
                                let numix = g.edge_obj(ex as u32).0;
                                let ncellsx = g.edge_obj(ex as u32).1;
                                if vx != v
                                    && ncellsx >= MIN_RATIO_KILL * ncells
                                    && numix >= MIN_RATIO_KILL * numi
                                {
                                    kills.push((v, w));
                                }
                            }
                        }
                    }
                }
            }
        }
        if g.n_to(v) > MAX_PARTNERS {
            badones.push(v);
        }
    }
    kills.sort();
    // presumably badly inefficient
    let mut to_delete = vec![false; tig_bc.len()];
    for i in 0..tig_bc.len() {
        for j1 in 0..tig_bc[i].len() {
            if tig_bc[i][j1].left {
                continue;
            }
            let x1 = &tig_bc[i][j1];
            let m1 = (x1.seq.clone(), x1.left, x1.cdr3_aa.clone(), x1.v_ref_id);
            let p1 = bin_position(&seqs, &m1) as usize;
            for j2 in 0..tig_bc[i].len() {
                if !tig_bc[i][j2].left {
                    continue;
                }
                let x2 = &tig_bc[i][j2];
                let m2 = (x2.seq.clone(), x2.left, x2.cdr3_aa.clone(), x2.v_ref_id);
                let p2 = bin_position(&seqs, &m2) as usize;
                if bin_member(&kills, &(p1, p2)) {
                    to_delete[i] = true;
                }
            }
        }
        if tig_bc[i].len() == 1 {
            let x0 = &tig_bc[i][0];
            let m0 = (x0.seq.clone(), x0.left, x0.cdr3_aa.clone(), x0.v_ref_id);
            let p = bin_position(&seqs, &m0) as usize;
            if bin_member(&badones, &p) {
                to_delete[i] = true;
            }
        }
    }
    erase_if(&mut tig_bc, &to_delete);
    if graph {
        fwriteln!(log, "");
        print!("{}", strme(&log));
    }

    // Kill weak branches from heavy to light chains.

    let mut log = Vec::<u8>::new();
    fwriteln!(log, "BRANCHING FROM HEAVY CHAINS");
    const MIN_RATIO_KILL_HEAVY: usize = 8;
    const MAX_KILL_HEAVY: usize = 6;
    const MAX_KILL_HEAVY_CELLS: usize = 1;
    let mut kills = Vec::<(usize, usize)>::new();
    for v in 0..g.node_count() {
        if g.n_from(v) > 1 {
            let mut stats = Vec::<((usize, usize), usize)>::new();
            fwriteln!(log, "\nheavy chain {} = {}", v, seqs[v].2);
            for i in 0..g.n_from(v) {
                let w = g.v_from(v, i);
                let e = g.e_from(v, i);
                let weight = g.edge_obj(e as u32);
                fwriteln!(
                    log,
                    "• light chain {} = {}, weight = {}/{}",
                    w,
                    seqs[w].2,
                    weight.1,
                    weight.0
                );
                stats.push((*weight, w));
            }
            reverse_sort(&mut stats);
            for i in 1..stats.len() {
                if (stats[0].0).0 >= MIN_RATIO_KILL_HEAVY * max(1, (stats[i].0).0)
                    && (stats[i].0).0 <= MAX_KILL_HEAVY
                    && (stats[i].0).1 <= MAX_KILL_HEAVY_CELLS
                {
                    kills.push((v, stats[i].1));
                }
            }
        }
    }
    kills.sort();
    // presumably badly inefficient
    let mut to_delete = vec![false; tig_bc.len()];
    for i in 0..tig_bc.len() {
        for j1 in 0..tig_bc[i].len() {
            if !tig_bc[i][j1].left {
                continue;
            }
            let x1 = &tig_bc[i][j1];
            let m1 = (x1.seq.clone(), x1.left, x1.cdr3_aa.clone(), x1.v_ref_id);
            let p1 = bin_position(&seqs, &m1) as usize;
            for j2 in 0..tig_bc[i].len() {
                if tig_bc[i][j2].left {
                    continue;
                }
                let x2 = &tig_bc[i][j2];
                let m2 = (x2.seq.clone(), x2.left, x2.cdr3_aa.clone(), x2.v_ref_id);
                let p2 = bin_position(&seqs, &m2) as usize;
                if bin_member(&kills, &(p1, p2)) {
                    to_delete[i] = true;
                }
            }
        }
    }
    erase_if(&mut tig_bc, &to_delete);
    if graph {
        fwriteln!(log, "");
        print!("{}", strme(&log));
    }
}
