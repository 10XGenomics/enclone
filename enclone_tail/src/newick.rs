// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Convert a rooted directed tree into Newick format,
// see https://en.wikipedia.org/wiki/Newick_format.
//
// This is not an efficient implementation.
//
// The input graph must be acyclic and connected.  Otherwise bad things may happen.  Not checked.
//
// Input data:
// 1. vertex names = vnames
// 2. index of root vertex = r
// 3. edges (v, w, edge-name) = edges, where v and w are zero-based indices of vertices.
//
// Upon entry, the edge-names should be string representations of weights.

use itertools::Itertools;
use std::cmp::max;

pub fn newick(vnames: &Vec<String>, r: usize, edges: &Vec<(usize, usize, String)>) -> String {
    // Set up.

    let mut edges = edges.clone();
    let mut n = 0;
    for i in 0..edges.len() {
        n = max(n, edges[i].0 + 1);
        n = max(n, edges[i].1 + 1);
    }
    assert!(r < n);
    let mut index = vec![Vec::<usize>::new(); n];
    for i in 0..edges.len() {
        index[edges[i].0].push(i);
        index[edges[i].1].push(i);
    }
    assert_eq!(n, vnames.len());

    // Incorporate the vertex names into the weights.

    for i in 0..edges.len() {
        edges[i].2 = format!("{}:{}", vnames[edges[i].1], edges[i].2);
    }

    // Gradually chew back the edges and roll up the labels as we go.

    let mut used = vec![false; n];
    let mut nused = 0;
    // while nused < n - 1 {
    loop {
        let nused0 = nused;
        for v in 0..n {
            if !used[v] && index[v].len() > 1 {
                let mut subterminal = true;
                for i in index[v].iter() {
                    let e = &edges[*i];
                    if e.0 == v {
                        let w = e.1;
                        if index[w].len() > 1 {
                            subterminal = false;
                            break;
                        }
                    }
                }
                if subterminal {
                    let mut labels = Vec::<String>::new();
                    let mut j = 0;
                    for i in index[v].iter() {
                        let e = &edges[*i];
                        if e.0 == v {
                            let w = e.1;
                            assert!(!used[w]);
                            used[w] = true;
                            nused += 1;
                            labels.push(e.2.clone());
                        } else {
                            j = *i;
                        }
                    }
                    index[v] = vec![j];
                    let label = labels.iter().format(",");
                    if v != r {
                        edges[j].2 = format!("({}){}", label, edges[j].2);
                    } else {
                        edges[j].2 = format!("({}){};", label, vnames[edges[j].0]);
                    }
                }
            }
        }
        if nused == nused0 {
            break;
        }
    }
    if nused == n - 1 {
        edges[index[r][0]].2.clone()
    } else {
        format!("({}){};", edges[index[r][0]].2, vnames[r])
    }
}
