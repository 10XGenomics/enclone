// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Text display of a rooted directed tree, in which each vertex has a name and each edge has a
// floating point length.  A width parameter sets the approximate page width in characters,
// with edge lengths scaled roughly to match this.
//
// The design is adapted from https://gitlab.com/Noughmad/ptree, by Miha Čančula.

use itertools::Itertools;
use std::cmp::max;
use vector_utils::sort_sync2;

// vnames: vertex names
// directed edges: {(v, w, weight)}
// r: index of the root vertex
// max_width: max page width

pub fn display_tree(
    vnames: &Vec<String>,
    edges: &Vec<(usize, usize, f64)>,
    r: usize,
    width: usize,
) -> String {
    // Test input data and create an index.

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

    // For each vertex, define a path, which is the sequence of edge indices from the root to it.

    let mut vpaths = vec![Vec::<usize>::new(); n];
    for v in 0..n {
        let mut w = v;
        while w != r {
            for j in index[w].iter() {
                if edges[*j].1 == w {
                    vpaths[v].push(*j);
                    w = edges[*j].0;
                }
            }
        }
        vpaths[v].reverse();
    }
    let mut vs = Vec::<usize>::new();
    for v in 0..n {
        vs.push(v);
    }

    // Sort.  The output lines now correspond to the entries of this vector.

    sort_sync2(&mut vpaths, &mut vs);

    // In the special case where every edge has length zero, change all to one.

    let mut max_e = 0 as f64;
    for e in edges.iter() {
        max_e = max_e.max(e.2);
    }
    if max_e == 0.0 {
        for i in 0..edges.len() {
            edges[i].2 = 1.0;
        }
    }

    // For each path, define its constant and variable length components.

    let mut clen = vec![0; n];
    let mut vlen = vec![Vec::<f64>::new(); n];
    for i in 0..vpaths.len() {
        clen[i] = vnames[i].chars().count() + 2 * vpaths[i].len();
        for j in vpaths[i].iter() {
            vlen[i].push(edges[*j].2.max(0.0));
        }
    }

    // Define length multiplier.  This is very inefficient.  Then scale the lengths.

    let mut mult = 1.0;
    {
        let fwidth = width as f64;
        let mut last_change = "".to_string();
        loop {
            let mut len = 0.0_f64;
            let mut max_w = 1.0_f64;
            for i in 0..vpaths.len() {
                let mut l = clen[i] as f64;
                for j in 0..vlen[i].len() {
                    let w = (vlen[i][j] * mult).round().max(1.0);
                    max_w = max_w.max(w);
                    l += w;
                }
                len = len.max(l);
            }
            if len <= fwidth && len >= 1.05 * fwidth {
                break;
            } else if max_w == 1.0 && len >= fwidth {
                break;
            } else if len > fwidth {
                mult *= 0.95;
                last_change = "minus".to_string();
            } else {
                if last_change == *"minus" {
                    break;
                }
                mult *= 1.05;
                last_change = "plus".to_string();
            }
        }
    }
    for i in 0..edges.len() {
        edges[i].2 = (edges[i].2 * mult).round().max(1.0);
    }

    // Generate the lines.

    let mut x = String::new();
    for i in 0..n {
        for j in 0..vpaths[i].len() {
            let e = vpaths[i][j];
            let t = edges[e].2 as usize;
            let hedge = format!("{}", vec!['═'; t].iter().format(""));
            let hnone = format!("{}", vec![' '; t].iter().format(""));
            let mut last_edge = true;
            for k in i + 1..n {
                if j >= vpaths[k].len() {
                    break;
                }
                if vpaths[k][j] != vpaths[i][j] {
                    last_edge = false;
                }
            }

            // alternatives (all variable length)

            // 1. "╠═══ "

            // 2. "║    "

            // 3. "╠═══ "

            // 4. "╚═══ "

            // 5. "     "

            if j >= vpaths[i - 1].len() {
                if !last_edge {
                    x += &format!("╠{} ", hedge);
                } else {
                    x += &format!("╚{} ", hedge);
                }
            } else if vpaths[i][j] != vpaths[i - 1][j] {
                if !last_edge {
                    x += &format!("╠{} ", hedge);
                } else {
                    x += &format!("╚{} ", hedge);
                }
            } else if !last_edge {
                x += &format!("║{} ", hnone);
            } else {
                x += &format!(" {} ", hnone);
            }
        }
        x += &format!("{}\n", vnames[vs[i]]);
    }
    x
}
