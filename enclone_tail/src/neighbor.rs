// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Implement the phylogenetic tree neighbor joining algorithm, with one tweak (see below).
//
// Saitou N., Nei M. (1987). The neighbor-joining method: a new method for reconstructing
// phylogenetic trees.  Molecular Biology and Evolution 4: 406–425. PMID 3447015.
//
// We follow https://en.wikipedia.org/wiki/Neighbor_joining.
//
// Tweak: negative edge lengths are replaced by zero as suggested by Kuhner and Felsenstein (1994).
// Kuhner M.K., Felsenstein J. (1994). A simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates. Molecular Biology and Evolution 11(3): 459-468. PMID 8015439.
//
// The single input argument should be a symmetric n x n matrix, n >= 1.
// The output is a vector of 2n-3 edges, represented as (v, w, distance).
//
// Note that this algorithm is O(n^3).

pub fn neighbor_joining(d: &Vec<Vec<f64>>) -> Vec<(usize, usize, f64)> {
    let (mut d, mut d2) = (d.clone(), d.clone());
    let n0 = d.len();
    assert!(n0 >= 1);
    for i in 0..n0 {
        assert_eq!(d[i].len(), n0);
    }
    for i in 0..n0 {
        for j in i + 1..n0 {
            assert_eq!(d[i][j], d[j][i]);
        }
    }
    if n0 == 1 {
        return Vec::new();
    } else if n0 == 2 {
        return vec![(0, 1, d[0][1])];
    }
    let mut verts = vec![0; n0];
    for i in 0..n0 {
        verts[i as usize] = i;
    }
    let mut edges = vec![(0, 0, 0.0); 2 * n0 - 3];
    let mut q = vec![vec![0.0; n0]; n0];
    for n in (3..=n0).rev() {
        for i in 0..n {
            for j in i + 1..n {
                q[i][j] = (n - 2) as f64 * d[i][j];
                for k in 0..n {
                    q[i][j] -= d[i][k] + d[j][k];
                }
                q[j][i] = q[i][j];
            }
        }
        let (mut f, mut g) = (0, 1);
        let mut m = q[0][1];
        for i in 0..n {
            for j in i + 1..n {
                if q[i][j] < m {
                    f = i;
                    g = j;
                    m = q[i][j];
                }
            }
        }
        let mut df = (n - 2) as f64 * d[f][g];
        for k in 0..n {
            df += d[f][k] - d[g][k];
        }
        df /= (2 * (n - 2)) as f64;
        let dg = d[f][g] - df;
        let vnew = n0 + (n0 - n);
        edges[2 * (n0 - n)] = (verts[f], vnew, df);
        edges[2 * (n0 - n) + 1] = (verts[g], vnew, dg);
        verts[f] = vnew;
        for k in g..n - 1 as usize {
            verts[k] = verts[k + 1];
        }
        for i in 0..n {
            for j in 0..n as usize {
                d2[i][j] = d[i][j];
            }
        }
        for k in 0..n {
            if k != f {
                d2[f][k] = (d[f][k] + d[g][k] - d[f][g]) / 2.0;
                d2[k][f] = d2[f][k];
            }
        }
        for i in 0..n {
            if i != g {
                for j in 0..n {
                    if j != g {
                        let (mut ip, mut jp) = (i, j);
                        if ip > g {
                            ip -= 1;
                        }
                        if jp > g {
                            jp -= 1;
                        }
                        d[ip][jp] = d2[i][j];
                    }
                }
            }
        }
        if n == 3 {
            edges[2 * n0 - 4] = (verts[0], verts[1], d[0][1]);
        }
    }
    for i in 0..edges.len() {
        edges[i].2 = edges[i].2.max(0.0);
    }
    edges
}
