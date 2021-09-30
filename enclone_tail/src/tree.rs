// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Generate experimental tree output (options NEWICK0 and TREE).

use crate::display_tree::display_tree;
use crate::neighbor::neighbor_joining;
use crate::newick::newick;
use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use enclone_proto::types::DonorReferenceItem;
use io_utils::{fwrite, fwriteln};
use std::cmp::max;
use std::collections::HashMap;
use std::io::Write;
use std::mem::swap;
use vdj_ann::refx::RefData;

pub fn print_tree(
    oo: usize,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    out_datas: &Vec<Vec<HashMap<String, String>>>,
    logx: &mut Vec<u8>,
) {
    if ctl.gen_opt.newick || ctl.gen_opt.tree_on {
        // Compute the n x n distance matrix for the exact subclonotypes.

        let n = exacts[oo].len();
        let cols = rsi[oo].mat.len();
        let mut dist = vec![vec![0; n]; n];
        for i1 in 0..n {
            for i2 in 0..n {
                let ex1 = &exact_clonotypes[exacts[oo][i1]];
                let ex2 = &exact_clonotypes[exacts[oo][i2]];
                for m in 0..cols {
                    if rsi[oo].mat[m][i1].is_some() && rsi[oo].mat[m][i2].is_some() {
                        let r1 = rsi[oo].mat[m][i1].unwrap();
                        let r2 = rsi[oo].mat[m][i2].unwrap();
                        let seq1 = &ex1.share[r1].seq_del_amino;
                        let seq2 = &ex2.share[r2].seq_del_amino;
                        for j in 0..seq1.len() {
                            if seq1[j] != seq2[j] {
                                dist[i1][i2] += 1;
                            }
                        }
                    }
                }
            }
        }

        // Add a zeroeth entry for a "root subclonotype" which is defined to have the
        // donor reference away from the recombination region, and is undefined within it.
        // Define its distance to actual exact subclonotypes by only computing away from
        // the recombination region.  This yields an (n+1) x (n+1) matrix.

        let mut droot = vec![0; n];
        for i in 0..n {
            let ex = &exact_clonotypes[exacts[oo][i]];
            for m in 0..cols {
                if rsi[oo].mat[m][i].is_some() {
                    let r = rsi[oo].mat[m][i].unwrap();
                    let seq = &ex.share[r].seq_del_amino;
                    let mut vref = refdata.refs[rsi[oo].vids[m]].to_ascii_vec();
                    if rsi[oo].vpids[m].is_some() {
                        vref = dref[rsi[oo].vpids[m].unwrap()].nt_sequence.clone();
                    }
                    let jref = refdata.refs[rsi[oo].jids[m]].to_ascii_vec();
                    let z = seq.len();
                    for p in 0..z {
                        let b = seq[p];
                        if p < vref.len() - ctl.heur.ref_v_trim && b != vref[p] {
                            droot[i] += 1;
                        }
                        if p >= z - (jref.len() - ctl.heur.ref_j_trim)
                            && b != jref[jref.len() - (z - p)]
                        {
                            droot[i] += 1;
                        }
                    }
                }
            }
        }
        let mut distp = vec![vec![0.0; n + 1]; n + 1];
        for i1 in 0..n {
            for i2 in 0..n {
                distp[i1 + 1][i2 + 1] = dist[i1][i2] as f64;
            }
        }
        for i in 0..n {
            distp[i + 1][0] = droot[i] as f64;
            distp[0][i + 1] = droot[i] as f64;
        }

        // Generate the neighborhood joining tree associated to these data.

        let mut tree = neighbor_joining(&distp);
        let mut nvert = 0;
        for i in 0..tree.len() {
            nvert = max(nvert, tree[i].0 + 1);
            nvert = max(nvert, tree[i].1 + 1);
        }

        // Use the root to direct the edges.

        let r = 0;
        let mut index = vec![Vec::<usize>::new(); nvert];
        for i in 0..tree.len() {
            index[tree[i].0].push(i);
            index[tree[i].1].push(i);
        }
        let mut rooted = vec![false; nvert];
        rooted[r] = true;
        let mut roots = vec![r];
        for i in 0..nvert {
            let v = roots[i];
            for j in index[v].iter() {
                let e = &mut tree[*j];

                if e.1 == v && !rooted[e.0] {
                    swap(&mut e.0, &mut e.1);
                }
                if e.0 == v && !rooted[e.1] {
                    rooted[e.1] = true;
                    roots.push(e.1);
                }
            }
        }

        // Output in Newick format.

        if ctl.gen_opt.newick {
            let mut vnames = Vec::<String>::new();
            for i in 0..=n {
                vnames.push(format!("{}", i));
            }
            let mut edges = Vec::<(usize, usize, String)>::new();
            for i in 0..tree.len() {
                edges.push((tree[i].0, tree[i].1, format!("{:.2}", tree[i].2)));
            }
            for i in n + 1..nvert {
                vnames.push(format!("I{}", i - n));
            }
            let nw = newick(&vnames, 0, &edges);
            fwriteln!(logx, "\n{}", nw);
        }

        // Output as visual tree.

        if ctl.gen_opt.tree_on {
            let mut edges = Vec::<(usize, usize, f64)>::new();
            let mut nvert = 0;
            for i in 0..tree.len() {
                edges.push((tree[i].0, tree[i].1, tree[i].2));
                nvert = max(nvert, tree[i].0 + 1);
                nvert = max(nvert, tree[i].1 + 1);
            }

            // Make edge names.

            let mut vnames = Vec::<String>::new();
            for i in 0..nvert {
                let mut len = 0.0;
                for j in 0..edges.len() {
                    if edges[j].1 == i {
                        len = edges[j].2;
                    }
                }
                let mut c = String::new();
                if i > 0 && i <= n && !ctl.gen_opt.tree.is_empty() {
                    let x = &out_datas[oo][i - 1];
                    for w in ctl.gen_opt.tree.iter() {
                        if x.contains_key(&*w) {
                            c += &format!(",{}={}", w, x[&*w]);
                        }
                    }
                }
                if i == 0 {
                    vnames.push("•".to_string());
                } else if i <= n {
                    if ctl.pretty {
                        vnames.push(format!("[01m[31m{}[0m [{:.2}{}]", i, len, c));
                    } else {
                        vnames.push(format!("{} [{:.2}{}]", i, len, c));
                    }
                } else {
                    vnames.push(format!("• [{:.2}{}]", len, c));
                }
            }

            // Display the tree.

            let nw = display_tree(&vnames, &edges, 0, 100);
            fwrite!(logx, "\n{}", nw);
        }
    }
}
