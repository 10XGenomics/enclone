// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Group and print clonotypes.  For now, limited grouping functionality.
//
// To keep compilation time down, this crate should not reach into the enclone crate.

use crate::clustal::*;
use crate::display_tree::*;
use crate::fasta::*;
use crate::grouper::*;
use crate::neighbor::*;
use crate::newick::*;
use crate::phylip::*;
use crate::plot_points::*;
use crate::print_stats::*;
use crate::requirements::*;
use ansi_escape::ansi_to_html::*;
use ansi_escape::*;
use enclone_core::defs::*;
use enclone_core::mammalian_fixed_len::*;
use enclone_core::print_tools::*;
use enclone_proto::types::*;
use io_utils::*;
use itertools::*;
use std::cmp::max;
use std::collections::HashMap;
use std::env;
use std::fs::remove_file;
use std::fs::File;
use std::io::Write;
use std::io::*;
use std::mem::swap;
use std::path::Path;
use std::time::Instant;
use string_utils::*;
use tables::*;
use tar::Builder;
use vdj_ann::refx::*;
use vector_utils::*;

pub fn group_and_print_clonotypes(
    tall: &Instant,
    refdata: &RefData,
    pics: &Vec<String>,
    exacts: &Vec<Vec<usize>>,
    in_center: &Vec<bool>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    join_info: &Vec<(usize, usize, bool, Vec<u8>)>,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    fate: &Vec<HashMap<String, String>>,
    dref: &Vec<DonorReferenceItem>,
) {
    // Build index to join info.

    let mut to_join_info = vec![Vec::<usize>::new(); exact_clonotypes.len()];
    for i in 0..join_info.len() {
        to_join_info[join_info[i].0].push(i);
        to_join_info[join_info[i].1].push(i);
    }

    // Set up for parseable output.

    let mut parseable_fields = Vec::<String>::new();
    set_speakers(&ctl, &mut parseable_fields);
    #[allow(bare_trait_objects)]
    let mut pout = match ctl.parseable_opt.pout.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        "stdouth" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.parseable_opt.pout);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
    };
    let mut pcols = ctl.parseable_opt.pcols.clone();
    for i in 0..pcols.len() {
        pcols[i] = pcols[i].replace("_Σ", "_sum");
        pcols[i] = pcols[i].replace("_μ", "_mean");
    }
    if pcols.is_empty() {
        pcols = parseable_fields.clone();
    }
    let mut pcols2 = Vec::<String>::new();
    for i in 0..pcols.len() {
        if pcols[i].contains(":") {
            pcols2.push(pcols[i].before(":").to_string());
        } else {
            pcols2.push(pcols[i].clone());
        }
    }
    pcols = pcols2;
    if !ctl.parseable_opt.pout.is_empty()
        && ctl.parseable_opt.pout != "stdout".to_string()
        && ctl.parseable_opt.pout != "stdouth".to_string()
    {
        fwriteln!(pout, "{}", pcols.iter().format(","));
    }

    // Set up for fasta output.

    #[allow(bare_trait_objects)]
    let mut fout = match ctl.gen_opt.fasta_filename.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.gen_opt.fasta_filename);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
    };
    #[allow(bare_trait_objects)]
    let mut faaout = match ctl.gen_opt.fasta_aa_filename.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.gen_opt.fasta_aa_filename);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
    };

    // Set up for clustal output.

    let (mut clustal_aa, mut clustal_dna) = (None, None);
    if ctl.gen_opt.clustal_aa.len() > 0 && ctl.gen_opt.clustal_aa != "stdout".to_string() {
        let file = File::create(&ctl.gen_opt.clustal_aa).unwrap();
        clustal_aa = Some(Builder::new(file));
    }
    if ctl.gen_opt.clustal_dna.len() > 0 && ctl.gen_opt.clustal_dna != "stdout".to_string() {
        let file = File::create(&ctl.gen_opt.clustal_dna).unwrap();
        clustal_dna = Some(Builder::new(file));
    }

    // Set up for phylip output.

    let (mut phylip_aa, mut phylip_dna) = (None, None);
    if ctl.gen_opt.phylip_aa.len() > 0 && ctl.gen_opt.phylip_aa != "stdout".to_string() {
        let file = File::create(&ctl.gen_opt.phylip_aa).unwrap();
        phylip_aa = Some(Builder::new(file));
    }
    if ctl.gen_opt.phylip_dna.len() > 0 && ctl.gen_opt.phylip_dna != "stdout".to_string() {
        let file = File::create(&ctl.gen_opt.phylip_dna).unwrap();
        phylip_dna = Some(Builder::new(file));
    }

    // Set up for peer group output.

    #[allow(bare_trait_objects)]
    let mut pgout = match ctl.gen_opt.peer_group_filename.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.gen_opt.peer_group_filename);
            Box::new(File::create(&path).unwrap()) as Box<Write>
        }
    };
    if ctl.gen_opt.peer_group_filename.len() > 0 {
        if ctl.gen_opt.peer_group_filename != "stdout".to_string() {
            if !ctl.gen_opt.peer_group_readable {
                fwriteln!(pgout, "group,clonotype,chain,pos,amino_acid,count");
            } else {
                fwriteln!(pgout, "group,clonotype,chain,pos,distribution");
            }
        }
    }

    // Group clonotypes.

    let groups = grouper(&refdata, &exacts, &in_center, &exact_clonotypes, &ctl);

    // Echo command.

    let mut last_width = 0;
    let mut logx = Vec::<u8>::new();
    if ctl.gen_opt.echo {
        let args: Vec<String> = env::args().collect();
        fwriteln!(logx, "\n{}", args.iter().format(" "));
        if ctl.gen_opt.html {
            fwriteln!(logx, "");
        }
    }

    // Now print clonotypes.

    let mut plot_xy_vals = Vec::<(f32, f32)>::new();
    for i in 0..groups.len() {
        let mut o = Vec::<i32>::new();
        for j in 0..groups[i].len() {
            o.push(groups[i][j].0);
        }
        let mut n = 0;
        for j in 0..o.len() {
            let x = o[j] as usize;
            let s = &exacts[x];
            for k in 0..s.len() {
                n += exact_clonotypes[s[k]].clones.len();
            }
        }

        // Generate human readable output.  Getting the newlines right is tricky, so
        // they're marked.

        if !ctl.gen_opt.noprint {
            if !ctl.gen_opt.html && !ctl.gen_opt.ngroup {
                fwriteln!(logx, ""); // NEWLINE 1
            }

            // If we just printed a clonotype box, output a bar.

            if last_width > 0 {
                if ctl.gen_opt.ngroup || ctl.gen_opt.html {
                    fwriteln!(logx, ""); // NEWLINE 2
                }
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_eight_bit_color_escape(&mut log, 44);
                    fwrite!(logx, "{}", strme(&log));
                }
                fwrite!(logx, "╺{}╸", "━".repeat(last_width - 2));
                if !ctl.gen_opt.ngroup {
                    fwriteln!(logx, ""); // NEWLINE 3
                }
                fwriteln!(logx, ""); // NEWLINE 4
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_end_escape(&mut log);
                    fwrite!(logx, "{}", strme(&log));
                }
            }

            // If NGROUP is not on, output a GROUP line, including a newline at the end.

            if !ctl.gen_opt.ngroup {
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_bold_escape(&mut log);
                    emit_eight_bit_color_escape(&mut log, 27);
                    fwrite!(logx, "{}", strme(&log));
                }
                fwrite!(
                    logx,
                    "[{}] GROUP = {} CLONOTYPES = {} CELLS",
                    i + 1,
                    o.len(),
                    n
                );
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_end_escape(&mut log);
                    fwrite!(logx, "{}", strme(&log));
                }
                fwriteln!(logx, ""); // NEWLINE 5
            }
        }
        let mut group_ncells = 0;
        for j in 0..o.len() {
            let oo = o[j] as usize;
            for l in 0..exacts[oo].len() {
                group_ncells += exact_clonotypes[exacts[oo][l]].ncells();
            }
        }
        for j in 0..o.len() {
            let oo = o[j] as usize;

            // Generate data for PLOT_XY.

            if ctl.plot_opt.plot_xy_filename.len() > 0 {
                let xvar = &ctl.plot_opt.plot_xy_xvar;
                let yvar = &ctl.plot_opt.plot_xy_yvar;
                for i in 0..out_datas[oo].len() {
                    let p = &out_datas[oo][i];
                    if p.contains_key(&xvar.clone()) {
                        let x = &p[&xvar.clone()];
                        if x.parse::<f64>().is_ok() {
                            let mut x = x.force_f64();
                            if ctl.plot_opt.plot_xy_x_log10 {
                                if x <= 0.0 {
                                    continue;
                                }
                                x = x.log10();
                            }
                            if p.contains_key(&yvar.clone()) {
                                let y = &p[&yvar.clone()];
                                if y.parse::<f64>().is_ok() {
                                    let mut y = y.force_f64();
                                    if ctl.plot_opt.plot_xy_y_log10 {
                                        if y <= 0.0 {
                                            continue;
                                        }
                                        y = y.log10();
                                    }
                                    plot_xy_vals.push((x as f32, y as f32));
                                }
                            }
                        }
                    }
                }
            }

            // Proceed.

            if !ctl.gen_opt.noprint {
                if i > 0 || j > 0 || !(ctl.gen_opt.html && ctl.gen_opt.ngroup) {
                    fwrite!(logx, "\n"); // NEWLINE 6
                }
                let mut s = format!("[{}.{}] {}", i + 1, j + 1, pics[oo]);
                if groups[i][j].1.len() > 0 {
                    s = format!("{}\n{}\n{}", s.before("\n"), groups[i][j].1, s.after("\n"));
                }
                if ctl.gen_opt.svg {
                    const FONT_SIZE: usize = 15;

                    // Generate svg.  This does not generate the shortest possible string.  One
                    // thing that could be done is to use only one text tag and instead use
                    // relative positions in the tspan tags to avoid repeating the font family,
                    // etc.  But there are probably other economizations.
                    //
                    // The other thing is that the aspect ratio is just a little bit off.

                    fwrite!(
                        logx,
                        "{}",
                        convert_text_with_ansi_escapes_to_svg(&s, "Menlo", FONT_SIZE)
                    );
                } else {
                    fwrite!(logx, "{}", s);
                }
            }
            let x = &pics[oo];
            let mut y = Vec::<char>::new();
            for c in x.chars() {
                y.push(c);
            }
            y.reverse();
            let mut m = 2;
            while m < y.len() {
                if y[m] == '\n' {
                    break;
                }
                m += 1;
            }
            last_width = m - 1;

            // Print join info.

            let mut ji = Vec::<usize>::new();
            for u in exacts[oo].iter() {
                ji.append(&mut to_join_info[*u].clone());
            }
            unique_sort(&mut ji);
            for i in 0..ji.len() {
                fwriteln!(logx, "{}", strme(&join_info[ji[i]].3));
            }

            // Generate clustal and phylip output.

            print_clustal(
                i,
                j,
                oo,
                &exacts,
                &rsi,
                &exact_clonotypes,
                &ctl,
                &mut logx,
                &mut clustal_aa,
                &mut clustal_dna,
            );
            print_phylip(
                i,
                j,
                oo,
                &exacts,
                &rsi,
                &exact_clonotypes,
                &ctl,
                &mut logx,
                &mut phylip_aa,
                &mut phylip_dna,
            );

            // Generate experimental tree output (options NEWICK0 and TREE).

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
                        if i > 0 && i <= n && ctl.gen_opt.tree.len() > 0 {
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

            // Generate peer group output.

            if ctl.gen_opt.peer_group_filename.len() > 0 {
                let pg = mammalian_fixed_len_peer_groups(&refdata);
                if !ctl.gen_opt.peer_group_readable {
                    if ctl.gen_opt.peer_group_filename == "stdout".to_string() {
                        fwriteln!(logx, "group,clonotype,chain,pos,amino_acid,count");
                    }
                } else {
                    if ctl.gen_opt.peer_group_filename == "stdout".to_string() {
                        fwriteln!(logx, "group,clonotype,chain,pos,distribution");
                    }
                }
                let chain_types = ["IGH", "IGK", "IGL", "TRA", "TRB"];
                for q in 0..rsi[oo].mat.len() {
                    let id = rsi[oo].vids[q];
                    if ctl.gen_opt.peer_group_filename != "stdout".to_string() {
                        if !ctl.gen_opt.peer_group_readable {
                            for y in pg[id].iter() {
                                fwriteln!(
                                    pgout,
                                    "{},{},{},{},{},{},{}",
                                    i + 1,
                                    j + 1,
                                    q + 1,
                                    chain_types[refdata.rtype[id] as usize],
                                    y.0,
                                    y.1 as char,
                                    y.2
                                );
                            }
                        } else {
                            let mut k = 0;
                            while k < pg[id].len() {
                                let l = next_diff1_3(&pg[id], k as i32) as usize;
                                let mut s = Vec::<String>::new();
                                for m in k..l {
                                    s.push(format!("{}={}", pg[id][m].1 as char, pg[id][m].2));
                                }
                                fwriteln!(
                                    pgout,
                                    "{},{},{},{},{},{}",
                                    i + 1,
                                    j + 1,
                                    q + 1,
                                    chain_types[refdata.rtype[id] as usize],
                                    pg[id][k].0,
                                    s.iter().format(":")
                                );
                                k = l;
                            }
                        }
                    } else {
                        if !ctl.gen_opt.peer_group_readable {
                            for y in pg[id].iter() {
                                fwriteln!(
                                    logx,
                                    "{},{},{},{},{},{},{}",
                                    i + 1,
                                    j + 1,
                                    q + 1,
                                    chain_types[refdata.rtype[id] as usize],
                                    y.0,
                                    y.1 as char,
                                    y.2
                                );
                            }
                        } else {
                            let mut k = 0;
                            while k < pg[id].len() {
                                let l = next_diff1_3(&pg[id], k as i32) as usize;
                                let mut s = Vec::<String>::new();
                                for m in k..l {
                                    s.push(format!("{}={}", pg[id][m].1 as char, pg[id][m].2));
                                }
                                fwriteln!(
                                    logx,
                                    "{},{},{},{},{},{}",
                                    i + 1,
                                    j + 1,
                                    q + 1,
                                    chain_types[refdata.rtype[id] as usize],
                                    pg[id][k].0,
                                    s.iter().format(":")
                                );
                                k = l;
                            }
                        }
                    }
                }
            }

            // Generate FASTA output.

            generate_fasta(
                i,
                j,
                oo,
                &exacts,
                &rsi,
                &exact_clonotypes,
                &ctl,
                &refdata,
                &mut logx,
                &mut fout,
                &mut faaout,
            );

            // Generate parseable output.

            if ctl.parseable_opt.pout.len() > 0 {
                let mut rows = Vec::<Vec<String>>::new();
                for m in 0..out_datas[oo].len() {
                    out_datas[oo][m].insert("group_id".to_string(), format!("{}", i + 1));
                    out_datas[oo][m]
                        .insert("group_ncells".to_string(), format!("{}", group_ncells));
                    out_datas[oo][m].insert("clonotype_id".to_string(), format!("{}", j + 1));
                }
                if ctl.parseable_opt.pout == "stdout".to_string() {
                    if !ctl.gen_opt.noprint || (i == 0 && j == 0) {
                        fwriteln!(logx, "{}", pcols.iter().format(","));
                    }
                }
                if ctl.parseable_opt.pout == "stdouth".to_string() {
                    rows.push(pcols.clone());
                }
                let x = &out_datas[oo];
                for (u, y) in x.iter().enumerate() {
                    if !ctl.parseable_opt.pbarcode {
                        if ctl.parseable_opt.pout != "stdouth".to_string() {
                            for (i, c) in pcols.iter().enumerate() {
                                if i > 0 {
                                    if ctl.parseable_opt.pout != "stdout".to_string() {
                                        fwrite!(pout, ",");
                                    } else {
                                        fwrite!(logx, ",");
                                    }
                                }
                                if y.contains_key(c) {
                                    let val = &y[c];
                                    let val = val.replace(POUT_SEP, ",");
                                    if !val.contains(',') {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, "{}", val);
                                        } else {
                                            fwrite!(logx, "{}", val);
                                        }
                                    } else {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, "\"{}\"", val);
                                        } else {
                                            fwrite!(logx, "\"{}\"", val);
                                        }
                                    }
                                } else {
                                    if ctl.parseable_opt.pout != "stdout".to_string() {
                                        fwrite!(pout, "");
                                    } else {
                                        fwrite!(logx, "");
                                    }
                                }
                            }
                            if ctl.parseable_opt.pout != "stdout".to_string() {
                                fwriteln!(pout, "");
                            } else {
                                fwriteln!(logx, "");
                            }
                        } else {
                            let mut row = Vec::<String>::new();
                            for c in pcols.iter() {
                                if y.contains_key(c) {
                                    let val = &y[c];
                                    row.push(val.clone());
                                } else {
                                    row.push("".to_string());
                                }
                            }
                            rows.push(row);
                        }
                    } else {
                        let ex = &exact_clonotypes[exacts[oo][u]];
                        let n = ex.ncells();
                        if ctl.parseable_opt.pout != "stdouth".to_string() {
                            // Traverse the cells in the exact subclonotype.
                            for m in 0..n {
                                // Traverse the parseable fields to be displayed.
                                for (i, c) in pcols.iter().enumerate() {
                                    if i > 0 {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, ",");
                                        } else {
                                            fwrite!(logx, ",");
                                        }
                                    }
                                    // Test for whether the out_data contain the field.
                                    if y.contains_key(c) {
                                        let mut id = 0;
                                        let vals = y[c].split(POUT_SEP).collect::<Vec<&str>>();
                                        if vals.len() > 1 {
                                            id = m;
                                        }
                                        if id >= vals.len() {
                                            panic!(
                                                "id >= vals.len() where id = {} and vals.len() \
                                                = {},\nparseable variable = {}, barcodes include \
                                                {}, n = {}, y[c] = {}",
                                                id,
                                                vals.len(),
                                                c,
                                                ex.clones[0][0].barcode,
                                                n,
                                                y[c],
                                            );
                                        }
                                        let val = vals[id];
                                        if !val.contains(',') {
                                            if ctl.parseable_opt.pout != "stdout".to_string() {
                                                fwrite!(pout, "{}", val);
                                            } else {
                                                fwrite!(logx, "{}", val);
                                            }
                                        } else {
                                            if ctl.parseable_opt.pout != "stdout".to_string() {
                                                fwrite!(pout, "\"{}\"", val);
                                            } else {
                                                fwrite!(logx, "\"{}\"", val);
                                            }
                                        }
                                    } else {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, "");
                                        } else {
                                            fwrite!(logx, "");
                                        }
                                    }
                                }
                                if ctl.parseable_opt.pout != "stdout".to_string() {
                                    fwriteln!(pout, "");
                                } else {
                                    fwriteln!(logx, "");
                                }
                            }
                        } else {
                            for m in 0..n {
                                let mut row = Vec::<String>::new();
                                for c in pcols.iter() {
                                    if y.contains_key(c) {
                                        let mut id = 0;
                                        let vals = y[c].split(POUT_SEP).collect::<Vec<&str>>();
                                        if vals.len() > 1 {
                                            id = m;
                                        }
                                        let val = vals[id];
                                        row.push(val.to_string());
                                    } else {
                                        row.push("".to_string());
                                    }
                                }
                                rows.push(row);
                            }
                        }
                    }
                }
                if ctl.parseable_opt.pout == "stdouth".to_string() {
                    let mut log = Vec::<u8>::new();
                    let mut justify = Vec::<u8>::new();
                    for x in rows[0].iter() {
                        justify.push(justification(&x));
                    }
                    print_tabular(&mut log, &rows, 2, Some(justify));
                    fwrite!(logx, "{}", strme(&log));
                }
            }
        }
    }

    // Execute PLOT_XY.

    if ctl.plot_opt.plot_xy_filename.len() > 0 {
        let mut xvar = ctl.plot_opt.plot_xy_xvar.clone();
        if ctl.plot_opt.plot_xy_x_log10 {
            xvar = format!("log10({})", xvar);
        }
        let mut yvar = ctl.plot_opt.plot_xy_yvar.clone();
        if ctl.plot_opt.plot_xy_y_log10 {
            yvar = format!("log10({})", yvar);
        }
        let filename = ctl.plot_opt.plot_xy_filename.clone();
        plot_points(&plot_xy_vals, &xvar, &yvar, &filename);
        if filename == "stdout" {
            let f = open_for_read!["stdout"];
            for line in f.lines() {
                let s = line.unwrap();
                println!("{}", s);
            }
            remove_file("stdout").unwrap();
        }
    }

    // Finish CLUSTAL.

    if clustal_aa.is_some() {
        clustal_aa.unwrap().finish().unwrap();
    }
    if clustal_dna.is_some() {
        clustal_dna.unwrap().finish().unwrap();
    }

    // Print stats.

    let mut nclono2 = 0;
    let mut two_chain = 0;
    print_stats(
        &tall,
        &exacts,
        &rsi,
        &exact_clonotypes,
        &ctl,
        &gex_info,
        &vdj_cells,
        &fate,
        &mut logx,
        &mut nclono2,
        &mut two_chain,
    );

    // Print to stdout.

    if !ctl.gen_opt.html {
        print!("{}", compress_ansi_escapes(&strme(&logx)));
    } else {
        // Note that we do not link to the css file, because it is less fragile then including
        // the font face information directly.  In particular, the css file could be accidentally
        // deleted or renamed, which would break previously generated user html files.  This
        // actually happened!
        let s = convert_text_with_ansi_escapes_to_html(
            strme(&logx),
            "", // source
            &ctl.gen_opt.html_title,
            &format!("<style type=\"text/css\">\n{}</style>", font_face_in_css()),
            "DejaVuSansMono",
            14,
        );
        print!("{}", s);
    }

    // Test requirements.

    test_requirements(&pics, &exacts, &exact_clonotypes, &ctl, nclono2, two_chain);
}
