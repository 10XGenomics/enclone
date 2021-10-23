// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Group and print clonotypes.  For now, limited grouping functionality.
// (NO LONGER DOES THE GROUPING.)
//
// To keep compilation time down, this crate should not reach into the enclone crate.

use crate::align_n::align_n;
use crate::clustal::print_clustal;
use crate::fasta::generate_fasta;
use crate::phylip::print_phylip;
use crate::plot::plot_clonotypes;
use crate::plot_points::plot_points;
use crate::print_stats::print_stats;
use crate::requirements::test_requirements;
use crate::sim_mat_plot::sim_mat_plot;
use crate::tree::print_tree;
use ansi_escape::ansi_to_html::{
    compress_ansi_escapes, convert_text_with_ansi_escapes_to_html,
    convert_text_with_ansi_escapes_to_svg,
};
use ansi_escape::{emit_bold_escape, emit_eight_bit_color_escape, emit_end_escape};
use enclone_core::combine_group_pics::combine_group_pics;
use enclone_core::defs::{
    justification, ColInfo, EncloneControl, ExactClonotype, GexInfo, POUT_SEP,
};
use enclone_core::mammalian_fixed_len::mammalian_fixed_len_peer_groups;
use enclone_core::print_tools::font_face_in_css;
use enclone_core::set_speakers::set_speakers;
use enclone_proto::types::DonorReferenceItem;
use io_utils::{fwrite, fwriteln, open_for_write_new};
use itertools::Itertools;
use std::cmp::max;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::stdout;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use string_utils::{stringme, strme, TextUtils};
use tables::print_tabular;
use tar::Builder;
use vdj_ann::refx::RefData;
use vector_utils::{next_diff1_3, unique_sort};

pub fn group_and_print_clonotypes(
    tall: &Instant,
    refdata: &RefData,
    pics: &Vec<String>,
    group_pics: &mut Vec<String>,
    last_widths: &mut Vec<u32>,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    join_info: &Vec<(usize, usize, bool, Vec<u8>)>,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    fate: &Vec<HashMap<String, String>>,
    dref: &Vec<DonorReferenceItem>,
    groups: &Vec<Vec<(i32, String)>>,
    opt_d_val: &Vec<(usize, Vec<Vec<Vec<usize>>>)>,
    svgs: &mut Vec<String>,
    summary: &mut String,
) -> Result<(), String> {
    // Build index to join info.

    let t = Instant::now();
    let mut to_join_info = vec![Vec::<usize>::new(); exact_clonotypes.len()];
    for i in 0..join_info.len() {
        to_join_info[join_info[i].0].push(i);
        to_join_info[join_info[i].1].push(i);
    }

    // Set up for parseable output.

    let mut parseable_fields = Vec::<String>::new();
    let mut max_chains = 0;
    for i in 0..rsi.len() {
        max_chains = max(max_chains, rsi[i].mat.len());
    }
    set_speakers(ctl, &mut parseable_fields, max_chains);
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
        if pcols[i].contains(':') {
            pcols2.push(pcols[i].before(":").to_string());
        } else {
            pcols2.push(pcols[i].clone());
        }
    }
    pcols = pcols2;
    if !ctl.parseable_opt.pout.is_empty()
        && ctl.parseable_opt.pout != *"stdout"
        && ctl.parseable_opt.pout != *"stdouth"
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
    if !ctl.gen_opt.clustal_aa.is_empty() && ctl.gen_opt.clustal_aa != *"stdout" {
        let file = File::create(&ctl.gen_opt.clustal_aa).unwrap();
        clustal_aa = Some(Builder::new(file));
    }
    if !ctl.gen_opt.clustal_dna.is_empty() && ctl.gen_opt.clustal_dna != *"stdout" {
        let file = File::create(&ctl.gen_opt.clustal_dna).unwrap();
        clustal_dna = Some(Builder::new(file));
    }

    // Set up for phylip output.

    let (mut phylip_aa, mut phylip_dna) = (None, None);
    if !ctl.gen_opt.phylip_aa.is_empty() && ctl.gen_opt.phylip_aa != *"stdout" {
        let file = File::create(&ctl.gen_opt.phylip_aa).unwrap();
        phylip_aa = Some(Builder::new(file));
    }
    if !ctl.gen_opt.phylip_dna.is_empty() && ctl.gen_opt.phylip_dna != *"stdout" {
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
    if !ctl.gen_opt.peer_group_filename.is_empty() && ctl.gen_opt.peer_group_filename != *"stdout" {
        if !ctl.gen_opt.peer_group_readable {
            fwriteln!(pgout, "group,clonotype,chain,pos,amino_acid,count");
        } else {
            fwriteln!(pgout, "group,clonotype,chain,pos,distribution");
        }
    }
    ctl.perf_stats(&t, "in group code 1");

    // Echo command.

    let t = Instant::now();
    let mut last_width = 0;
    let mut logx = Vec::<u8>::new();
    if ctl.gen_opt.echo {
        let args: Vec<String> = env::args().collect();
        fwriteln!(logx, "\n{}", args.iter().format(" "));
        if ctl.gen_opt.html {
            fwriteln!(logx, "");
        }
    }

    // Parallized precompute for ALIGN<n>.

    let width = 100;
    let align_out = align_n(
        refdata,
        exacts,
        rsi,
        exact_clonotypes,
        ctl,
        dref,
        groups,
        width,
        false,
    );
    let jun_align_out = align_n(
        refdata,
        exacts,
        rsi,
        exact_clonotypes,
        ctl,
        dref,
        groups,
        width,
        true,
    );

    // Test for consistency if requested.

    let mut align_out_test = HashMap::<(usize, usize), Vec<u8>>::new();
    let mut jun_align_out_test = HashMap::<(usize, usize), Vec<u8>>::new();
    if ctl.gen_opt.align_jun_align_consistency {
        let width = 1000;
        align_out_test = align_n(
            refdata,
            exacts,
            rsi,
            exact_clonotypes,
            ctl,
            dref,
            groups,
            width,
            false,
        );
        jun_align_out_test = align_n(
            refdata,
            exacts,
            rsi,
            exact_clonotypes,
            ctl,
            dref,
            groups,
            width,
            true,
        );
    }

    // Now print clonotypes.

    let mut plot_xy_vals = Vec::<(f32, f32)>::new();
    let mut plot_xy_comments = Vec::<String>::new();
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

        // Generate human readable output.

        let mut glog = Vec::<u8>::new();
        if !ctl.gen_opt.noprint {
            // If NGROUP is not on, output a GROUP line, including a newline at the end.

            if !ctl.clono_group_opt.ngroup {
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_bold_escape(&mut log);
                    emit_eight_bit_color_escape(&mut log, 27);
                    fwrite!(glog, "{}", strme(&log));
                }
                fwrite!(
                    glog,
                    "[{}] GROUP = {} CLONOTYPES = {} CELLS",
                    i + 1,
                    o.len(),
                    n
                );
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_end_escape(&mut log);
                    fwrite!(glog, "{}", strme(&log));
                }
                fwriteln!(glog, ""); // NEWLINE 5
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

            if !ctl.plot_opt.plot_xy_filename.is_empty() {
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
                                    let group_id = i;
                                    let clonotype_id = j;
                                    let com = format!(
                                        "data-tooltip='{{\"group_id\":\"{}\",\"clonotype_id\":\"{}\"}}'",
                                        group_id + 1,
                                        clonotype_id + 1,
                                    );
                                    plot_xy_comments.push(com);
                                }
                            }
                        }
                    }
                }
            }

            // Proceed.

            if !ctl.gen_opt.noprint {
                if i > 0 || j > 0 || !(ctl.gen_opt.html && ctl.clono_group_opt.ngroup) {
                    fwrite!(glog, "\n"); // NEWLINE 6
                }
                let mut s = format!("[{}.{}] {}", i + 1, j + 1, pics[oo]);
                if !groups[i][j].1.is_empty() {
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
                        glog,
                        "{}",
                        convert_text_with_ansi_escapes_to_svg(&s, "Menlo", FONT_SIZE)
                    );
                } else {
                    fwrite!(glog, "{}", s);
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
                fwriteln!(glog, "{}", strme(&join_info[ji[i]].3));
            }

            // Implement ALIGN<n> and JALIGN<n>.

            if !ctl.gen_opt.noprint {
                glog.append(&mut align_out[&(i, j)].clone());
                glog.append(&mut jun_align_out[&(i, j)].clone());
            }
            if ctl.gen_opt.align_jun_align_consistency {
                let x = stringme(&align_out_test[&(i, j)]);
                let y = stringme(&jun_align_out_test[&(i, j)]);
                let mut xlines = Vec::<String>::new();
                for line in x.lines() {
                    xlines.push(line.to_string());
                }
                let mut ylines = Vec::<String>::new();
                for line in y.lines() {
                    ylines.push(line.to_string());
                }
                let mut ok = true;
                for n in 1..=3 {
                    if !xlines[xlines.len() - n].contains(&ylines[ylines.len() - n]) {
                        let err = format!(
                            "\nERROR\n{}\ndoes not contain\n{}\n",
                            xlines[xlines.len() - n],
                            ylines[ylines.len() - n],
                        );
                        glog.append(&mut err.as_bytes().to_vec());
                        ok = false;
                        break;
                    }
                }
                if !ok {
                    glog.append(&mut b"consistency test failed\n".to_vec());
                    return Err(strme(&glog).to_string());
                }
            }

            // Generate clustal and phylip output.

            print_clustal(
                i,
                j,
                oo,
                exacts,
                rsi,
                exact_clonotypes,
                ctl,
                &mut glog,
                &mut clustal_aa,
                &mut clustal_dna,
            );
            print_phylip(
                i,
                j,
                oo,
                exacts,
                rsi,
                exact_clonotypes,
                ctl,
                &mut glog,
                &mut phylip_aa,
                &mut phylip_dna,
            );

            // Generate experimental tree output (options NEWICK0 and TREE).

            print_tree(
                oo,
                exacts,
                rsi,
                exact_clonotypes,
                ctl,
                refdata,
                dref,
                out_datas,
                &mut glog,
            );

            // Generate peer group output.

            if !ctl.gen_opt.peer_group_filename.is_empty() {
                let pg = mammalian_fixed_len_peer_groups(refdata);
                if !ctl.gen_opt.peer_group_readable {
                    if ctl.gen_opt.peer_group_filename == *"stdout" {
                        fwriteln!(glog, "group,clonotype,chain,pos,amino_acid,count");
                    }
                } else if ctl.gen_opt.peer_group_filename == *"stdout" {
                    fwriteln!(glog, "group,clonotype,chain,pos,distribution");
                }
                let chain_types = ["IGH", "IGK", "IGL", "TRA", "TRB"];
                for q in 0..rsi[oo].mat.len() {
                    let id = rsi[oo].vids[q];
                    if ctl.gen_opt.peer_group_filename != *"stdout" {
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
                    } else if !ctl.gen_opt.peer_group_readable {
                        for y in pg[id].iter() {
                            fwriteln!(
                                glog,
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
                                glog,
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

            // Generate FASTA output.

            generate_fasta(
                i,
                j,
                oo,
                exacts,
                rsi,
                exact_clonotypes,
                ctl,
                refdata,
                &mut glog,
                &mut fout,
                &mut faaout,
            );

            // Generate parseable output.

            if !ctl.parseable_opt.pout.is_empty() {
                let mut rows = Vec::<Vec<String>>::new();
                for m in 0..out_datas[oo].len() {
                    out_datas[oo][m].insert("group_id".to_string(), format!("{}", i + 1));
                    out_datas[oo][m]
                        .insert("group_ncells".to_string(), format!("{}", group_ncells));
                    out_datas[oo][m].insert("clonotype_id".to_string(), format!("{}", j + 1));
                }
                if ctl.parseable_opt.pout == *"stdout"
                    && (!ctl.gen_opt.noprint || (i == 0 && j == 0))
                {
                    fwriteln!(glog, "{}", pcols.iter().format(","));
                }
                if ctl.parseable_opt.pout == *"stdouth" {
                    rows.push(pcols.clone());
                }
                let x = &out_datas[oo];
                for (u, y) in x.iter().enumerate() {
                    if !ctl.parseable_opt.pbarcode {
                        if ctl.parseable_opt.pout != *"stdouth" {
                            for (i, c) in pcols.iter().enumerate() {
                                if i > 0 {
                                    if ctl.parseable_opt.pout != *"stdout" {
                                        fwrite!(pout, ",");
                                    } else {
                                        fwrite!(glog, ",");
                                    }
                                }
                                if y.contains_key(c) {
                                    let val = &y[c];
                                    let val = val.replace(POUT_SEP, ",");
                                    if !val.contains(',') {
                                        if ctl.parseable_opt.pout != *"stdout" {
                                            fwrite!(pout, "{}", val);
                                        } else {
                                            fwrite!(glog, "{}", val);
                                        }
                                    } else if ctl.parseable_opt.pout != *"stdout" {
                                        fwrite!(pout, "\"{}\"", val);
                                    } else {
                                        fwrite!(glog, "\"{}\"", val);
                                    }
                                } else if ctl.parseable_opt.pout != *"stdout" {
                                    fwrite!(pout, "");
                                } else {
                                    fwrite!(glog, "");
                                }
                            }
                            if ctl.parseable_opt.pout != *"stdout" {
                                fwriteln!(pout, "");
                            } else {
                                fwriteln!(glog, "");
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
                        if ctl.parseable_opt.pout != *"stdouth" {
                            // Traverse the cells in the exact subclonotype.
                            for m in 0..n {
                                // Traverse the parseable fields to be displayed.
                                for (i, c) in pcols.iter().enumerate() {
                                    if i > 0 {
                                        if ctl.parseable_opt.pout != *"stdout" {
                                            fwrite!(pout, ",");
                                        } else {
                                            fwrite!(glog, ",");
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
                                            if ctl.parseable_opt.pout != *"stdout" {
                                                fwrite!(pout, "{}", val);
                                            } else {
                                                fwrite!(glog, "{}", val);
                                            }
                                        } else if ctl.parseable_opt.pout != *"stdout" {
                                            fwrite!(pout, "\"{}\"", val);
                                        } else {
                                            fwrite!(glog, "\"{}\"", val);
                                        }
                                    } else if ctl.parseable_opt.pout != *"stdout" {
                                        fwrite!(pout, "");
                                    } else {
                                        fwrite!(glog, "");
                                    }
                                }
                                if ctl.parseable_opt.pout != *"stdout" {
                                    fwriteln!(pout, "");
                                } else {
                                    fwriteln!(glog, "");
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
                if ctl.parseable_opt.pout == *"stdouth" {
                    if ctl.gen_opt.noprint {
                        for (k, r) in rows.iter().enumerate() {
                            if k > 0 || (i == 0 && j == 0) {
                                let s = format!("{}\n", r.iter().format("\t"));
                                glog.append(&mut s.as_bytes().to_vec());
                            }
                        }
                    } else {
                        let mut log = Vec::<u8>::new();
                        let mut justify = Vec::<u8>::new();
                        for x in rows[0].iter() {
                            justify.push(justification(x));
                        }
                        print_tabular(&mut log, &rows, 2, Some(justify));
                        if ctl.gen_opt.noprint && (i > 0 || j > 0) {
                            let mut x = String::new();
                            for (k, line) in strme(&log).lines().enumerate() {
                                if k > 0 {
                                    x += &mut format!("{}\n", line);
                                }
                            }
                            log = x.as_bytes().to_vec();
                        }
                        fwrite!(glog, "{}", strme(&log));
                    }
                }
            }
        }
        group_pics.push(stringme(&glog));
        last_widths.push(last_width as u32);
    }
    if ctl.gen_opt.group_post_filter.as_ref().is_some() {
        let x = &ctl.gen_opt.group_post_filter.as_ref().unwrap();
        if !x.is_empty() {
            if !x.is_empty() && x[x.len() - 1] > group_pics.len() {
                return Err(
                    "\nArgument to G= references a group id that exceeds the number of groups.\n"
                        .to_string(),
                );
            }
            let mut group_pics2 = Vec::<String>::new();
            let mut last_widths2 = Vec::<u32>::new();
            for i in 0..x.len() {
                group_pics2.push(group_pics[x[i] - 1].clone());
                last_widths2.push(last_widths[x[i] - 1]);
            }
            *group_pics = group_pics2;
            *last_widths = last_widths2;
        }
    }
    if !ctl.gen_opt.noprintx {
        logx.append(
            &mut combine_group_pics(
                group_pics,
                last_widths,
                ctl.parseable_opt.pout == "stdouth",
                ctl.gen_opt.noprint,
                ctl.gen_opt.noprintx,
                ctl.gen_opt.html,
                ctl.clono_group_opt.ngroup,
                ctl.pretty,
            )
            .as_bytes()
            .to_vec(),
        );
    }

    // Execute SIM_MAT_PLOT.

    sim_mat_plot(ctl, groups, out_datas, svgs);

    // Execute PLOT_XY.

    if !ctl.plot_opt.plot_xy_filename.is_empty() {
        let mut xvar = ctl.plot_opt.plot_xy_xvar.clone();
        if ctl.plot_opt.plot_xy_x_log10 {
            xvar = format!("log10({})", xvar);
        }
        let mut yvar = ctl.plot_opt.plot_xy_yvar.clone();
        if ctl.plot_opt.plot_xy_y_log10 {
            yvar = format!("log10({})", yvar);
        }
        let filename = ctl.plot_opt.plot_xy_filename.clone();
        let mut svg = String::new();
        plot_points(
            &plot_xy_vals,
            &xvar,
            &yvar,
            &mut svg,
            ctl.plot_opt.plot_xy_sym,
        )?;
        if filename == "stdout" || filename == "gui_stdout" {
            for line in svg.lines() {
                println!("{}", line);
            }
        } else if filename == "gui" {
            // Add tooltip notes.

            let mut svg2 = String::new();
            let mut count = 0;
            for line in svg.lines() {
                let mut s = line.to_string();
                if s.starts_with("<circle ") {
                    s = format!("<circle {}{}", plot_xy_comments[count], s.after("<circle"));
                    count += 1;
                }
                svg2 += &mut format!("{}\n", s);
            }

            // Save.

            svgs.push(svg2);
        } else {
            let mut f = open_for_write_new![&filename];
            fwrite!(f, "{}", svg);
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
    let mut three_chain = 0;
    let mut four_chain = 0;
    let mut slog = Vec::<u8>::new();
    print_stats(
        tall,
        exacts,
        rsi,
        exact_clonotypes,
        groups,
        ctl,
        gex_info,
        vdj_cells,
        fate,
        &mut slog,
        &mut nclono2,
        &mut two_chain,
        &mut three_chain,
        &mut four_chain,
        opt_d_val,
    );
    *summary = stringme(&slog);
    logx.append(&mut slog);

    // Print to stdout.

    if !ctl.gen_opt.html {
        if !ctl.visual_mode {
            print!("{}", compress_ansi_escapes(strme(&logx)));
        }
    } else {
        // Remove initial newline if present.
        loop {
            if !logx.is_empty() && logx[0] == b'\n' {
                logx = logx[1..].to_vec();
            } else {
                break;
            }
        }

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
    ctl.perf_stats(&t, "in group code, before plotting clonotypes");

    // Plot clonotypes.

    let mut svg = String::new();
    let plot_opt = ctl.plot_opt.clone();
    plot_clonotypes(
        ctl,
        &plot_opt,
        refdata,
        exacts,
        exact_clonotypes,
        out_datas,
        groups,
        &mut svg,
    )?;

    // Output clonotype plot (if it was generated and directed to stdout).

    let t = Instant::now();
    if ctl.plot_opt.plot_file == "stdout" || ctl.plot_opt.plot_file == "gui_stdout" {
        print!("{}", svg);
        if !ctl.gen_opt.noprint {
            println!();
        }
    } else if ctl.plot_opt.plot_file == "gui" {
        svgs.push(svg);
    }

    // Test requirements.

    test_requirements(
        pics,
        exacts,
        exact_clonotypes,
        ctl,
        nclono2,
        two_chain,
        three_chain,
        four_chain,
    )?;
    ctl.perf_stats(&t, "in group code 2");
    Ok(())
}
