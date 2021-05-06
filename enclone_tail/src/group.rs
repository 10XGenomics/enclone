// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Group and print clonotypes.  For now, limited grouping functionality.
// (NO LONGER DOES THE GROUPING.)
//
// To keep compilation time down, this crate should not reach into the enclone crate.

use crate::align_n::*;
use crate::clustal::*;
use crate::fasta::*;
use crate::phylip::*;
use crate::plot::*;
use crate::plot_points::*;
use crate::print_stats::*;
use crate::requirements::*;
use crate::tree::*;
use ansi_escape::ansi_to_html::*;
use ansi_escape::*;
use enclone_core::defs::*;
use enclone_core::mammalian_fixed_len::*;
use enclone_core::print_tools::*;
use enclone_proto::types::*;
use io_utils::*;
use itertools::*;
use std::collections::HashMap;
use std::env;
use std::fs::remove_file;
use std::fs::File;
use std::io::Write;
use std::io::*;
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
) {
    // Build index to join info.

    let t = Instant::now();
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
        &refdata,
        &exacts,
        &rsi,
        &exact_clonotypes,
        &ctl,
        &dref,
        &groups,
        width,
        false,
    );
    let jun_align_out = align_n(
        &refdata,
        &exacts,
        &rsi,
        &exact_clonotypes,
        &ctl,
        &dref,
        &groups,
        width,
        true,
    );

    // Test for consistency if requested.

    let mut align_out_test = HashMap::<(usize, usize), Vec<u8>>::new();
    let mut jun_align_out_test = HashMap::<(usize, usize), Vec<u8>>::new();
    if ctl.gen_opt.align_jun_align_consistency {
        let width = 1000;
        align_out_test = align_n(
            &refdata,
            &exacts,
            &rsi,
            &exact_clonotypes,
            &ctl,
            &dref,
            &groups,
            width,
            false,
        );
        jun_align_out_test = align_n(
            &refdata,
            &exacts,
            &rsi,
            &exact_clonotypes,
            &ctl,
            &dref,
            &groups,
            width,
            true,
        );
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
            if !ctl.gen_opt.html && !ctl.clono_group_opt.ngroup {
                fwriteln!(logx, ""); // NEWLINE 1
            }

            // If we just printed a clonotype box, output a bar.

            if last_width > 0 {
                if ctl.clono_group_opt.ngroup || ctl.gen_opt.html {
                    fwriteln!(logx, ""); // NEWLINE 2
                }
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_eight_bit_color_escape(&mut log, 44);
                    fwrite!(logx, "{}", strme(&log));
                }
                fwrite!(logx, "╺{}╸", "━".repeat(last_width - 2));
                if !ctl.clono_group_opt.ngroup {
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

            if !ctl.clono_group_opt.ngroup {
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
                if i > 0 || j > 0 || !(ctl.gen_opt.html && ctl.clono_group_opt.ngroup) {
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

            // Implement ALIGN<n> and JALIGN<n>.

            if !ctl.gen_opt.noprint {
                logx.append(&mut align_out[&(i, j)].clone());
                logx.append(&mut jun_align_out[&(i, j)].clone());
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
                        logx.append(&mut err.as_bytes().to_vec());
                        ok = false;
                        break;
                    }
                }
                if !ok {
                    logx.append(&mut b"consistency test failed\n".to_vec());
                    println!("{}", strme(&logx));
                    std::process::exit(1);
                }
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

            print_tree(
                oo,
                &exacts,
                &rsi,
                &exact_clonotypes,
                &ctl,
                &refdata,
                &dref,
                &out_datas,
                &mut logx,
            );

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
        &opt_d_val,
    );

    // Print to stdout.

    if !ctl.gen_opt.html {
        print!("{}", compress_ansi_escapes(&strme(&logx)));
    } else {
        // Remove initial newline if present.
        loop {
            if logx.len() > 0 && logx[0] == b'\n' {
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

    // Plot clonotypes.

    let mut svg = String::new();
    let plot_opt = ctl.plot_opt.clone();
    plot_clonotypes(
        &ctl,
        &plot_opt,
        &refdata,
        &exacts,
        &exact_clonotypes,
        &groups,
        &mut svg,
    );

    // Output clonotype plot (if it was generated and directed to stdout).

    if ctl.plot_opt.plot_file == "stdout".to_string() {
        print!("{}", svg);
        if !ctl.gen_opt.noprint {
            println!("");
        }
    }

    // Test requirements.

    test_requirements(&pics, &exacts, &exact_clonotypes, &ctl, nclono2, two_chain);
    ctl.perf_stats(&t, "in group code 2");
}
