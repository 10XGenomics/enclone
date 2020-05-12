// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Group and print clonotypes.  For now, limited grouping functionality.

use crate::print_utils5::*;
use amino::*;
use ansi_escape::ansi_to_html::*;
use ansi_escape::*;
use enclone_core::defs::*;
use equiv::EquivRel;
use io_utils::*;
use itertools::*;
use perf_stats::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::*;
use std::path::Path;
use std::time::Instant;
use string_utils::*;
use tables::*;
use vdj_ann::refx::*;
use vector_utils::*;

pub fn group_and_print_clonotypes(
    tall: &Instant,
    refdata: &RefData,
    pics: &Vec<String>,
    exacts: &Vec<Vec<usize>>,
    mat: &Vec<Vec<Vec<Option<usize>>>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    parseable_fields: &Vec<String>,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    join_info: &Vec<(usize, usize, bool, Vec<u8>)>,
) {
    // Build index to join info.

    let mut to_join_info = vec![Vec::<usize>::new(); exact_clonotypes.len()];
    for i in 0..join_info.len() {
        to_join_info[join_info[i].0].push(i);
        to_join_info[join_info[i].1].push(i);
    }

    // Set up for parseable output.

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

    // Group clonotypes and make output.

    let mut last_width = 0;
    let mut e: EquivRel = EquivRel::new(pics.len() as i32);
    if ctl.clono_group_opt.heavy_cdr3_aa {
        let mut all = Vec::<(String, usize)>::new();
        for i in 0..pics.len() {
            for x in exacts[i].iter() {
                for m in 0..exact_clonotypes[*x].share.len() {
                    let y = &exact_clonotypes[*x].share[m];
                    if y.left {
                        all.push((y.cdr3_aa.clone(), i));
                    }
                }
            }
        }
        all.sort();
        let mut i = 0;
        while i < all.len() {
            let j = next_diff1_2(&all, i as i32) as usize;
            for k in i + 1..j {
                e.join(all[i].1 as i32, all[k].1 as i32);
            }
            i = j;
        }
    }
    if ctl.clono_group_opt.vj_refname {
        let mut all = Vec::<(Vec<String>, usize)>::new();
        for i in 0..pics.len() {
            let ex = &exact_clonotypes[exacts[i][0]];
            let mut s = Vec::<String>::new();
            for j in 0..ex.share.len() {
                s.push(refdata.name[ex.share[j].v_ref_id].clone());
                s.push(refdata.name[ex.share[j].j_ref_id].clone());
            }
            s.sort();
            all.push((s, i));
        }
        // Note duplication with above code.
        all.sort();
        let mut i = 0;
        while i < all.len() {
            let j = next_diff1_2(&all, i as i32) as usize;
            for k in i + 1..j {
                e.join(all[i].1 as i32, all[k].1 as i32);
            }
            i = j;
        }
    }
    let mut groups = 0;
    let mut greps = Vec::<i32>::new();
    e.orbit_reps(&mut greps);

    // Sort so that larger groups (as measured by cells) come first.

    let mut grepsn = Vec::<(usize, usize)>::new();
    for i in 0..greps.len() {
        let mut o = Vec::<i32>::new();
        e.orbit(greps[i], &mut o);
        if o.len() < ctl.clono_group_opt.min_group {
            continue;
        }
        let mut n = 0;
        for j in 0..o.len() {
            let x = o[j] as usize;
            let s = &exacts[x];
            for k in 0..s.len() {
                n += exact_clonotypes[s[k]].clones.len();
            }
        }
        grepsn.push((n, i));
    }
    reverse_sort(&mut grepsn);

    // Now print.

    let mut logx = Vec::<u8>::new();
    for z in 0..grepsn.len() {
        let i = grepsn[z].1;
        let n = grepsn[z].0;
        let mut o = Vec::<i32>::new();
        e.orbit(greps[i], &mut o);
        groups += 1;

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
                    groups,
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
            if !ctl.gen_opt.noprint {
                if z > 0 || j > 0 || !(ctl.gen_opt.html && ctl.gen_opt.ngroup) {
                    fwrite!(logx, "\n"); // NEWLINE 6
                }
                if ctl.gen_opt.svg {
                    const FONT_SIZE: usize = 15;
                    let s = format!("[{}.{}] {}", groups, j + 1, pics[oo]);

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
                    fwrite!(logx, "[{}.{}] {}", groups, j + 1, pics[oo]);
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

            // Generate fasta output.

            if ctl.gen_opt.fasta_filename.len() > 0 {
                for (k, u) in exacts[oo].iter().enumerate() {
                    for m in 0..mat[oo].len() {
                        if mat[oo][m][k].is_some() {
                            let r = mat[oo][m][k].unwrap();
                            let ex = &exact_clonotypes[*u];
                            if ctl.gen_opt.fasta_filename != "stdout".to_string() {
                                fwriteln!(
                                    fout,
                                    ">group{}.clonotype{}.exact{}.chain{}",
                                    groups,
                                    j + 1,
                                    k + 1,
                                    m + 1
                                );
                            } else {
                                fwriteln!(
                                    logx,
                                    ">group{}.clonotype{}.exact{}.chain{}",
                                    groups,
                                    j + 1,
                                    k + 1,
                                    m + 1
                                );
                            }
                            let mut seq = ex.share[r].seq.clone();
                            let mut cid = ex.share[r].c_ref_id;
                            if cid.is_none() {
                                for l in 0..exacts[oo].len() {
                                    if mat[oo][m][l].is_some() {
                                        let r2 = mat[oo][m][l].unwrap();
                                        let ex2 = &exact_clonotypes[exacts[oo][l]];
                                        let cid2 = ex2.share[r2].c_ref_id;
                                        if cid2.is_some() {
                                            cid = cid2;
                                            break;
                                        }
                                    }
                                }
                            }
                            if cid.is_some() {
                                let mut cseq = refdata.refs[cid.unwrap()].to_ascii_vec();
                                seq.append(&mut cseq);
                                if ctl.gen_opt.fasta_filename != "stdout".to_string() {
                                    fwriteln!(fout, "{}", strme(&seq));
                                } else {
                                    fwriteln!(logx, "{}", strme(&seq));
                                }
                            }
                        }
                    }
                }
            }

            // Generate fasta amino acid output.

            if ctl.gen_opt.fasta_aa_filename.len() > 0 {
                for (k, u) in exacts[oo].iter().enumerate() {
                    for m in 0..mat[oo].len() {
                        if mat[oo][m][k].is_some() {
                            let r = mat[oo][m][k].unwrap();
                            let ex = &exact_clonotypes[*u];
                            if ctl.gen_opt.fasta_aa_filename != "stdout".to_string() {
                                fwriteln!(
                                    faaout,
                                    ">group{}.clonotype{}.exact{}.chain{}",
                                    groups,
                                    j + 1,
                                    k + 1,
                                    m + 1
                                );
                            } else {
                                fwriteln!(
                                    logx,
                                    ">group{}.clonotype{}.exact{}.chain{}",
                                    groups,
                                    j + 1,
                                    k + 1,
                                    m + 1
                                );
                            }
                            let mut seq = ex.share[r].seq.clone();
                            let mut cid = ex.share[r].c_ref_id;
                            if cid.is_none() {
                                for l in 0..exacts[oo].len() {
                                    if mat[oo][m][l].is_some() {
                                        let r2 = mat[oo][m][l].unwrap();
                                        let ex2 = &exact_clonotypes[exacts[oo][l]];
                                        let cid2 = ex2.share[r2].c_ref_id;
                                        if cid2.is_some() {
                                            cid = cid2;
                                            break;
                                        }
                                    }
                                }
                            }
                            if cid.is_some() {
                                let mut cseq = refdata.refs[cid.unwrap()].to_ascii_vec();
                                seq.append(&mut cseq);
                                if ctl.gen_opt.fasta_aa_filename != "stdout".to_string() {
                                    fwriteln!(faaout, "{}", strme(&aa_seq(&seq, 0)));
                                } else {
                                    fwriteln!(logx, "{}", strme(&aa_seq(&seq, 0)));
                                }
                            }
                        }
                    }
                }
            }

            // Generate parseable output.

            if ctl.parseable_opt.pout.len() > 0
                && (!ctl.gen_opt.noprint
                    || (ctl.parseable_opt.pout != "stdout".to_string()
                        && ctl.parseable_opt.pout != "stdouth".to_string()))
            {
                let mut rows = Vec::<Vec<String>>::new();
                for m in 0..out_datas[oo].len() {
                    out_datas[oo][m].insert("group_id".to_string(), format!("{}", groups));
                    out_datas[oo][m]
                        .insert("group_ncells".to_string(), format!("{}", group_ncells));
                    out_datas[oo][m].insert("clonotype_id".to_string(), format!("{}", j + 1));
                }
                if ctl.parseable_opt.pout == "stdout".to_string() {
                    fwriteln!(logx, "{}", pcols.iter().format(","));
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
                            for m in 0..n {
                                for (i, c) in pcols.iter().enumerate() {
                                    if i > 0 {
                                        if ctl.parseable_opt.pout != "stdout".to_string() {
                                            fwrite!(pout, ",");
                                        } else {
                                            fwrite!(logx, ",");
                                        }
                                    }
                                    if y.contains_key(c) {
                                        let mut id = 0;
                                        let vals = y[c].split(';').collect::<Vec<&str>>();
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
                                        let vals = y[c].split(';').collect::<Vec<&str>>();
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

    // Compute two umi stats.

    let nclono = exacts.len();
    let mut umish = Vec::<usize>::new();
    let mut umisl = Vec::<usize>::new();
    for i in 0..nclono {
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            for k in 0..ex.clones.len() {
                for l in 0..ex.share.len() {
                    if ex.share[l].left {
                        umish.push(ex.clones[k][l].umi_count);
                    } else {
                        umisl.push(ex.clones[k][l].umi_count);
                    }
                }
            }
        }
    }
    umish.sort();
    umisl.sort();
    let (mut middleh, mut denomh) = (0, 0);
    for j in umish.len() / 3..(2 * umish.len()) / 3 {
        middleh += umish[j];
        denomh += 1;
    }
    let mut middle_mean_umish = 0.0;
    if denomh > 0 {
        middle_mean_umish = (middleh as f64) / (denomh as f64);
    }
    let (mut middlel, mut denoml) = (0, 0);
    for j in umisl.len() / 3..(2 * umisl.len()) / 3 {
        middlel += umisl[j];
        denoml += 1;
    }
    let mut middle_mean_umisl = 0.0;
    if denoml > 0 {
        middle_mean_umisl = (middlel as f64) / (denoml as f64);
    }

    // Compute n23.

    let mut n23 = 0;
    for i in 0..nclono {
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            if ex.share.len() == 2 || ex.share.len() == 3 {
                n23 += ex.ncells();
            }
        }
    }

    // Print summary stats.

    if ctl.gen_opt.summary {
        fwriteln!(logx, "\nSUMMARY STATISTICS");
        fwriteln!(logx, "1. overall");
        let mut nclono2 = 0;
        let mut ncells = 0;
        let mut nchains = Vec::<usize>::new();
        let mut sd = Vec::<(Option<usize>, Option<usize>)>::new();
        for i in 0..nclono {
            let mut n = 0;
            for j in 0..exacts[i].len() {
                let ex = &exact_clonotypes[exacts[i][j]];
                n += ex.ncells();
                for k in 0..ex.clones.len() {
                    let x = &ex.clones[k][0];
                    sd.push((x.sample_index, x.donor_index));
                }
            }
            if n >= 2 {
                nclono2 += 1;
            }
            ncells += n;
            nchains.push(mat[i].len());
        }
        sd.sort();
        let mut sdx = Vec::<(Option<usize>, Option<usize>, usize)>::new();
        let mut i = 0;
        while i < sd.len() {
            let j = next_diff(&sd, i);
            sdx.push((sd[i].0, sd[i].1, j - i));
            i = j;
        }
        fwriteln!(logx, "   • number of datasets = {}", ctl.sample_info.n());
        fwriteln!(logx, "   • number of donors = {}", ctl.sample_info.donors);

        // Print mean reads per cell if known.

        let mut known = true;
        for i in 0..ctl.sample_info.n() {
            if ctl.sample_info.cells_cellranger[i].is_none() {
                known = false;
            } else if ctl.sample_info.mean_read_pairs_per_cell_cellranger[i].is_none() {
                known = false;
            }
        }
        if known {
            let (mut cells, mut read_pairs) = (0, 0);
            for i in 0..ctl.sample_info.n() {
                let c = ctl.sample_info.cells_cellranger[i].unwrap();
                let rpc = ctl.sample_info.mean_read_pairs_per_cell_cellranger[i].unwrap();
                cells += c;
                read_pairs += cells * rpc;
            }
            let rpc = ((read_pairs as f64) / (cells as f64)).round();
            fwriteln!(logx, "   • cells (from cellranger) = {}", cells);
            fwriteln!(logx, "   • read pairs per cell (from cellranger) = {}", rpc);
        }

        // Print computational performance stats.

        if !ctl.gen_opt.summary_clean {
            fwriteln!(
                logx,
                "   • total elapsed time = {:.1} seconds",
                elapsed(&tall)
            );
            #[cfg(not(target_os = "macos"))]
            fwriteln!(logx, "   • peak memory = {:.1} GB", peak_mem_usage_gb());
        }
        fwriteln!(logx, "2. for the selected clonotypes");
        fwriteln!(logx, "   • number of clonotypes = {}", nclono);
        fwriteln!(
            logx,
            "   • number of clonotypes having at least two cells = {}",
            nclono2
        );
        fwriteln!(logx, "   • number of cells = {}", ncells);
        fwriteln!(logx, "   • number of cells having 2 or 3 chains = {}", n23);
        nchains.sort();
        let mut i = 0;
        while i < nchains.len() {
            let j = next_diff(&nchains, i);
            fwriteln!(
                logx,
                "   • number of clonotypes having {} chains = {}",
                nchains[i],
                j - i
            );
            i = j;
        }
        fwriteln!(
            logx,
            "   • mean over middle third of contig UMI counts (heavy chain/ TRB) = {:.2}",
            middle_mean_umish,
        );
        fwriteln!(
            logx,
            "   • mean over middle third of contig UMI counts (light chain/ TRA) = {:.2}",
            middle_mean_umisl,
        );
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec!["sample".to_string(), "donor".to_string(), "n".to_string()];
        rows.push(row);
        let row = vec!["\\hline".to_string(); 3];
        rows.push(row);
        for i in 0..sdx.len() {
            let mut row = Vec::<String>::new();
            if sdx[i].0.is_some() {
                row.push(format!(
                    "{}",
                    ctl.sample_info.sample_list[sdx[i].0.unwrap()]
                ));
            } else {
                row.push("?".to_string());
            }
            if sdx[i].1.is_some() {
                row.push(format!("{}", ctl.sample_info.donor_list[sdx[i].1.unwrap()]));
            } else {
                row.push("?".to_string());
            }
            row.push(format!("{}", sdx[i].2));
            rows.push(row);
        }
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"llr".to_vec(), false, false);
        log = log.replace("\n", "\n   ");
        fwrite!(logx, "   {}", log);
    }

    // Print summary csv stats.

    if ctl.gen_opt.summary_csv {
        println!("\nmiddle_mean_umis_heavy,middle_mean_umis_light,n_twothreesie");
        println!("{:.2},{:.2},{}", middle_mean_umish, middle_mean_umisl, n23);
    }

    // Print to stdout.

    if !ctl.gen_opt.html {
        print!("{}", compress_ansi_escapes(&strme(&logx)));
    } else {
        let s = convert_text_with_ansi_escapes_to_html(
            strme(&logx),
            "", // source
            &ctl.gen_opt.html_title,
            "<link href='https://10xgenomics.github.io/enclone/pages/enclone.css' \
             rel='stylesheet' type='text/css'>",
            "DejaVuSansMono",
            14,
        );
        print!("{}", s);
    }

    // Test for required number of false positives.

    if ctl.gen_opt.required_fps.is_some() {
        let mut fps = 0;
        for i in 0..pics.len() {
            if pics[i].contains("WARNING:") {
                fps += 1;
            }
        }
        if fps != ctl.gen_opt.required_fps.unwrap() {
            eprintln!(
                "\nA \"false positive\" is a clonotype that contains cells from multiple\n\
                 donors.  You invoked enclone with the argument REQUIRED_FPS={}, but we found\n\
                 {} false positives, so the requirement is not met.\n",
                ctl.gen_opt.required_fps.unwrap(),
                fps
            );
            std::process::exit(1);
        }
    }
}
