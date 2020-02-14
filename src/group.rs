// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Group and print clonotypes.  For now, limited grouping functionality.

use crate::defs::*;
use ansi_escape::*;
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
) {
    // Set up for parseable output.

    #[allow(bare_trait_objects)]
    let mut pout = match ctl.parseable_opt.pout.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.parseable_opt.pout);
            (Box::new(File::create(&path).unwrap()) as Box<Write>)
        }
    };
    let mut pcols = ctl.parseable_opt.pcols.clone();
    if pcols.is_empty() {
        pcols = parseable_fields.clone();
    }
    if !ctl.parseable_opt.pout.is_empty() && ctl.parseable_opt.pout != "stdout".to_string() {
        fwriteln!(pout, "{}", pcols.iter().format(","));
    }

    // Set up for fasta output.

    #[allow(bare_trait_objects)]
    let mut fout = match ctl.gen_opt.fasta_filename.as_str() {
        "" => (Box::new(stdout()) as Box<Write>),
        "stdout" => (Box::new(stdout()) as Box<Write>),
        _ => {
            let path = Path::new(&ctl.gen_opt.fasta_filename);
            (Box::new(File::create(&path).unwrap()) as Box<Write>)
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
    if ctl.clono_group_opt.vj2 {
        let mut all = Vec::<((String, String, String, String), usize)>::new();
        for i in 0..pics.len() {
            for x in exacts[i].iter() {
                let ex = &exact_clonotypes[*x];
                for j1 in 0..ex.share.len() {
                    for j2 in 0..ex.share.len() {
                        if j1 == j2 {
                            continue;
                        }
                        let y1 = &ex.share[j1];
                        let y2 = &ex.share[j2];
                        all.push((
                            (
                                refdata.name[y1.v_ref_id].clone(),
                                refdata.name[y1.j_ref_id].clone(),
                                refdata.name[y2.v_ref_id].clone(),
                                refdata.name[y2.j_ref_id].clone(),
                            ),
                            i,
                        ));
                    }
                }
            }
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
        groups += 1;

        // Generate human readable output.

        if !ctl.gen_opt.noprint {
            println!("");
            if last_width > 0 {
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_eight_bit_color_escape(&mut log, 44);
                    print!("{}", strme(&log));
                }
                print!("╺");
                for _ in 0..last_width - 2 {
                    print!("━");
                }
                print!("╸");
                println!("\n");
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_end_escape(&mut log);
                    print!("{}", strme(&log));
                }
            }
            if ctl.pretty {
                let mut log = Vec::<u8>::new();
                emit_bold_escape(&mut log);
                emit_eight_bit_color_escape(&mut log, 27);
                print!("{}", strme(&log));
            }
            print!("[{}] GROUP = {} CLONOTYPES = {} CELLS", groups, o.len(), n);
            if ctl.pretty {
                let mut log = Vec::<u8>::new();
                emit_end_escape(&mut log);
                print!("{}", strme(&log));
            }
            println!("");
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
                print!("\n[{}.{}] {}", groups, j + 1, pics[oo]);
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

            // Generate fasta output.

            if ctl.gen_opt.fasta_filename.len() > 0 {
                for (k, u) in exacts[oo].iter().enumerate() {
                    for m in 0..mat[oo].len() {
                        if mat[oo][m][k].is_some() {
                            let r = mat[oo][m][k].unwrap();
                            let ex = &exact_clonotypes[*u];
                            fwriteln!(
                                fout,
                                ">group{}.clonotype{}.exact{}.chain{}",
                                groups,
                                j + 1,
                                k + 1,
                                m + 1
                            );
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
                                fwriteln!(fout, "{}", strme(&seq));
                            }
                        }
                    }
                }
            }

            // Generate parseable output.

            if ctl.parseable_opt.pout.len() > 0 {
                for m in 0..out_datas[oo].len() {
                    out_datas[oo][m].insert("group_id".to_string(), format!("{}", groups));
                    out_datas[oo][m]
                        .insert("group_ncells".to_string(), format!("{}", group_ncells));
                    out_datas[oo][m].insert("clonotype_id".to_string(), format!("{}", j + 1));
                }
                if ctl.parseable_opt.pout == "stdout".to_string() {
                    fwriteln!(pout, "{}", pcols.iter().format(","));
                }
                let x = &out_datas[oo];
                for y in x.iter() {
                    for (i, c) in pcols.iter().enumerate() {
                        if i > 0 {
                            fwrite!(pout, ",");
                        }
                        if y.contains_key(c) {
                            let val = &y[c];
                            if !val.contains(',') {
                                fwrite!(pout, "{}", val);
                            } else {
                                fwrite!(pout, "\"{}\"", val);
                            }
                        } else {
                            fwrite!(pout, "");
                        }
                    }
                    fwriteln!(pout, "");
                }
            }
        }
    }

    // Print summary stats.

    if ctl.gen_opt.summary {
        println!("\nSUMMARY STATISTICS");
        println!("1. overall");
        let nclono = exacts.len();
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
        println!("   • number of datasets = {}", ctl.sample_info.n());
        println!("   • number of donors = {}", ctl.sample_info.donors);
        println!("   • total elapsed time = {:.2} seconds", elapsed(&tall));
        println!("   • peak memory = {:.2} GB", peak_mem_usage_gb());
        println!("2. for the selected clonotypes");
        println!("   • number of clonotypes = {}", nclono);
        println!(
            "   • number of clonotypes having at least two cells = {}",
            nclono2
        );
        println!("   • number of cells in clonotypes = {}", ncells);
        nchains.sort();
        let mut i = 0;
        while i < nchains.len() {
            let j = next_diff(&nchains, i);
            println!(
                "   • number of clonotypes having {} chains = {}",
                nchains[i],
                j - i
            );
            i = j;
        }
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "sample".to_string(),
            "donor".to_string(),
            "ncells".to_string(),
        ];
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
        print_tabular_vbox(&mut log, &rows, 2, &b"llr".to_vec(), false);
        log = log.replace("\n", "\n   ");
        print!("   {}", log);
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
