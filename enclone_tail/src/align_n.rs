// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Parallized precompute for ALIGN<n> and JALIGN<n>.

use align_tools::*;
use ansi_escape::*;
use enclone_core::align_to_vdj_ref::*;
use enclone_core::defs::*;
use enclone_core::opt_d::*;
use enclone_proto::types::*;
use io_utils::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn print_vis_align(
    seq: &[u8],
    concat: &[u8],
    col: usize,
    k: usize,
    vref: &[u8],
    dref: &[u8],
    d2ref: &[u8],
    jref: &[u8],
    drefname: &str,
    left: bool, // will be ex.share[r].left
    ctl: &EncloneControl,
    logx: &mut Vec<u8>,
    width: usize,
    add: &str,
) {
    // Make alignment.

    let (ops, _score) = align_to_vdj_ref(&seq, &vref, &dref, &d2ref, &jref, &drefname, left, &ctl);

    // Make visual alignment.

    let mut vis = vis_align(&seq, &concat, &ops, width);

    // Colorize it.

    let mut vdj;
    if left {
        vdj = "VDJ".to_string();
    } else {
        vdj = "VJ".to_string();
    }
    if ctl.pretty {
        let mut vis_new = String::new();
        let mut pos = 0;
        let vcolor = 3;
        let dcolor = 1;
        let d2color = 0;
        let jcolor = 2;
        let mut vdj_bytes = Vec::<u8>::new();
        print_color(vcolor, &mut vdj_bytes);
        vdj_bytes.push(b'V');
        if left {
            print_color(dcolor, &mut vdj_bytes);
            vdj_bytes.push(b'D');
        }
        if d2ref.len() > 0 {
            print_color(d2color, &mut vdj_bytes);
            vdj_bytes.push(b'D');
        }
        print_color(jcolor, &mut vdj_bytes);
        vdj_bytes.push(b'J');
        emit_end_escape(&mut vdj_bytes);
        vdj = stringme(&vdj_bytes);
        for (i, line) in vis.lines().enumerate() {
            if i % 4 != 2 {
                vis_new += &line.clone();
                vis_new += "\n";
            } else {
                let mut log = Vec::<u8>::new();
                let line = line.as_bytes();
                if pos < vref.len() {
                    print_color(vcolor, &mut log);
                } else if d2ref.len() > 0
                    && pos >= vref.len() + dref.len()
                    && pos < vref.len() + dref.len() + d2ref.len()
                {
                    print_color(d2color, &mut log);
                } else if left && pos < vref.len() + dref.len() + d2ref.len() {
                    print_color(dcolor, &mut log);
                } else {
                    print_color(jcolor, &mut log);
                }
                for j in 0..line.len() {
                    log.push(line[j]);
                    if line[j] != b' ' {
                        pos += 1;
                        if j != line.len() - 1 {
                            if left {
                                if d2ref.len() > 0 && pos == vref.len() + dref.len() {
                                    print_color(d2color, &mut log);
                                } else if pos == vref.len() + dref.len() + d2ref.len() {
                                    print_color(jcolor, &mut log);
                                } else if pos == vref.len() {
                                    print_color(dcolor, &mut log);
                                }
                            } else {
                                if pos == vref.len() {
                                    print_color(jcolor, &mut log);
                                }
                            }
                        }
                    }
                }
                emit_end_escape(&mut log);
                log.push(b'\n');
                vis_new += &strme(&log);
            }
        }
        vis = vis_new;
    }

    // Print alignment.

    fwrite!(
        logx,
        "\nCHAIN {} OF #{}  •  CONCATENATED {} REFERENCE{}\n{}",
        col,
        k + 1,
        vdj,
        add,
        vis,
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn align_n(
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    dref: &Vec<DonorReferenceItem>,
    groups: &Vec<Vec<(i32, String)>>,
    width: usize,
    jun: bool,
) -> HashMap<(usize, usize), Vec<u8>> {
    let mut align_out = HashMap::<(usize, usize), Vec<u8>>::new();
    let mut results = Vec::<(usize, Vec<Vec<u8>>)>::new();
    for i in 0..groups.len() {
        results.push((i, Vec::new()));
    }
    let (cols1, cols2);
    if !jun {
        cols1 = &ctl.gen_opt.chains_to_align;
        cols2 = &ctl.gen_opt.chains_to_align2;
    } else {
        cols1 = &ctl.gen_opt.chains_to_jun_align;
        cols2 = &ctl.gen_opt.chains_to_jun_align2;
    }
    let mut all_cols = cols1.clone();
    all_cols.append(&mut cols2.clone());
    unique_sort(&mut all_cols);
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut o = Vec::<i32>::new();
        for j in 0..groups[i].len() {
            o.push(groups[i][j].0);
        }
        for j in 0..o.len() {
            let mut logx = Vec::<u8>::new();
            let oo = o[j] as usize;
            for col in all_cols.iter() {
                let m = col - 1;
                for k in 0..exacts[oo].len() {
                    let ex = &exact_clonotypes[exacts[oo][k]];
                    if m < rsi[oo].mat.len() && rsi[oo].mat[m][k].is_some() {
                        let r = rsi[oo].mat[m][k].unwrap();
                        for pass in 1..=2 {
                            let mut seq = ex.share[r].seq_del.clone();
                            let mut vref = refdata.refs[rsi[oo].vids[m]].to_ascii_vec();
                            if rsi[oo].vpids[m].is_none() {
                            } else {
                                vref = dref[rsi[oo].vpids[m].unwrap()].nt_sequence.clone();
                            }
                            let vstart = vref.len() - vflank(&seq, &vref);
                            if jun {
                                vref = vref[vstart..vref.len()].to_vec();
                            }
                            if pass == 1 && !cols1.contains(&*col) {
                                continue;
                            }
                            if pass == 2 && !cols2.contains(&*col) {
                                continue;
                            }
                            let mut concat = vref.clone();
                            let mut drefx = Vec::<u8>::new();
                            let mut d2ref = Vec::<u8>::new();
                            let mut drefname = String::new();
                            if ex.share[r].left {
                                let mut opt = Vec::new();
                                let mut opt2 = Vec::new();
                                let mut delta = 0.0;
                                opt_d(
                                    &ex, m, k, &rsi[oo], &refdata, &dref, &mut opt, &mut opt2,
                                    &mut delta, &ctl,
                                );
                                if pass == 1 {
                                    for j in 0..opt.len() {
                                        let d = opt[j];
                                        if j == 0 {
                                            drefx = refdata.refs[d].to_ascii_vec();
                                        } else {
                                            d2ref = refdata.refs[d].to_ascii_vec();
                                        }
                                        if j > 0 {
                                            drefname += ":";
                                        }
                                        drefname += &mut refdata.name[d].clone();
                                    }
                                } else {
                                    for j in 0..opt2.len() {
                                        let d = opt2[j];
                                        if j == 0 {
                                            drefx = refdata.refs[d].to_ascii_vec();
                                        } else {
                                            d2ref = refdata.refs[d].to_ascii_vec();
                                        }
                                        if j > 0 {
                                            drefname += ":";
                                        }
                                        drefname += &mut refdata.name[d].clone();
                                    }
                                }
                                concat.append(&mut drefx.clone());
                                concat.append(&mut d2ref.clone());
                            }
                            let mut jref = refdata.refs[rsi[oo].jids[m]].to_ascii_vec();
                            if jun {
                                let jend = jflank(&seq, &jref);
                                let mut seq_start = vstart as isize;
                                // probably not exactly right
                                if ex.share[r].annv.len() > 1 {
                                    let q1 = ex.share[r].annv[0].0 + ex.share[r].annv[0].1;
                                    let q2 = ex.share[r].annv[1].0;
                                    seq_start += q2 as isize - q1 as isize;
                                }
                                let seq_end = seq.len() - (jref.len() - jend);
                                // very flaky bug workaround
                                // asserted on BCR=180030 CDR3=CARERDLIWFGPW JALIGN1
                                if seq_start as usize > seq_end {
                                    seq_start = vstart as isize;
                                }
                                seq = seq[seq_start as usize..seq_end].to_vec();
                                jref = jref[0..jend].to_vec();
                            }
                            concat.append(&mut jref.clone());

                            // Make and print alignment.

                            let mut add = String::new();
                            if ex.share[r].left {
                                let mut drefname = drefname.clone();
                                if drefname == "".to_string() {
                                    drefname = "none".to_string();
                                }
                                let rank;
                                if pass == 1 {
                                    rank = "1ST";
                                } else {
                                    rank = "2ND";
                                }
                                add = format!("  •  D = {} = {}", rank, drefname);
                            }
                            let widthx;
                            if !jun {
                                widthx = width;
                            } else {
                                widthx = 1000;
                            }
                            print_vis_align(
                                &seq,
                                &concat,
                                *col,
                                k,
                                &vref,
                                &drefx,
                                &d2ref,
                                &jref,
                                &drefname,
                                ex.share[r].left,
                                &ctl,
                                &mut logx,
                                widthx,
                                &add,
                            );
                        }
                    }
                }
            }
            res.1.push(logx);
        }
    });
    for i in 0..groups.len() {
        for j in 0..groups[i].len() {
            align_out.insert((i, j), results[i].1[j].clone());
        }
    }
    align_out
}
