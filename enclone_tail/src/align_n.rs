// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Parallized precompute for ALIGN<n> and JUN_ALIGN<n>.

use align_tools::*;
use ansi_escape::*;
use enclone_core::align_to_vdj_ref::*;
use enclone_core::defs::*;
use enclone_core::opt_d::*;
use enclone_proto::types::*;
use io_utils::*;
use rayon::prelude::*;
use std::cmp::min;
use std::collections::HashMap;
use std::io::Write;
use string_utils::*;
use vdj_ann::refx::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn print_vis_align(
    seq: &[u8],
    concat: &[u8],
    col: usize,
    k: usize,
    vref: &[u8],
    dref: &[u8],
    jref: &[u8],
    left: bool, // will be ex.share[r].left
    ctl: &EncloneControl,
    logx: &mut Vec<u8>,
    width: usize,
) {
    // Make alignment.

    let al = align_to_vdj_ref(&seq, &vref, &dref, &jref);

    // Make visual alignment.

    let mut vis = vis_align(&seq, &concat, &al, width);

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
        let jcolor = 5;
        let mut vdj_bytes = Vec::<u8>::new();
        print_color(vcolor, &mut vdj_bytes);
        vdj_bytes.push(b'V');
        if left {
            print_color(dcolor, &mut vdj_bytes);
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
                } else if left && pos < vref.len() + dref.len() {
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
                                if pos == vref.len() + dref.len() {
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
        "\nALIGNMENT OF CHAIN {} FOR EXACT SUBCLONOTYPE {} TO CONCATENATED {} REFERENCE\n{}",
        col,
        k + 1,
        vdj,
        vis,
    );
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// There is terrible duplication between the functions below.

pub fn align_n(
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    dref: &Vec<DonorReferenceItem>,
    groups: &Vec<Vec<(i32, String)>>,
    width: usize,
) -> HashMap<(usize, usize), Vec<u8>> {
    let mut align_out = HashMap::<(usize, usize), Vec<u8>>::new();
    let mut results = Vec::<(usize, Vec<Vec<u8>>)>::new();
    for i in 0..groups.len() {
        results.push((i, Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut o = Vec::<i32>::new();
        for j in 0..groups[i].len() {
            o.push(groups[i][j].0);
        }
        for j in 0..o.len() {
            let mut logx = Vec::<u8>::new();
            let oo = o[j] as usize;
            for col in ctl.gen_opt.chains_to_align.iter() {
                let m = col - 1;
                for k in 0..exacts[oo].len() {
                    let ex = &exact_clonotypes[exacts[oo][k]];
                    if m < rsi[oo].mat.len() && rsi[oo].mat[m][k].is_some() {
                        let r = rsi[oo].mat[m][k].unwrap();
                        let seq = &ex.share[r].seq;
                        let mut concat = Vec::<u8>::new();
                        let mut vref = refdata.refs[rsi[oo].vids[m]].to_ascii_vec();
                        if rsi[oo].vpids[m].is_none() {
                        } else {
                            vref = dref[rsi[oo].vpids[m].unwrap()].nt_sequence.clone();
                        }
                        concat.append(&mut vref.clone());
                        let mut drefx = Vec::<u8>::new();
                        if ex.share[r].left {
                            let mut opt = None;
                            let mut opt2 = None;
                            let mut delta = 0;
                            opt_d(
                                &ex, m, k, &rsi[oo], &refdata, &dref, &mut opt, &mut opt2,
                                &mut delta,
                            );
                            if opt.is_some() {
                                let opt = opt.unwrap();
                                drefx = refdata.refs[opt].to_ascii_vec();
                                concat.append(&mut drefx.clone());
                            }
                        }
                        let jref = refdata.refs[rsi[oo].jids[m]].to_ascii_vec();
                        concat.append(&mut jref.clone());

                        // Make and print alignment.

                        print_vis_align(
                            &seq,
                            &concat,
                            *col,
                            k,
                            &vref,
                            &drefx,
                            &jref,
                            ex.share[r].left,
                            &ctl,
                            &mut logx,
                            width,
                        );
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn align2_n(
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    dref: &Vec<DonorReferenceItem>,
    groups: &Vec<Vec<(i32, String)>>,
    width: usize,
) -> HashMap<(usize, usize), Vec<u8>> {
    let mut align_out = HashMap::<(usize, usize), Vec<u8>>::new();
    let mut results = Vec::<(usize, Vec<Vec<u8>>)>::new();
    for i in 0..groups.len() {
        results.push((i, Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut o = Vec::<i32>::new();
        for j in 0..groups[i].len() {
            o.push(groups[i][j].0);
        }
        for j in 0..o.len() {
            let mut logx = Vec::<u8>::new();
            let oo = o[j] as usize;
            for col in ctl.gen_opt.chains_to_align2.iter() {
                let m = col - 1;
                for k in 0..exacts[oo].len() {
                    let ex = &exact_clonotypes[exacts[oo][k]];
                    if m < rsi[oo].mat.len() && rsi[oo].mat[m][k].is_some() {
                        let r = rsi[oo].mat[m][k].unwrap();
                        let seq = &ex.share[r].seq;
                        let mut concat = Vec::<u8>::new();
                        let mut vref = refdata.refs[rsi[oo].vids[m]].to_ascii_vec();
                        if rsi[oo].vpids[m].is_none() {
                        } else {
                            vref = dref[rsi[oo].vpids[m].unwrap()].nt_sequence.clone();
                        }
                        concat.append(&mut vref.clone());
                        let mut drefx = Vec::<u8>::new();
                        if ex.share[r].left {
                            let mut opt = None;
                            let mut opt2 = None;
                            let mut delta = 0;
                            opt_d(
                                &ex, m, k, &rsi[oo], &refdata, &dref, &mut opt, &mut opt2,
                                &mut delta,
                            );
                            if opt2.is_some() {
                                let opt2 = opt2.unwrap();
                                drefx = refdata.refs[opt2].to_ascii_vec();
                                concat.append(&mut drefx.clone());
                            }
                        }
                        let jref = refdata.refs[rsi[oo].jids[m]].to_ascii_vec();
                        concat.append(&mut jref.clone());

                        // Make and print alignment.

                        print_vis_align(
                            &seq,
                            &concat,
                            *col,
                            k,
                            &vref,
                            &drefx,
                            &jref,
                            ex.share[r].left,
                            &ctl,
                            &mut logx,
                            width,
                        );
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn jun_align_n(
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    dref: &Vec<DonorReferenceItem>,
    groups: &Vec<Vec<(i32, String)>>,
    width: usize,
) -> HashMap<(usize, usize), Vec<u8>> {
    let mut align_out = HashMap::<(usize, usize), Vec<u8>>::new();
    let mut results = Vec::<(usize, Vec<Vec<u8>>)>::new();
    for i in 0..groups.len() {
        results.push((i, Vec::new()));
    }
    const LFLANK: usize = 15;
    const RFLANK: usize = 35;
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut o = Vec::<i32>::new();
        for j in 0..groups[i].len() {
            o.push(groups[i][j].0);
        }
        for j in 0..o.len() {
            let mut logx = Vec::<u8>::new();
            let oo = o[j] as usize;
            // This line differs from the previous function.
            for col in ctl.gen_opt.chains_to_jun_align.iter() {
                let m = col - 1;
                for k in 0..exacts[oo].len() {
                    let ex = &exact_clonotypes[exacts[oo][k]];
                    if m < rsi[oo].mat.len() && rsi[oo].mat[m][k].is_some() {
                        let r = rsi[oo].mat[m][k].unwrap();
                        let seq = &ex.share[r].seq;
                        let mut concat = Vec::<u8>::new();
                        let mut vref = refdata.refs[rsi[oo].vids[m]].to_ascii_vec();
                        if rsi[oo].vpids[m].is_none() {
                        } else {
                            vref = dref[rsi[oo].vpids[m].unwrap()].nt_sequence.clone();
                        }

                        // Different from previous function:

                        let mut vstart = 0;
                        if vref.len() >= LFLANK {
                            vstart = vref.len() - LFLANK;
                        }
                        vref = vref[vstart..vref.len()].to_vec();

                        concat.append(&mut vref.clone());
                        let mut drefx = Vec::<u8>::new();
                        if ex.share[r].left {
                            let mut opt = None;
                            let mut opt2 = None;
                            let mut delta = 0;
                            opt_d(
                                &ex, m, k, &rsi[oo], &refdata, &dref, &mut opt, &mut opt2,
                                &mut delta,
                            );
                            if opt.is_some() {
                                let opt = opt.unwrap();
                                drefx = refdata.refs[opt].to_ascii_vec();
                                concat.append(&mut drefx.clone());
                            }
                        }

                        // Different from previous function:

                        let jref = refdata.refs[rsi[oo].jids[m]].to_ascii_vec();
                        let jend = min(RFLANK, jref.len());
                        let mut seq_start = vstart as isize;
                        // probably not exactly right
                        if ex.share[r].annv.len() > 1 {
                            let q1 = ex.share[r].annv[0].0 + ex.share[r].annv[0].1;
                            let q2 = ex.share[r].annv[1].0;

                            seq_start += q2 as isize - q1 as isize;
                        }
                        let seq_end = seq.len() - (jref.len() - jend);
                        let seq = seq[seq_start as usize..seq_end].to_vec();
                        let jref = jref[0..jend].to_vec();
                        concat.append(&mut jref.clone());

                        // Make and print alignment.

                        print_vis_align(
                            &seq,
                            &concat,
                            *col,
                            k,
                            &vref,
                            &drefx,
                            &jref,
                            ex.share[r].left,
                            &ctl,
                            &mut logx,
                            width,
                        );
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn jun_align2_n(
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    dref: &Vec<DonorReferenceItem>,
    groups: &Vec<Vec<(i32, String)>>,
    width: usize,
) -> HashMap<(usize, usize), Vec<u8>> {
    let mut align_out = HashMap::<(usize, usize), Vec<u8>>::new();
    let mut results = Vec::<(usize, Vec<Vec<u8>>)>::new();
    for i in 0..groups.len() {
        results.push((i, Vec::new()));
    }
    const LFLANK: usize = 15;
    const RFLANK: usize = 35;
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let mut o = Vec::<i32>::new();
        for j in 0..groups[i].len() {
            o.push(groups[i][j].0);
        }
        for j in 0..o.len() {
            let mut logx = Vec::<u8>::new();
            let oo = o[j] as usize;
            // This line differs from the previous function.
            for col in ctl.gen_opt.chains_to_jun_align2.iter() {
                let m = col - 1;
                for k in 0..exacts[oo].len() {
                    let ex = &exact_clonotypes[exacts[oo][k]];
                    if m < rsi[oo].mat.len() && rsi[oo].mat[m][k].is_some() {
                        let r = rsi[oo].mat[m][k].unwrap();
                        let seq = &ex.share[r].seq;
                        let mut concat = Vec::<u8>::new();
                        let mut vref = refdata.refs[rsi[oo].vids[m]].to_ascii_vec();
                        if rsi[oo].vpids[m].is_none() {
                        } else {
                            vref = dref[rsi[oo].vpids[m].unwrap()].nt_sequence.clone();
                        }

                        // Different from previous function:

                        let mut vstart = 0;
                        if vref.len() >= LFLANK {
                            vstart = vref.len() - LFLANK;
                        }
                        vref = vref[vstart..vref.len()].to_vec();

                        concat.append(&mut vref.clone());
                        let mut drefx = Vec::<u8>::new();
                        if ex.share[r].left {
                            let mut opt = None;
                            let mut opt2 = None;
                            let mut delta = 0;
                            opt_d(
                                &ex, m, k, &rsi[oo], &refdata, &dref, &mut opt, &mut opt2,
                                &mut delta,
                            );
                            if opt2.is_some() {
                                let opt2 = opt2.unwrap();
                                drefx = refdata.refs[opt2].to_ascii_vec();
                                concat.append(&mut drefx.clone());
                            }
                        }

                        // Different from previous function:

                        let jref = refdata.refs[rsi[oo].jids[m]].to_ascii_vec();
                        let jend = min(RFLANK, jref.len());
                        let mut seq_start = vstart as isize;
                        // probably not exactly right
                        if ex.share[r].annv.len() > 1 {
                            let q1 = ex.share[r].annv[0].0 + ex.share[r].annv[0].1;
                            let q2 = ex.share[r].annv[1].0;

                            seq_start += q2 as isize - q1 as isize;
                        }
                        let seq_end = seq.len() - (jref.len() - jend);
                        let seq = seq[seq_start as usize..seq_end].to_vec();
                        let jref = jref[0..jend].to_vec();
                        concat.append(&mut jref.clone());

                        // Make and print alignment.

                        print_vis_align(
                            &seq,
                            &concat,
                            *col,
                            k,
                            &vref,
                            &drefx,
                            &jref,
                            ex.share[r].left,
                            &ctl,
                            &mut logx,
                            width,
                        );
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
