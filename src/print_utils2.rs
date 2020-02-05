// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::defs::*;
use crate::print_utils1::*;
use crate::types::*;
use amino::*;
use ansi_escape::*;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use itertools::*;
use ndarray::s;
use stats_utils::*;
use std::cmp::{max, min};
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// The following code creates a row in the enclone output table for a clonotype.  Simultaneously
// it generates a row of parseable output.  And it does some other things that are not described
// here.
//
// Awful interface, should work to improve.

pub fn row_fill(
    u: usize,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    mults: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    gex_info: &GexInfo,
    refdata: &RefData,
    varmat: &Vec<Vec<Vec<u8>>>,
    vars_amino: &Vec<Vec<usize>>,
    show_aa: &Vec<Vec<usize>>,
    bads: &mut Vec<bool>,
    gex_low: &mut usize,
    row: &mut Vec<String>,                       // row of human-readable output
    out_data: &mut Vec<HashMap<String, String>>, // row of parseable output
    cx: &mut Vec<Vec<String>>,
    d_all: &mut Vec<Vec<u32>>,
    ind_all: &mut Vec<Vec<u32>>,
    rsi: &ColInfo,
    dref: &Vec<DonorReferenceItem>,
    groups: &HashMap<usize, Vec<usize>>,
) {
    // Redefine some things to reduce dependencies.

    let mat = &rsi.mat;
    let lvars = &ctl.clono_print_opt.lvars;
    macro_rules! speak {
        ($u:expr, $var:expr, $val:expr) => {
            if ctl.parseable_opt.pout.len() > 0 {
                if ctl.parseable_opt.pcols.is_empty()
                    || bin_member(&ctl.parseable_opt.pcols, &$var.to_string())
                {
                    out_data[$u].insert($var.to_string(), $val);
                }
            }
        };
    }
    let pcols_sort = &ctl.parseable_opt.pcols_sort;
    macro_rules! speakc {
        ($u:expr, $col:expr, $var:expr, $val:expr) => {
            if ctl.parseable_opt.pout.len() > 0 && $col + 1 <= ctl.parseable_opt.pchains {
                let varc = format!("{}{}", $var, $col + 1);
                if pcols_sort.is_empty() || bin_member(&pcols_sort, &varc) {
                    out_data[$u].insert(varc, format!("{}", $val));
                }
            }
        };
    }
    let cols = varmat[0].len();

    // Precompute for near and far.

    let mut fp = vec![Vec::<usize>::new(); varmat.len()]; // footprints
    for i in 0..varmat.len() {
        for j in 0..varmat[i].len() {
            if varmat[i][j] != vec![b'-'] {
                fp[i].push(j);
            }
        }
    }

    // Set up lead variable macro.  This is the mechanism for generating
    // both human-readable and parseable output for lead variables.

    macro_rules! lvar {
        ($var:expr, $val:expr) => {
            row.push($val);
            speak!(u, $var.to_string(), $val);
        };
    }

    // Compute dataset indices, gex_med, gex_max, n_gex.

    let mut dataset_indices = Vec::<usize>::new();
    let clonotype_id = exacts[u];
    let ex = &exact_clonotypes[clonotype_id];
    for l in 0..ex.clones.len() {
        dataset_indices.push(ex.clones[l][0].dataset_index);
    }
    unique_sort(&mut dataset_indices);
    let mut lenas = Vec::<String>::new();
    for l in dataset_indices.iter() {
        lenas.push(ctl.sample_info.dataset_id[*l].clone());
    }
    row.push("".to_string()); // row number (#), filled in below
    let mut counts = Vec::<usize>::new();
    let mut n_gex = 0;
    for l in 0..ex.clones.len() {
        let li = ex.clones[l][0].dataset_index;
        let bc = ex.clones[l][0].barcode.clone();
        if !gex_info.gex_barcodes.is_empty() {
            if bin_member(&gex_info.gex_cell_barcodes[li], &bc) {
                n_gex += 1;
            }
            let mut count = 0;
            let p = bin_position(&gex_info.gex_barcodes[li], &bc);
            if p >= 0 {
                let mut raw_count = 0;
                if !ctl.gen_opt.h5 {
                    let row = gex_info.gex_matrices[li].row(p as usize);
                    for j in 0..row.len() {
                        let f = row[j].0;
                        let n = row[j].1;
                        if gex_info.is_gex[li][f] {
                            raw_count += n;
                        }
                    }
                } else {
                    let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                    let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize; // is p+1 OK??
                    let d: Vec<u32> = gex_info.h5_data[li]
                        .as_ref()
                        .unwrap()
                        .as_reader()
                        .read_slice(&s![z1..z2])
                        .unwrap()
                        .to_vec();
                    d_all[l] = d.clone();
                    let ind: Vec<u32> = gex_info.h5_indices[li]
                        .as_ref()
                        .unwrap()
                        .as_reader()
                        .read_slice(&s![z1..z2])
                        .unwrap()
                        .to_vec();
                    ind_all[l] = ind.clone();
                    for j in 0..d.len() {
                        if gex_info.is_gex[li][ind[j] as usize] {
                            raw_count += d[j] as usize;
                        }
                    }
                }
                count = (raw_count as f64 * gex_info.gex_mults[li]).round() as usize;
            }
            counts.push(count);
        }
    }
    counts.sort();
    for n in counts.iter() {
        if *n < 100 {
            *gex_low += 1;
        }
    }
    let (mut gex_median, mut gex_max) = (0, 0);
    if counts.len() > 0 {
        gex_median = counts[counts.len() / 2];
        gex_max = counts[counts.len() - 1];
    }

    // Output lead variable columns.

    for i in 0..lvars.len() {
        if lvars[i].starts_with('g') && lvars[i].after("g").parse::<usize>().is_ok() {
            let d = lvars[i].after("g").force_usize();
            lvar![lvars[i], format!("{}", groups[&d][u] + 1)];
        } else if lvars[i] == "datasets".to_string() {
            lvar![lvars[i], format!("{}", lenas.iter().format(","))];
        } else if lvars[i] == "donors".to_string() {
            let mut donors = Vec::<String>::new();
            for lena in dataset_indices.iter() {
                donors.push(ctl.sample_info.donor_id[ctl.sample_info.donor_index[*lena]].clone());
            }
            unique_sort(&mut donors);
            lvar![lvars[i], format!("{}", donors.iter().format(","))];
        } else if lvars[i] == "ncells".to_string() {
            lvar![lvars[i], format!("{}", mults[u])];
        } else if lvars[i].starts_with("n_") && lvars[i] != "n_gex".to_string() {
            let name = lvars[i].after("n_");
            let indices = ctl.sample_info.name_list[name].clone();
            let mut count = 0;
            for j in 0..ex.clones.len() {
                if bin_member(&indices, &ex.clones[j][0].dataset_index) {
                    count += 1;
                }
            }
            lvar![lvars[i], format!("{}", count)];
        } else if lvars[i] == "near".to_string() {
            let mut dist = 1_000_000;
            for i2 in 0..varmat.len() {
                if i2 == u || fp[i2] != fp[u] {
                    continue;
                }
                let mut d = 0;
                for c in fp[u].iter() {
                    for j in 0..varmat[u][*c].len() {
                        if varmat[u][*c][j] != varmat[i2][*c][j] {
                            d += 1;
                        }
                    }
                }
                dist = min(dist, d);
            }
            if dist == 1_000_000 {
                lvar![lvars[i], "-".to_string()];
            } else {
                lvar![lvars[i], format!("{}", dist)];
            }
        } else if lvars[i] == "far".to_string() {
            let mut dist = -1 as isize;
            for i2 in 0..varmat.len() {
                if i2 == u || fp[i2] != fp[u] {
                    continue;
                }
                let mut d = 0 as isize;
                for c in fp[u].iter() {
                    for j in 0..varmat[u][*c].len() {
                        if varmat[u][*c][j] != varmat[i2][*c][j] {
                            d += 1;
                        }
                    }
                }
                dist = max(dist, d);
            }
            if dist == -1 as isize {
                lvar![lvars[i], "-".to_string()];
            } else {
                lvar![lvars[i], format!("{}", dist)];
            }
        } else if lvars[i] == "gex_med".to_string() {
            lvar![lvars[i], format!("{}", gex_median)];
        } else if lvars[i] == "n_gex".to_string() {
            lvar![lvars[i], format!("{}", n_gex)];
        } else if lvars[i] == "gex_max".to_string() {
            lvar![lvars[i], format!("{}", gex_max)];
        } else if lvars[i] == "ext".to_string() {
            let mut exts = Vec::<String>::new();
            for l in 0..ex.clones.len() {
                let li = ctl.sample_info.dataset_id[ex.clones[l][0].dataset_index].clone();
                let bc = ex.clones[l][0].barcode.clone();
                if ctl.gen_opt.extc.contains_key(&(li.clone(), bc.clone())) {
                    exts.push(ctl.gen_opt.extc[&(li, bc)].clone());
                }
            }
            exts.sort();
            let mut s = String::new();
            let mut j = 0;
            while j < exts.len() {
                let k = next_diff(&exts, j);
                if j > 0 {
                    s += ",";
                }
                s += &format!(
                    "{}[{}/{}]",
                    exts[j],
                    k - j,
                    ctl.gen_opt.extn[&exts[j].clone()]
                );
                j = k;
            }
            lvar![lvars[i], s.clone()];
        } else {
            let mut count = 0;
            for l in 0..ex.clones.len() {
                let li = ex.clones[l][0].dataset_index;
                let bc = ex.clones[l][0].barcode.clone();
                let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                if p >= 0 {
                    if gex_info.feature_id[li].contains_key(&lvars[i]) {
                        let fid = gex_info.feature_id[li][&lvars[i]];
                        let mut raw_count = 0 as f64;
                        if !ctl.gen_opt.h5 {
                            raw_count = gex_info.gex_matrices[li].value(p as usize, fid) as f64;
                        } else {
                            for j in 0..d_all[l].len() {
                                if ind_all[l][j] == fid as u32 {
                                    raw_count = d_all[l][j] as f64;
                                    break;
                                }
                            }
                        }
                        let mult: f64;
                        if lvars[i].ends_with("_g") {
                            mult = gex_info.gex_mults[li];
                        } else {
                            mult = gex_info.fb_mults[li];
                        }
                        count += (raw_count * mult).round() as usize;
                    }
                }
            }
            lvar![lvars[i], format!("{}", count / ex.ncells())];
        }
    }

    // Get the relevant barcodes.

    let mut bli = Vec::<(String, usize, usize)>::new();
    for l in 0..ex.clones.len() {
        bli.push((
            ex.clones[l][0].barcode.clone(),
            ex.clones[l][0].dataset_index,
            l,
        ));
    }
    bli.sort();

    // Traverse the chains.

    for col in 0..cols {
        let mid = mat[col][u];
        if mid.is_none() {
            continue;
        }
        let mid = mid.unwrap();
        let ex = &exact_clonotypes[clonotype_id];
        let seq_amino = rsi.seqss_amino[col][u].clone();

        // Get UMI and read stats.

        let mut numis = Vec::<usize>::new();
        let mut nreads = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
            nreads.push(ex.clones[j][mid].read_count);
        }
        numis.sort();
        let median_numis = numis[numis.len() / 2];
        let utot: usize = numis.iter().sum();
        let umax = *numis.iter().max().unwrap();
        nreads.sort();
        let median_nreads = nreads[nreads.len() / 2];

        // Set up chain variable macro.  This is the mechanism for generating
        // both human-readable and parseable output for chain variables.

        macro_rules! cvar {
            ($i: expr, $var:expr, $val:expr) => {
                cx[col][$i] = $val.clone();
                speakc!(u, col, $var, $val);
            };
        }

        // Speak quality score column entries.

        if ctl.parseable_opt.pout.len() > 0 && col + 1 <= ctl.parseable_opt.pchains {
            for i in 0..pcols_sort.len() {
                if pcols_sort[i].starts_with('q')
                    && pcols_sort[i].ends_with(&format!("_{}", col + 1))
                {
                    let n = pcols_sort[i].after("q").rev_before("_").force_usize();
                    if n < ex.share[mid].seq.len() {
                        let mut quals = Vec::<u8>::new();
                        for j in 0..ex.clones.len() {
                            quals.push(ex.clones[j][mid].quals[n]);
                        }
                        let q = format!("{}", quals.iter().format(","));
                        out_data[u].insert(pcols_sort[i].clone(), q);
                    }
                }
            }
        }

        // Speak some other column entries.

        speakc!(u, col, "vj_seq", stringme(&ex.share[mid].seq));
        speakc!(u, col, "seq", stringme(&ex.share[mid].full_seq));
        speakc!(u, col, "v_start", ex.share[mid].v_start);
        let cid = ex.share[mid].c_ref_id;
        if cid.is_some() {
            let cid = cid.unwrap();
            speakc!(u, col, "const_id", refdata.id[cid]);
        }
        let uid = ex.share[mid].u_ref_id;
        if uid.is_some() {
            let uid = uid.unwrap();
            speakc!(u, col, "utr_id", refdata.id[uid]);
            speakc!(u, col, "utr_name", refdata.name[uid]);
        }
        speakc!(u, col, "cdr3_start", ex.share[mid].cdr3_start);
        speakc!(u, col, "cdr3_aa", ex.share[mid].cdr3_aa);
        let mut vv = Vec::<usize>::new();
        for x in vars_amino[col].iter() {
            vv.push(*x / 3);
        }
        unique_sort(&mut vv);
        let mut varaa = Vec::<u8>::new();
        for p in vv.iter() {
            // what does it mean if this fails?
            if 3 * p + 3 <= seq_amino.len() {
                if seq_amino[3 * p..3 * p + 3].to_vec() == b"---".to_vec() {
                    varaa.push(b'-');
                } else {
                    varaa.push(codon_to_aa(&seq_amino[3 * p..3 * p + 3]));
                }
            }
        }
        speakc!(u, col, "var_aa", strme(&varaa));

        // Create column entry.

        for j in 0..rsi.cvars[col].len() {
            if rsi.cvars[col][j] == "amino".to_string() {
                let cs = rsi.cdr3_starts[col] / 3;
                let n = rsi.cdr3_lens[col];
                for k in 0..show_aa[col].len() {
                    let p = show_aa[col][k];
                    if k > 0 && p == cs {
                        cx[col][j] += " ";
                    }
                    if 3 * p + 3 <= seq_amino.len()
                        && seq_amino[3 * p..3 * p + 3].to_vec() == b"---".to_vec()
                    {
                        cx[col][j] += "-";
                    } else if 3 * p + 3 > seq_amino.len()
                        || seq_amino[3 * p..3 * p + 3].contains(&b'-')
                    {
                        cx[col][j] += "*";
                    } else {
                        let mut log = Vec::<u8>::new();
                        emit_codon_color_escape(&seq_amino[3 * p..3 * p + 3], &mut log);
                        let aa = codon_to_aa(&seq_amino[3 * p..3 * p + 3]);
                        log.push(aa);
                        emit_end_escape(&mut log);
                        cx[col][j] += strme(&log);
                    }
                    if k < show_aa[col].len() - 1 && p == cs + n - 1 {
                        cx[col][j] += " ";
                    }
                }
            } else if rsi.cvars[col][j] == "comp".to_string() {
                let mut comp = 1000000;
                let td = &ex.share[mid];
                let tig = &td.seq;
                let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
                let mut aligner = Aligner::new(-6, -1, &score);
                let mut z = 1;
                if ex.share[mid].left {
                    z = refdata.refs.len();
                }

                // Go through every reference record.  Nuts.

                for d in 0..z {
                    if ex.share[mid].left && !refdata.is_d(d) {
                        continue;
                    }
                    let mut concat = Vec::<u8>::new();
                    let mut x = refdata.refs[rsi.vids[col]].to_ascii_vec();
                    if rsi.vpids[col].is_none() {
                        concat.append(&mut x);
                    } else {
                        let mut y = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
                        concat.append(&mut y);
                    }
                    if ex.share[mid].left {
                        let mut x = refdata.refs[d].to_ascii_vec();
                        concat.append(&mut x);
                    }
                    let mut x = refdata.refs[rsi.jids[col]].to_ascii_vec();
                    concat.append(&mut x);
                    let al = aligner.semiglobal(&tig, &concat);
                    let mut m = 0;
                    let mut pos = al.xstart;
                    let mut count = 0;
                    let start = td.cdr3_start;
                    let stop = td.cdr3_start + 3 * td.cdr3_aa.len();
                    while m < al.operations.len() {
                        let n = next_diff(&al.operations, m);
                        match al.operations[m] {
                            Match => {
                                pos += 1;
                            }
                            Subst => {
                                if pos >= start && pos < stop {
                                    count += 1;
                                }
                                pos += 1;
                            }
                            Del => {
                                if pos >= start && pos < stop {
                                    count += 1;
                                }
                                pos += n - m;
                                m = n - 1;
                            }
                            Ins => {
                                if pos >= start && pos < stop {
                                    count += 1;
                                }
                                m = n - 1;
                            }
                            _ => {}
                        };
                        m += 1;
                    }
                    comp = min(comp, count);
                }
                cvar![j, rsi.cvars[col][j], format!("{}", comp)];
            } else if rsi.cvars[col][j] == "cdr3_dna".to_string() {
                cvar![j, rsi.cvars[col][j], ex.share[mid].cdr3_dna.clone()];
            } else if rsi.cvars[col][j] == "ulen".to_string() {
                cvar![j, rsi.cvars[col][j], format!("{}", ex.share[mid].v_start)];
            } else if rsi.cvars[col][j] == "clen".to_string() {
                cvar![
                    j,
                    rsi.cvars[col][j],
                    format!("{}", ex.share[mid].full_seq.len() - ex.share[mid].j_stop)
                ];
            } else if rsi.cvars[col][j].starts_with("ndiff") {
                let u0 = rsi.cvars[col][j].after("ndiff").force_usize() - 1;
                if u0 < exacts.len() && mat[col][u0].is_some() && mat[col][u].is_some() {
                    let m0 = mat[col][u0].unwrap();
                    let m = mat[col][u].unwrap();
                    let mut ndiff = 0;
                    let ex0 = &exact_clonotypes[exacts[u0]];
                    let ex = &exact_clonotypes[exacts[u]];
                    for p in 0..ex0.share[m0].seq_del.len() {
                        if ex0.share[m0].seq_del[p] != ex.share[m].seq_del[p] {
                            ndiff += 1;
                        }
                    }
                    cvar![j, rsi.cvars[col][j], format!("{}", ndiff)];
                } else {
                    cvar![j, rsi.cvars[col][j], "_".to_string()];
                }
            } else if rsi.cvars[col][j] == "cdiff".to_string() {
                let cstart = ex.share[mid].j_stop;
                let clen = ex.share[mid].full_seq.len() - cstart;
                let cid = ex.share[mid].c_ref_id;
                let mut cdiff = String::new();
                let mut ndiffs = 0;
                if cid.is_some() {
                    let r = &refdata.refs[cid.unwrap()];
                    let mut extra = 0;
                    if clen > r.len() {
                        extra = clen - r.len();
                    }
                    for i in 0..min(clen, r.len()) {
                        let tb = ex.share[mid].full_seq[cstart + i];
                        let rb = r.to_ascii_vec()[i];
                        if tb != rb {
                            ndiffs += 1;
                            if ndiffs <= 5 {
                                cdiff += &format!("{}{}", i, tb as char);
                            }
                        }
                    }
                    if ndiffs > 5 {
                        cdiff += "...";
                    }
                    if extra > 0 {
                        cdiff += &format!("%{}", extra);
                    }
                } else if clen > 0 {
                    cdiff = format!("%{}", clen);
                }
                cvar![j, rsi.cvars[col][j], cdiff];
            } else if rsi.cvars[col][j] == "udiff".to_string() {
                let ulen = ex.share[mid].v_start;
                let uid = ex.share[mid].u_ref_id;
                let mut udiff = String::new();
                let mut ndiffs = 0;
                if uid.is_some() {
                    let r = &refdata.refs[uid.unwrap()];
                    let mut extra = 0;
                    if ulen > r.len() {
                        extra = ulen - r.len();
                    }
                    for i in 0..ulen {
                        let mut rpos = i;
                        if ulen < r.len() {
                            rpos += r.len() - ulen;
                        } else {
                            if i + r.len() < ulen {
                                continue;
                            }
                            rpos -= ulen - r.len();
                        }
                        let tb = ex.share[mid].full_seq[i];
                        let rb = r.to_ascii_vec()[rpos];
                        if tb != rb {
                            ndiffs += 1;
                            if ndiffs <= 5 {
                                udiff += &format!("{}{}", rpos, tb as char);
                            }
                        }
                    }
                    if ndiffs > 5 {
                        udiff += "...";
                    }
                    if extra > 0 {
                        udiff += &format!("%{}", extra);
                    }
                } else if ulen > 0 {
                    udiff = format!("%{}", ulen);
                }
                cvar![j, rsi.cvars[col][j], udiff];
            } else if rsi.cvars[col][j] == "notes".to_string() {
                cvar![j, rsi.cvars[col][j], ex.share[mid].vs_notesx.clone()];
            } else if rsi.cvars[col][j] == "var".to_string() {
                cvar![j, rsi.cvars[col][j], stringme(&varmat[u][col])];
            } else if rsi.cvars[col][j] == "umed".to_string() {
                cvar![j, rsi.cvars[col][j], format!("{}", median_numis)];
            } else if rsi.cvars[col][j] == "umax".to_string() {
                cvar![j, rsi.cvars[col][j], format!("{}", umax)];
            } else if rsi.cvars[col][j] == "utot".to_string() {
                cvar![j, rsi.cvars[col][j], format!("{}", utot)];
            } else if rsi.cvars[col][j] == "rmed".to_string() {
                cvar![j, rsi.cvars[col][j], format!("{}", median_nreads)];
            } else if rsi.cvars[col][j] == "const".to_string() {
                let mut constx = Vec::<String>::new();
                let cid = ex.share[mid].c_ref_id;
                if cid.is_some() {
                    constx.push(refdata.name[cid.unwrap()].clone());
                } else {
                    constx.push("?".to_string());
                }
                unique_sort(&mut constx);
                // This is overcomplicated because there is now at most one
                // const entry per exact subclonotype.
                cvar![
                    j,
                    rsi.cvars[col][j],
                    format!("{}", constx.iter().format(","))
                ];

            // Compute potential whitelist contamination percent.  And filter.
            // This is an undocumented option.
            } else if rsi.cvars[col][j] == "white".to_string() || ctl.clono_filt_opt.whitef {
                let mut bch = vec![Vec::<(usize, String, usize, usize)>::new(); 2];
                for l in 0..ex.clones.len() {
                    let li = ex.clones[l][0].dataset_index;
                    let bc = &ex.clones[l][0].barcode;
                    let mut numi = 0;
                    for j in 0..ex.clones[l].len() {
                        numi += ex.clones[l][j].umi_count;
                    }
                    bch[0].push((li, bc[0..8].to_string(), numi, l));
                    bch[1].push((li, bc[8..16].to_string(), numi, l));
                }
                let mut junk = 0;
                let mut bad = vec![false; ex.clones.len()];
                for l in 0..2 {
                    bch[l].sort();
                    let mut m = 0;
                    while m < bch[l].len() {
                        let n = next_diff12_4(&bch[l], m as i32) as usize;
                        for u1 in m..n {
                            for u2 in m..n {
                                if bch[l][u1].2 >= 10 * bch[l][u2].2 {
                                    bad[bch[l][u2].3] = true;
                                }
                            }
                        }
                        m = n;
                    }
                }
                for u in 0..bad.len() {
                    if bad[u] {
                        junk += 1;
                    }
                }
                // Don't look at very large clones because of course they
                // show overlap.
                /* // BROKEN AND WAS UGLY ANYWAY
                const MAX_WHITELIST_CLONE: usize = 100;
                if ex.clones.len() <= MAX_WHITELIST_CLONE {
                    res.3 += junk;
                    res.4 += ex.clones.len();
                }
                */
                let junk_rate = percent_ratio(junk, ex.clones.len());
                if rsi.cvars[col][j] == "white".to_string() {
                    cx[col][j] = format!("{:.1}", junk_rate);
                }
                // WRONG!  THIS IS SUPPOSED TO BE EXECUTED ON PASS 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if ctl.clono_filt_opt.whitef && junk_rate == 0.0
                /* && pass == 1 */
                {
                    bads[u] = true;
                }
            }
        }
    }
}
