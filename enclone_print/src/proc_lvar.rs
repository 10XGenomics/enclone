// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_lvar,
// plus a small helper function get_gex_matrix_entry.

use amino::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use itertools::*;
use regex::Regex;
use std::cmp::{max, min};
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

fn median_f64(x: &[f64]) -> f64 {
    let h = x.len() / 2;
    if x.len() % 2 == 1 {
        x[h]
    } else {
        (x[h - 1] + x[h]) / 2.0
    }
}

pub fn get_gex_matrix_entry(
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    fid: usize,
    d_all: &Vec<Vec<u32>>,
    ind_all: &Vec<Vec<u32>>,
    li: usize,
    l: usize,
    p: usize,
    y: &str,
) -> f64 {
    let mut raw_count = 0 as f64;
    if gex_info.gex_matrices[li].initialized() {
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
    if y.ends_with("_g") {
        mult = gex_info.gex_mults[li];
    } else {
        mult = gex_info.fb_mults[li];
    }
    if !ctl.gen_opt.full_counts {
        raw_count *= mult;
    }
    raw_count
}

pub fn proc_lvar(
    i: usize,
    x: &String,
    pass: usize,
    u: usize,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    mults: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    gex_info: &GexInfo,
    refdata: &RefData,
    varmat: &Vec<Vec<Vec<u8>>>,
    fp: &Vec<Vec<usize>>,
    row: &mut Vec<String>,
    out_data: &mut Vec<HashMap<String, String>>,
    d_all: &mut Vec<Vec<u32>>,
    ind_all: &mut Vec<Vec<u32>>,
    rsi: &ColInfo,
    dref: &Vec<DonorReferenceItem>,
    groups: &HashMap<usize, Vec<usize>>,
    stats: &mut Vec<(String, Vec<f64>)>,
    vdj_cells: &Vec<Vec<String>>,
    n_vdj_gex: &Vec<usize>,
    nd_fields: &Vec<String>,
    lvars: &Vec<String>,
    lenas: &Vec<String>,
    alt_bcs: &Vec<String>,
    n_gex: usize,
    n_gexs: &Vec<usize>,
    gex_min: usize,
    gex_max: usize,
    gex_mean: f64,
    gex_sum: f64,
    gex_median: usize,
    count_unsorted: &Vec<usize>,
    entropy: f64,
    entropies_unsorted: &Vec<f64>,
    fcounts: &Vec<f64>,
) {
    let mat = &rsi.mat;
    let clonotype_id = exacts[u];
    let ex = &exact_clonotypes[clonotype_id];
    let cols = varmat[0].len();

    // Set up speak macro.

    macro_rules! speak {
        ($u:expr, $var:expr, $val:expr) => {
            if pass == 2 && ctl.parseable_opt.pout.len() > 0 {
                let mut v = $var.to_string();
                v = v.replace("_Σ", "_sum");
                v = v.replace("_μ", "_mean");
                if ctl.parseable_opt.pcols.is_empty()
                    || bin_member(&ctl.parseable_opt.pcols_sortx, &v)
                {
                    out_data[$u].insert(v, $val);
                }
            }
        };
    }

    // Set up lead variable macro.  This is the mechanism for generating
    // both human-readable and parseable output for lead variables.

    macro_rules! lvar {
        ($i: expr, $var:expr, $val:expr) => {
            if $i < lvars.len() {
                row.push($val)
            }
            if pass == 2 {
                speak!(u, $var.to_string(), $val);
            }
        };
    }

    // Proceed.

    if x.starts_with('g') && x.after("g").parse::<usize>().is_ok() {
        let d = x.after("g").force_usize();
        lvar![i, x, format!("{}", groups[&d][u] + 1)];
    } else if x == "origins" {
        let mut origins = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            if ex.clones[j][0].origin_index.is_some() {
                origins.push(
                    ctl.origin_info.origin_list[ex.clones[j][0].origin_index.unwrap()].clone(),
                );
            } else {
                origins.push("?".to_string());
            }
        }
        unique_sort(&mut origins);
        lvar![i, x, format!("{}", origins.iter().format(","))];
    } else if x == "datasets" {
        lvar![i, x, format!("{}", lenas.iter().format(","))];
    } else if x == "donors" {
        let mut donors = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            if ex.clones[j][0].donor_index.is_some() {
                donors
                    .push(ctl.origin_info.donor_list[ex.clones[j][0].donor_index.unwrap()].clone());
            } else {
                donors.push("?".to_string());
            }
        }
        unique_sort(&mut donors);
        lvar![i, x, format!("{}", donors.iter().format(","))];
    } else if x == "n" {
        lvar![i, x, format!("{}", mults[u])];
        let counts = vec![1.0; mults[u]];
        stats.push((x.to_string(), counts));
    } else if x == "clust" {
        let mut clust = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            let mut cid = 0;
            let bc = &ex.clones[j][0].barcode;
            let li = ex.clones[j][0].dataset_index;
            if gex_info.cluster[li].contains_key(&bc.clone()) {
                cid = gex_info.cluster[li][&bc.clone()];
            }
            clust.push(cid);
        }
        clust.sort();
        lvar![i, x, format!("{}", abbrev_list(&clust))];
    } else if x == "n_other" {
        let mut n = 0;
        for j in 0..ex.clones.len() {
            let mut found = false;
            let di = ex.clones[j][0].dataset_index;
            let f = format!("n_{}", ctl.origin_info.dataset_id[di]);
            for i in 0..nd_fields.len() {
                if f == nd_fields[i] {
                    found = true;
                }
            }
            if !found {
                n += 1;
            }
        }
        lvar![i, x, format!("{}", n)];
    } else if x == "n_b" {
        let mut n_b = 0;
        for j in 0..ex.clones.len() {
            let bc = &ex.clones[j][0].barcode;
            let li = ex.clones[j][0].dataset_index;
            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                    n_b += 1;
                }
            }
        }
        lvar![i, x, format!("{}", n_b)];
    } else if x == "type" {
        let mut cell_types = Vec::<String>::new();
        /*
        for j in 0..ex.clones.len() {
            let mut cell_type = "".to_string();
            let bc = &ex.clones[j][0].barcode;
            let li = ex.clones[j][0].dataset_index;
            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                cell_type = gex_info.cell_type[li][&bc.clone()].clone();
            }
            cell_types.push(cell_type);
        }
        */
        cell_types.sort();
        lvar![i, x, format!("{}", abbrev_list(&cell_types))];
    } else if x == "filter" {
        lvar![i, x, String::new()];
    } else if x == "mark" {
        let mut n = 0;
        for j in 0..ex.clones.len() {
            if ex.clones[j][0].marked {
                n += 1;
            }
        }
        lvar![i, x, format!("{}", n)];
    } else if x == "inkt" {
        let mut s = String::new();
        let alpha_g = ex.share[0].inkt_alpha_chain_gene_match;
        let alpha_j = ex.share[0].inkt_alpha_chain_junction_match;
        let beta_g = ex.share[0].inkt_beta_chain_gene_match;
        let beta_j = ex.share[0].inkt_beta_chain_junction_match;
        if alpha_g || alpha_j {
            s += "𝝰";
            if alpha_g {
                s += "g";
            }
            if alpha_j {
                s += "j";
            }
        }
        if beta_g || beta_j {
            s += "𝝱";
            if beta_g {
                s += "g";
            }
            if beta_j {
                s += "j";
            }
        }
        lvar![i, x, s.clone()];
    } else if x == "mait" {
        let mut s = String::new();
        let alpha_g = ex.share[0].mait_alpha_chain_gene_match;
        let alpha_j = ex.share[0].mait_alpha_chain_junction_match;
        let beta_g = ex.share[0].mait_beta_chain_gene_match;
        let beta_j = ex.share[0].mait_beta_chain_junction_match;
        if alpha_g || alpha_j {
            s += "𝝰";
            if alpha_g {
                s += "g";
            }
            if alpha_j {
                s += "j";
            }
        }
        if beta_g || beta_j {
            s += "𝝱";
            if beta_g {
                s += "g";
            }
            if beta_j {
                s += "j";
            }
        }
        lvar![i, x, s.clone()];
    } else if x.starts_with("pe") {
        lvar![i, x, format!("")];
    } else if x.starts_with("npe") {
        lvar![i, x, format!("")];
    } else if x.starts_with("ppe") {
        lvar![i, x, format!("")];
    } else if x == "cred" || x == "cred_cell" {
        let mut credsx = Vec::<f64>::new();
        for l in 0..ex.clones.len() {
            let bc = &ex.clones[l][0].barcode;
            let li = ex.clones[l][0].dataset_index;
            if gex_info.pca[li].contains_key(&bc.clone()) {
                let mut creds = 0;
                let mut z = Vec::<(f64, String)>::new();
                let x = &gex_info.pca[li][&bc.clone()];
                for y in gex_info.pca[li].iter() {
                    let mut dist2 = 0.0;
                    for m in 0..x.len() {
                        dist2 += (y.1[m] - x[m]) * (y.1[m] - x[m]);
                    }
                    z.push((dist2, y.0.clone()));
                }
                z.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let top = n_vdj_gex[li];
                for i in 0..top {
                    if bin_member(&vdj_cells[li], &z[i].1) {
                        creds += 1;
                    }
                }
                let pc = 100.0 * creds as f64 / top as f64;
                credsx.push(pc);
            } else {
                credsx.push(0.0);
            }
        }
        credsx.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if x == "cred" {
            if credsx.is_empty() {
                lvar![i, x, format!("")];
            } else {
                lvar![i, x, format!("{:.1}", credsx[credsx.len() / 2])];
            }
        } else {
            if pass == 2 {
                let mut r = Vec::<String>::new();
                for j in 0..credsx.len() {
                    r.push(format!("{:.1}", credsx[j]));
                }
                speak!(u, x, format!("{}", r.iter().format(POUT_SEP)));
            }
        }
    } else if bin_member(&alt_bcs, x) {
        lvar![i, x, format!("")];
        if pass == 2 {
            let mut r = Vec::<String>::new();
            for l in 0..ex.clones.len() {
                let li = ex.clones[l][0].dataset_index;
                let bc = ex.clones[l][0].barcode.clone();
                let mut val = String::new();
                let alt = &ctl.origin_info.alt_bc_fields[li];
                for j in 0..alt.len() {
                    if alt[j].0 == *x {
                        if alt[j].1.contains_key(&bc.clone()) {
                            val = alt[j].1[&bc.clone()].clone();
                        }
                    }
                }
                r.push(val);
            }
            speak!(u, x, format!("{}", r.iter().format(POUT_SEP)));
        }
    } else if x.starts_with("n_") && !x.starts_with("n_gex") {
        let name = x.after("n_");
        let mut count = 0;
        let mut counts = Vec::<f64>::new();
        for j in 0..ex.clones.len() {
            let x = &ex.clones[j][0];
            if ctl.origin_info.dataset_id[x.dataset_index] == name {
                count += 1;
                counts.push(1.0);
            } else if x.origin_index.is_some()
                && ctl.origin_info.origin_list[x.origin_index.unwrap()] == name
            {
                count += 1;
                counts.push(1.0);
            } else if x.donor_index.is_some()
                && ctl.origin_info.donor_list[x.donor_index.unwrap()] == name
            {
                count += 1;
                counts.push(1.0);
            } else if x.tag_index.is_some()
                && ctl.origin_info.tag_list[x.tag_index.unwrap()] == name
            {
                count += 1;
                counts.push(1.0);
            }
        }
        lvar![i, x, format!("{}", count)];
        stats.push((x.to_string(), counts));
    } else if x == "sec" && ctl.gen_opt.using_secmem {
        let mut n = 0;
        let mut y = Vec::<f64>::new();
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = &ex.clones[l][0].barcode;
            let mut count = 0;
            if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                count = ctl.origin_info.secmem[li][&bc.clone()].0;
                n += count;
            }
            y.push(count as f64);
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x.to_string(), y));
    } else if x == "mem" && ctl.gen_opt.using_secmem {
        let mut n = 0;
        let mut y = Vec::<f64>::new();
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = &ex.clones[l][0].barcode;
            let mut count = 0;
            if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                count = ctl.origin_info.secmem[li][&bc.clone()].1;
                n += count;
            }
            y.push(count as f64);
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x.to_string(), y));
    } else if x == "dref" {
        let mut diffs = 0;
        for m in 0..cols {
            if mat[m][u].is_some() {
                let r = mat[m][u].unwrap();
                let seq = &ex.share[r].seq_del_amino;
                let mut vref = refdata.refs[rsi.vids[m]].to_ascii_vec();
                if rsi.vpids[m].is_some() {
                    vref = dref[rsi.vpids[m].unwrap()].nt_sequence.clone();
                }
                let jref = refdata.refs[rsi.jids[m]].to_ascii_vec();
                let z = seq.len();
                for p in 0..z {
                    let b = seq[p];
                    if p < vref.len() - ctl.heur.ref_v_trim && b != vref[p] {
                        diffs += 1;
                    }
                    if p >= z - (jref.len() - ctl.heur.ref_j_trim)
                        && b != jref[jref.len() - (z - p)]
                    {
                        diffs += 1;
                    }
                }
            }
        }
        lvar![i, x, format!("{}", diffs)];
    } else if x == "dref_aa" {
        let mut diffs = 0;
        for m in 0..cols {
            if mat[m][u].is_some() {
                let r = mat[m][u].unwrap();
                let aa_seq = &ex.share[r].aa_mod_indel;
                let mut vref = refdata.refs[rsi.vids[m]].to_ascii_vec();
                if rsi.vpids[m].is_some() {
                    vref = dref[rsi.vpids[m].unwrap()].nt_sequence.clone();
                }
                let jref = refdata.refs[rsi.jids[m]].to_ascii_vec();
                let z = 3 * aa_seq.len() + 1;
                for p in 0..aa_seq.len() {
                    if aa_seq[p] == b'-' {
                        diffs += 1;
                        continue;
                    }
                    if 3 * p + 3 <= vref.len() - ctl.heur.ref_v_trim
                        && aa_seq[p] != codon_to_aa(&vref[3 * p..3 * p + 3])
                    {
                        diffs += 1;
                    }
                    if 3 * p > z - (jref.len() - ctl.heur.ref_j_trim) + 3
                        && aa_seq[p]
                            != codon_to_aa(
                                &jref[jref.len() - (z - 3 * p)..jref.len() - (z - 3 * p) + 3],
                            )
                    {
                        diffs += 1;
                    }
                }
            }
        }
        lvar![i, x, format!("{}", diffs)];
    } else if x == "near" {
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
            lvar![i, x, "".to_string()];
        } else {
            lvar![i, x, format!("{}", dist)];
        }
    } else if x == "far" {
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
            lvar![i, x, "".to_string()];
        } else {
            lvar![i, x, format!("{}", dist)];
        }
    } else if x.starts_with("count_cdr1_") || x.contains(":count_cdr1_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_cdr1_") {
            x = x.before(":count_cdr1_").to_string();
        }
        y = y.after("count_cdr1_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr1_start.is_some() && ex.share[j].fr2_start.is_some() {
                let cdr1 = ex.share[j].cdr1_start.unwrap();
                let fwr2 = ex.share[j].fr2_start.unwrap();
                if cdr1 < fwr2 {
                    let aa = aa_seq(&ex.share[j].seq[cdr1..fwr2], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_cdr2_") || x.contains(":count_cdr2_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_cdr2_") {
            x = x.before(":count_cdr2_").to_string();
        }
        y = y.after("count_cdr2_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr2_start.is_some() && ex.share[j].fr3_start.is_some() {
                let cdr2 = ex.share[j].cdr2_start.unwrap();
                let fwr3 = ex.share[j].fr3_start.unwrap();
                if cdr2 < fwr3 {
                    let aa = aa_seq(&ex.share[j].seq[cdr2..fwr3], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_cdr3_") || x.contains(":count_cdr3_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_cdr3_") {
            x = x.before(":count_cdr3_").to_string();
        }
        y = y.after("count_cdr3_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            let cdr3 = ex.share[j].cdr3_start;
            let fwr4 = cdr3 + 3 * ex.share[j].cdr3_aa.len();
            let aa = aa_seq(&ex.share[j].seq[cdr3..fwr4], 0);
            n += reg.find_iter(&strme(&aa)).count();
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_fwr1_") || x.contains(":count_fwr1_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr1_") {
            x = x.before(":count_fwr1_").to_string();
        }
        y = y.after("count_fwr1_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr1_start.is_some() {
                let fwr1 = ex.share[j].fr1_start;
                let cdr1 = ex.share[j].cdr1_start.unwrap();
                if fwr1 < cdr1 {
                    let aa = aa_seq(&ex.share[j].seq[fwr1..cdr1], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_fwr2_") || x.contains(":count_fwr2_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr2_") {
            x = x.before(":count_fwr2_").to_string();
        }
        y = y.after("count_fwr2_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].fr2_start.is_some() && ex.share[j].cdr2_start.is_some() {
                let fwr2 = ex.share[j].fr2_start.unwrap();
                let cdr2 = ex.share[j].cdr2_start.unwrap();
                if fwr2 < cdr2 {
                    let aa = aa_seq(&ex.share[j].seq[fwr2..cdr2], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_fwr3_") || x.contains(":count_fwr3_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr3_") {
            x = x.before(":count_fwr3_").to_string();
        }
        y = y.after("count_fwr3_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].fr3_start.is_some() {
                let fwr3 = ex.share[j].fr3_start.unwrap();
                let cdr3 = ex.share[j].cdr3_start;
                if fwr3 < cdr3 {
                    let aa = aa_seq(&ex.share[j].seq[fwr3..cdr3], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_fwr4_") || x.contains(":count_fwr4_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr4_") {
            x = x.before(":count_fwr4_").to_string();
        }
        y = y.after("count_fwr4_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            let fwr4 = ex.share[j].cdr3_start + 3 * ex.share[j].cdr3_aa.len();
            let aa = aa_seq(&ex.share[j].seq[fwr4..], 0);
            n += reg.find_iter(&strme(&aa)).count();
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_cdr_") || x.contains(":count_cdr_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_cdr_") {
            x = x.before(":count_cdr_").to_string();
        }
        y = y.after("count_cdr_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr1_start.is_some() && ex.share[j].fr2_start.is_some() {
                let cdr1 = ex.share[j].cdr1_start.unwrap();
                let fwr2 = ex.share[j].fr2_start.unwrap();
                if cdr1 < fwr2 {
                    let aa = aa_seq(&ex.share[j].seq[cdr1..fwr2], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
            if ex.share[j].cdr2_start.is_some() && ex.share[j].fr3_start.is_some() {
                let cdr2 = ex.share[j].cdr2_start.unwrap();
                let fwr3 = ex.share[j].fr3_start.unwrap();
                if cdr2 < fwr3 {
                    let aa = aa_seq(&ex.share[j].seq[cdr2..fwr3], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
            let cdr3 = ex.share[j].cdr3_start;
            let fwr4 = cdr3 + 3 * ex.share[j].cdr3_aa.len();
            let aa = aa_seq(&ex.share[j].seq[cdr3..fwr4], 0);
            n += reg.find_iter(&strme(&aa)).count();
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_fwr_") || x.contains(":count_fwr_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_fwr_") {
            x = x.before(":count_fwr_").to_string();
        }
        y = y.after("count_fwr_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            if ex.share[j].cdr1_start.is_some() {
                let fwr1 = ex.share[j].fr1_start;
                let cdr1 = ex.share[j].cdr1_start.unwrap();
                if fwr1 < cdr1 {
                    let aa = aa_seq(&ex.share[j].seq[fwr1..cdr1], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
            if ex.share[j].fr2_start.is_some() && ex.share[j].cdr2_start.is_some() {
                let fwr2 = ex.share[j].fr2_start.unwrap();
                let cdr2 = ex.share[j].cdr2_start.unwrap();
                if fwr2 < cdr2 {
                    let aa = aa_seq(&ex.share[j].seq[fwr2..cdr2], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
            if ex.share[j].fr3_start.is_some() {
                let fwr3 = ex.share[j].fr3_start.unwrap();
                let cdr3 = ex.share[j].cdr3_start;
                if fwr3 < cdr3 {
                    let aa = aa_seq(&ex.share[j].seq[fwr3..cdr3], 0);
                    n += reg.find_iter(&strme(&aa)).count();
                }
            }
            let fwr4 = ex.share[j].cdr3_start + 3 * ex.share[j].cdr3_aa.len();
            let aa = aa_seq(&ex.share[j].seq[fwr4..], 0);
            n += reg.find_iter(&strme(&aa)).count();
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x.starts_with("count_") || x.contains(":count_") {
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(":count_") {
            x = x.before(":count_").to_string();
        }
        y = y.after("count_").to_string();
        let reg = Regex::new(&y).unwrap(); // seems inefficient
        let mut n = 0;
        for j in 0..ex.share.len() {
            let aa = aa_seq(&ex.share[j].seq, 0); // seems inefficient
            n += reg.find_iter(&strme(&aa)).count();
        }
        lvar![i, x, format!("{}", n)];
        stats.push((x, vec![n as f64; ex.ncells()]));
    } else if x == "gex" {
        lvar![i, x, format!("{}", gex_median)];
    } else if x == "gex_cell" {
        if pass == 2 {
            speak!(u, x, format!("{}", count_unsorted.iter().format(POUT_SEP)));
        }
    } else if x == "n_gex" {
        lvar![i, x, format!("{}", n_gex)];
    } else if x == "n_gex_cell" {
        if i < lvars.len() {
            row.push("".to_string());
        }
        if pass == 2 {
            speak!(
                u,
                "n_gex_cell".to_string(),
                format!("{}", n_gexs.iter().format(POUT_SEP))
            );
        }
    } else if x == "entropy" {
        lvar![i, x, format!("{:.2}", entropy)];
    } else if x == "entropy_cell" {
        let mut e = Vec::<String>::new();
        for x in entropies_unsorted.iter() {
            e.push(format!("{:.2}", x));
        }
        speak!(u, x, format!("{}", e.iter().format(POUT_SEP)));
    } else if x == "gex_min" {
        lvar![i, x, format!("{}", gex_min)];
    } else if x == "gex_max" {
        lvar![i, x, format!("{}", gex_max)];
    } else if x == "gex_μ" {
        lvar![i, x, format!("{}", gex_mean.round() as usize)];
    } else if x == "gex_Σ" {
        lvar![i, x, format!("{}", gex_sum.round() as usize)];
    } else if x == "ext" {
        let mut exts = Vec::<String>::new();
        for l in 0..ex.clones.len() {
            let li = ctl.origin_info.dataset_id[ex.clones[l][0].dataset_index].clone();
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
        lvar![i, x, s.clone()];
    } else {
        let (mut counts_sub, mut fcounts_sub) = (Vec::<f64>::new(), Vec::<f64>::new());
        let xorig = x.clone();
        let (mut x, mut y) = (x.to_string(), x.to_string());
        if x.contains(':') {
            x = x.before(":").to_string();
        }
        if y.contains(':') {
            y = y.after(":").to_string();
        }
        let y0 = y.clone();
        for _ in 1..=2 {
            let suffixes = ["_min", "_max", "_μ", "_Σ", "_cell", "_%"];
            for s in suffixes.iter() {
                if y.ends_with(s) {
                    y = y.rev_before(&s).to_string();
                    break;
                }
            }
        }
        let mut computed = false;
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = ex.clones[l][0].barcode.clone();
            let mut ux = Vec::<usize>::new();
            if ctl.clono_print_opt.regex_match[li].contains_key(&y) {
                ux = ctl.clono_print_opt.regex_match[li][&y].clone();
            }
            if ux.len() > 0 {
                let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                if p >= 0 {
                    computed = true;
                    let mut raw_count = 0.0;
                    for fid in ux.iter() {
                        let raw_counti = get_gex_matrix_entry(
                            &ctl, &gex_info, *fid, &d_all, &ind_all, li, l, p as usize, &y,
                        );
                        raw_count += raw_counti;
                    }
                    counts_sub.push(raw_count.round() as f64);
                    fcounts_sub.push(raw_count);
                }
            } else {
                if gex_info.feature_id[li].contains_key(&y) {
                    computed = true;
                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                    if p >= 0 {
                        let fid = gex_info.feature_id[li][&y];
                        let raw_count = get_gex_matrix_entry(
                            &ctl, &gex_info, fid, &d_all, &ind_all, li, l, p as usize, &y,
                        );
                        counts_sub.push(raw_count.round() as f64);
                        fcounts_sub.push(raw_count);
                    }
                }
            }
        }
        if computed {
            if !y0.ends_with("_%") {
                stats.push((x.clone(), fcounts_sub.clone()));
            } else {
                let mut f = Vec::<f64>::new();
                for i in 0..fcounts_sub.len() {
                    let mut x = 0.0;
                    if gex_mean > 0.0 {
                        x = 100.0 * fcounts_sub[i] / gex_mean;
                    }
                    f.push(x);
                }
                stats.push((x.clone(), f));
            }
            let mut counts_sub_sorted = counts_sub.clone();
            counts_sub_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let sum = fcounts_sub.iter().sum::<f64>();
            let mean = sum / counts_sub.len() as f64;

            if xorig.ends_with("_%_cell") {
                if pass == 2 {
                    let mut c = Vec::<String>::new();
                    for j in 0..counts_sub.len() {
                        c.push(format!("{:.2}", 100.0 * counts_sub[j] as f64 / fcounts[j]));
                    }
                    let val = format!("{}", c.iter().format(POUT_SEP));
                    speak!(u, x, val);
                }
            } else if xorig.ends_with("_cell") {
                if pass == 2 {
                    let val = format!("{}", counts_sub.iter().format(POUT_SEP));
                    speak!(u, x, val);
                }
            } else {
                if y0.ends_with("_min") {
                    lvar![i, x, format!("{}", counts_sub_sorted[0].round())];
                } else if y0.ends_with("_max") {
                    lvar![
                        i,
                        x,
                        format!("{}", counts_sub_sorted[counts_sub.len() - 1].round())
                    ];
                } else if y0.ends_with("_μ") {
                    lvar![i, x, format!("{}", mean.round())];
                } else if y0.ends_with("_Σ") {
                    lvar![i, x, format!("{}", sum.round())];
                } else if y0.ends_with("_%") {
                    lvar![i, x, format!("{:.2}", (100.0 * sum) / gex_sum)];
                } else {
                    let mut median = 0.0;
                    if counts_sub_sorted.len() > 0 {
                        median = median_f64(&counts_sub_sorted);
                    }
                    lvar![i, x, format!("{}", median)];
                }
            }
        } else if i < lvars.len() {
            lvar![i, x, "".to_string()];
        }
    }
}
