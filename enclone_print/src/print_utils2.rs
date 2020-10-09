// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This file contains the single function row_fill,
// plus a small helper function get_gex_matrix_entry.

use crate::print_utils1::*;
use amino::*;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use itertools::*;
use ndarray::s;
use stats_utils::*;
use std::cmp::{max, min};
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

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

// The following code creates a row in the enclone output table for a clonotype.  Simultaneously
// it generates a row of parseable output.  And it does some other things that are not described
// here.
//
// TODO: Awful interface, should work to improve.

pub fn row_fill(
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
    vars_amino: &Vec<Vec<usize>>,
    show_aa: &Vec<Vec<usize>>,
    field_types: &Vec<Vec<u8>>,
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
    d_readers: &Vec<Option<hdf5::Reader>>,
    ind_readers: &Vec<Option<hdf5::Reader>>,
    h5_data: &Vec<(usize, Vec<u32>, Vec<u32>)>,
    stats: &mut Vec<(String, Vec<f64>)>,
    vdj_cells: &Vec<Vec<String>>,
    n_vdj_gex: &Vec<usize>,
    lvarsc: &Vec<String>,
    nd_fields: &Vec<String>,
    peer_groups: &Vec<Vec<(usize, u8, u32)>>,
) {
    // Redefine some things to reduce dependencies.

    let mat = &rsi.mat;
    let cvars = &ctl.clono_print_opt.cvars;
    let lvars = lvarsc.clone();
    let clonotype_id = exacts[u];
    let ex = &exact_clonotypes[clonotype_id];
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
    let mut pcols_sort = ctl.parseable_opt.pcols_sort.clone();
    for i in 0..pcols_sort.len() {
        pcols_sort[i] = pcols_sort[i].replace("_Σ", "_sum");
        pcols_sort[i] = pcols_sort[i].replace("_μ", "_mean");
    }
    pcols_sort.sort();
    macro_rules! speakc {
        ($u:expr, $col:expr, $var:expr, $val:expr) => {
            if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && $col + 1 <= ctl.parseable_opt.pchains
            {
                let mut v = $var.clone();
                v = v.replace("_Σ", "_sum");
                v = v.replace("_μ", "_mean");

                // Strip escape character sequences from val.  Can happen in notes, maybe
                // other places.

                let mut val_clean = String::new();
                let mut chars = Vec::<char>::new();
                let valx = format!("{}", $val);
                for c in valx.chars() {
                    chars.push(c);
                }
                let mut escaped = false;
                for l in 0..chars.len() {
                    if chars[l] == '' {
                        escaped = true;
                    }
                    if escaped {
                        if chars[l] == 'm' {
                            escaped = false;
                        }
                        continue;
                    }
                    val_clean.push(chars[l]);
                }

                // Proceed.

                let varc = format!("{}{}", v, $col + 1);
                if pcols_sort.is_empty() || bin_member(&pcols_sort, &varc) {
                    out_data[$u].insert(varc, val_clean);
                }
            }
        };
    }
    let cols = varmat[0].len();

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

    // Compute dataset indices, gex, gex_min, gex_max, gex_mean, gex_sum,
    // n_gex_cell, n_gex, entropy.

    let mut dataset_indices = Vec::<usize>::new();
    for l in 0..ex.clones.len() {
        dataset_indices.push(ex.clones[l][0].dataset_index);
    }
    unique_sort(&mut dataset_indices);
    let mut lenas = Vec::<String>::new();
    for l in dataset_indices.iter() {
        lenas.push(ctl.origin_info.dataset_id[*l].clone());
    }
    row.push("".to_string()); // row number (#), filled in below
    let mut counts = Vec::<usize>::new();
    let mut fcounts = Vec::<f64>::new();
    let mut n_gex = 0;
    let mut n_gexs = Vec::<usize>::new();
    let mut total_counts = Vec::<usize>::new();
    // It might be possible to speed this up a lot by pulling part of the "let d" and
    // "let ind" constructs out of the loop.
    if lvars.contains(&"entropy".to_string()) {
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = ex.clones[l][0].barcode.clone();
            if !gex_info.gex_barcodes.is_empty() {
                let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                if p >= 0 {
                    let mut raw_count = 0;
                    if gex_info.gex_matrices[li].initialized() {
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
                        let d: Vec<u32>;
                        let ind: Vec<u32>;
                        if ctl.gen_opt.h5_pre {
                            d = h5_data[li].1[z1..z2].to_vec();
                            ind = h5_data[li].2[z1..z2].to_vec();
                        } else {
                            d = d_readers[li]
                                .as_ref()
                                .unwrap()
                                .read_slice(s![z1..z2])
                                .unwrap()
                                .to_vec();
                            ind = ind_readers[li]
                                .as_ref()
                                .unwrap()
                                .read_slice(s![z1..z2])
                                .unwrap()
                                .to_vec();
                        }
                        for j in 0..d.len() {
                            if gex_info.is_gex[li][ind[j] as usize] {
                                raw_count += d[j] as usize;
                            }
                        }
                        d_all[l] = d;
                        ind_all[l] = ind;
                    }
                    total_counts.push(raw_count);
                }
            }
        }
    }
    let mut entropies = Vec::<f64>::new();
    for l in 0..ex.clones.len() {
        let li = ex.clones[l][0].dataset_index;
        let bc = ex.clones[l][0].barcode.clone();
        if !gex_info.gex_barcodes.is_empty() {
            if bin_member(&gex_info.gex_cell_barcodes[li], &bc) {
                n_gex += 1;
                n_gexs.push(1);
            } else {
                n_gexs.push(0);
            }
            let mut count = 0;
            let mut fcount = 0.0;
            let mut entropy = 0.0;
            let p = bin_position(&gex_info.gex_barcodes[li], &bc);
            if p >= 0 {
                let mut raw_count = 0;
                if gex_info.gex_matrices[li].initialized() {
                    let row = gex_info.gex_matrices[li].row(p as usize);
                    for j in 0..row.len() {
                        let f = row[j].0;
                        let n = row[j].1;
                        if gex_info.is_gex[li][f] {
                            if lvars.contains(&"entropy".to_string()) {
                                let q = n as f64 / total_counts[l] as f64;
                                entropy -= q * q.log2();
                            }
                            raw_count += n;
                        }
                    }
                } else {
                    let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                    let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize; // is p+1 OK??
                    let d: Vec<u32>;
                    let ind: Vec<u32>;
                    if ctl.gen_opt.h5_pre {
                        d = h5_data[li].1[z1..z2].to_vec();
                        ind = h5_data[li].2[z1..z2].to_vec();
                    } else {
                        d = d_readers[li]
                            .as_ref()
                            .unwrap()
                            .read_slice(s![z1..z2])
                            .unwrap()
                            .to_vec();
                        ind = ind_readers[li]
                            .as_ref()
                            .unwrap()
                            .read_slice(s![z1..z2])
                            .unwrap()
                            .to_vec();
                    }
                    for j in 0..d.len() {
                        if gex_info.is_gex[li][ind[j] as usize] {
                            let n = d[j] as usize;
                            if lvars.contains(&"entropy".to_string()) {
                                let q = n as f64 / total_counts[l] as f64;
                                entropy -= q * q.log2();
                            }
                            raw_count += n;
                        }
                    }
                    d_all[l] = d;
                    ind_all[l] = ind;
                }
                if !ctl.gen_opt.full_counts {
                    count = (raw_count as f64 * gex_info.gex_mults[li]).round() as usize;
                    fcount = raw_count as f64 * gex_info.gex_mults[li];
                } else {
                    count = (raw_count as f64).round() as usize;
                    fcount = raw_count as f64;
                }
            }
            counts.push(count);
            fcounts.push(fcount);
            entropies.push(entropy);
        }
    }
    let count_unsorted = counts.clone();
    counts.sort();
    for n in counts.iter() {
        if *n < 100 {
            *gex_low += 1;
        }
    }
    let (mut gex_median, mut gex_min, mut gex_max, mut gex_mean, mut gex_sum) = (0, 0, 0, 0.0, 0.0);
    if counts.len() > 0 {
        gex_median = counts[counts.len() / 2];
        gex_min = counts[0];
        gex_max = counts[counts.len() - 1];
        gex_sum = fcounts.iter().sum::<f64>();
        gex_mean = gex_sum / fcounts.len() as f64;
    }
    let entropies_unsorted = entropies.clone();
    entropies.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut entropy = 0.0;
    if entropies.len() > 0 {
        entropy = entropies[entropies.len() / 2];
    }

    // Output lead variable columns.
    // WARNING!  If you add lead variables, you may need to add them to the function
    // LinearCondition::require_valid_variables.

    let mut all_lvars = lvars.clone();
    if ctl.parseable_opt.pout.len() == 0 {
    } else if ctl.parseable_opt.pcols.is_empty() {
        for i in 0..LVARS_ALLOWED.len() {
            if !lvars.contains(&LVARS_ALLOWED[i].to_string()) {
                all_lvars.push(LVARS_ALLOWED[i].to_string());
            }
        }
    } else {
        for i in 0..ctl.parseable_opt.pcols.len() {
            if !lvars.contains(&ctl.parseable_opt.pcols[i].to_string()) {
                all_lvars.push(ctl.parseable_opt.pcols[i].to_string());
            }
        }
    }
    let mut alt_bcs = Vec::<String>::new();
    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.clone());
        }
    }
    unique_sort(&mut alt_bcs);
    for i in 0..all_lvars.len() {
        let x = &all_lvars[i];
        if x.starts_with('g') && x.after("g").parse::<usize>().is_ok() {
            let d = x.after("g").force_usize();
            lvar![i, x, format!("{}", groups[&d][u] + 1)];
        } else if x == "origins" {
            let mut origins = Vec::<String>::new();
            for j in 0..ex.clones.len() {
                if ex.clones[j][0].origin_index.is_some() {
                    origins.push(
                        ctl.origin_info.origin_id[ex.clones[j][0].origin_index.unwrap()].clone(),
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
                    donors.push(
                        ctl.origin_info.donor_list[ex.clones[j][0].donor_index.unwrap()].clone(),
                    );
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
                    speak!(u, x, format!("{}", r.iter().format(";")));
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
                speak!(u, x, format!("{}", r.iter().format(";")));
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
        } else if x == "gex" {
            lvar![i, x, format!("{}", gex_median)];
        } else if x == "gex_cell" {
            if pass == 2 {
                speak!(u, x, format!("{}", count_unsorted.iter().format(";")));
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
                    format!("{}", n_gexs.iter().format(";"))
                );
            }
        } else if x == "entropy" {
            lvar![i, x, format!("{:.2}", entropy)];
        } else if x == "entropy_cell" {
            let mut e = Vec::<String>::new();
            for x in entropies_unsorted.iter() {
                e.push(format!("{:.2}", x));
            }
            speak!(u, x, format!("{}", e.iter().format(";")));
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
                        let val = format!("{}", c.iter().format(";"));
                        speak!(u, x, val);
                    }
                } else if xorig.ends_with("_cell") {
                    if pass == 2 {
                        let val = format!("{}", counts_sub.iter().format(";"));
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
                            median = counts_sub_sorted[counts_sub_sorted.len() / 2].round();
                        }
                        lvar![i, x, format!("{}", median)];
                    }
                }
            }
        }
    }

    // Sanity check.  It's here because if it fails and that failure was not detected, something
    // exceptionally cryptic would happen downstream.

    assert_eq!(row.len(), lvars.len() + 1);

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
        let u_mean = (utot as f64 / numis.len() as f64).round() as usize;
        let u_min = *numis.iter().min().unwrap();
        let u_max = *numis.iter().max().unwrap();
        nreads.sort();
        let rtot: usize = nreads.iter().sum();
        let r_mean = (rtot as f64 / nreads.len() as f64).round() as usize;
        let r_min = *nreads.iter().min().unwrap();
        let r_max = *nreads.iter().max().unwrap();
        let median_nreads = nreads[nreads.len() / 2];

        // Set up chain variable macro.  This is the mechanism for generating
        // both human-readable and parseable output for chain variables.

        macro_rules! cvar {
            ($i: expr, $var:expr, $val:expr) => {
                if $i < rsi.cvars[col].len() && cvars.contains(&$var) {
                    cx[col][$i] = $val.clone();
                }
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

        let xm = &ex.share[mid];
        speakc!(u, col, "vj_aa".to_string(), stringme(&aa_seq(&xm.seq, 0)));
        speakc!(
            u,
            col,
            "vj_aa_nl".to_string(),
            stringme(&aa_seq(&xm.seq, xm.fr1_start))
        );
        speakc!(u, col, "vj_seq".to_string(), stringme(&xm.seq));
        speakc!(
            u,
            col,
            "vj_seq_nl".to_string(),
            stringme(&xm.seq[xm.fr1_start..])
        );
        speakc!(u, col, "seq".to_string(), stringme(&xm.full_seq));
        speakc!(u, col, "v_start".to_string(), xm.v_start);
        let (mut d_start, mut d_frame) = (String::new(), String::new());
        if xm.d_start.is_some() {
            d_start = format!("{}", xm.d_start.unwrap());
            d_frame = format!("{}", (xm.d_start.unwrap() - xm.v_start) % 3);
        }
        speakc!(u, col, "d_start".to_string(), d_start);
        speakc!(u, col, "d_frame".to_string(), d_frame);
        let cid = xm.c_ref_id;
        if cid.is_some() {
            let cid = cid.unwrap();
            speakc!(u, col, "const_id".to_string(), refdata.id[cid]);
        }
        let uid = ex.share[mid].u_ref_id;
        if uid.is_some() {
            let uid = uid.unwrap();
            speakc!(u, col, "utr_id".to_string(), refdata.id[uid]);
            speakc!(u, col, "utr_name".to_string(), refdata.name[uid]);
        }
        speakc!(u, col, "cdr3_start".to_string(), xm.cdr3_start);
        speakc!(u, col, "cdr3_aa".to_string(), xm.cdr3_aa);
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
        speakc!(u, col, "var_aa".to_string(), strme(&varaa));

        // Create column entry.

        let rsi_vars = &rsi.cvars[col];
        let mut all_vars = rsi_vars.clone();
        for j in 0..CVARS_ALLOWED.len() {
            let var = &CVARS_ALLOWED[j];
            if !rsi_vars.contains(&var.to_string()) {
                all_vars.push(var.to_string());
            }
        }
        for j in 0..CVARS_ALLOWED_PCELL.len() {
            let var = &CVARS_ALLOWED_PCELL[j];
            if !rsi_vars.contains(&var.to_string()) {
                all_vars.push(var.to_string());
            }
        }
        let mut somelist = vec![false; all_vars.len()];
        for j in 0..all_vars.len() {
            let var = &all_vars[j];
            let varc = format!("{}{}", var, col + 1);
            if j < rsi.cvars[col].len() && cvars.contains(&var) {
                somelist[j] = true;
            } else if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && col + 1 <= ctl.parseable_opt.pchains
                && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc))
            {
                somelist[j] = true;
            }
        }
        for j in 0..all_vars.len() {
            // Decide if there is nothing to compute.  This is almost certainly not optimal.

            let mut needed = false;
            let var = &all_vars[j];
            let varc = format!("{}{}", var, col + 1);
            if j < rsi.cvars[col].len() && cvars.contains(&var) {
                needed = true;
            } else if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && col + 1 <= ctl.parseable_opt.pchains
                && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc))
            {
                needed = true;
            } else if *var == "amino".to_string() {
                needed = true;
            } else if *var == "u_cell".to_string() || *var == "r_cell".to_string() {
                needed = true;
            } else if *var == "white".to_string() || ctl.clono_filt_opt.whitef {
                needed = true;
            }
            if !needed {
                continue;
            }
            let col_var = j < rsi_vars.len();
            if !col_var && ctl.parseable_opt.pout.len() == 0 {
                continue;
            }

            // Compute.

            if *var == "amino".to_string() && col_var {
                let mut last_color = "black".to_string();
                for k in 0..show_aa[col].len() {
                    let p = show_aa[col][k];
                    if k > 0 && field_types[col][k] != field_types[col][k - 1] {
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
                        let x = &peer_groups[rsi.vids[col]];
                        let last = k == show_aa[col].len() - 1;
                        let log = color_codon(&ctl, &seq_amino, &x, p, &mut last_color, last);
                        cx[col][j] += strme(&log);
                    }
                }
            } else if *var == "comp".to_string() || *var == "edit".to_string() {
                let mut comp = 1000000;
                let mut edit = String::new();
                let td = &ex.share[mid];
                let tig = &td.seq;
                let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
                let mut aligner = Aligner::new(-6, -1, &score);

                // Go through passes.  If IGH/TRB, we go through every D segment.  Otherwise
                // there is just one pass.

                let mut z = 1;
                if ex.share[mid].left {
                    z = refdata.ds.len();
                }
                for di in 0..z {
                    let mut d = 0;
                    if ex.share[mid].left {
                        d = refdata.ds[di];
                    }

                    // Start to build reference concatenation.  First append the V segment.

                    let mut concat = Vec::<u8>::new();
                    let mut vref = refdata.refs[rsi.vids[col]].to_ascii_vec();
                    if rsi.vpids[col].is_none() {
                    } else {
                        vref = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
                    }
                    concat.append(&mut vref.clone());

                    // Append the D segment if IGH/TRB.

                    if ex.share[mid].left {
                        let mut x = refdata.refs[d].to_ascii_vec();
                        concat.append(&mut x);
                    }

                    // Append the J segment.

                    let mut x = refdata.refs[rsi.jids[col]].to_ascii_vec();
                    concat.append(&mut x);

                    // Align the V..J sequence on the contig to the reference concatenation.

                    let al = aligner.semiglobal(&tig, &concat);
                    let mut m = 0;
                    let mut pos = al.xstart;
                    let mut rpos = (al.ystart as isize) - (vref.len() as isize);
                    let mut count = 0;
                    let start = td.cdr3_start;
                    let stop = td.j_stop - td.v_start;
                    let mut edits = Vec::<String>::new();
                    while m < al.operations.len() {
                        let n = next_diff(&al.operations, m);
                        match al.operations[m] {
                            Match => {
                                pos += 1;
                                rpos += 1;
                            }
                            Subst => {
                                if pos >= start && pos < stop {
                                    count += 1;
                                    edits.push(format!("S{}", rpos));
                                }
                                pos += 1;
                                rpos += 1;
                            }
                            Del => {
                                if pos >= start && pos < stop {
                                    count += 1;
                                    edits.push(format!("D{}:{}", rpos, n - m));
                                }
                                pos += n - m;
                                m = n - 1;
                            }
                            Ins => {
                                if pos >= start && pos < stop {
                                    count += 1;
                                    edits.push(format!("I{}:{}", rpos, n - m));
                                }
                                rpos += (n - m) as isize;
                                m = n - 1;
                            }
                            _ => {}
                        };
                        m += 1;
                    }
                    if count < comp {
                        comp = count;
                        edit = format!("{}", edits.iter().format("•"));
                    }
                }
                if *var == "comp".to_string() {
                    cvar![j, var, format!("{}", comp)];
                } else {
                    cvar![j, var, format!("{}", edit)];
                }
            } else if *var == "cdr1_dna".to_string()
                || *var == "cdr1_aa".to_string()
                || *var == "cdr1_len".to_string()
            {
                let x = &ex.share[mid];
                let mut y = "unknown".to_string();
                if x.cdr1_start.is_some()
                    && x.fr2_start.is_some()
                    && x.cdr1_start.unwrap() <= x.fr2_start.unwrap()
                {
                    let mut dna = Vec::<u8>::new();
                    for p in x.cdr1_start.unwrap()..x.fr2_start.unwrap() {
                        for j in 0..x.ins.len() {
                            if x.ins[j].0 == p {
                                let mut z = x.ins[j].1.clone();
                                dna.append(&mut z);
                            }
                        }
                        if x.seq_del_amino[p] != b'-' {
                            dna.push(x.seq_del_amino[p]);
                        }
                    }

                    // Test for internal error.

                    let mut found = false;
                    for i in 0..x.seq.len() {
                        if x.seq[i..].starts_with(&dna) {
                            found = true;
                        }
                    }
                    if !found {
                        eprintln!(
                            "\nInternal error, failed to find {}, CDR3 = {}.\n",
                            strme(&dna),
                            x.cdr3_aa
                        );
                        std::process::exit(1);
                    }
                    if *var == "cdr1_dna".to_string() {
                        y = stringme(&dna);
                    } else if *var == "cdr1_aa".to_string() {
                        y = stringme(&aa_seq(&dna, 0));
                    } else {
                        y = format!("{}", dna.len() / 3);
                    }
                }
                cvar![j, var, y];
            } else if *var == "cdr2_dna".to_string()
                || *var == "cdr2_aa".to_string()
                || *var == "cdr2_len".to_string()
            {
                let x = &ex.share[mid];
                let mut y = "unknown".to_string();
                if x.cdr2_start.is_some()
                    && x.fr3_start.is_some()
                    && x.cdr2_start.unwrap() <= x.fr3_start.unwrap()
                {
                    let mut dna = Vec::<u8>::new();
                    for p in x.cdr2_start.unwrap()..x.fr3_start.unwrap() {
                        for j in 0..x.ins.len() {
                            if x.ins[j].0 == p {
                                let mut z = x.ins[j].1.clone();
                                dna.append(&mut z);
                            }
                        }
                        if x.seq_del_amino[p] != b'-' {
                            dna.push(x.seq_del_amino[p]);
                        }
                    }

                    // Test for internal error.

                    let mut found = false;
                    for i in 0..x.seq.len() {
                        if x.seq[i..].starts_with(&dna) {
                            found = true;
                        }
                    }
                    if !found {
                        eprintln!(
                            "\nInternal error, failed to find {}, CDR3 = {}.\n",
                            strme(&dna),
                            x.cdr3_aa
                        );
                        std::process::exit(1);
                    }
                    if *var == "cdr2_dna".to_string() {
                        y = stringme(&dna);
                    } else if *var == "cdr2_aa".to_string() {
                        y = stringme(&aa_seq(&dna, 0));
                    } else {
                        y = format!("{}", dna.len() / 3);
                    }
                }
                cvar![j, var, y];
            } else if *var == "cdr3_aa".to_string() {
                cvar![j, var, ex.share[mid].cdr3_aa.clone()];
            } else if *var == "cdr3_dna".to_string() {
                cvar![j, var, ex.share[mid].cdr3_dna.clone()];
            } else if *var == "cdr3_len".to_string() {
                cvar![j, var, ex.share[mid].cdr3_aa.len().to_string()];
            } else if *var == "fwr1_dna".to_string()
                || *var == "fwr1_aa".to_string()
                || *var == "fwr1_len".to_string()
            {
                let x = &ex.share[mid];
                let mut y = "unknown".to_string();
                if x.cdr1_start.is_some() && x.fr1_start <= x.cdr1_start.unwrap() {
                    let mut dna = Vec::<u8>::new();
                    for p in x.fr1_start..x.cdr1_start.unwrap() {
                        for j in 0..x.ins.len() {
                            if x.ins[j].0 == p {
                                let mut z = x.ins[j].1.clone();
                                dna.append(&mut z);
                            }
                        }
                        if x.seq_del_amino[p] != b'-' {
                            dna.push(x.seq_del_amino[p]);
                        }
                    }

                    // Test for internal error.

                    let mut found = false;
                    for i in 0..x.seq.len() {
                        if x.seq[i..].starts_with(&dna) {
                            found = true;
                        }
                    }
                    if !found {
                        eprintln!(
                            "\nInternal error, failed to find {}, CDR3 = {}.\n",
                            strme(&dna),
                            x.cdr3_aa
                        );
                        std::process::exit(1);
                    }
                    if *var == "fwr1_dna".to_string() {
                        y = stringme(&dna);
                    } else if *var == "fwr1_aa".to_string() {
                        y = stringme(&aa_seq(&dna, 0));
                    } else {
                        y = format!("{}", dna.len() / 3);
                    }
                }
                cvar![j, var, y];
            } else if *var == "fwr2_dna".to_string()
                || *var == "fwr2_aa".to_string()
                || *var == "fwr2_len".to_string()
            {
                let x = &ex.share[mid];
                let mut y = "unknown".to_string();
                if x.fr2_start.unwrap() <= x.cdr2_start.unwrap() {
                    let mut dna = Vec::<u8>::new();
                    for p in x.fr2_start.unwrap()..x.cdr2_start.unwrap() {
                        for j in 0..x.ins.len() {
                            if x.ins[j].0 == p {
                                let mut z = x.ins[j].1.clone();
                                dna.append(&mut z);
                            }
                        }
                        if x.seq_del_amino[p] != b'-' {
                            dna.push(x.seq_del_amino[p]);
                        }
                    }

                    // Test for internal error.

                    let mut found = false;
                    for i in 0..x.seq.len() {
                        if x.seq[i..].starts_with(&dna) {
                            found = true;
                        }
                    }
                    if !found {
                        eprintln!(
                            "\nInternal error, failed to find {}, CDR3 = {}.\n",
                            strme(&dna),
                            x.cdr3_aa
                        );
                        std::process::exit(1);
                    }
                    if *var == "fwr2_dna".to_string() {
                        y = stringme(&dna);
                    } else if *var == "fwr2_aa".to_string() {
                        y = stringme(&aa_seq(&dna, 0));
                    } else {
                        y = format!("{}", dna.len() / 3);
                    }
                }
                cvar![j, var, y];
            } else if *var == "fwr3_dna".to_string()
                || *var == "fwr3_aa".to_string()
                || *var == "fwr3_len".to_string()
            {
                let x = &ex.share[mid];
                let mut y = "unknown".to_string();
                if x.fr3_start.is_some() && x.fr3_start.unwrap() <= x.cdr3_start {
                    let mut dna = Vec::<u8>::new();
                    for p in x.fr3_start.unwrap()..x.cdr3_start {
                        for j in 0..x.ins.len() {
                            if x.ins[j].0 == p {
                                let mut z = x.ins[j].1.clone();
                                dna.append(&mut z);
                            }
                        }
                        if x.seq_del_amino[p] != b'-' {
                            dna.push(x.seq_del_amino[p]);
                        }
                    }

                    // Test for internal error.

                    let mut found = false;
                    for i in 0..x.seq.len() {
                        if x.seq[i..].starts_with(&dna) {
                            found = true;
                        }
                    }
                    if !found {
                        eprintln!(
                            "\nInternal error, failed to find {}, CDR3 = {}.\n",
                            strme(&dna),
                            x.cdr3_aa
                        );
                        std::process::exit(1);
                    }
                    if *var == "fwr3_dna".to_string() {
                        y = stringme(&dna);
                    } else if *var == "fwr3_aa".to_string() {
                        y = stringme(&aa_seq(&dna, 0));
                    } else {
                        y = format!("{}", dna.len() / 3);
                    }
                }
                cvar![j, var, y];
            } else if *var == "fwr4_dna".to_string()
                || *var == "fwr4_aa".to_string()
                || *var == "fwr4_len".to_string()
            {
                let x = &ex.share[mid];
                let start = x.cdr3_start + 3 * x.cdr3_aa.len();
                let stop = x.seq_del_amino.len();
                let dna = &x.seq_del_amino[start..stop];
                let y;
                if *var == "fwr4_dna".to_string() {
                    y = stringme(&dna);
                } else if *var == "fwr4_aa".to_string() {
                    y = stringme(&aa_seq(&dna.to_vec(), 0));
                } else {
                    y = format!("{}", dna.len() / 3);
                }
                cvar![j, var, y];
            } else if *var == "ulen".to_string() {
                cvar![j, *var, format!("{}", ex.share[mid].v_start)];
            } else if *var == "clen".to_string() {
                cvar![
                    j,
                    var,
                    format!("{}", ex.share[mid].full_seq.len() - ex.share[mid].j_stop)
                ];
            } else if *var == "aa%".to_string() {
                let xm = &ex.share[mid];
                let mut diffs = 0;
                let mut denom = 0;
                let aa_seq = &xm.aa_mod_indel;
                let mut vref = refdata.refs[xm.v_ref_id].to_ascii_vec();
                if xm.v_ref_id_donor_alt_id.is_some() {
                    vref = dref[xm.v_ref_id_donor.unwrap()].nt_sequence.clone();
                }
                let jref = refdata.refs[xm.j_ref_id].to_ascii_vec();
                let z = 3 * aa_seq.len() + 1;
                for p in 0..aa_seq.len() {
                    if aa_seq[p] == b'-' {
                        diffs += 1;
                        denom += 1;
                        continue;
                    }
                    if 3 * p + 3 <= vref.len() - ctl.heur.ref_v_trim {
                        denom += 1;
                        if aa_seq[p] != codon_to_aa(&vref[3 * p..3 * p + 3]) {
                            diffs += 1;
                        }
                    }
                    if 3 * p > z - (jref.len() - ctl.heur.ref_j_trim) + 3 {
                        denom += 1;
                        if aa_seq[p]
                            != codon_to_aa(
                                &jref[jref.len() - (z - 3 * p)..jref.len() - (z - 3 * p) + 3],
                            )
                        {
                            diffs += 1;
                        }
                    }
                }
                cvar![
                    j,
                    *var,
                    format!("{:.1}", percent_ratio(denom - diffs, denom))
                ];
            } else if *var == "dna%".to_string() {
                let xm = &ex.share[mid];
                let mut diffs = 0;
                let mut denom = 0;
                let seq = &xm.seq_del_amino;
                let mut vref = refdata.refs[xm.v_ref_id].to_ascii_vec();
                if xm.v_ref_id_donor_alt_id.is_some() {
                    vref = dref[xm.v_ref_id_donor.unwrap()].nt_sequence.clone();
                }
                let jref = refdata.refs[xm.j_ref_id].to_ascii_vec();
                let z = seq.len();
                for p in 0..z {
                    let b = seq[p];
                    if b == b'-' {
                        diffs += 1;
                        denom += 1;
                        continue;
                    }
                    if p < vref.len() - ctl.heur.ref_v_trim {
                        denom += 1;
                        if b != vref[p] {
                            diffs += 1;
                        }
                    }
                    if p >= z - (jref.len() - ctl.heur.ref_j_trim) {
                        denom += 1;
                        if b != jref[jref.len() - (z - p)] {
                            diffs += 1;
                        }
                    }
                }
                cvar![
                    j,
                    *var,
                    format!("{:.1}", percent_ratio(denom - diffs, denom))
                ];
            } else if *var == "vjlen".to_string() {
                cvar![
                    j,
                    var,
                    format!("{}", ex.share[mid].j_stop - ex.share[mid].v_start)
                ];
            } else if var.starts_with("ndiff") {
                let u0 = var.between("ndiff", "vj").force_usize() - 1;
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
                    cvar![j, *var, format!("{}", ndiff)];
                } else {
                    cvar![j, *var, "_".to_string()];
                }
            } else if *var == "cdiff".to_string() {
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
                        cdiff += &format!("+{}", extra);
                    }
                } else if clen > 0 {
                    cdiff = format!("+{}", clen);
                }
                cvar![j, var, cdiff];
            } else if *var == "udiff".to_string() {
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
                        udiff += &format!("+{}", extra);
                    }
                } else if ulen > 0 {
                    udiff = format!("+{}", ulen);
                }
                cvar![j, var, udiff];
            } else if *var == "d_univ".to_string() {
                let vid = ex.share[mid].v_ref_id;
                let vref = &refdata.refs[vid].to_ascii_vec();
                let jid = ex.share[mid].j_ref_id;
                let jref = &refdata.refs[jid].to_ascii_vec();
                let tig = &ex.share[mid].seq_del;
                let n = tig.len();
                let mut diffs = 0;
                for p in 0..n {
                    if tig[p] == b'-' {
                        continue;
                    }
                    if p < vref.len() - ctl.heur.ref_v_trim && tig[p] != vref[p] {
                        diffs += 1;
                    } else if p >= n - (jref.len() - ctl.heur.ref_j_trim)
                        && tig[p] != jref[jref.len() - (n - p)]
                    {
                        diffs += 1;
                    }
                }
                cvar![j, var, format!("{}", diffs)];
            } else if *var == "d_donor".to_string() {
                let vid = ex.share[mid].v_ref_id;
                let mut vref = refdata.refs[vid].to_ascii_vec();
                if rsi.vpids[col].is_some() {
                    vref = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
                }
                let jid = ex.share[mid].j_ref_id;
                let jref = &refdata.refs[jid].to_ascii_vec();
                let tig = &ex.share[mid].seq_del;
                let n = tig.len();
                let mut diffs = 0;
                for p in 0..n {
                    if tig[p] == b'-' {
                        continue;
                    }
                    if p < vref.len() - ctl.heur.ref_v_trim && tig[p] != vref[p] {
                        diffs += 1;
                    } else if p >= n - (jref.len() - ctl.heur.ref_j_trim)
                        && tig[p] != jref[jref.len() - (n - p)]
                    {
                        diffs += 1;
                    }
                }
                cvar![j, var, format!("{}", diffs)];
            } else if *var == "notes".to_string() {
                cvar![j, var, ex.share[mid].vs_notesx.clone()];
            } else if *var == "var".to_string() {
                cvar![j, var, stringme(&varmat[u][col])];
            } else if *var == "u".to_string() {
                cvar![j, var, format!("{}", median_numis)];
            } else if *var == "u_cell".to_string() {
                let var = var.clone();
                if col + 1 <= ctl.parseable_opt.pchains {
                    let varc = format!("{}{}", var, col + 1);
                    if pcols_sort.is_empty() || bin_member(&pcols_sort, &varc) {
                        let mut vals = String::new();
                        for k in 0..ex.ncells() {
                            if k > 0 {
                                vals += ";";
                            }
                            vals += &format!("{}", ex.clones[k][mid].umi_count);
                        }
                        out_data[u].insert(varc, format!("{}", vals));
                    }
                }
            } else if *var == "u_min".to_string() {
                cvar![j, var, format!("{}", u_min)];
            } else if *var == "u_max".to_string() {
                cvar![j, var, format!("{}", u_max)];
            } else if *var == "u_μ".to_string() {
                cvar![j, var, format!("{}", u_mean)];
            } else if *var == "u_Σ".to_string() {
                cvar![j, var, format!("{}", utot)];
            } else if *var == "r".to_string() {
                cvar![j, var, format!("{}", median_nreads)];
            } else if *var == "r_min".to_string() {
                cvar![j, var, format!("{}", r_min)];
            } else if *var == "r_max".to_string() {
                cvar![j, var, format!("{}", r_max)];
            } else if *var == "r_μ".to_string() {
                cvar![j, var, format!("{}", r_mean)];
            } else if *var == "r_Σ".to_string() {
                cvar![j, var, format!("{}", rtot)];
            } else if *var == "r_cell".to_string() {
                let var = var.clone();
                if col + 1 <= ctl.parseable_opt.pchains {
                    let varc = format!("{}{}", var, col + 1);
                    if pcols_sort.is_empty() || bin_member(&pcols_sort, &varc) {
                        let mut vals = String::new();
                        for k in 0..ex.ncells() {
                            if k > 0 {
                                vals += ";";
                            }
                            vals += &format!("{}", ex.clones[k][mid].read_count);
                        }
                        out_data[u].insert(varc, format!("{}", vals));
                    }
                }
            } else if *var == "const".to_string() {
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
                cvar![j, var, format!("{}", constx.iter().format(","))];

            // Compute potential whitelist contamination percent and filter.
            // This is an undocumented option.
            } else if *var == "white".to_string() || ctl.clono_filt_opt.whitef {
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
                if *var == "white".to_string() && col_var {
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
