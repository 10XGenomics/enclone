// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function proc_lvar,
// plus a small helper function get_gex_matrix_entry.

use enclone_core::defs::*;
use enclone_core::median::*;
use enclone_proto::types::*;
use itertools::*;
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

pub fn proc_lvar1(
    i: usize,
    x: &String,
    pass: usize,
    u: usize,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    mults: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    gex_info: &GexInfo,
    _refdata: &RefData,
    _varmat: &Vec<Vec<Vec<u8>>>,
    _fp: &Vec<Vec<usize>>,
    row: &mut Vec<String>,
    out_data: &mut Vec<HashMap<String, String>>,
    _d_all: &mut Vec<Vec<u32>>,
    _ind_all: &mut Vec<Vec<u32>>,
    rsi: &ColInfo,
    _dref: &Vec<DonorReferenceItem>,
    groups: &HashMap<usize, Vec<usize>>,
    stats: &mut Vec<(String, Vec<String>)>,
    vdj_cells: &Vec<Vec<String>>,
    n_vdj_gex: &Vec<usize>,
    nd_fields: &Vec<String>,
    lvars: &Vec<String>,
    lenas: &Vec<String>,
    alt_bcs: &Vec<String>,
    _n_gex: usize,
    _n_gexs: &Vec<usize>,
    _gex_min: usize,
    _gex_max: usize,
    _gex_mean: f64,
    _gex_sum: f64,
    _gex_median: usize,
    _count_unsorted: &Vec<usize>,
    _entropy: f64,
    _entropies_unsorted: &Vec<f64>,
    _fcounts: &Vec<f64>,
    extra_args: &Vec<String>,
    fate: &Vec<HashMap<String, String>>,
) -> bool {
    let clonotype_id = exacts[u];
    let ex = &exact_clonotypes[clonotype_id];
    let verbose = ctl.gen_opt.row_fill_verbose;

    // Set up speak macro.

    macro_rules! speak {
        ($u:expr, $var:expr, $val:expr) => {
            if pass == 2 && (ctl.parseable_opt.pout.len() > 0 || extra_args.len() > 0) {
                let mut v = $var.to_string();
                v = v.replace("_Σ", "_sum");
                v = v.replace("_μ", "_mean");
                if ctl.parseable_opt.pcols.is_empty()
                    || bin_member(&ctl.parseable_opt.pcols_sortx, &v)
                    || bin_member(&extra_args, &v)
                {
                    out_data[$u].insert(v, $val);
                }
            }
        };
    }

    // Set up lead variable macros.  This is the mechanism for generating
    // both human-readable and parseable output for lead variables.

    macro_rules! lvar {
        ($i: expr, $var:expr, $val:expr) => {
            if verbose {
                eprint!("lvar {} ==> {}; ", $var, $val);
                eprintln!("$i = {}, lvars.len() = {}", $i, lvars.len());
            }
            if $i < lvars.len() {
                row.push($val)
            }
            if pass == 2 {
                speak!(u, $var.to_string(), $val);
            }
        };
    }
    macro_rules! lvar_stats1 {
        ($i: expr, $var:expr, $val:expr) => {
            if verbose {
                eprint!("lvar {} ==> {}; ", $var, $val);
                eprintln!("$i = {}, lvars.len() = {}", $i, lvars.len());
            }
            if $i < lvars.len() {
                row.push($val)
            }
            if pass == 2 {
                speak!(u, $var.to_string(), $val);
            }
            stats.push(($var.to_string(), vec![$val; ex.ncells()]));
        };
    }
    macro_rules! lvar_stats {
        ($i: expr, $var:expr, $val:expr, $stats: expr) => {
            if verbose {
                eprint!("lvar {} ==> {}; ", $var, $val);
                eprintln!("$i = {}, lvars.len() = {}", $i, lvars.len());
            }
            if $i < lvars.len() {
                row.push($val)
            }
            if pass == 2 {
                speak!(u, $var.to_string(), $val);
            }
            stats.push(($var.to_string(), $stats.clone()));
        };
    }

    // Check INFO.

    if ctl.gen_opt.info.is_some() {
        for q in 0..ctl.gen_opt.info_fields.len() {
            if *x == ctl.gen_opt.info_fields[q] {
                let mut found = false;
                let mut lvarred = false;
                if ex.share.len() == 2 && ex.share[0].left != ex.share[1].left {
                    let mut tag = String::new();
                    for j in 0..ex.share.len() {
                        if ex.share[j].left {
                            tag += &strme(&ex.share[j].seq);
                        }
                    }
                    tag += "_";
                    for j in 0..ex.share.len() {
                        if !ex.share[j].left {
                            tag += &strme(&ex.share[j].seq);
                        }
                    }
                    if ctl.gen_opt.info_data.contains_key(&tag) {
                        let val = &ctl.gen_opt.info_data[&tag][q];
                        lvar_stats1![i, x, val.clone()];
                        lvarred = true;
                        found = true;
                    }
                }
                if !lvarred {
                    lvar_stats1![i, x, String::new()];
                }
                if !found {
                    stats.push((x.to_string(), vec![]));
                }
                return true;
            }
        }
    }

    // Proceed.

    if x.starts_with('g') && x.after("g").parse::<usize>().is_ok() {
        let d = x.after("g").force_usize();
        if groups.contains_key(&d) {
            lvar_stats1![i, x, format!("{}", groups[&d][u] + 1)];
            return true;
        }
    }

    if x == "origins" {
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
        let origins_unsorted = origins.clone();
        unique_sort(&mut origins);
        lvar_stats![
            i,
            x,
            format!("{}", origins.iter().format(",")),
            origins_unsorted
        ];
    } else if x == "origins_cell" {
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
        if pass == 2 {
            speak!(u, x, format!("{}", origins.iter().format(POUT_SEP)));
        }
    } else if x == "clonotype_ncells" {
        let mut n = 0;
        for u in exacts.iter() {
            n += exact_clonotypes[*u].ncells();
        }
        lvar_stats1![i, x, format!("{}", n)];
    } else if x == "nchains" {
        lvar_stats1![i, x, format!("{}", rsi.mat.len())];
    } else if x == "nchains_present" {
        lvar_stats1![i, x, format!("{}", exact_clonotypes[exacts[u]].share.len())];
    } else if x == "datasets" {
        lvar_stats1![i, x, format!("{}", lenas.iter().format(","))];
    } else if x == "datasets_cell" {
        let mut datasets = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            datasets.push(ctl.origin_info.dataset_id[ex.clones[j][0].dataset_index].clone());
        }
        if pass == 2 {
            speak!(u, x, format!("{}", datasets.iter().format(POUT_SEP)));
        }
        stats.push((x.to_string(), datasets));
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
        let donors_unsorted = donors.clone();
        unique_sort(&mut donors);
        lvar_stats![
            i,
            x,
            format!("{}", donors.iter().format(",")),
            donors_unsorted
        ];
    } else if x == "donors_cell" {
        let mut donors = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            if ex.clones[j][0].donor_index.is_some() {
                donors
                    .push(ctl.origin_info.donor_list[ex.clones[j][0].donor_index.unwrap()].clone());
            } else {
                donors.push("?".to_string());
            }
        }
        if pass == 2 {
            speak!(u, x, format!("{}", donors.iter().format(POUT_SEP)));
        }
    } else if x == "n" {
        let counts = vec!["1.0".to_string(); mults[u]];
        lvar_stats![i, x, format!("{}", mults[u]), counts];
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
        let mut clustf = Vec::<String>::new();
        for x in clust.iter() {
            clustf.push(format!("{}", x));
        }
        clust.sort();
        lvar_stats![i, x, format!("{}", abbrev_list(&clust)), clustf];
    } else if x == "n_other" {
        let mut n = 0;
        let mut ns = Vec::<String>::new();
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
                ns.push("1.0".to_string());
            } else {
                ns.push("0.0".to_string());
            }
        }
        lvar_stats![i, x, format!("{}", n), ns];
    } else if x == "n_b" {
        let mut n_b = 0;
        let mut ns = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            let bc = &ex.clones[j][0].barcode;
            let li = ex.clones[j][0].dataset_index;
            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                    n_b += 1;
                    ns.push("1.0".to_string());
                } else {
                    ns.push("0.0".to_string());
                }
            }
        }
        lvar_stats![i, x, format!("{}", n_b), ns];
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
        let mut fates = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            let mut f = String::new();
            let bc = &ex.clones[j][0].barcode;
            let li = ex.clones[j][0].dataset_index;
            if fate[li].contains_key(&bc.clone()) {
                f = fate[li][&bc.clone()].clone();
                f = f.between(" ", " ").to_string();
            }
            fates.push(f);
        }
        lvar_stats![i, x, String::new(), fates];
        if pass == 2 {
            speak!(u, x, format!("{}", fates.iter().format(POUT_SEP)));
        }
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
        lvar_stats1![i, x, s.clone()];
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
        lvar_stats1![i, x, s.clone()];
    } else if x.starts_with("pe") {
        lvar_stats1![i, x, format!("")];
    } else if x.starts_with("npe") {
        lvar_stats1![i, x, format!("")];
    } else if x.starts_with("ppe") {
        lvar_stats1![i, x, format!("")];
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
        let credsx_unsorted = credsx.clone();
        credsx.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if x == "cred" {
            if credsx.is_empty() {
                lvar![i, x, format!("")];
            } else {
                lvar_stats1![i, x, format!("{:.1}", median_f64(&credsx))];
            }
        } else {
            let mut r = Vec::<String>::new();
            for j in 0..credsx_unsorted.len() {
                r.push(format!("{}", credsx_unsorted[j]));
            }
            stats.push((x.to_string(), r));
            if pass == 2 {
                let mut r = Vec::<String>::new();
                for j in 0..credsx_unsorted.len() {
                    r.push(format!("{:.1}", credsx_unsorted[j]));
                }
                speak!(u, x, format!("{}", r.iter().format(POUT_SEP)));
            }
        }
    } else if bin_member(&alt_bcs, x) {
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
        lvar_stats![i, x, format!(""), r];
        if pass == 2 {
            speak!(u, x, format!("{}", r.iter().format(POUT_SEP)));
        }
    } else if x.starts_with("n_") && !x.starts_with("n_gex") {
        let name = x.after("n_");
        let mut count = 0;
        let mut counts = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            let x = &ex.clones[j][0];
            if ctl.origin_info.dataset_id[x.dataset_index] == name {
                count += 1;
                counts.push("1.0".to_string());
            } else if x.origin_index.is_some()
                && ctl.origin_info.origin_list[x.origin_index.unwrap()] == name
            {
                count += 1;
                counts.push("1.0".to_string());
            } else if x.donor_index.is_some()
                && ctl.origin_info.donor_list[x.donor_index.unwrap()] == name
            {
                count += 1;
                counts.push("1.0".to_string());
            } else if x.tag_index.is_some()
                && ctl.origin_info.tag_list[x.tag_index.unwrap()] == name
            {
                count += 1;
                counts.push("1.0".to_string());
            }
        }
        lvar_stats![i, x, format!("{}", count), counts];
    } else {
        return false;
    }
    return true;
}
