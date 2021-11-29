// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
// This file is auto-generated by the crate enclone_vars, please do not edit.

use amino::*;
// use crate::print_utils1::*;
// use crate::print_utils3::*;
// use enclone_core::align_to_vdj_ref::*;
use enclone_core::defs::*;
use enclone_core::median::*;
// use enclone_core::opt_d::*;
use enclone_proto::types::*;
use itertools::Itertools;
// use stats_utils::*;
use std::cmp::{max, min};
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::RefData;
use vector_utils::*;

pub fn proc_lvar_auto(
    i: usize,
    pass: usize,
    var: &String,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    u: usize,
    rsi: &ColInfo,
    refdata: &RefData,
    ctl: &EncloneControl,
    extra_args: &Vec<String>,
    out_data: &mut Vec<HashMap<String, String>>,
    stats: &mut Vec<(String, Vec<String>)>,
    lvars: &Vec<String>,
    row: &mut Vec<String>,
    fate: &Vec<HashMap<String, String>>,
    dref: &Vec<DonorReferenceItem>,
    varmat: &Vec<Vec<Vec<u8>>>,
    fp: &Vec<Vec<usize>>,
    n_vdj_gex: &Vec<usize>,
    vdj_cells: &Vec<Vec<String>>,
    gex_info: &GexInfo,
    groups: &HashMap<usize, Vec<usize>>,
    mults: &Vec<usize>,
    nd_fields: &Vec<String>,
) -> Result<bool, String> {
    let clonotype_id = exacts[u];
    let ex = &exact_clonotypes[clonotype_id];
    let mat = &rsi.mat;
    let cols = varmat[0].len();
    let verbose = ctl.gen_opt.row_fill_verbose;

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

    let val = if false {
        (String::new(), Vec::<String>::new(), String::new())
    } else if var == "clonotype_ncells" {
        let mut n = 0;
        for u in exacts.iter() {
            n += exact_clonotypes[*u].ncells();
        }

        (format!("{}", n), Vec::new(), "clono".to_string())
    } else if var == "clust" {
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
        clust.sort_unstable();

        (abbrev_list(&clust), clustf, "cell-exect".to_string())
    } else if var == "cred" {
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
        let mut r = Vec::<String>::new();
        for j in 0..credsx_unsorted.len() {
            r.push(format!("{:.1}", credsx_unsorted[j]));
        }

        (
            format!("{:.1}", median_f64(&credsx)),
            r,
            "cell-exact".to_string(),
        )
    } else if var == "cred_cell" {
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
        let mut r = Vec::<String>::new();
        for j in 0..credsx_unsorted.len() {
            r.push(format!("{:.1}", credsx_unsorted[j]));
        }

        let _exact = format!("{:.1}", median_f64(&credsx));
        (String::new(), r, "cell-exact".to_string())
    } else if var == "datasets" {
        let mut datasets = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            datasets.push(ctl.origin_info.dataset_id[ex.clones[j][0].dataset_index].clone());
        }
        let mut datasets_unique = datasets.clone();
        unique_sort(&mut datasets_unique);

        (
            format!("{}", datasets_unique.iter().format(",")),
            datasets,
            "cell-exact".to_string(),
        )
    } else if var == "datasets_cell" {
        let mut datasets = Vec::<String>::new();
        for j in 0..ex.clones.len() {
            datasets.push(ctl.origin_info.dataset_id[ex.clones[j][0].dataset_index].clone());
        }
        let mut datasets_unique = datasets.clone();
        unique_sort(&mut datasets_unique);

        let _exact = format!("{}", datasets_unique.iter().format(","));
        (String::new(), datasets, "cell-exact".to_string())
    } else if var == "donors" {
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

        (
            format!("{}", donors.iter().format(",")),
            donors_unsorted,
            "cell-exact".to_string(),
        )
    } else if var == "donors_cell" {
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

        let _exact = format!("{}", donors.iter().format(","));
        (String::new(), donors_unsorted, "cell-exact".to_string())
    } else if var == "dref" {
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

        (format!("{}", diffs), Vec::new(), "exact".to_string())
    } else if var == "dref_aa" {
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

        (format!("{}", diffs), Vec::new(), "exact".to_string())
    } else if var == "far" {
        let mut dist = -1_isize;
        for i2 in 0..varmat.len() {
            if i2 == u || fp[i2] != fp[u] {
                continue;
            }
            let mut d = 0_isize;
            for c in fp[u].iter() {
                for j in 0..varmat[u][*c].len() {
                    if varmat[u][*c][j] != varmat[i2][*c][j] {
                        d += 1;
                    }
                }
            }
            dist = max(dist, d);
        }
        let d;
        if dist == -1_isize {
            d = "".to_string();
        } else {
            d = format!("{}", dist);
        }

        (d, Vec::new(), "exact".to_string())
    } else if var.starts_with("fb")
        && var.ends_with("")
        && var.between2("fb", "").parse::<i64>().is_ok()
        && var.between2("fb", "").force_i64() >= 1
    {
        let arg1 = var.between2("fb", "").force_i64();
        let ncols = gex_info.fb_top_matrices[0].ncols();
        let n = (arg1 - 1) as usize;
        let fb = if n < ncols {
            gex_info.fb_top_matrices[0].col_label(n)
        } else {
            String::new()
        };

        ((*fb).to_string(), Vec::new(), "exact".to_string())
    } else if var.starts_with("fb")
        && var.ends_with("_n")
        && var.between2("fb", "_n").parse::<i64>().is_ok()
        && var.between2("fb", "_n").force_i64() >= 1
    {
        let arg1 = var.between2("fb", "_n").force_i64();
        let ncols = gex_info.fb_top_matrices[0].ncols();
        let n = (arg1 - 1) as usize;
        let median;
        let mut counts;
        if n >= ncols {
            median = 0;
            counts = vec!["0".to_string(); ex.ncells()];
        } else {
            counts = Vec::<String>::new();
            let mut counts_sorted = Vec::<usize>::new();
            for l in 0..ex.clones.len() {
                let bc = ex.clones[l][0].barcode.clone();
                let p = bin_position(&gex_info.fb_top_barcodes[0], &bc);
                if p < 0 {
                    counts.push("0".to_string());
                    counts_sorted.push(0);
                } else {
                    let x = gex_info.fb_top_matrices[0].value(p as usize, n);
                    counts.push(format!("{}", x));
                    counts_sorted.push(x);
                }
            }
            counts_sorted.sort_unstable();
            median = rounded_median(&counts_sorted);
        }

        (format!("{}", median), counts, "cell-exact".to_string())
    } else if var.starts_with("fb")
        && var.ends_with("_n_cell")
        && var.between2("fb", "_n_cell").parse::<i64>().is_ok()
        && var.between2("fb", "_n_cell").force_i64() >= 1
    {
        let arg1 = var.between2("fb", "_n_cell").force_i64();
        let ncols = gex_info.fb_top_matrices[0].ncols();
        let n = (arg1 - 1) as usize;
        let median;
        let mut counts;
        if n >= ncols {
            median = 0;
            counts = vec!["0".to_string(); ex.ncells()];
        } else {
            counts = Vec::<String>::new();
            let mut counts_sorted = Vec::<usize>::new();
            for l in 0..ex.clones.len() {
                let bc = ex.clones[l][0].barcode.clone();
                let p = bin_position(&gex_info.fb_top_barcodes[0], &bc);
                if p < 0 {
                    counts.push("0".to_string());
                    counts_sorted.push(0);
                } else {
                    let x = gex_info.fb_top_matrices[0].value(p as usize, n);
                    counts.push(format!("{}", x));
                    counts_sorted.push(x);
                }
            }
            counts_sorted.sort_unstable();
            median = rounded_median(&counts_sorted);
        }

        let _exact = format!("{}", median);
        (String::new(), counts, "cell-exact".to_string())
    } else if var == "filter" {
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

        (String::new(), fates, "cell".to_string())
    } else if var.starts_with("g")
        && var.ends_with("")
        && var.between2("g", "").parse::<i64>().is_ok()
        && var.between2("g", "").force_i64() >= 0
    {
        let arg1 = var.between2("g", "").force_i64();
        let d = arg1 as usize;
        let answer = if groups.contains_key(&d) {
            format!("{}", groups[&d][u] + 1)
        } else {
            String::new()
        };

        (answer, Vec::new(), "exact".to_string())
    } else if var == "inkt" {
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

        (s, Vec::new(), "exact".to_string())
    } else if var == "mait" {
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

        (s, Vec::new(), "exact".to_string())
    } else if var == "mem" {
        let mut n = 0;
        let mut y = Vec::<String>::new();
        if ctl.gen_opt.using_secmem {
            for l in 0..ex.clones.len() {
                let li = ex.clones[l][0].dataset_index;
                let bc = &ex.clones[l][0].barcode;
                let mut count = 0;
                if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                    count = ctl.origin_info.secmem[li][&bc.clone()].1;
                    n += count;
                }
                y.push(format!("{}", count));
            }
        }

        (format!("{}", n), y, "cell-exact".to_string())
    } else if var == "mem_cell" {
        let mut n = 0;
        let mut y = Vec::<String>::new();
        if ctl.gen_opt.using_secmem {
            for l in 0..ex.clones.len() {
                let li = ex.clones[l][0].dataset_index;
                let bc = &ex.clones[l][0].barcode;
                let mut count = 0;
                if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                    count = ctl.origin_info.secmem[li][&bc.clone()].1;
                    n += count;
                }
                y.push(format!("{}", count));
            }
        }

        let _exact = format!("{}", n);
        (String::new(), y, "cell-exact".to_string())
    } else if var == "n" {
        let counts = vec!["1.0".to_string(); mults[u]];

        (format!("{}", mults[u]), counts, "cell-exact".to_string())
    } else if var == "n_cell" {
        let counts = vec!["1.0".to_string(); mults[u]];

        let _exact = format!("{}", mults[u]);
        (String::new(), counts, "cell-exact".to_string())
    } else if var == "n_b" {
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

        (format!("{}", n_b), ns, "cell-exact".to_string())
    } else if var == "n_b_cell" {
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

        let _exact = format!("{}", n_b);
        (String::new(), ns, "cell-exact".to_string())
    } else if var == "n_other" {
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

        (format!("{}", n), ns, "cell-exact".to_string())
    } else if var == "n_other_cell" {
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

        let _exact = format!("{}", n);
        (String::new(), ns, "cell-exact".to_string())
    } else if var == "nchains" {
        (
            format!("{}", rsi.mat.len()),
            Vec::new(),
            "clono".to_string(),
        )
    } else if var == "nchains_present" {
        (
            format!("{}", exact_clonotypes[exacts[u]].share.len()),
            Vec::new(),
            "exact".to_string(),
        )
    } else if var == "near" {
        let near;
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
            near = "".to_string()
        } else {
            near = format!("{}", dist)
        }

        (near, Vec::new(), "exact".to_string())
    } else if var == "origins" {
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

        (
            format!("{}", origins.iter().format(",")),
            origins_unsorted,
            "cell-exact".to_string(),
        )
    } else if var == "origins_cell" {
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

        let _exact = format!("{}", origins.iter().format(","));
        (String::new(), origins_unsorted, "cell-exact".to_string())
    } else if var == "sec" {
        let mut n = 0;
        let mut y = Vec::<String>::new();
        if ctl.gen_opt.using_secmem {
            for l in 0..ex.clones.len() {
                let li = ex.clones[l][0].dataset_index;
                let bc = &ex.clones[l][0].barcode;
                let mut count = 0;
                if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                    count = ctl.origin_info.secmem[li][&bc.clone()].0;
                    n += count;
                }
                y.push(format!("{}", count));
            }
        }

        (format!("{}", n), y, "cell-exact".to_string())
    } else if var == "sec_cell" {
        let mut n = 0;
        let mut y = Vec::<String>::new();
        if ctl.gen_opt.using_secmem {
            for l in 0..ex.clones.len() {
                let li = ex.clones[l][0].dataset_index;
                let bc = &ex.clones[l][0].barcode;
                let mut count = 0;
                if ctl.origin_info.secmem[li].contains_key(&bc.clone()) {
                    count = ctl.origin_info.secmem[li][&bc.clone()].0;
                    n += count;
                }
                y.push(format!("{}", count));
            }
        }

        let _exact = format!("{}", n);
        (String::new(), y, "cell-exact".to_string())
    } else if var == "type" {
        let mut cell_types = Vec::<String>::new(); /*
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

        (abbrev_list(&cell_types), Vec::new(), "exact".to_string())
    } else {
        (
            "$UNDEFINED".to_string(),
            Vec::<String>::new(),
            String::new(),
        )
    };
    if val.0 == "$UNDEFINED" {
        return Ok(false);
    } else {
        let (exact, cell, level) = &val;
        if level == "cell" && !var.ends_with("_cell") {
            lvar_stats![i, var, String::new(), cell];
            if pass == 2 {
                speak!(u, var, format!("{}", cell.iter().format(POUT_SEP)));
            }
        } else if (exact.len() > 0 && !var.ends_with("_cell")) || cell.len() == 0 {
            if verbose {
                eprint!("lvar {} ==> {}; ", var, exact);
                eprintln!("i = {}, lvars.len() = {}", i, lvars.len());
            }
            if i < lvars.len() {
                row.push(exact.clone())
            }
            if pass == 2 {
                speak!(u, var.to_string(), exact.to_string());
            }
            if cell.len() == 0 {
                stats.push((var.to_string(), vec![exact.to_string(); ex.ncells()]));
            } else {
                stats.push((var.to_string(), cell.to_vec()));
            }
        } else if cell.len() > 0 {
            if pass == 2 {
                speak!(u, var, format!("{}", cell.iter().format(POUT_SEP)));
            }
            stats.push((var.to_string(), cell.to_vec()));
        }
        return Ok(true);
    }
}
