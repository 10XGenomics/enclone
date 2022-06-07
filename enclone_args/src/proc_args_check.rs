// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Check lvars, cvars, and pcols.

use enclone_core::allowed_vars::{
    CVARS_ALLOWED, CVARS_ALLOWED_PCELL, GVARS_ALLOWED, LVARS_ALLOWED, PCVARS_ALLOWED,
    PLVARS_ALLOWED,
};
use enclone_core::defs::{EncloneControl, GexInfo};
use itertools::Itertools;
use rayon::prelude::*;
use regex::Regex;
use std::time::Instant;
use string_utils::{strme, TextUtils};
use vector_utils::{bin_member, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Get known features.  This code is inefficient.

pub fn get_known_features(gex_info: &GexInfo) -> Result<Vec<String>, String> {
    let mut known_features = Vec::<String>::new();
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let suffixes_g = ["", "_min", "_max", "_μ", "_Σ", "_%"];
    let mut results = Vec::<(usize, Vec<String>, String)>::new();
    for i in 0..gex_info.gex_features.len() {
        results.push((i, Vec::<String>::new(), String::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        for j in 0..gex_info.gex_features[i].len() {
            let f = &gex_info.gex_features[i][j];
            let ff = f.split('\t').collect::<Vec<&str>>();
            if ff.len() != 3 {
                res.2 = format!(
                    "\nUnexpected structure of features file, at this line\n{}\n\
                    Giving up.\n",
                    f
                );
                return;
            }
            for z in 0..2 {
                if ff[2].starts_with("Antibody") {
                    for s in suffixes.iter() {
                        res.1.push(format!("{}_ab{}", ff[z], s));
                    }
                } else if ff[2].starts_with("CRISPR") {
                    for s in suffixes.iter() {
                        res.1.push(format!("{}_cr{}", ff[z], s));
                    }
                } else if ff[2].starts_with("CUSTOM") {
                    for s in suffixes.iter() {
                        res.1.push(format!("{}_cu{}", ff[z], s));
                    }
                } else if ff[2].starts_with("Antigen") {
                    for s in suffixes.iter() {
                        res.1.push(format!("{}_ag{}", ff[z], s));
                    }
                } else {
                    for s in suffixes_g.iter() {
                        res.1.push(format!("{}_g{}", ff[z], s));
                    }
                }
            }
        }
    });
    for i in 0..results.len() {
        if !results[i].2.is_empty() {
            return Err(results[i].2.clone());
        }
    }
    for i in 0..results.len() {
        known_features.append(&mut results[i].1.clone());
    }
    known_features.par_sort();
    known_features.dedup();
    Ok(known_features)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn involves_gex_fb(x: &String) -> bool {
    let ends0 = [
        "_g", "_ab", "_ag", "_cr", "_cu", "_g_μ", "_ab_μ", "_ag_μ", "_cr_μ", "_cu_μ", "_g_%",
    ];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let mut ends = Vec::<String>::new();
    for z in ends0.iter() {
        for y in suffixes.iter() {
            ends.push(format!("{}{}", z, y));
        }
    }
    let mut x = x.clone();
    if x.contains(':') {
        x = x.rev_after(":").to_string();
    }
    if x.ends_with("_cell") {
        x = x.rev_before("_cell").to_string();
    }
    for y in ends.iter() {
        if x.ends_with(y) {
            return true;
        }
    }
    x == "gex"
        || x.starts_with("gex_")
        || x == "n_gex"
        || x == "clust"
        || x == "type"
        || x == "entropy"
        || x == "cred"
        || x == "cred_cell"
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn is_pattern(x: &String, parseable: bool) -> bool {
    let ends0 = [
        "_g", "_ab", "_ag", "_cr", "_cu", "_g_μ", "_ab_μ", "_ag_μ", "_cr_μ", "_cu_μ", "_g_%",
    ];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let mut ends = Vec::<String>::new();
    for z in ends0.iter() {
        for y in suffixes.iter() {
            ends.push(format!("{}{}", z, y));
        }
    }
    let mut x = x.clone();
    if x.contains(':') {
        x = x.rev_after(":").to_string();
    }
    if parseable && x.ends_with("_cell") {
        x = x.rev_before("_cell").to_string();
    }
    let mut pat = false;
    for y in ends.iter() {
        if x.ends_with(y) {
            let p = x.rev_before(y);
            if !p.is_empty() && Regex::new(p).is_ok() {
                let mut ok = true;
                let mut special = false;
                let p = p.as_bytes();
                for i in 0..p.len() {
                    if !((p[i] >= b'A' && p[i] <= b'Z')
                        || (p[i] >= b'a' && p[i] <= b'z')
                        || (p[i] >= b'0' && p[i] <= b'9')
                        || b".-_[]()|*".contains(&p[i]))
                    {
                        ok = false;
                        break;
                    }
                    if b"[]()|*".contains(&p[i]) {
                        special = true;
                    }
                }
                if ok && special {
                    pat = true;
                    break;
                }
            }
        }
    }
    pat
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn check_gene_fb(
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    to_check: &Vec<String>,
    category: &str,
) -> Result<(), String> {
    let g_ends0 = ["_g"];
    let fb_ends0 = ["_ab", "_cr", "_cu", "_ag"];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let suffixes_g = ["", "_min", "_max", "_μ", "_Σ", "_%"];
    let (mut g_ends, mut fb_ends) = (Vec::<String>::new(), Vec::<String>::new());
    for x in g_ends0.iter() {
        for y in suffixes_g.iter() {
            g_ends.push(format!("{}{}", x, y));
        }
    }
    for x in fb_ends0.iter() {
        for y in suffixes.iter() {
            fb_ends.push(format!("{}{}", x, y));
        }
    }
    for x in to_check.iter() {
        let mut x = x.to_string();
        if x.contains(':') {
            x = x.after(":").to_string();
        }
        if !gex_info.have_gex
            && !gex_info.have_fb
            && (*x == "n_gex".to_string() || *x == "n_gex_cell".to_string())
        {
            if category == "parseable" {
                return Err(format!(
                    "\nParseable field {} does not make sense because neither gene expression \
                     nor feature barcode data\nwere provided as input.\n",
                    x
                ));
            } else {
                return Err(format!(
                    "\nLead variable {} does not make sense because neither gene expression \
                     not feature barcode data\nwere provided as input.\n",
                    x
                ));
            }
        }
        if !gex_info.have_gex {
            let mut problem = false;
            for y in g_ends.iter() {
                if x.ends_with(y) {
                    problem = true;
                }
            }
            if problem
                || *x == "gex".to_string()
                || x.starts_with("gex_")
                || *x == "clust".to_string()
                || *x == "type".to_string()
                || *x == "entropy".to_string()
                || *x == "cred".to_string()
                || *x == "cred_cell".to_string()
            {
                if category == "parseable" {
                    return Err(format!(
                        "\nParseable field {} does not make sense because gene expression \
                         data\nwere not provided as input.\n",
                        x
                    ));
                } else {
                    return Err(format!(
                        "\nLead variable {} does not make sense because gene expression \
                         data\nwere not provided as input.\n",
                        x
                    ));
                }
            }
        }
        if !gex_info.have_fb {
            for y in fb_ends.iter() {
                if x.ends_with(y) {
                    if category == "parseable" {
                        return Err(format!(
                            "\nParseable field {} does not make sense because feature \
                             barcode data\nwere not provided as input.\n",
                            x
                        ));
                    } else {
                        return Err(format!(
                            "\nLead variable {} does not make sense because feature barcode \
                             data\nwere not provided as input.\n",
                            x
                        ));
                    }
                }
            }
        }
    }

    // Get known features.  This code is inefficient.

    let known_features = get_known_features(gex_info)?;

    // Do the check.

    for i in 0..to_check.len() {
        let mut x = to_check[i].clone();
        if x.contains(':') {
            x = x.after(":").to_string();
        }
        let mut y = x.clone();
        if category == "parseable" && y.ends_with("_cell") {
            y = y.before("_cell").to_string();
        }
        if !bin_member(&known_features, &y) {
            let mut n_var = false;
            if x.starts_with("n_") {
                n_var = true;
                let mut is_dataset_name = false;
                let mut is_origin_name = false;
                let mut is_donor_name = false;
                let mut is_tag_name = false;
                let name = x.after("n_").to_string();
                let s = ctl.origin_info.n();
                for j in 0..s {
                    if ctl.origin_info.dataset_id[j] == name {
                        is_dataset_name = true;
                    }
                }
                for j in 0..ctl.origin_info.origin_list.len() {
                    if ctl.origin_info.origin_list[j] == name {
                        is_origin_name = true;
                    }
                }
                for j in 0..ctl.origin_info.donor_list.len() {
                    if ctl.origin_info.donor_list[j] == name {
                        is_donor_name = true;
                    }
                }
                for j in 0..ctl.origin_info.tag_list.len() {
                    if ctl.origin_info.tag_list[j] == name {
                        is_tag_name = true;
                    }
                }
                let msg = "\nSuggested reading: \"enclone help input\" and \
                           \"enclone help glossary\".\n";
                if !is_dataset_name && !is_origin_name && !is_donor_name && !is_tag_name {
                    return Err(format!(
                        "\nYou've used the {} variable {}, and yet {} \
                         does not name a dataset, nor an origin,\nnor a donor, nor a tag.\n{}",
                        category, x, name, msg
                    ));
                }
                let mut types = 0;
                if is_dataset_name {
                    types += 1;
                }
                if is_origin_name {
                    types += 1;
                }
                if is_donor_name {
                    types += 1;
                }
                if is_tag_name {
                    types += 1;
                }
                if is_dataset_name && is_origin_name && is_donor_name {
                    return Err(format!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a dataset, an origin, and a donor.  That's ambiguous.\n{}",
                        category, x, name, msg
                    ));
                }
                if is_dataset_name && is_origin_name {
                    return Err(format!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a dataset and an origin.  That's ambiguous.\n{}",
                        category, x, name, msg
                    ));
                }
                if is_dataset_name && is_donor_name {
                    return Err(format!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a dataset and a donor.  That's ambiguous.\n{}",
                        category, x, name, msg
                    ));
                }
                if is_origin_name && is_donor_name {
                    return Err(format!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names an origin and a donor.  That's ambiguous.\n{}",
                        category, x, name, msg
                    ));
                }
                if types != 1 {
                    return Err(format!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a tag and also a dataset, origin or donor.\n\
                         That's ambiguous.\n{}",
                        category, x, name, msg
                    ));
                }
            }
            if !n_var {
                let mut alts = Vec::<String>::new();
                let mut xl = x.clone();
                xl.make_ascii_lowercase();
                for y in known_features.iter() {
                    let mut yl = y.to_string();
                    yl.make_ascii_lowercase();
                    if xl == yl {
                        alts.push(y.to_string());
                    }
                }
                if category == "lead" {
                    if x.is_empty() {
                        continue;
                    }
                    if !alts.is_empty() {
                        return Err(format!(
                            "\nThe variable {} for LVARS is unrecognized.  Might you have \
                            meant {}?\nPlease type \"enclone help lvars\".\n",
                            x,
                            alts.iter().format(" or "),
                        ));
                    }
                    return Err(format!(
                        "\nThe variable {} for LVARS is unrecognized.  Please type \
                         \"enclone help lvars\".\n",
                        x
                    ));
                } else {
                    if !alts.is_empty() {
                        return Err(format!(
                            "\nUnrecognized parseable variable {}.  Might you have meant {}?\n\
                            Please type \
                             \"enclone help parseable\".\nIf the variable is a chain variable \
                            (cvar), please make sure it is suffixed with the chain index.\n",
                            x,
                            alts.iter().format(" or "),
                        ));
                    }
                    return Err(format!(
                        "\nUnrecognized parseable variable {}.  Please type \
                         \"enclone help parseable\".\nIf the variable is a chain variable (cvar), \
                        please make sure it is suffixed with the chain index.\n",
                        x
                    ));
                }
            }
        }
    }
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check pcols args.

pub fn check_pcols(
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    cols: &Vec<String>,
    allow_cell: bool,
) -> Result<(), String> {
    let mut alt_bcs = Vec::<String>::new();
    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.clone());
        }
    }
    unique_sort(&mut alt_bcs);
    let mut to_check = Vec::<String>::new();
    let pchains = &ctl.parseable_opt.pchains;
    let ends = build_ends();
    let mut nd_used = false;
    for x in cols.iter() {
        let mut x = x.to_string();
        if x.contains(':') {
            x = x.after(":").to_string();
        }
        let mut ok = false;
        // Note that the following test is probably redundant with some of the testing below.
        if check_one_lvar(&*x, ctl, gex_info, &mut nd_used, &ends, false)? {
            ok = true;
        }
        for i in 0..ctl.gen_opt.info_fields.len() {
            if *x == ctl.gen_opt.info_fields[i] {
                ok = true;
            }
        }
        if bin_member(&alt_bcs, &x) {
            ok = true;
        }
        for y in ctl.clono_print_opt.lvars.iter() {
            if y.contains(':') {
                let y = y.before(":");
                if x == y {
                    ok = true;
                }
            }
        }
        for y in PLVARS_ALLOWED.iter() {
            if x == *y {
                ok = true;
            }
        }
        for y in ctl.origin_info.dataset_list.iter() {
            if *x == format!("{}_barcodes", y) {
                ok = true;
            }
        }
        if ctl.parseable_opt.pbarcode {
            if x == "barcode" {
                ok = true;
            }
            for y in ctl.origin_info.dataset_list.iter() {
                if *x == format!("{}_barcode", y) {
                    ok = true;
                }
            }
        }
        let gpvar = x.starts_with('g') && x.after("g").parse::<usize>().is_ok();

        if !gex_info.have_gex && !gex_info.have_fb && x.starts_with("n_gex") {
            return Err(format!(
                "\nCan't use parseable variable {} without having gene \
                 expression or feature barcode data.\n",
                x
            ));
        }
        if !gex_info.have_gex && (x.starts_with("gex") || x == "clust") || x == "type" {
            return Err(format!(
                "\nCan't use parseable variable {} without having gene \
                 expression data.\n",
                x
            ));
        }
        if LVARS_ALLOWED.contains(&x.as_str()) || gpvar {
            ok = true;
        } else if is_pattern(&x, true) {
            ok = true;
        } else {
            let mut y = Vec::<u8>::new();
            for c in x.chars().rev() {
                if c.is_ascii_digit() {
                    y.push(c as u8);
                } else {
                    break;
                }
            }
            y.reverse();
            let ps = strme(&y);
            if !ps.is_empty()
                && (pchains == "max"
                    || (ps.force_usize() > 0 && ps.force_usize() <= pchains.force_usize()))
            {
                let y = x.rev_before(ps);
                if CVARS_ALLOWED.contains(&y) || (allow_cell && CVARS_ALLOWED_PCELL.contains(&y)) {
                    ok = true;
                } else if PCVARS_ALLOWED.contains(&y) {
                    ok = true;
                } else if y.starts_with("ndiff")
                    && y.ends_with("vj")
                    && y.between("ndiff", "vj").parse::<usize>().is_ok()
                    && y.between("ndiff", "vj").force_usize() >= 1
                {
                    ok = true;
                } else if (y.starts_with("cdr1_aa_")
                    || y.starts_with("cdr2_aa_")
                    || y.starts_with("cdr3_aa_"))
                    && y.after("aa_").contains('_')
                    && y.between("aa_", "_").parse::<isize>().is_ok()
                    && y.after("aa_").after("_").ends_with("_ext")
                    && y.after("aa_").between("_", "_ext").parse::<isize>().is_ok()
                {
                    ok = true;
                }
            }
        }
        if !ok {
            to_check.push(x.to_string());
        }
    }
    if !to_check.is_empty() {
        check_gene_fb(ctl, gex_info, &to_check, "parseable")?;
    }
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check cvars args.

pub fn check_cvars(ctl: &EncloneControl) -> Result<(), String> {
    for x in ctl.clono_print_opt.cvars.iter() {
        let mut x = x.to_string();
        if x.contains(':') {
            x = x.after(":").to_string();
        }
        let mut ok = CVARS_ALLOWED.contains(&x.as_str());
        if x.starts_with("ndiff")
            && x.ends_with("vj")
            && x.between("ndiff", "vj").parse::<usize>().is_ok()
            && x.between("ndiff", "vj").force_usize() >= 1
        {
            ok = true;
        }
        if (x.starts_with("cdr1_aa_") || x.starts_with("cdr2_aa_") || x.starts_with("cdr3_aa_"))
            && x.after("aa_").contains('_')
            && x.between("aa_", "_").parse::<usize>().is_ok()
            && x.after("aa_").after("_").ends_with("_ext")
            && x.after("aa_").between("_", "_ext").parse::<usize>().is_ok()
        {
            ok = true;
        } else if x.starts_with('q')
            && x.ends_with('_')
            && x.after("q").rev_before("_").parse::<usize>().is_ok()
        {
            ok = true;
        }
        if !ok {
            return Err(format!(
                "\nUnrecognized variable {} for CVARS or CVARSP.  \
                 Please type \"enclone help cvars\".\n",
                x
            ));
        }
    }
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn check_one_lvar(
    x: &str,
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    nd_used: &mut bool,
    ends: &Vec<String>,
    is_lvar: bool,
) -> Result<bool, String> {
    for i in 0..ctl.gen_opt.info_fields.len() {
        if *x == ctl.gen_opt.info_fields[i] {
            return Ok(true);
        }
    }
    let mut x = x.to_string();
    if x.contains(':') {
        x = x.after(":").to_string();
    }

    // See if type is ok.

    if x == "type" {
        let mut specified = false;
        for i in 0..gex_info.cell_type_specified.len() {
            if gex_info.cell_type_specified[i] {
                specified = true;
            }
        }
        if !ctl.gen_opt.internal_run && !x.is_empty() {
            return Err(format!(
                "\nUnrecognized variable {} for LVARS or PCOLS.  Please type \
                 \"enclone help lvars\".\n",
                x
            ));
        }
        if !specified {
            return Err(
                "\nYou've used the lead or parseable variable \"type\", but the file \
                cell_types.csv was not found.\n\
                This could be because you're using a GEX pipestance that was \
                run using too old a version of Cell Ranger.\n\
                Or it might have been generated using the CS pipeline.\n\
                Or you might have copied the pipestance outs but not included \
                that file.\n"
                    .to_string(),
            );
        }
    }

    // Check alt_bc_fields.

    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            if ctl.origin_info.alt_bc_fields[li][i].0 == x {
                return Ok(true);
            }
        }
    }

    // Check names defined by VAR_DEF.

    for i in 0..ctl.gen_opt.var_def.len() {
        if x == ctl.gen_opt.var_def[i].0 {
            return Ok(true);
        }
    }

    // Check for fb<n> and fb<n>_n, and _cell versions.

    if x.starts_with("fb") {
        let mut y = x.after("fb").to_string();
        if y.ends_with("_cell") {
            y = y.rev_before("_cell").to_string();
        }
        if y.ends_with("_n") {
            y = y.rev_before("_n").to_string();
        }
        if y.parse::<usize>().is_ok() && y.force_usize() >= 1 {
            if ctl.origin_info.n() != 1 {
                return Err(
                    "\nThe variables fb<n> and fb<n>_n can only be used if there is just one \
                        dataset.\n"
                        .to_string(),
                );
            }
            if !gex_info.fb_top_matrices[0].initialized() {
                return Err(
                    "\nThe variables fb<n> and fb<n>_n can only be used if the file \
                        feature_barcode_matrix_top.bin was generated.\n"
                        .to_string(),
                );
            }
            return Ok(true);
        }
    }

    // Check for nd<k>.

    if x.starts_with("nd")
        && x.after("nd").parse::<usize>().is_ok()
        && x.after("nd").force_usize() >= 1
    {
        if *nd_used {
            return Err("\nOnly one instance of the lead variable nd<k> is allowed.\n".to_string());
        }
        *nd_used = true;
        return Ok(true);
    }

    // Check for [abbr:]count_<regex> and similar.

    if x.starts_with("count_") || x.contains(":count_") {
        let mut z = x.to_string();
        if x.contains(":count_") {
            z = x.after(":").to_string();
        }
        let mut class = "count_".to_string();
        if z.starts_with("count_cdr1_")
            || z.starts_with("count_cdr2_")
            || z.starts_with("count_cdr3_")
            || z.starts_with("count_fwr1_")
            || z.starts_with("count_fwr2_")
            || z.starts_with("count_fwr3_")
            || z.starts_with("count_fwr4_")
            || z.starts_with("count_cdr_")
            || z.starts_with("count_fwr_")
        {
            class = format!("count_{}_", z.between("_", "_"));
        }
        let y = z.after(&class);
        let reg = Regex::new(y);
        if reg.is_err() || y.contains('_') {
            return Err(format!(
                "\nThe string after {} in your lead or parseable variable {} is not a valid \
                regular expression for amino acids.\n",
                class, x
            ));
        }
        return Ok(true);
    }

    // Check for pe<n> and npe<n> and ppe<n>.

    if x.starts_with("pe") && x.after("pe").parse::<usize>().is_ok() {
        return Ok(true);
    }
    if x.starts_with("npe") && x.after("npe").parse::<usize>().is_ok() {
        return Ok(true);
    }
    if x.starts_with("ppe") && x.after("ppe").parse::<usize>().is_ok() {
        return Ok(true);
    }

    // Check for patterns.

    if is_pattern(&x, false) {
        return Ok(true);
    }

    // The rest.

    if !gex_info.have_gex && !gex_info.have_fb && x.starts_with("n_gex") {
        return Err(format!(
            "\nCan't use LVARS or LVARSP or PCOLS variable {} without having gene \
             expression or feature barcode data.\n",
            x
        ));
    }
    if !gex_info.have_gex && (x.starts_with("gex") || x == "clust" || x == "type") {
        return Err(format!(
            "\nCan't use LVARS or LVARSP or PCOLS variable {} without having gene \
             expression data.\n",
            x
        ));
    }
    let gpvar = x.starts_with('g') && x.after("g").parse::<usize>().is_ok();
    if gpvar {
        return Ok(true);
    }
    if !LVARS_ALLOWED.contains(&x.as_str()) {
        let mut end_ok = false;
        for i in 0..ends.len() {
            if x.ends_with(&ends[i]) {
                end_ok = true;
            }
        }
        if end_ok {
            return Ok(false);
        }
        if is_lvar && !x.starts_with("n_") && !x.is_empty() {
            return Err(format!(
                "\nUnrecognized variable {} for LVARS.  Please type \
                 \"enclone help lvars\".\n",
                x
            ));
        } else {
            return Ok(false);
        }
    }
    Ok(true)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn build_ends() -> Vec<String> {
    let mut ends = Vec::<String>::new();
    let ends0 = [
        "_g", "_ab", "_ag", "_cr", "_cu", "_g_μ", "_ab_μ", "_ag_μ", "_cr_μ", "_cu_μ", "_g_%",
    ];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    for x in ends0.iter() {
        for y in suffixes.iter() {
            ends.push(format!("{}{}", x, y));
        }
    }
    ends
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check lvars args.

pub fn check_lvars(ctl: &EncloneControl, gex_info: &GexInfo) -> Result<(), String> {
    let t = Instant::now();
    let mut to_check = Vec::<String>::new();
    let ends = build_ends();
    let mut nd_used = false;
    for x in ctl.clono_print_opt.lvars.iter() {
        if x.ends_with("_cell") {
            return Err(
                "\nFields ending with _cell cannot be used in LVARS or LVARSP.\n".to_string(),
            );
        }
        if !check_one_lvar(&*x, ctl, gex_info, &mut nd_used, &ends, true)? {
            to_check.push(x.clone());
        }
    }
    ctl.perf_stats(&t, "checking lvars top");
    let t = Instant::now();
    if !to_check.is_empty() {
        check_gene_fb(ctl, gex_info, &to_check, "lead")?;
    }
    ctl.perf_stats(&t, "checking gene");
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check gvars args.

pub fn check_gvars(ctl: &EncloneControl) -> Result<(), String> {
    for x in ctl.gen_opt.gvars.iter() {
        if !GVARS_ALLOWED.contains(&x.as_str()) {
            return Err(format!("\nUnknown global variable {}.\n", x));
        }
    }
    Ok(())
}
