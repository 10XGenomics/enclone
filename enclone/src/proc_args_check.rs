// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Check lvars, cvars, and pcols.

use enclone_core::defs::*;
use regex::Regex;
use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn is_pattern(x: &String, parseable: bool) -> bool {
    let ends0 = [
        "_g", "_ab", "_cr", "_cu", "_g_μ", "_ab_μ", "_cr_μ", "_cu_μ", "_g_%",
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
            if !p.is_empty() && Regex::new(&p).is_ok() {
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

fn check_gene_fb(ctl: &EncloneControl, gex_info: &GexInfo, to_check: &Vec<String>, category: &str) {
    let g_ends0 = ["_g"];
    let fb_ends0 = ["_ab", "_cr", "_cu"];
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
        if !gex_info.have_gex && !gex_info.have_fb {
            if *x == "n_gex".to_string() || *x == "n_gex_cell".to_string() {
                if category == "parseable" {
                    eprintln!(
                        "\nParseable field {} does not make sense because neither gene expression \
                         nor feature barcode data\nwere provided as input.\n",
                        x
                    );
                } else {
                    eprintln!(
                        "\nLead variable {} does not make sense because neither gene expression \
                         not feature barcode data\nwere provided as input.\n",
                        x
                    );
                }
                std::process::exit(1);
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
                    eprintln!(
                        "\nParseable field {} does not make sense because gene expression \
                         data\nwere not provided as input.\n",
                        x
                    );
                } else {
                    eprintln!(
                        "\nLead variable {} does not make sense because gene expression \
                         data\nwere not provided as input.\n",
                        x
                    );
                }
                std::process::exit(1);
            }
        }
        if !gex_info.have_fb {
            for y in fb_ends.iter() {
                if x.ends_with(y) {
                    if category == "parseable" {
                        eprintln!(
                            "\nParseable field {} does not make sense because feature \
                             barcode data\nwere not provided as input.\n",
                            x
                        );
                    } else {
                        eprintln!(
                            "\nLead variable {} does not make sense because feature barcode \
                             data\nwere not provided as input.\n",
                            x
                        );
                    }
                    std::process::exit(1);
                }
            }
        }
    }
    let mut known_features = Vec::<String>::new();
    for i in 0..gex_info.gex_features.len() {
        for j in 0..gex_info.gex_features[i].len() {
            let f = &gex_info.gex_features[i][j];
            let ff = f.split('\t').collect::<Vec<&str>>();
            if ff.len() != 3 {
                eprintln!("Unexpected structure of features file, at this line\n{}", f);
                eprintln!("Giving up.\n");
                std::process::exit(1);
            }
            for z in 0..2 {
                if ff[2].starts_with("Antibody") {
                    for s in suffixes.iter() {
                        known_features.push(format!("{}_ab{}", ff[z], s));
                    }
                } else if ff[2].starts_with("CRISPR") {
                    for s in suffixes.iter() {
                        known_features.push(format!("{}_cr{}", ff[z], s));
                    }
                } else if ff[2].starts_with("CUSTOM") {
                    for s in suffixes.iter() {
                        known_features.push(format!("{}_cu{}", ff[z], s));
                    }
                } else {
                    for s in suffixes_g.iter() {
                        known_features.push(format!("{}_g{}", ff[z], s));
                    }
                }
            }
        }
    }
    unique_sort(&mut known_features);
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
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         does not name a dataset, nor an origin,\nnor a donor, nor a tag.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
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
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a dataset, an origin, and a donor.  That's ambiguous.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
                }
                if is_dataset_name && is_origin_name {
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a dataset and an origin.  That's ambiguous.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
                }
                if is_dataset_name && is_donor_name {
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a dataset and a donor.  That's ambiguous.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
                }
                if is_origin_name && is_donor_name {
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names an origin and a donor.  That's ambiguous.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
                }
                if types != 1 {
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a tag and also a dataset, origin or donor.\n\
                         That's ambiguous.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
                }
            }
            if !n_var {
                if category == "lead" {
                    eprintln!(
                        "\nThe variable {} for LVARS is unrecognized.  Please type \
                         \"enclone help lvars\".\n",
                        x
                    );
                } else {
                    eprintln!(
                        "\nUnrecognized parseable variable {}.  Please type \
                         \"enclone help parseable\".\n",
                        x
                    );
                }
                std::process::exit(1);
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check pcols args.

pub fn check_pcols(ctl: &EncloneControl, gex_info: &GexInfo) {
    let mut alt_bcs = Vec::<String>::new();
    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.clone());
        }
    }
    unique_sort(&mut alt_bcs);
    let mut to_check = Vec::<String>::new();
    let pchains = ctl.parseable_opt.pchains;
    for x in ctl.parseable_opt.pcols.iter() {
        let mut ok = false;
        if bin_member(&alt_bcs, x) {
            ok = true;
        }
        for y in PLVARS_ALLOWED.iter() {
            if *x == *y {
                ok = true;
            }
        }
        for y in ctl.origin_info.dataset_list.iter() {
            if *x == format!("{}_barcodes", y) {
                ok = true;
            }
        }
        if ctl.parseable_opt.pbarcode {
            if *x == "barcode" {
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
            eprintln!(
                "\nCan't use parseable variable {} without having gene \
                 expression or feature barcode data.\n",
                x
            );
            std::process::exit(1);
        }
        if !gex_info.have_gex && (x.starts_with("gex") || x == "clust") || x == "type" {
            eprintln!(
                "\nCan't use parseable variable {} without having gene \
                 expression data.\n",
                x
            );
            std::process::exit(1);
        }
        if LVARS_ALLOWED.contains(&x.as_str()) || gpvar {
            ok = true;
        } else if is_pattern(&x, true) {
            ok = true;
        } else {
            for p in 1..=pchains {
                let ps = format!("{}", p);
                if x.ends_with(&ps) {
                    let y = x.rev_before(&ps);
                    if CVARS_ALLOWED.contains(&y)
                        || (ctl.parseable_opt.pbarcode && CVARS_ALLOWED_PCELL.contains(&y))
                    {
                        ok = true;
                    } else if PCVARS_ALLOWED.contains(&y) {
                        ok = true;
                    } else if y.starts_with('q')
                        && y.ends_with('_')
                        && y.between("q", "_").parse::<usize>().is_ok()
                    {
                        ok = true;
                    } else if y.starts_with("ndiff")
                        && y.ends_with("vj")
                        && y.between("ndiff", "vj").parse::<usize>().is_ok()
                        && y.between("ndiff", "vj").force_usize() >= 1
                    {
                        ok = true;
                        break;
                    }
                }
            }
        }
        if !ok {
            to_check.push(x.clone());
        }
    }
    if !to_check.is_empty() {
        check_gene_fb(&ctl, &gex_info, &to_check, "parseable");
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check cvars args.

pub fn check_cvars(ctl: &EncloneControl) {
    for x in ctl.clono_print_opt.cvars.iter() {
        let mut ok = CVARS_ALLOWED.contains(&(*x).as_str());
        if x.starts_with("ndiff")
            && x.ends_with("vj")
            && x.between("ndiff", "vj").parse::<usize>().is_ok()
            && x.between("ndiff", "vj").force_usize() >= 1
        {
            ok = true;
        }
        if !ok {
            eprintln!(
                "\nUnrecognized variable {} for CVARS or CVARSP.  \
                 Please type \"enclone help cvars\".\n",
                x
            );
            std::process::exit(1);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check lvars args.

pub fn check_lvars(ctl: &EncloneControl, gex_info: &GexInfo) {
    let mut to_check = Vec::<String>::new();
    let ends0 = [
        "_g", "_ab", "_cr", "_cu", "_g_μ", "_ab_μ", "_cr_μ", "_cu_μ", "_g_%",
    ];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let mut ends = Vec::<String>::new();
    for x in ends0.iter() {
        for y in suffixes.iter() {
            ends.push(format!("{}{}", x, y));
        }
    }
    let mut nd_used = false;
    'main_loop: for x in ctl.clono_print_opt.lvars.iter() {
        // See if type is ok.

        if *x == "type" {
            let mut specified = false;
            for i in 0..gex_info.cell_type_specified.len() {
                if gex_info.cell_type_specified[i] {
                    specified = true;
                }
            }
            if !ctl.gen_opt.internal_run {
                eprintln!(
                    "\nUnrecognized variable {} for LVARS.  Please type \
                     \"enclone help lvars\".\n",
                    x
                );
                std::process::exit(1);
            }
            if !specified {
                eprintln!(
                    "\nYou've used the lead variable \"type\", but the file \
                    cell_types.csv was not found.\n\
                    This could be because you're using a GEX pipestance that was \
                    run using too old a version of Cell Ranger.\n\
                    Or it might have been generated using the CS pipeline.\n\
                    Or you might have copied the pipestance outs but not included \
                    that file.\n"
                );
                std::process::exit(1);
            }
        }

        // Check alt_bc_fields.

        for li in 0..ctl.origin_info.alt_bc_fields.len() {
            for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
                if ctl.origin_info.alt_bc_fields[li][i].0 == *x {
                    continue 'main_loop;
                }
            }
        }

        // Check for nd<k>.

        if x.starts_with("nd")
            && x.after("nd").parse::<usize>().is_ok()
            && x.after("nd").force_usize() >= 1
        {
            if nd_used {
                eprintln!("\nOnly one instance of the lead variable nd<k> is allowed.\n");
                std::process::exit(1);
            }
            nd_used = true;
            continue;
        }

        // Check for [abbr:]count_<regex>.

        if x.starts_with("count_") || x.contains(":count_") {
            let y = x.after("count_");
            let reg = Regex::new(&y);
            if !reg.is_ok() {
                eprintln!(
                    "\nThe string after count_ in your lead variable {} is not a valid \
                    regular expression.\n",
                    x
                );
                std::process::exit(1);
            }
            continue;
        }

        // Check for pe<n> and npe<n> and ppe<n>.

        if x.starts_with("pe") && x.after("pe").parse::<usize>().is_ok() {
            continue;
        }
        if x.starts_with("npe") && x.after("npe").parse::<usize>().is_ok() {
            continue;
        }
        if x.starts_with("ppe") && x.after("ppe").parse::<usize>().is_ok() {
            continue;
        }

        // Check for patterns.

        if is_pattern(&x, false) {
            continue;
        }

        // The rest.

        if !gex_info.have_gex && !gex_info.have_fb && x.starts_with("n_gex") {
            eprintln!(
                "\nCan't use LVARS or LVARSP variable {} without having gene \
                 expression or feature barcode data.\n",
                x
            );
            std::process::exit(1);
        }
        if !gex_info.have_gex && (x.starts_with("gex") || x == "clust" || x == "type") {
            eprintln!(
                "\nCan't use LVARS or LVARSP variable {} without having gene \
                 expression data.\n",
                x
            );
            std::process::exit(1);
        }
        if x.ends_with("_cell") {
            eprintln!("\nFields ending with _cell cannot be used in LVARS or LVARSP.\n");
            std::process::exit(1);
        }
        let gpvar = x.starts_with('g') && x.after("g").parse::<usize>().is_ok();
        if !(LVARS_ALLOWED.contains(&x.as_str()) || gpvar) {
            let mut end_ok = false;
            for i in 0..ends.len() {
                if x.ends_with(&ends[i]) {
                    end_ok = true;
                }
            }
            if !end_ok && !x.starts_with("n_") {
                eprintln!(
                    "\nUnrecognized variable {} for LVARS.  Please type \
                     \"enclone help lvars\".\n",
                    x
                );
                std::process::exit(1);
            } else {
                to_check.push(x.clone());
            }
        }
    }
    if !to_check.is_empty() {
        check_gene_fb(&ctl, &gex_info, &to_check, "lead");
    }
}
