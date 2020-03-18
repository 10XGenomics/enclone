// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Check lvars, cvars, and pcols.  Sadly repetitive.

use crate::defs::*;
use itertools::*;
use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn check_gene_fb(
    ctl: &EncloneControl, 
    gex_features: &Vec<Vec<String>>,
    to_check: &Vec<String>,
    category: &str,
) {
    let mut known_features = Vec::<String>::new();
    for i in 0..gex_features.len() {
        for j in 0..gex_features[i].len() {
            let f = &gex_features[i][j];
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
                } else if ff[2].starts_with("Antigen") {
                    for s in suffixes.iter() {
                        known_features.push(format!("{}_ag{}", ff[z], s));
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
                    for s in suffixes.iter() {
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
        if !bin_member(&known_features, &x) {
            let mut n_var = false;
            if x.starts_with("n_") {
                n_var = true;
                let mut is_dataset_name = false;
                let mut is_sample_name = false;
                let mut is_donor_name = false;
                let mut is_tag_name = false;
                let name = x.after("n_").to_string();
                let s = ctl.sample_info.n();
                for j in 0..s {
                    if ctl.sample_info.dataset_id[j] == name {
                        is_dataset_name = true;
                    }
                }
                for j in 0..ctl.sample_info.sample_list.len() {
                    if ctl.sample_info.sample_list[j] == name {
                        is_sample_name = true;
                    }
                }
                for j in 0..ctl.sample_info.donor_list.len() {
                    if ctl.sample_info.donor_list[j] == name {
                        is_donor_name = true;
                    }
                }
                for j in 0..ctl.sample_info.tag_list.len() {
                    if ctl.sample_info.tag_list[j] == name {
                        is_tag_name = true;
                    }
                }
                let msg = "\nSuggested reading: \"enclone help input\" and \
                           \"enclone help glossary\".\n";
                if !is_dataset_name && !is_sample_name && !is_donor_name && !is_tag_name {
                    eprintln!(
                        // The following line is cryptic as a user message:
                        "\ntags = {}\n\
                         You've used the {} variable {}, and yet {} \
                         does not name a dataset, nor a sample,\nnor a donor, nor a tag.\n{}",
                        category,
                        ctl.sample_info.tag_list.iter().format(","),
                        x,
                        name,
                        msg
                    );
                    std::process::exit(1);
                }
                let mut types = 0;
                if is_dataset_name {
                    types += 1;
                }
                if is_sample_name {
                    types += 1;
                }
                if is_donor_name {
                    types += 1;
                }
                if is_tag_name {
                    types += 1;
                }
                if is_dataset_name && is_sample_name && is_donor_name {
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a dataset, a sample, and a donor.  That's ambiguous.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
                }
                if is_dataset_name && is_sample_name {
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a dataset and a sample.  That's ambiguous.\n{}",
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
                if is_sample_name && is_donor_name {
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a sample and a donor.  That's ambiguous.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
                }
                if types != 1 {
                    eprintln!(
                        "\nYou've used the {} variable {}, and yet {} \
                         names a tag and also a dataset, sample or donor.\n\
                         That's ambiguous.\n{}",
                        category, x, name, msg
                    );
                    std::process::exit(1);
                }
            }
            if !n_var {
                if category == "lead" {
                    eprintln!(
                        "\nUnrecognized variable {} for LVARS.  Please type \
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

pub fn check_pcols(ctl: &EncloneControl) {
    let mut to_check = Vec::<String>::new();
    let pchains = ctl.parseable_opt.pchains;
    let ends0 = [
        "_g", "_ab", "_ag", "_cr", "_cu", "_g_μ", "_ab_μ", "_ag_μ", "_cr_μ", "_cu_μ",
    ];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let mut ends = Vec::<String>::new();
    for x in ends0.iter() {
        for y in suffixes.iter() {
            ends.push(format!("{}{}", x, y));
        }
    }
    let mut have_gex = false;
    for i in 0..ctl.sample_info.gex_path.len() {
        if ctl.sample_info.gex_path[i].len() > 0 {
            have_gex = true;
        }
    }
    for x in ctl.parseable_opt.pcols.iter() {
        if !have_gex {
            if *x == "gex".to_string()
                || x.starts_with("gex_")
                || x.ends_with("_g")
                || x.ends_with("_g_μ")
                || *x == "n_gex_cell".to_string()
                || *x == "n_gex".to_string()
                || *x == "entropy".to_string()
            {
                eprint("\nParseable field {} does not make sense because gene expression \
                    data were not provided as input.\n");
                std::process::exit(1);
            }
        }

        *** what about fb data ***

        let mut ok = false;






    for col in 0..ctl.parseable_opt.pchains {
        for x in CVARS_ALLOWED.iter() {
            speakerc!(col, x);
        }
        if ctl.parseable_opt.pbarcode {
            for x in CVARS_ALLOWED_PCELL.iter() {
                speakerc!(col, x);
            }
        }
        for x in &["v_name", "d_name", "j_name", "v_id", "d_id", "j_id"] {
            speakerc!(col, x);
        }
        for x in &[
            "var_indices_dna",
            "var_indices_aa",
            "share_indices_dna",
            "share_indices_aa",
        ] {
            speakerc!(col, x);
        }
        for x in &[
            "v_start",
            "const_id",
            "utr_id",
            "utr_name",
            "cdr3_start",
            "cdr3_aa",
        ] {
            speakerc!(col, x);
        }
        for x in &["seq", "vj_seq", "var_aa"] {
            speakerc!(col, x);
        }
        for i in 0..pcols_sort.len() {
            if pcols_sort[i].starts_with('q') && pcols_sort[i].ends_with(&format!("_{}", col + 1)) {
                let x = pcols_sort[i].after("q").rev_before("_");
                if x.parse::<usize>().is_ok() {
                    parseable_fields.push(pcols_sort[i].clone());
                }
            }
        }
    }
    speaker!("group_id");
    speaker!("group_ncells");
    speaker!("clonotype_id");
    speaker!("clonotype_ncells");
    speaker!("nchains");
    speaker!("exact_subclonotype_id");
    speaker!("barcodes");
    for x in ctl.sample_info.dataset_list.iter() {
        speaker!(&format!("{}_barcodes", x));
    }
    if ctl.parseable_opt.pbarcode {
        speaker!("barcode");
        for x in ctl.sample_info.dataset_list.iter() {
            speaker!(&format!("{}_barcode", x));
        }
    }
    let mut pfsort = parseable_fields.clone();
    unique_sort(&mut pfsort);
    for x in pcols_sort.iter() {
        let y = x.clone();
        let mut y_alt = y.clone();
        y_alt = y_alt.replace("_Σ", "_sum");
        y_alt = y_alt.replace("_μ", "_mean");
        if !bin_member(&pfsort, &y) && !bin_member(&pfsort, &y_alt) {
            eprintln!("\nUnknown parseable output field: {}.\n", x);
            eprintln!(
                "Note that the allowed fields depend on your specification for the \
                 LVARS or LVARSP,\nand CVARS or CVARSP options, and that for the latter two,\n\
                 suffixing by the column number is required.  Please see \
                 \"enclone help parseable\".\n"
            );
            std::process::exit(1);
        }
    }
}





        let gpvar = x.starts_with('g') && x.after("g").parse::<usize>().is_ok();
        if LVARS_ALLOWED.contains(&x.as_str()) || gpvar {
            ok = true;
        } else {
            for p in 1..=pchains {
                let ps = format!( "{}", p );
                if x.ends_with(&ps) {
                    let y = x.rev_before(&ps);
                    if CVARS_ALLOWED.contains(&y) {
                        ok = true;
                    } else if y.starts_with("ndiff") {
                        && y.after("ndiff").parse::<usize>().is_ok()
                        && y.after("ndiff").force_usize() >= 1
                    {
                        ok = true;
                        break;
                    }
                }
            }
        }
        if !ok {
            let mut end_ok = false;
            for i in 0..ends.len() {
                if x.ends_with(&ends[i]) {
                    end_ok = true;
                }
            }
            if !end_ok && !x.starts_with("n_") {
                eprintln!(
                    "\nUnrecognized variable {} for PCOLS.  Please type \
                     \"enclone help parseable\".\n",
                    x
                );
                std::process::exit(1);
            } else {
                to_check.push(x.clone());
            }
        }
    }
    if !to_check.is_empty() {
        check_gene_fb(&ctl, &gex_features, &to_check, "parseable");
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check cvars args.

pub fn check_cvars(ctl: &EncloneControl) {
    for x in ctl.clono_print_opt.cvars.iter() {
        let mut ok = CVARS_ALLOWED.contains(&(*x).as_str());
        if x.starts_with("ndiff")
            && x.after("ndiff").parse::<usize>().is_ok()
            && x.after("ndiff").force_usize() >= 1
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

pub fn check_lvars(ctl: &EncloneControl, gex_features: &Vec<Vec<String>>) {
    let mut to_check = Vec::<String>::new();
    let ends0 = [
        "_g", "_ab", "_ag", "_cr", "_cu", "_g_μ", "_ab_μ", "_ag_μ", "_cr_μ", "_cu_μ",
    ];
    let suffixes = ["", "_min", "_max", "_μ", "_Σ"];
    let mut ends = Vec::<String>::new();
    for x in ends0.iter() {
        for y in suffixes.iter() {
            ends.push(format!("{}{}", x, y));
        }
    }
    for x in ctl.clono_print_opt.lvars.iter() {
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
        check_gene_fb(&ctl, &gex_features, &to_check, "lead");
    }
}
