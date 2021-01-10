// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the two functions proc_xcr and proc_meta.

use enclone_core::defs::*;
use io_utils::*;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::thread;
use std::time;
use std::time::Instant;
use string_utils::*;
use tilde_expand::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn expand_integer_ranges(x: &str) -> String {
    let mut tokens = Vec::<String>::new();
    let mut token = String::new();
    for c in x.chars() {
        if c == ',' || c == ':' || c == ';' {
            if token.len() > 0 {
                tokens.push(token.clone());
                token.clear();
            }
            tokens.push(c.to_string());
        } else {
            token.push(c);
        }
    }
    if token.len() > 0 {
        tokens.push(token);
    }
    let mut tokens2 = Vec::<String>::new();
    for i in 0..tokens.len() {
        if tokens[i].contains("-")
            && tokens[i].before("-").parse::<usize>().is_ok()
            && tokens[i].after("-").parse::<usize>().is_ok()
        {
            let n1 = tokens[i].before("-").force_usize();
            let n2 = tokens[i].after("-").force_usize();
            if n1 <= n2 {
                for n in n1..=n2 {
                    if n > n1 {
                        tokens2.push(",".to_string());
                    }
                    tokens2.push(format!("{}", n));
                }
                continue;
            }
        }
        tokens2.push(tokens[i].clone());
    }
    let mut y = String::new();
    for i in 0..tokens2.len() {
        y += &tokens2[i];
    }
    y
}

fn expand_analysis_sets(x: &str) -> String {
    let mut tokens = Vec::<String>::new();
    let mut token = String::new();
    for c in x.chars() {
        if c == ',' || c == ':' || c == ';' {
            if token.len() > 0 {
                tokens.push(token.clone());
                token.clear();
            }
            tokens.push(c.to_string());
        } else {
            token.push(c);
        }
    }
    if token.len() > 0 {
        tokens.push(token);
    }
    let mut tokens2 = Vec::<String>::new();
    for i in 0..tokens.len() {
        if tokens[i].starts_with('S') {
            let setid = tokens[i].after("S");
            // do not use xena.txgmesh.net, does not work from inside enclone
            let url = format!("https://xena.fuzzplex.com/api/analysis_sets/{}", setid);
            let o = Command::new("curl")
                .arg(url)
                .output()
                .expect("failed to execute xena http");
            let m = String::from_utf8(o.stdout).unwrap();
            if m.contains("502 Bad Gateway") {
                // do not use xena.txgmesh.net, does not work from inside enclone
                eprintln!(
                    "\nWell, this is sad.  The URL \
                    http://xena.fuzzplex.com/api/analysis_sets/{} returned a 502 Bad Gateway \
                    message.  Please try again later or ask someone for help.\n\n",
                    setid
                );
                std::process::exit(1);
            }
            // printme!(m);
            if m.contains("\"analysis_ids\":[") {
                let mut ids = m.between("\"analysis_ids\":[", "]").to_string();
                ids = ids.replace(" ", "");
                ids = ids.replace("\n", "");
                let ids = ids.split(',').collect::<Vec<&str>>();
                let mut ids2 = Vec::<String>::new();

                // Remove wiped analysis IDs.

                for j in 0..ids.len() {
                    // do not use xena.txgmesh.net, does not work from inside enclone
                    let url = format!("https://xena.fuzzplex.com/api/analyses/{}", ids[j]);
                    let o = Command::new("curl")
                        .arg(url)
                        .output()
                        .expect("failed to execute xena http");
                    let m = String::from_utf8(o.stdout).unwrap();
                    if m.contains("502 Bad Gateway") {
                        // do not use xena.txgmesh.net, does not work from inside enclone
                        eprintln!(
                            "\nWell, this is sad.  The URL \
                            http://xena.fuzzplex.com/api/analyses/{} returned a 502 Bad Gateway \
                            message.  Please try again later or ask someone for help.\n",
                            ids[j]
                        );
                        std::process::exit(1);
                    }
                    if !m.contains("\"wiped\"") {
                        ids2.push(ids[j].to_string());
                    }
                }

                // Proceed.

                for j in 0..ids2.len() {
                    if j > 0 {
                        tokens2.push(",".to_string());
                    }
                    tokens2.push(ids2[j].to_string());
                }
                continue;
            } else {
                eprintln!(
                    "\nIt looks like you've provided an incorrect analysis set ID {}.\n",
                    setid
                );
                std::process::exit(1);
            }
        }
        tokens2.push(tokens[i].clone());
    }
    let mut y = String::new();
    for i in 0..tokens2.len() {
        y += &tokens2[i];
    }
    y
}

// Functions to find the path to data.

pub fn get_path_fail(p: &str, ctl: &EncloneControl, source: &str) -> String {
    for x in ctl.gen_opt.pre.iter() {
        let pp = format!("{}/{}", x, p);
        if path_exists(&pp) {
            return pp;
        }
    }
    if !path_exists(&p) {
        if ctl.gen_opt.pre.is_empty() {
            let path = std::env::current_dir().unwrap();
            eprintln!(
                "\nIn directory {}, unable to find the path {}.  This came from the {} argument.\n",
                path.display(),
                p,
                source
            );
        } else {
            let path = std::env::current_dir().unwrap();
            eprintln!(
                "\nIn directory {}, unable to find the\npath {},\n\
                even if prepended by any of the directories \
                in\nPRE={}.\nThis came from the {} argument.\n",
                path.display(),
                p,
                ctl.gen_opt.pre.iter().format(","),
                source
            );
        }
        if ctl.gen_opt.cellranger {
            panic!("This should not have happened, so a traceback follows.\n");
        }
        std::process::exit(1);
    }
    p.to_string()
}

fn get_path(p: &str, ctl: &EncloneControl, ok: &mut bool) -> String {
    *ok = false;
    for x in ctl.gen_opt.pre.iter() {
        let mut pp = format!("{}/{}", x, p);
        if pp.starts_with("~") {
            pp = stringme(&tilde_expand(&pp.as_bytes()));
        }
        if path_exists(&pp) {
            *ok = true;
            return pp;
        }
    }
    let mut pp = p.to_string();
    if pp.starts_with("~") {
        pp = stringme(&tilde_expand(&pp.as_bytes()));
    }
    *ok = path_exists(&pp);
    pp
}

fn get_path_or_internal_id(
    p: &str,
    ctl: &EncloneControl,
    source: &str,
    _current_ref: &mut bool,
) -> String {
    let mut ok = false;
    let mut pp = get_path(&p, &ctl, &mut ok);
    if !ok {
        if !ctl.gen_opt.internal_run {
            get_path_fail(&pp, &ctl, source);
        } else {
            // For internal runs, try much harder.  This is so that internal users can
            // just type an internal numerical id for a dataset and have it always
            // work.  The code that's used here should be placed somewhere else.

            if p.parse::<usize>().is_ok() {
                // do not use xena.txgmesh.net, does not work from inside enclone
                let url = format!("https://xena.fuzzplex.com/api/analyses/{}", p);
                let o = Command::new("curl")
                    .arg(url.clone())
                    .output()
                    .expect("failed to execute xena http");
                let m = String::from_utf8(o.stdout).unwrap();
                if o.status.code() != Some(0) {
                    let merr = String::from_utf8(o.stderr).unwrap();
                    eprintln!(
                        "\nSomething went wrong. the URL \
                        \nhttp://xena.fuzzplex.com/api/analyses/{}\n\
                        failed with the following stderr:\n{}\n\
                        and the following stdout:\n{}\n",
                        p, merr, m
                    );
                    std::process::exit(1);
                }
                if m.contains("502 Bad Gateway") {
                    // do not use xena.txgmesh.net, does not work from inside enclone
                    eprintln!(
                        "\nWell this is sad.  The URL \
                        http://xena.fuzzplex.com/api/analyses/{} yielded a 502 Bad Gateway \
                        message.  Please try again later or ask someone for help.\n",
                        p
                    );
                    std::process::exit(1);
                }
                if m.contains("\"path\":\"") {
                    let path = m.between("\"path\":\"", "\"").to_string();
                    pp = format!("{}/outs", path);
                    if !path_exists(&pp) {
                        thread::sleep(time::Duration::from_millis(100));
                        if path_exists(&pp) {
                            eprintln!(
                                "\nYou are experiencing unstable filesystem access: \
                                100 milliseconds ago, \
                                the path\n\
                                {}\nwas not visible, but now it is.  You might consider posting \
                                this problem on \
                                the slack channel #seqops-bespin.\nOr retry again.  enclone is \
                                giving up because \
                                if filesystem access blinks in and out of existence,\n\
                                other more cryptic events are likely to occur.\n",
                                pp
                            );
                        } else {
                            eprintln!(
                                "\nIt looks like you've provided an analysis ID for \
                                which the pipeline outs folder\n{}\nhas not yet been generated.\n\
                                This path did not exist:\n{}\n",
                                p, pp
                            );
                            eprintln!("Here is the stdout:\n{}\n", m);
                        }
                        std::process::exit(1);
                    }
                } else {
                    eprintln!(
                        "\nIt looks like you've provided either an incorrect \
                        analysis ID {} or else one for which\n\
                        the pipeline outs folder has not yet been generated.\n\
                        This URL\n{}\ndid not provide a path.\n",
                        p, url
                    );
                    std::process::exit(1);
                }
            } else {
                eprintln!(
                    "\nAfter searching high and low, your path\n{}\nfor {} \
                    cannot be found.\nPlease check its value and also the value \
                    for PRE if you provided that.\n",
                    p, source
                );
                std::process::exit(1);
            }
        }
    }
    if !pp.ends_with("/outs") && path_exists(&format!("{}/outs", pp)) {
        pp = format!("{}/outs", pp);
    }
    pp
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Parse barcode-level information file.

fn parse_bc(mut bc: String, ctl: &mut EncloneControl, call_type: &str) {
    let mut origin_for_bc = HashMap::<String, String>::new();
    let mut donor_for_bc = HashMap::<String, String>::new();
    let mut tag = HashMap::<String, String>::new();
    let mut barcode_color = HashMap::<String, String>::new();
    let mut alt_bc_fields = Vec::<(String, HashMap<String, String>)>::new();
    if bc != "".to_string() {
        bc = get_path_fail(&bc, &ctl, call_type);
        let f = open_for_read![&bc];
        let mut first = true;
        let mut fieldnames = Vec::<String>::new();
        let mut barcode_pos = 0;
        let (mut origin_pos, mut donor_pos, mut tag_pos, mut color_pos) = (None, None, None, None);
        let mut to_alt = Vec::<isize>::new();
        for line in f.lines() {
            let s = line.unwrap();
            if first {
                let fields = s.split(',').collect::<Vec<&str>>();
                to_alt = vec![-1 as isize; fields.len()];
                if !fields.contains(&"barcode") {
                    let mut origin = "from the bc field used in META";
                    if call_type == "BC" {
                        origin = "from the BC argument";
                    }
                    if bc.ends_with(".tsv") {
                        eprintln!(
                            "\nThe file\n{}\n{}\nappears to be a TSV file rather than a \
                            CSV file, as required.\n",
                            bc, origin,
                        );
                    } else {
                        eprintln!(
                            "\nThe file\n{}\n{}\nis missing the barcode field.\n",
                            bc, origin,
                        );
                    }
                    std::process::exit(1);
                }
                for x in fields.iter() {
                    fieldnames.push(x.to_string());
                }
                for i in 0..fields.len() {
                    if fields[i] == "barcode" {
                        barcode_pos = i;
                    } else if fields[i] == "origin" {
                        origin_pos = Some(i);
                    } else if fields[i] == "donor" {
                        donor_pos = Some(i);
                    } else if fields[i] == "tag" {
                        tag_pos = Some(i);
                    } else if fields[i] == "color" {
                        color_pos = Some(i);
                    } else {
                        to_alt[i] = alt_bc_fields.len() as isize;
                        alt_bc_fields
                            .push((fields[i].to_string(), HashMap::<String, String>::new()));
                    }
                }
                first = false;
            } else {
                let fields = s.split(',').collect::<Vec<&str>>();
                if fields.len() != fieldnames.len() {
                    let mut origin = "bc in META";
                    if call_type == "BC" {
                        origin = "BC";
                    }
                    eprintln!(
                        "\nThere is a line\n{}\nin a CSV file defined by {}\n\
                         that has {} fields, which isn't right, because the header line \
                         has {} fields.  This is for the file\n{}.\n",
                        s,
                        origin,
                        fields.len(),
                        fieldnames.len(),
                        bc
                    );
                    std::process::exit(1);
                }
                for i in 0..fields.len() {
                    if to_alt[i] >= 0 {
                        alt_bc_fields[to_alt[i] as usize]
                            .1
                            .insert(fields[barcode_pos].to_string(), fields[i].to_string());
                    }
                }
                if !fields[barcode_pos].contains('-') {
                    let mut origin = "bc in META";
                    if call_type == "BC" {
                        origin = "BC";
                    }
                    eprintln!(
                        "\nThe barcode \"{}\" appears in the file\n{}\ndefined \
                         by {}.  That doesn't make sense because a barcode\n\
                         should include a hyphen.\n",
                        fields[barcode_pos], bc, origin
                    );
                    std::process::exit(1);
                }
                if origin_pos.is_some() {
                    origin_for_bc.insert(
                        fields[barcode_pos].to_string(),
                        fields[origin_pos.unwrap()].to_string(),
                    );
                }
                if donor_pos.is_some() {
                    donor_for_bc.insert(
                        fields[barcode_pos].to_string(),
                        fields[donor_pos.unwrap()].to_string(),
                    );
                }
                if tag_pos.is_some() {
                    let tag_pos = tag_pos.unwrap();
                    tag.insert(fields[barcode_pos].to_string(), fields[tag_pos].to_string());
                }
                if color_pos.is_some() {
                    let color_pos = color_pos.unwrap();
                    barcode_color.insert(
                        fields[barcode_pos].to_string(),
                        fields[color_pos].to_string(),
                    );
                }
            }
        }
    }
    ctl.origin_info.origin_for_bc.push(origin_for_bc);
    ctl.origin_info.donor_for_bc.push(donor_for_bc);
    ctl.origin_info.tag.push(tag);
    ctl.origin_info.barcode_color.push(barcode_color);
    ctl.origin_info.alt_bc_fields.push(alt_bc_fields);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_xcr(f: &str, gex: &str, bc: &str, have_gex: bool, mut ctl: &mut EncloneControl) {
    ctl.origin_info = OriginInfo::default();
    if (ctl.gen_opt.tcr && f.starts_with("BCR=")) || (ctl.gen_opt.bcr && f.starts_with("TCR=")) {
        eprintln!("\nOnly one of TCR or BCR can be specified.\n");
        std::process::exit(1);
    }
    let t = Instant::now();
    ctl.gen_opt.tcr = f.starts_with("TCR=");
    ctl.gen_opt.bcr = f.starts_with("BCR=");
    let mut val: String;
    if ctl.gen_opt.tcr {
        val = f.after("TCR=").to_string();
    } else if ctl.gen_opt.bcr {
        val = f.after("BCR=").to_string();
    } else {
        val = f.to_string();
    }
    if val == "".to_string() {
        eprintln!(
            "\nYou can't write {} with no value on the right hand side.",
            f
        );
        eprintln!("Perhaps you need to remove some white space from your command line.\n");
        std::process::exit(1);
    }
    val = expand_integer_ranges(&val);
    if ctl.gen_opt.internal_run {
        val = expand_analysis_sets(&val);
    }
    let donor_groups;
    if ctl.gen_opt.cellranger {
        donor_groups = vec![&val[..]];
    } else {
        donor_groups = val.split(';').collect::<Vec<&str>>();
    }
    let mut gex2 = expand_integer_ranges(&gex);
    if ctl.gen_opt.internal_run {
        gex2 = expand_analysis_sets(&gex2);
    }
    let donor_groups_gex;
    if ctl.gen_opt.cellranger {
        donor_groups_gex = vec![&gex2[..]];
    } else {
        donor_groups_gex = gex2.split(';').collect::<Vec<&str>>();
    }
    let donor_groups_bc = bc.split(';').collect::<Vec<&str>>();
    let mut xcr = "TCR".to_string();
    if ctl.gen_opt.bcr {
        xcr = "BCR".to_string();
    }
    if have_gex && donor_groups_gex.len() != donor_groups.len() {
        eprintln!(
            "\nThere are {} {} donor groups and {} GEX donor groups, so \
             the {} and GEX arguments do not exactly mirror each \
             other's structure.\n",
            xcr,
            donor_groups.len(),
            donor_groups_gex.len(),
            xcr
        );
        std::process::exit(1);
    }
    if !bc.is_empty() && donor_groups_bc.len() != donor_groups.len() {
        eprintln!(
            "\nThe {} and BC arguments do not exactly mirror each \
             other's structure.\n",
            xcr
        );
        std::process::exit(1);
    }
    ctl.perf_stats(&t, "in proc_xcr 1");
    let t = Instant::now();
    for (id, d) in donor_groups.iter().enumerate() {
        let origin_groups;
        if ctl.gen_opt.cellranger {
            origin_groups = vec![&d[..]];
        } else {
            origin_groups = (*d).split(':').collect::<Vec<&str>>();
        }
        let mut origin_groups_gex = Vec::<&str>::new();
        if have_gex {
            if ctl.gen_opt.cellranger {
                origin_groups_gex = vec![&donor_groups_gex[id][..]];
            } else {
                origin_groups_gex = donor_groups_gex[id].split(':').collect::<Vec<&str>>();
            }
            if origin_groups_gex.len() != origin_groups.len() {
                eprintln!(
                    "\nFor donor {}, there are {} {} origin groups and {} GEX origin groups, so \
                     the {} and GEX arguments do not exactly mirror each \
                     other's structure.\n",
                    id + 1,
                    xcr,
                    origin_groups.len(),
                    origin_groups_gex.len(),
                    xcr
                );
                std::process::exit(1);
            }
        }
        let mut origin_groups_bc = Vec::<&str>::new();
        if !bc.is_empty() {
            origin_groups_bc = donor_groups_bc[id].split(':').collect::<Vec<&str>>();
            if origin_groups_bc.len() != origin_groups.len() {
                eprintln!(
                    "\nThe {} and BC arguments do not exactly mirror each \
                     other's structure.\n",
                    xcr
                );
                std::process::exit(1);
            }
        }
        for (is, s) in origin_groups.iter().enumerate() {
            let datasets;
            if ctl.gen_opt.cellranger {
                datasets = vec![&s[..]];
            } else {
                datasets = (*s).split(',').collect::<Vec<&str>>();
            }
            let datasets_gex: Vec<&str>;
            let mut datasets_bc = Vec::<&str>::new();
            if have_gex {
                if ctl.gen_opt.cellranger {
                    datasets_gex = vec![&origin_groups_gex[is][..]];
                } else {
                    datasets_gex = origin_groups_gex[is].split(',').collect::<Vec<&str>>();
                }
                if datasets_gex.len() != datasets.len() {
                    eprintln!(
                        "\nSee {} {} datasets and {} GEX datasets, so \
                         the {} and GEX arguments do not exactly mirror each \
                         other's structure.\n",
                        xcr,
                        datasets.len(),
                        datasets_gex.len(),
                        xcr
                    );
                    std::process::exit(1);
                }
            }
            if !bc.is_empty() {
                datasets_bc = origin_groups_bc[is].split(',').collect::<Vec<&str>>();
                if datasets_bc.len() != datasets.len() {
                    eprintln!(
                        "\nThe {} and BC arguments do not exactly mirror each \
                         other's structure.\n",
                        xcr
                    );
                    std::process::exit(1);
                }
            }
            for (ix, x) in datasets.iter().enumerate() {
                ctl.origin_info.color.push("".to_string());
                ctl.origin_info.tag.push(HashMap::<String, String>::new());
                let donor_name = format!("d{}", id + 1);
                let origin_name = format!("s{}", is + 1);
                ctl.origin_info.donor_id.push(donor_name);
                ctl.origin_info.origin_id.push(origin_name);
                let mut dataset_name = (*x).to_string();
                if dataset_name.contains('/') {
                    dataset_name = dataset_name.rev_after("/").to_string();
                }
                ctl.origin_info.descrips.push(dataset_name.clone());
                ctl.origin_info.dataset_id.push(dataset_name.clone());

                // Now work on the BC path.

                let mut bcx = String::new();
                if !bc.is_empty() {
                    bcx = datasets_bc[ix].to_string();
                }
                parse_bc(bcx, &mut ctl, "BC");
            }
        }
    }

    // Get paths.  This will need to change when cellranger switches to multi.  This code is
    // parallelized because this code can indirectly make many calls to path_exists, and the wall
    // clock time for these can add up.  There should be a way to do this that does not involve
    // multithreading.

    let mut source = f.clone();
    if f.contains('=') {
        source = f.before("=");
    }
    let mut results = Vec::<(String, String, bool)>::new();
    for (id, d) in donor_groups.iter().enumerate() {
        let origin_groups = (*d).split(':').collect::<Vec<&str>>();
        let mut origin_groups_gex = Vec::<&str>::new();
        if have_gex {
            origin_groups_gex = donor_groups_gex[id].split(':').collect::<Vec<&str>>();
        }
        for (is, s) in origin_groups.iter().enumerate() {
            let datasets = (*s).split(',').collect::<Vec<&str>>();
            let mut datasets_gex = Vec::<&str>::new();
            if have_gex {
                datasets_gex = origin_groups_gex[is].split(',').collect::<Vec<&str>>();
            }
            for (ix, x) in datasets.iter().enumerate() {
                let p = (*x).to_string();
                let mut pg = String::new();
                if have_gex {
                    pg = datasets_gex[ix].to_string();
                }
                results.push((p, pg, false));
            }
        }
    }
    results.par_iter_mut().for_each(|res| {
        let (p, pg) = (&mut res.0, &mut res.1);
        let mut current_ref = &mut res.2;
        *p = get_path_or_internal_id(&p, &ctl, source, &mut current_ref);
        if ctl.gen_opt.bcr && path_exists(&format!("{}/vdj_b", p)) {
            *p = format!("{}/vdj_b", p);
        }
        if ctl.gen_opt.tcr && path_exists(&format!("{}/vdj_t", p)) {
            *p = format!("{}/vdj_t", p);
        }
        if have_gex {
            *pg = get_path_or_internal_id(&pg, &ctl, "GEX", &mut current_ref);
            if path_exists(&format!("{}/count", pg)) {
                *pg = format!("{}/count", pg);
            }
        }
    });
    for i in 0..results.len() {
        ctl.origin_info.dataset_path.push(results[i].0.clone());
        ctl.origin_info.gex_path.push(results[i].1.clone());
    }
    ctl.perf_stats(&t, "in proc_xcr 2");
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_meta_core(lines: &Vec<String>, mut ctl: &mut EncloneControl) {
    let mut fields = Vec::<String>::new();
    let mut donors = Vec::<String>::new();
    for (count, s) in lines.iter().enumerate() {
        if count == 0 {
            let x = s.split(',').collect::<Vec<&str>>();
            for i in 0..x.len() {
                fields.push(x[i].to_string());
            }
            let mut fields_sorted = fields.clone();
            unique_sort(&mut fields_sorted);
            if fields_sorted.len() < fields.len() {
                eprintln!(
                    "\nThe CSV file that you specified using the META or METAX argument \
                     has duplicate field names\nin its first line.\n"
                );
                std::process::exit(1);
            }
            let allowed_fields = vec![
                "bc".to_string(),
                "bcr".to_string(),
                "donor".to_string(),
                "gex".to_string(),
                "origin".to_string(),
                "tcr".to_string(),
                "color".to_string(),
            ];
            for x in fields.iter() {
                if !allowed_fields.contains(&x) {
                    eprintln!(
                        "\nThe CSV file that you specified using the META or METAX argument \
                         has an illegal field name ({}) in its first line.\n",
                        x
                    );
                    std::process::exit(1);
                }
            }
            ctl.gen_opt.tcr = fields.contains(&"tcr".to_string());
            ctl.gen_opt.bcr = fields.contains(&"bcr".to_string());
            if !ctl.gen_opt.tcr && !ctl.gen_opt.bcr {
                eprintln!(
                    "\nThe CSV file that you specified using the META or METAX argument \
                     has neither the field tcr or bcr in its first line.\n"
                );
                std::process::exit(1);
            }
            if ctl.gen_opt.tcr && ctl.gen_opt.bcr {
                eprintln!(
                    "\nThe CSV file that you specified using the META or METAX argument \
                     has both the fields tcr and bcr in its first line.\n"
                );
                std::process::exit(1);
            }
        } else if !s.starts_with('#') {
            let val = s.split(',').collect::<Vec<&str>>();
            if val.len() != fields.len() {
                eprintln!(
                    "\nMETA or METAX file line {} has a different number of fields than the \
                     first line of the file.\n",
                    count + 1
                );
                std::process::exit(1);
            }
            let mut path = String::new();
            let mut abbr = String::new();
            let mut gpath = String::new();
            let mut origin = "s1".to_string();
            let mut donor = "d1".to_string();
            let mut color = "".to_string();
            let mut bc = "".to_string();
            for i in 0..fields.len() {
                let x = &fields[i];
                let mut y = val[i].to_string();
                if y.starts_with('"') && y.ends_with('"') {
                    y = y.after("\"").rev_before("\"").to_string();
                }
                if *x == "tcr" || *x == "bcr" {
                    if y.contains(':') {
                        path = y.after(":").to_string();
                        abbr = y.before(":").to_string();
                    } else {
                        path = y.to_string();
                        if path.contains("/") {
                            abbr = path.rev_after("/").to_string();
                        } else {
                            abbr = path.clone();
                        }
                    }
                } else if *x == "gex" {
                    gpath = y.to_string();
                } else if *x == "origin" {
                    origin = y.to_string();
                } else if *x == "donor" {
                    donor = y.to_string();
                } else if *x == "color" {
                    color = y.to_string();
                } else if *x == "bc" && y.len() > 0 {
                    bc = y.to_string();
                }
            }

            // Parse bc and finish up.

            parse_bc(bc.clone(), &mut ctl, "META");
            let mut current_ref = false;
            path = get_path_or_internal_id(&path, &ctl, "META", &mut current_ref);
            if ctl.gen_opt.bcr && path_exists(&format!("{}/vdj_b", path)) {
                path = format!("{}/vdj_b", path);
            }
            if ctl.gen_opt.tcr && path_exists(&format!("{}/vdj_t", path)) {
                path = format!("{}/vdj_t", path);
            }
            if gpath.len() > 0 {
                gpath = get_path_or_internal_id(&gpath, &mut ctl, "META", &mut current_ref);
                if path_exists(&format!("{}/count", gpath)) {
                    gpath = format!("{}/count", gpath);
                }
            }
            if current_ref {
                ctl.gen_opt.current_ref = true;
            }
            let mut dp = None;
            for j in 0..donors.len() {
                if donor == donors[j] {
                    dp = Some(j);
                    break;
                }
            }
            if dp.is_none() {
                donors.push(donor.clone());
            }
            ctl.origin_info.descrips.push(abbr.clone());
            ctl.origin_info.dataset_path.push(path);
            ctl.origin_info.gex_path.push(gpath);
            ctl.origin_info.dataset_id.push(abbr);
            ctl.origin_info.donor_id.push(donor);
            ctl.origin_info.origin_id.push(origin);
            ctl.origin_info.color.push(color);
        }
    }
}

pub fn proc_meta(f: &str, mut ctl: &mut EncloneControl) {
    if !path_exists(&f) {
        eprintln!("\nCan't find the file referenced by your META argument.\n");
        std::process::exit(1);
    }
    let fx = File::open(&f);
    if fx.is_err() {
        eprintln!(
            "\nProblem with META: unable to read from the file\n\
             \"{}\".\nPlease check that that path makes sense and that you have read \
             permission for it.\n",
            f
        );
        std::process::exit(1);
    }
    let f = BufReader::new(fx.unwrap());
    let mut lines = Vec::<String>::new();
    for line in f.lines() {
        let s = line.unwrap();
        lines.push(s);
    }
    proc_meta_core(&lines, &mut ctl);
}
