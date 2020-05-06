// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This file contains the two functions proc_xcr and proc_meta.

use enclone_core::defs::*;
use io_utils::*;
use itertools::Itertools;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::Command;
use string_utils::*;
use tilde_expand::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Functions to find the path to data.

fn get_path(p: &str, ctl: &EncloneControl) -> String {
    for x in ctl.gen_opt.pre.iter() {
        let mut pp = format!("{}/{}", x, p);
        if pp.starts_with("~") {
            pp = stringme(&tilde_expand(&pp.as_bytes()));
        }
        if path_exists(&pp) {
            return pp;
        }
    }
    let mut pp = p.to_string();
    if pp.starts_with("~") {
        pp = stringme(&tilde_expand(&pp.as_bytes()));
    }
    pp
}

pub fn get_path_fail(p: &str, ctl: &EncloneControl, source: &str) -> String {
    for x in ctl.gen_opt.pre.iter() {
        let pp = format!("{}/{}", x, p);
        if path_exists(&pp) {
            return pp;
        }
    }
    if !path_exists(&p) {
        if ctl.gen_opt.pre.is_empty() {
            eprintln!(
                "\nUnable to find the path {}.  This came from the {} argument.\n",
                p, source
            );
        } else {
            eprintln!(
                "\nUnable to find the path {}, even if prepended by any of the directories \
                in\nPRE={}.\nThis came from the {} argument.\n",
                p,
                ctl.gen_opt.pre.iter().format(","),
                source
            );
        }
        std::process::exit(1);
    }
    p.to_string()
}

fn get_path_or_internal_id(p: &str, ctl: &mut EncloneControl, source: &str) -> String {
    let mut pp = get_path(&p, &ctl);
    if !path_exists(&pp) {
        if !ctl.gen_opt.internal_run {
            get_path_fail(&pp, &ctl, source);
        } else {
            // For internal runs, try much harder.  This is so that internal users can
            // just type an internal numerical id for a dataset and have it always
            // work.  The code that's used here should be placed somewhere else.

            if p.parse::<usize>().is_ok() {
                let url = format!("https://xena.fuzzplex.com/api/analyses/{}", p);
                let o = Command::new("curl")
                    .arg(url)
                    .output()
                    .expect("failed to execute xena http");
                let m = String::from_utf8(o.stdout).unwrap();
                if m.contains("502 Bad Gateway") {
                    eprintln!(
                        "\nWell this is sad.  The URL \
                        http://xena/api/analyses/{} yielded a 502 Bad Geteway \
                        message.  Either try again later or ask someone for help.\n",
                        p
                    );
                    std::process::exit(1);
                }
                if m.contains("\"path\":\"") {
                    let path = m.between("\"path\":\"", "\"").to_string();
                    ctl.gen_opt.current_ref = true;
                    pp = format!("{}/outs", path);
                    if !path_exists(&pp) {
                        eprintln!(
                            "\nIt looks like you've provided a xena analysis id for \
                            which the pipeline outs folder\n{}\nhas not yet been generated.\n\n",
                            p
                        );
                        std::process::exit(1);
                    }
                } else {
                    eprintln!(
                        "\nIt looks like you've provided either an incorrect \
                        xena id {} or else one for which\n\
                        the pipeline outs folder has not yet been generated.\n",
                        p
                    );
                    std::process::exit(1);
                }
            } else {
                eprintln!(
                    "\nAfter searching high and low, your path for {} \
                    cannot be found.\nPlease check its value and also the value \
                    for PRE if you provided that.\n",
                    source
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
    let mut sample_for_bc = HashMap::<String, String>::new();
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
        let (mut sample_pos, mut donor_pos, mut tag_pos, mut color_pos) = (None, None, None, None);
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
                    eprintln!(
                        "\nThe file\n{}\n{}\nis missing the barcode field.\n",
                        bc, origin,
                    );
                    std::process::exit(1);
                }
                for x in fields.iter() {
                    fieldnames.push(x.to_string());
                }
                for i in 0..fields.len() {
                    if fields[i] == "barcode" {
                        barcode_pos = i;
                    } else if fields[i] == "sample" {
                        sample_pos = Some(i);
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
                if sample_pos.is_some() {
                    sample_for_bc.insert(
                        fields[barcode_pos].to_string(),
                        fields[sample_pos.unwrap()].to_string(),
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
    ctl.sample_info.sample_for_bc.push(sample_for_bc);
    ctl.sample_info.donor_for_bc.push(donor_for_bc);
    ctl.sample_info.tag.push(tag);
    ctl.sample_info.barcode_color.push(barcode_color);
    ctl.sample_info.alt_bc_fields.push(alt_bc_fields);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_xcr(f: &str, gex: &str, bc: &str, have_gex: bool, mut ctl: &mut EncloneControl) {
    ctl.sample_info = SampleInfo::default();
    if (ctl.gen_opt.tcr && f.starts_with("BCR=")) || (ctl.gen_opt.bcr && f.starts_with("TCR=")) {
        eprintln!("\nOnly one of TCR or BCR can be specified.\n");
        std::process::exit(1);
    }
    ctl.gen_opt.tcr = f.starts_with("TCR=");
    ctl.gen_opt.bcr = f.starts_with("BCR=");
    let val: String;
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
    let donor_groups = val.split(';').collect::<Vec<&str>>();
    let donor_groups_gex = gex.split(';').collect::<Vec<&str>>();
    let donor_groups_bc = bc.split(';').collect::<Vec<&str>>();
    let mut xcr = "TCR".to_string();
    if ctl.gen_opt.bcr {
        xcr = "BCR".to_string();
    }
    if have_gex && donor_groups_gex.len() != donor_groups.len() {
        eprintln!(
            "\nThe {} and GEX arguments do not exactly mirror each \
             other's structure.\n",
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
    for (id, d) in donor_groups.iter().enumerate() {
        let sample_groups = (*d).split(':').collect::<Vec<&str>>();
        let mut sample_groups_gex = Vec::<&str>::new();
        if have_gex {
            sample_groups_gex = donor_groups_gex[id].split(':').collect::<Vec<&str>>();
            if sample_groups_gex.len() != sample_groups.len() {
                eprintln!(
                    "\nThe {} and GEX arguments do not exactly mirror each \
                     other's structure.\n",
                    xcr
                );
                std::process::exit(1);
            }
        }
        let mut sample_groups_bc = Vec::<&str>::new();
        if !bc.is_empty() {
            sample_groups_bc = donor_groups_bc[id].split(':').collect::<Vec<&str>>();
            if sample_groups_bc.len() != sample_groups.len() {
                eprintln!(
                    "\nThe {} and BC arguments do not exactly mirror each \
                     other's structure.\n",
                    xcr
                );
                std::process::exit(1);
            }
        }
        for (is, s) in sample_groups.iter().enumerate() {
            let datasets = (*s).split(',').collect::<Vec<&str>>();
            let mut datasets_gex = Vec::<&str>::new();
            let mut datasets_bc = Vec::<&str>::new();
            let mut datasetsx = Vec::<String>::new();
            for i in 0..datasets.len() {
                if datasets[i].contains('-') {
                    let (f1, f2) = (datasets[i].before("-"), datasets[i].after("-"));
                    if f1.parse::<usize>().is_ok() && f2.parse::<usize>().is_ok() {
                        let (l1, l2) = (f1.force_usize(), f2.force_usize());
                        if l1 <= l2 {
                            for l in l1..=l2 {
                                datasetsx.push(format!("{}", l));
                            }
                        } else {
                            eprintln!("\nIllegal argument {}.\n", datasets[i]);
                            std::process::exit(1);
                        }
                    } else {
                        datasetsx.push(datasets[i].to_string());
                    }
                } else {
                    datasetsx.push(datasets[i].to_string());
                }
            }
            if have_gex {
                datasets_gex = sample_groups_gex[is].split(',').collect::<Vec<&str>>();
                if datasets_gex.len() != datasetsx.len() {
                    eprintln!(
                        "\nThe {} and GEX arguments do not exactly mirror each \
                         other's structure.\n",
                        xcr
                    );
                    std::process::exit(1);
                }
            }
            if !bc.is_empty() {
                datasets_bc = sample_groups_bc[is].split(',').collect::<Vec<&str>>();
                if datasets_bc.len() != datasetsx.len() {
                    eprintln!(
                        "\nThe {} and BC arguments do not exactly mirror each \
                         other's structure.\n",
                        xcr
                    );
                    std::process::exit(1);
                }
            }
            for (ix, x) in datasetsx.iter().enumerate() {
                let mut p = (*x).to_string();
                // ◼ In CR 4.0, the way we get to outs below will need to change.

                let mut source = f.clone();
                if f.contains('=') {
                    source = f.before("=");
                }
                p = get_path_or_internal_id(&p, &mut ctl, source);

                // Now work on the BC path.

                let mut bcx = String::new();
                if !bc.is_empty() {
                    bcx = datasets_bc[ix].to_string();
                }
                parse_bc(bcx, &mut ctl, "BC");

                // Now work on the GEX path.

                let mut pg = String::new();
                if have_gex {
                    pg = datasets_gex[ix].to_string();
                    pg = get_path_or_internal_id(&pg, &mut ctl, "GEX");
                }

                // OK everything worked, all set.

                let donor_name = format!("d{}", id + 1);
                let sample_name = format!("s{}", is + 1);
                let mut dataset_name = (*x).to_string();
                if dataset_name.contains('/') {
                    dataset_name = dataset_name.rev_after("/").to_string();
                }
                ctl.sample_info.descrips.push(dataset_name.clone());
                ctl.sample_info.dataset_path.push(p);
                ctl.sample_info.gex_path.push(pg);
                ctl.sample_info.dataset_id.push(dataset_name.clone());
                ctl.sample_info.donor_id.push(donor_name);
                ctl.sample_info.color.push("".to_string());
                ctl.sample_info.sample_id.push(sample_name);
                ctl.sample_info.tag.push(HashMap::<String, String>::new());
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

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
    let mut fields = Vec::<String>::new();
    let mut donors = Vec::<String>::new();
    for (count, line) in f.lines().enumerate() {
        let s = line.unwrap();
        if count == 0 {
            let x = s.split(',').collect::<Vec<&str>>();
            for i in 0..x.len() {
                fields.push(x[i].to_string());
            }
            let mut fields_sorted = fields.clone();
            unique_sort(&mut fields_sorted);
            if fields_sorted.len() < fields.len() {
                eprintln!(
                    "\nThe CSV file that you specified using the META argument \
                     has duplicate field names\nin its first line.\n"
                );
                std::process::exit(1);
            }
            let allowed_fields = vec![
                "bc".to_string(),
                "bcr".to_string(),
                "donor".to_string(),
                "gex".to_string(),
                "sample".to_string(),
                "tcr".to_string(),
                "color".to_string(),
            ];
            for x in fields.iter() {
                if !allowed_fields.contains(&x) {
                    eprintln!(
                        "\nThe CSV file that you specified using the META argument \
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
                    "\nThe CSV file that you specified using the META argument \
                     has neither the field tcr or bcr in its first line.\n"
                );
                std::process::exit(1);
            }
            if ctl.gen_opt.tcr && ctl.gen_opt.bcr {
                eprintln!(
                    "\nThe CSV file that you specified using the META argument \
                     has both the fields tcr and bcr in its first line.\n"
                );
                std::process::exit(1);
            }
        } else if !s.starts_with('#') {
            let val = s.split(',').collect::<Vec<&str>>();
            if val.len() != fields.len() {
                eprintln!(
                    "\nMETA file line {} has a different number of fields than the \
                     first line of the file.\n",
                    count + 1
                );
                std::process::exit(1);
            }
            let mut path = String::new();
            let mut abbr = String::new();
            let mut gpath = String::new();
            let mut sample = "s1".to_string();
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
                } else if *x == "sample" {
                    sample = y.to_string();
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
            path = get_path_or_internal_id(&path, &mut ctl, "META");
            if gpath.len() > 0 {
                gpath = get_path_or_internal_id(&gpath, &mut ctl, "META");
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
            ctl.sample_info.descrips.push(abbr.clone());
            ctl.sample_info.dataset_path.push(path);
            ctl.sample_info.gex_path.push(gpath);
            ctl.sample_info.dataset_id.push(abbr);
            ctl.sample_info.donor_id.push(donor);
            ctl.sample_info.sample_id.push(sample);
            ctl.sample_info.color.push(color);
        }
    }
}
