// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// This file contains the two functions proc_xcr and proc_meta.

use crate::defs::*;
use io_utils::*;
use marsoc::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::Command;
use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_xcr(f: &str, gex: &str, have_gex: bool, ctl: &mut EncloneControl) {
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
        for (is, s) in sample_groups.iter().enumerate() {
            let datasets = (*s).split(',').collect::<Vec<&str>>();
            let mut datasets_gex = Vec::<&str>::new();
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
            for (ix, x) in datasetsx.iter().enumerate() {
                let mut p = (*x).to_string();
                // ◼ In CR 4.0, the way we get to outs below will need to change.

                // Specify the "outs" path p.
                //
                // Use case 1.  PRE is specified.

                if ctl.gen_opt.pre != "" {
                    if !p.ends_with("/outs") {
                        p = format!("{}/{}/outs", ctl.gen_opt.pre, p);
                    } else {
                        p = format!("{}/{}", ctl.gen_opt.pre, p);
                    }

                // Use case 2.  It's an internal run, an id has been provided, and PRE
                // was not specified.  Then we look internally.
                } else if ctl.gen_opt.internal_run
                    && p.parse::<u32>().is_ok()
                    && ctl.gen_opt.pre == ""
                {
                    p = format!("{}", get_outs(&p));

                // Use case 3.  All else.
                } else if !p.ends_with("/outs") {
                    p = format!("{}/outs", p);
                }

                // For internal runs, try much harder.  This is so that internal users can just
                // type an internal numerical id for a dataset and have it always work.
                // The code that's used here should be placed somewhere else.

                if !path_exists(&p) && ctl.gen_opt.internal_run {
                    let url = format!("https://xena.fuzzplex.com/api/analyses/{}", x);
                    let o = Command::new("curl")
                        .arg(url)
                        .output()
                        .expect("failed to execute xena http");
                    let m = String::from_utf8(o.stdout).unwrap();
                    if m.contains("502 Bad Gateway") {
                        panic!("502 Bad Gateway from http://xena/api/analyses/{}", x);
                    }
                    let mut path = String::new();
                    if m.contains("\"path\":\"") {
                        path = m.between("\"path\":\"", "\"").to_string();
                        ctl.gen_opt.current_ref = true;
                    }
                    p = format!("{}/outs", path);
                }

                // Now, possibly, we should remove the /outs suffix.  We do this to allow for the
                // case where the customer has copied Cell Ranger output files but not preserved
                // the directory structure.  Or perhaps they appended /outs to their path.

                if !path_exists(&p) {
                    p = p.rev_before("/outs").to_string();
                }

                // If the path p doesn't exist, we have to give up.

                if !path_exists(&p) {
                    if !f.contains("=") {
                        eprintln!("\nCan't find the path {}.\n", p);
                        std::process::exit(1);
                    } else if ctl.gen_opt.pre != "".to_string() {
                        if !p.contains("/outs") {
                            eprintln!(
                                "\nThe value given for {} on the enclone command line \
                                 includes\n{}, which after prefixing by PRE yields\n\
                                 {},\n\
                                 and that path does not exist.\n",
                                f.before("="),
                                x,
                                p
                            );
                        } else {
                            eprintln!(
                                "\nThe value given for {} on the enclone command line \
                                 includes\n{}, which after prefixing by PRE yields\n\
                                 {},\n\
                                 and that path does not contain a subdirectory outs.\n",
                                f.before("="),
                                x,
                                p.rev_before("/outs")
                            );
                        }
                    } else {
                        eprintln!(
                            "\nThe value given for {} on the enclone command line \
                             includes\n{}, and that path does not contain a subdirectory outs.\n",
                            f.before("="),
                            x
                        );
                    }
                    std::process::exit(1);
                }

                // Now work on the GEX path.

                let mut pg = String::new();
                if have_gex {
                    pg = datasets_gex[ix].to_string();
                }
                if pg != "".to_string() {
                    let pg0 = pg.clone();
                    if ctl.gen_opt.internal_run
                        && ctl.gen_opt.h5
                        && ctl.gen_opt.pre != ""
                        && !path_exists(&format!(
                            "{}/{}/outs/raw_gene_bc_matrices_h5.h5",
                            ctl.gen_opt.pre, pg
                        ))
                        && !path_exists(&format!(
                            "{}/{}/outs/raw_feature_bc_matrix.h5",
                            ctl.gen_opt.pre, pg
                        ))
                        && pg.parse::<u32>().is_ok()
                    {
                        pg = format!("{}", get_outs(&pg));
                    } else if ctl.gen_opt.internal_run
                        && ctl.gen_opt.pre != ""
                        && !path_exists(&format!("{}/{}/outs", ctl.gen_opt.pre, pg))
                        && pg.parse::<u32>().is_ok()
                    {
                        pg = format!("{}", get_outs(&pg));
                    } else if ctl.gen_opt.pre != "" {
                        pg = format!("{}/{}/outs", ctl.gen_opt.pre, pg);
                    } else {
                        if pg.parse::<i32>().is_ok() && ctl.gen_opt.internal_run {
                            pg = format!("{}", get_outs(&pg));
                        } else {
                            pg = format!("{}/outs", pg);
                        }
                    }

                    // Now, possibly, we should remove the /outs suffix, see discussion above.

                    if !path_exists(&pg) {
                        pg = pg.rev_before("/outs").to_string();
                    }

                    // Check for nonexistent path

                    if !path_exists(&pg) {
                        if ctl.gen_opt.pre != "".to_string() {
                            let mut pg = pg.clone();
                            if pg.contains("/outs") {
                                pg = pg.rev_before("/outs").to_string();
                            }
                            eprintln!(
                                "\nThe value given for GEX on the enclone command line \
                                 includes\n{}, which after prefixing by PRE yields\n\
                                 {},\n\
                                 and that path does not contain a subdirectory outs.\n",
                                pg0, pg
                            );
                        } else {
                            eprintln!(
                                "\nThe value given for GEX on the enclone command line \
                                 includes\n{}, and that path does not contain a subdirectory \
                                 outs.\n",
                                pg.rev_before("/outs")
                            );
                        }
                        std::process::exit(1);
                    }
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
                ctl.sample_info
                    .sample_donor
                    .push(HashMap::<String, (String, String)>::new());
                ctl.sample_info.tag.push(HashMap::<String, String>::new());
                ctl.sample_info
                    .barcode_color
                    .push(HashMap::<String, String>::new());
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_meta(f: &str, ctl: &mut EncloneControl) {
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

            // Parse bc.

            let mut sample_donor = HashMap::<String, (String, String)>::new();
            let mut tag = HashMap::<String, String>::new();
            let mut barcode_color = HashMap::<String, String>::new();
            if bc != "".to_string() {
                if ctl.gen_opt.pre != "".to_string() {
                    bc = format!("{}/{}", ctl.gen_opt.pre, bc);
                }
                if !path_exists(&bc) {
                    eprintln!(
                        "\nIn your META file, a value for bc implies the existence of \
                         a file\n{}\nbut that file does not exist.\n",
                        bc
                    );
                    std::process::exit(1);
                }
                let f = open_for_read![&bc];

                // Parse one bc file.

                let mut first = true;
                let mut fieldnames = Vec::<String>::new();
                let mut barcode_pos = 0;
                let mut sample_pos = 0;
                let mut donor_pos = 0;
                let mut tag_pos = None;
                let mut color_pos = None;
                for line in f.lines() {
                    let s = line.unwrap();
                    if first {
                        let allowed_fields = vec![
                            "barcode".to_string(),
                            "sample".to_string(),
                            "donor".to_string(),
                            "tag".to_string(),
                            "color".to_string(),
                        ];
                        let fields = s.split(',').collect::<Vec<&str>>();
                        for x in fields.iter() {
                            if !allowed_fields.contains(&x.to_string()) {
                                eprintln!(
                                    "\nThe file\n\
                                     {}\n\
                                     from the bc field used in META\n\
                                     has an illegal field name ({}) in its first line.\n",
                                    bc, x
                                );
                                std::process::exit(1);
                            }
                        }
                        let required = vec!["barcode", "sample", "donor"];
                        for f in required.iter() {
                            if !fields.contains(f) {
                                eprintln!(
                                    "\nThe file\n\
                                     {}\n\
                                     from the bc field used in META\n\
                                     is missing the field {}.\n",
                                    bc, f
                                );
                                std::process::exit(1);
                            }
                        }
                        for x in fields.iter() {
                            fieldnames.push(x.to_string());
                        }
                        for i in 0..fields.len() {
                            if fields[i] == "barcode" {
                                barcode_pos = i;
                            } else if fields[i] == "sample" {
                                sample_pos = i;
                            } else if fields[i] == "donor" {
                                donor_pos = i;
                            } else if fields[i] == "tag" {
                                tag_pos = Some(i);
                            } else if fields[i] == "color" {
                                color_pos = Some(i);
                            }
                        }
                        first = false;
                    } else {
                        let fields = s.split(',').collect::<Vec<&str>>();
                        if fields.len() != fieldnames.len() {
                            eprintln!(
                                "\nThere is a line\n{}\n\
                                 in the CSV file defined by bc in META\n\
                                 that has {} fields, which isn't right, because the header line\n\
                                 has {} fields..  This is for the file\n{}\ndefined by bc.\n",
                                s,
                                fields.len(),
                                fieldnames.len(),
                                bc
                            );
                            std::process::exit(1);
                        }
                        if !fields[barcode_pos].contains('-') {
                            eprintln!(
                                "\nThe barcode \"{}\" appears in the file\n{}\ndefined \
                                 by bc in META.  That doesn't make sense because a barcode\n\
                                 should include a hyphen.\n",
                                fields[barcode_pos], bc
                            );
                            std::process::exit(1);
                        }
                        sample_donor.insert(
                            fields[barcode_pos].to_string(),
                            (
                                fields[sample_pos].to_string(),
                                fields[donor_pos].to_string(),
                            ),
                        );
                        if tag_pos.is_some() {
                            let tag_pos = tag_pos.unwrap();
                            tag.insert(
                                fields[barcode_pos].to_string(),
                                fields[tag_pos].to_string(),
                            );
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

            // Finish up.

            if ctl.gen_opt.pre != "".to_string() {
                path = format!("{}/{}/outs", ctl.gen_opt.pre, path);
                if gpath != "".to_string() {
                    gpath = format!("{}/{}/outs", ctl.gen_opt.pre, gpath);
                }
            } else {
                path = format!("{}/outs", path);
                if gpath != "".to_string() {
                    gpath = format!("{}/outs", gpath);
                }
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
            ctl.sample_info.sample_donor.push(sample_donor);
            ctl.sample_info.tag.push(tag);
            ctl.sample_info.barcode_color.push(barcode_color);
        }
    }
}
