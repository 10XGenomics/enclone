// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::defs::*;
use io_utils::*;
use marsoc::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_xcr(f: &str, gex: &str, have_gex: bool, internal_run: bool, ctl: &mut EncloneControl) {
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
                    p = format!("{}/{}/outs", ctl.gen_opt.pre, p);

                // Use case 2.  It's an internal run, a lena id has been provided, and PRE
                // was not specified.  Then we look on marsoc.
                } else if internal_run && p.parse::<u32>().is_ok() && ctl.gen_opt.pre == "" {
                    p = format!("{}", get_outs(&p));

                // Use case 3.  All else.
                } else {
                    p = format!("{}/outs", p);
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
                        eprintln!(
                            "\nThe value given for {} on the enclone command line \
                             includes\n{}, which after prefixing by PRE yields\n\
                             {},\n\
                             and that path does not contain a subdirectory outs.\n",
                            f.before("="),
                            x,
                            p.rev_before("/outs")
                        );
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
                    if internal_run
                        && ctl.gen_opt.pre != ""
                        && !path_exists(&format!("{}/{}/outs", ctl.gen_opt.pre, pg))
                        && pg.parse::<u32>().is_ok()
                    {
                        pg = format!("{}", get_outs(&pg));
                    } else if ctl.gen_opt.pre != "" {
                        pg = format!("{}/{}/outs", ctl.gen_opt.pre, pg);
                    } else {
                        if pg.parse::<i32>().is_ok() && internal_run {
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
                            eprintln!(
                                "\nThe value given for GEX on the enclone command line \
                                 includes\n{}, which after prefixing by PRE yields\n\
                                 {},\n\
                                 and that path does not contain a subdirectory outs.\n",
                                pg0,
                                pg.rev_before("/outs")
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
                ctl.sample_info.donor_index.push(id);
                ctl.sample_info.donor_id.push(donor_name);
                ctl.sample_info.sample_id.push(sample_name);
                ctl.sample_info
                    .sample_donor
                    .push(HashMap::<String, (String, String)>::new());
            }
        }
    }
    let mut i = 0;
    while i < ctl.sample_info.donor_index.len() {
        let j = next_diff(&ctl.sample_info.donor_index, i);
        let mut x = Vec::<usize>::new();
        for k in i..j {
            x.push(k);
        }
        i = j;
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
                } else if *x == "bc" && y.len() > 0 {
                    bc = y.to_string();
                }
            }
            if bc != "" && ( fields.contains(&"sample".to_string()) || fields.contains(&"donor".to_string()) ) {
                eprintln!(
                    "\nIf bc is specified in META, for a given dataset, it does not\n\
                     make sense to also specify sample or donor.\n"
                );
                std::process::exit(1);
            }
            let mut sample_donor = HashMap::<String, (String, String)>::new();
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
                let mut first = true;
                for line in f.lines() {
                    let s = line.unwrap();
                    if first {
                        if s != "barcode,sample,donor".to_string() {
                            eprintln!(
                                "\nThe first line in the CSV file defined by bc in META \
                                 must be\nbarcode,sample,donor\nbut it's not.  This is for the \
                                 file\n{}\ndefined by bc.\n",
                                bc
                            );
                            std::process::exit(1);
                        }
                        first = false;
                    } else {
                        let fields = s.split(',').collect::<Vec<&str>>();
                        if fields.len() != 3 {
                            eprintln!(
                                "\nThere is a line in the CSV file defined by bc in META\n\
                                 that does not have three fields.  That's wrong.  This is for the \
                                 file\n{}\ndefined by bc.\n",
                                bc
                            );
                            std::process::exit(1);
                        }
                        if !fields[0].contains('-') {
                            eprintln!(
                                "\nThe barcode \"{}\" appears in the file\n{}\ndefined \
                                 by bc in META.  That doesn't make sense because a barcode\n\
                                 should include a hyphen.\n",
                                fields[0], bc
                            );
                            std::process::exit(1);
                        }
                        sample_donor.insert(
                            fields[0].to_string(),
                            (fields[1].to_string(), fields[2].to_string()),
                        );
                    }
                }
            }
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
                    ctl.sample_info.donor_index.push(j);
                    break;
                }
            }
            if dp.is_none() {
                ctl.sample_info.donor_index.push(donors.len());
                donors.push(donor.clone());
            }
            ctl.sample_info.descrips.push(abbr.clone());
            ctl.sample_info.dataset_path.push(path);
            ctl.sample_info.gex_path.push(gpath);
            ctl.sample_info.dataset_id.push(abbr);
            ctl.sample_info.donor_id.push(donor);
            ctl.sample_info.sample_id.push(sample);
            ctl.sample_info.sample_donor.push(sample_donor);
        }
    }
}
