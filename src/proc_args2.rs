// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::defs::*;
use crate::help1::*;
use crate::help2::*;
use crate::help3::*;
use crate::help4::*;
use crate::help5::*;
use crate::misc1::*;
use crate::proc_args::*;
use io_utils::*;
use itertools::*;
use perf_stats::*;
use pretty_trace::*;
use std::{
    env,
    fs::File,
    io::{BufRead, BufReader},
    time::Instant,
};
use string_utils::*;
use vector_utils::*;

const VERSION_STRING: &'static str = env!("VERSION_STRING");

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Simple arguments.  We test for e.g. PLAIN or PLAIN=, the latter to allow for the case
// where the argument has been set by an environment variable.

pub fn is_simple_arg(arg: &str, x: &str) -> bool {
    if arg == x || arg == &format!("{}=", x) {
        return true;
    } else if arg.starts_with(&format!("{}=", x)) {
        eprintln!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}\".\n",
            arg, x
        );
        std::process::exit(1);
    }
    return false;
}

// Usize arguments.  We require that these are nonnegative integers.

pub fn is_usize_arg(arg: &str, x: &str) -> bool {
    if arg == x {
        eprintln!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}=n\", where n >= 0 is an integer.\n",
            arg, x
        );
        std::process::exit(1);
    } else if arg.starts_with(&format!("{}=", x)) {
        let val = arg.after(&format!("{}=", x)).parse::<usize>();
        if val.is_ok() {
            return true;
        } else {
            eprintln!(
                "\nYour command line includes \"{}\", which is not a valid argument.\n\
                 Perhaps you meant \"{}=n\", where n >= 0 is an integer.\n",
                arg, x
            );
            std::process::exit(1);
        }
    }
    return false;
}

pub fn is_f64_arg(arg: &str, x: &str) -> bool {
    if arg == x {
        eprintln!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}=n\", where n is a floating point number.\n",
            arg, x
        );
        std::process::exit(1);
    } else if arg.starts_with(&format!("{}=", x)) {
        let val = arg.after(&format!("{}=", x)).parse::<f64>();
        if val.is_ok() {
            return true;
        } else {
            eprintln!(
                "\nYour command line includes \"{}\", which is not a valid argument.\n\
                 Perhaps you meant \"{}=n\", where n is a floating point number.\n",
                arg, x
            );
            std::process::exit(1);
        }
    }
    return false;
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Check lvars args.

pub fn check_lvars(ctl: &mut EncloneControl, gex_features: &Vec<Vec<String>>) {
    let mut to_check = Vec::<String>::new();
    for x in ctl.clono_print_opt.lvars.iter() {
        let gpvar = x.starts_with('g') && x.after("g").parse::<usize>().is_ok();
        if !(*x == "datasets"
            || *x == "donors"
            || *x == "ncells"
            || *x == "gex_med"
            || *x == "gex_max"
            || *x == "n_gex"
            || *x == "entropy"
            || *x == "near"
            || *x == "far"
            || *x == "ext"
            || gpvar)
        {
            if !x.ends_with("_g")
                && !x.ends_with("_ab")
                && !x.starts_with("_ag")
                && !x.starts_with("_cr")
                && !x.starts_with("_cu")
                && !x.starts_with("n_")
            {
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
                        known_features.push(format!("{}_ab", ff[z]));
                    } else if ff[2].starts_with("Antigen") {
                        known_features.push(format!("{}_ag", ff[z]));
                    } else if ff[2].starts_with("CRISPR") {
                        known_features.push(format!("{}_cr", ff[z]));
                    } else if ff[2].starts_with("CUSTOM") {
                        known_features.push(format!("{}_cu", ff[z]));
                    } else {
                        known_features.push(format!("{}_g", ff[z]));
                    }
                }
            }
        }
        unique_sort(&mut known_features);
        for i in 0..to_check.len() {
            let x = to_check[i].clone();
            if !bin_member(&known_features, &x) {
                let mut n_var = false;
                if x.starts_with("n_") {
                    n_var = true;
                    let mut indices = Vec::<usize>::new();
                    let mut is_dataset_name = false;
                    let mut is_sample_name = false;
                    let mut is_donor_name = false;
                    let name = x.after("n_").to_string();
                    let s = ctl.sample_info.dataset_path.len();
                    for j in 0..s {
                        if ctl.sample_info.dataset_id[j] == name {
                            is_dataset_name = true;
                            indices.push(j);
                        }
                        if ctl.sample_info.sample_id[j] == name {
                            is_sample_name = true;
                            indices.push(j);
                        }
                        if ctl.sample_info.donor_id[j] == name {
                            is_donor_name = true;
                            indices.push(j);
                        }
                    }
                    let msg = "Suggested reading: \"enclone help input\" and \
                               \"enclone help glossary\".\n";
                    if !is_dataset_name && !is_sample_name && !is_donor_name {
                        eprintln!(
                            "\nYou've used the lead variable {}, and yet {} \
                             does not name a dataset, or a sample,\nor a donor.\n{}",
                            x, name, msg
                        );
                        std::process::exit(1);
                    }
                    if is_dataset_name && indices.len() > 1 {
                        eprintln!(
                            "\nYou've used the lead variable {}, and yet {} \
                             names more than one dataset.  That's ambiguous.\n{}",
                            x, name, msg
                        );
                        std::process::exit(1);
                    }
                    if is_dataset_name && is_sample_name && is_donor_name {
                        eprintln!(
                            "\nYou've used the lead variable {}, and yet {} \
                             names a dataset, a sample, and a donor.  That's ambiguous.\n{}",
                            x, name, msg
                        );
                        std::process::exit(1);
                    }
                    if is_dataset_name && is_sample_name {
                        eprintln!(
                            "\nYou've used the lead variable {}, and yet {} \
                             names a dataset and a sample.  That's ambiguous.\n{}",
                            x, name, msg
                        );
                        std::process::exit(1);
                    }
                    if is_dataset_name && is_donor_name {
                        eprintln!(
                            "\nYou've used the lead variable {}, and yet {} \
                             names a dataset and a donor.  That's ambiguous.\n{}",
                            x, name, msg
                        );
                        std::process::exit(1);
                    }
                    if is_sample_name && is_donor_name {
                        eprintln!(
                            "\nYou've used the lead variable {}, and yet {} \
                             names a sample and a donor.  That's ambiguous.\n{}",
                            x, name, msg
                        );
                        std::process::exit(1);
                    }
                    // could this get called twice on the same name, and what would that do?
                    ctl.sample_info.name_list.insert(name, indices);
                }
                if !n_var {
                    eprintln!(
                        "\nUnrecognized variable {} for LVARS.  Please type \
                         \"enclone help lvars\".\n",
                        x
                    );
                    std::process::exit(1);
                }
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn setup(mut ctl: &mut EncloneControl, args: &Vec<String>) {
    // Provide help if requested.

    help1(&args);
    help2(&args);
    help3(&args);
    help4(&args);
    help5(&args);

    // Pretest for some options.

    ctl.pretty = true;
    let mut nopretty = false;
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "PLAIN") {
            ctl.pretty = false;
        }
        if is_simple_arg(&args[i], "NOPRETTY") {
            nopretty = true;
        }
        if is_simple_arg(&args[i], "COMP") {
            ctl.comp = true;
        }
        if is_simple_arg(&args[i], "CELLRANGER") {
            ctl.gen_opt.cellranger = true;
        }
    }

    // Test for happening mode and turn on pretty trace.

    if !nopretty {
        let mut happening = 0;
        let mut ctrlc = false;
        for i in 1..args.len() {
            if args[i].starts_with("HAPS=") {
                // should actually test for usize
                happening = args[i].after("HAPS=").force_usize();
            }
            if is_simple_arg(&args[i], "CTRLC") {
                ctrlc = true;
            }
        }
        let thread_message = new_thread_message();
        if happening > 0 {
            PrettyTrace::new()
                .message(&thread_message)
                .profile(happening)
                .whitelist(&vec![
                    "amino",
                    "ansi_escape",
                    "binary_vec_io",
                    "enclone",
                    "equiv",
                    "graph_simple",
                    "io_utils",
                    "marsoc",
                    "mirror_sparse_matrix",
                    "perf_stats",
                    "stats_utils",
                    "stirling_numbers",
                    "string_utils",
                    "tables",
                    "vector_utils",
                ])
                .ctrlc()
                .on();
        } else if ctrlc {
            PrettyTrace::new().message(&thread_message).ctrlc().on();
        } else {
            let exit_message: String;
            if !ctl.gen_opt.cellranger {
                exit_message =
                    format!( "Something has gone badly wrong.  Please check to make \
                    sure that none\nof your input files are corrupted.  If they are all OK, then \
                    you have probably\n\
                    encountered an internal error in enclone.\n\
                    Please email us at enclone@10xgenomics.com, including the traceback shown\n\
                    above and also the following version information:\n\
                    {} = {}.\n\n\
                    Thank you and have a nice day!", 
                    env!("CARGO_PKG_VERSION"), VERSION_STRING );
            } else {
                exit_message = format!(
                    "Something has gone badly wrong.  You have probably \
                     encountered an internal error\nin cellranger.  \
                     Please email us at help@10xgenomics.com, including the traceback\nshown \
                     above."
                );
            }
            PrettyTrace::new().exit_message(&exit_message).on();
            let mut nopager = false;
            for i in 1..args.len() {
                if args[i] == "NOPAGER" {
                    nopager = true;
                }
            }
            setup_pager(!nopager);
        }
    }

    // Process args (and set defaults for them).

    proc_args(&mut ctl, &args);

    // Dump lenas.

    for i in 1..args.len() {
        if is_simple_arg(&args[i], "DUMP_LENAS") {
            let mut x = Vec::<usize>::new();
            for y in ctl.sample_info.dataset_id.iter() {
                x.push(y.force_usize());
            }
            x.sort();
            println!("\n{}\n", x.iter().format(","));
            std::process::exit(0);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_args_tail(ctl: &mut EncloneControl, args: &Vec<String>, internal_run: bool) {
    let tall = Instant::now();
    let mut lvars_specified = false;
    for i in 1..args.len() {
        if args[i].starts_with("LVARS=") {
            lvars_specified = true;
        }
    }
    if !ctl.clono_print_opt.amino.is_empty() {
        ctl.clono_print_opt.cvars.insert(0, "amino".to_string());
    }
    if ctl.gen_opt.mouse && ctl.gen_opt.refname.len() > 0 {
        eprintln!(
            "\nIf you specify REF, please do not also specify MOUSE.  It is enough to\n\
             set REF to a mouse reference sequence.\n"
        );
        std::process::exit(1);
    }

    // Remove "datasets" from lvars if there is only one dataset and LVARS not specified.

    if !lvars_specified && ctl.sample_info.dataset_path.len() == 1 {
        ctl.clono_print_opt.lvars.remove(0);
    }

    // Print command line arguments and dataset summary.

    if !ctl.silent {
        println!("");
        for i in 0..args.len() {
            let mut x = args[i].clone();
            if i == 0 && x.contains("/") {
                x = x.rev_after("/").to_string();
            }
            if i > 0 {
                print!(" ");
            }
            print!("{}", x);
        }
        println!("");
        println!(
            "\nThere are {} datasets from {} donors.",
            ctl.sample_info.dataset_path.len(),
            ctl.sample_info.donors
        );
    }

    // Check for duplicated directory paths.

    let mut dp = ctl.sample_info.dataset_path.clone();
    dp.sort();
    let mut i = 0;
    while i < dp.len() {
        let j = next_diff(&dp, i);
        if j - i > 1 {
            eprintln!("\nInput dataset path {} is duplicated.\n", dp[i]);
            std::process::exit(1);
        }
        i = j;
    }
    if !ctl.silent {
        println!("");
    }

    // Get sample descriptions.  Flaky and particularly flaky when lena args are paths,
    // since it will look in outs for the file.

    let tinv = Instant::now();
    if internal_run {
        ctl.sample_info.descrips.clear();
        for i in 0..ctl.sample_info.dataset_path.len() {
            let mut d = ctl.sample_info.dataset_id[i].clone();
            let dir = ctl.sample_info.dataset_path[i].rev_before("/outs");
            let invo = format!("{}/_invocation", dir);
            if path_exists(&invo) {
                let f = open_for_read![invo];
                for line in f.lines() {
                    let s = line.unwrap();
                    if s.contains("sample_desc ") {
                        d = s.between("\"", "\"").to_string();
                    }
                }
            }
            ctl.sample_info.descrips.push(d);
        }
        if ctl.gen_opt.descrip {
            println!("");
            for i in 0..ctl.sample_info.n() {
                println!(
                    "dataset {} ==> sample {} ==> donor {} ==> dataset descrip = {}",
                    ctl.sample_info.dataset_id[i],
                    ctl.sample_info.sample_id[i],
                    ctl.sample_info.donor_id[i],
                    ctl.sample_info.descrips[i]
                );
            }
        }
    }
    if ctl.comp {
        println!("-- used {:.2} seconds reading invocation", elapsed(&tinv));
    }

    // Set up control datastructure (EncloneControl).  This is stuff that is constant for a given
    // run of enclone.

    if ctl.comp {
        println!("used {:.2} seconds in proc_args_tail", elapsed(&tall));
    }
}
