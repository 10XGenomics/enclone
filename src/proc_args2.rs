// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::defs::*;
use crate::help1::*;
use crate::help2::*;
use crate::help3::*;
use crate::help4::*;
use crate::help5::*;
use crate::help_utils::*;
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

pub fn setup(mut ctl: &mut EncloneControl, args: &Vec<String>) {
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "STABLE_DOC") {
            ctl.gen_opt.stable_doc = true;
        }
    }

    // Provide help if requested.

    {
        if args.len() == 2 && (args[1] == "version" || args[1] == "--version") {
            println!("{} : {}", env!("CARGO_PKG_VERSION"), VERSION_STRING);
            std::process::exit(0);
        }
        let mut args = args.clone();
        let mut to_delete = vec![false; args.len()];
        let mut nopager = false;
        let mut plain = false;
        let mut long_help = false;
        for i in 1..args.len() {
            if args[i] == "NOPAGER" {
                nopager = true;
                to_delete[i] = true;
            } else if args[i] == "HTML" {
                ctl.gen_opt.html = true;
                to_delete[i] = true;
            } else if args[i] == "SVG" {
                ctl.gen_opt.svg = true;
                to_delete[i] = true;
            } else if args[i] == "FORCE_EXTERNAL" {
                to_delete[i] = true;
            } else if args[i] == "LONG_HELP" {
                long_help = true;
                to_delete[i] = true;
            } else if args[i].starts_with("MAX_CORES=") {
                to_delete[i] = true;
            } else if args[i].starts_with("PRE=") {
                to_delete[i] = true;
            } else if args[i] == "PLAIN" {
                to_delete[i] = true;
                plain = true;
                unsafe {
                    PLAIN = true;
                }
            }
        }
        if ctl.gen_opt.html && ctl.gen_opt.svg {
            eprintln!("\nBoth HTML and SVG cannot be used at the same time.\n");
            std::process::exit(1);
        }
        erase_if(&mut args, &to_delete);
        if args.len() == 1 || args.contains(&"help".to_string()) {
            PrettyTrace::new().on();
            setup_pager(!nopager);
        }
        let mut help_all = false;
        if args.len() >= 3 && args[1] == "help" && args[2] == "all" {
            unsafe {
                HELP_ALL = true;
            }
            help_all = true;
        }
        let mut h = HelpDesk::new(plain, help_all, long_help, ctl.gen_opt.html);
        help1(&args, &mut h);
        help2(&args, &ctl, &mut h);
        help3(&args, &mut h);
        help4(&args, &mut h);
        help5(&args, &ctl, &mut h);
    }

    // Pretest for some options.

    ctl.pretty = true;
    let mut nopretty = false;
    ctl.gen_opt.h5 = true;
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
        if is_simple_arg(&args[i], "COMP2") {
            ctl.comp = true;
            ctl.comp2 = true;
        }
        if is_simple_arg(&args[i], "CELLRANGER") {
            ctl.gen_opt.cellranger = true;
        }
        if is_simple_arg(&args[i], "NH5") {
            ctl.gen_opt.h5 = false;
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
                exit_message = format!(
                    "Something has gone badly wrong.  You have probably encountered an internal \
                    error in enclone.\n\n\
                    Please email us at enclone@10xgenomics.com, including the traceback shown\n\
                    above and also the following version information:\n\
                    {} : {}.\n\n\
                    Thank you and have a nice day!",
                    env!("CARGO_PKG_VERSION"),
                    VERSION_STRING
                );
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

pub fn proc_args_tail(ctl: &mut EncloneControl, args: &Vec<String>) {
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
    if ctl.gen_opt.internal_run {
        ctl.sample_info.descrips.clear();
        for i in 0..ctl.sample_info.dataset_path.len() {
            let mut d = ctl.sample_info.dataset_id[i].clone();
            let mut dir = ctl.sample_info.dataset_path[i].clone();
            if dir.ends_with("/outs") {
                dir = dir.rev_before("/outs").to_string();
            }
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
                if i > 0 {
                    println!("");
                }
                println!(
                    "dataset {} ==> sample {} ==> donor {} ==> dataset descrip = {}",
                    ctl.sample_info.dataset_id[i],
                    // sample_id and donor_id don't make sense if bc specified in META
                    ctl.sample_info.sample_id[i],
                    ctl.sample_info.donor_id[i],
                    ctl.sample_info.descrips[i]
                );
                println!("vdj path = {}", ctl.sample_info.dataset_path[i]);
                if !ctl.sample_info.gex_path.is_empty() {
                    println!("gex path = {}", ctl.sample_info.gex_path[i]);
                }
            }
        }
    }
    if ctl.comp2 {
        println!("-- used {:.2} seconds reading invocation", elapsed(&tinv));
    }

    // Set up control datastructure (EncloneControl).  This is stuff that is constant for a given
    // run of enclone.

    if ctl.comp {
        if !ctl.comp2 {
            println!("");
        }
        println!("used {:.2} seconds in proc_args_tail", elapsed(&tall));
    }
}
