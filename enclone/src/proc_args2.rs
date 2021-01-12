// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::*;
use io_utils::*;
use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    time::Instant,
};
use string_utils::*;
use vector_utils::*;

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

pub fn is_string_arg(arg: &str, x: &str) -> bool {
    if arg == x {
        eprintln!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}=s\" for some string s.\n",
            arg, x
        );
        std::process::exit(1);
    } else if arg.starts_with(&format!("{}=", x)) {
        return true;
    }
    return false;
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

    if !lvars_specified && ctl.origin_info.dataset_path.len() == 1 {
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
            ctl.origin_info.dataset_path.len(),
            ctl.origin_info.donors
        );
    }

    // Check for duplicated directory paths.

    let mut dp = ctl.origin_info.dataset_path.clone();
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

    // Get origin descriptions.  Flaky and particularly flaky when internal origin args are paths,
    // since it will look in outs for the file.

    if ctl.gen_opt.internal_run || ctl.gen_opt.descrip {
        ctl.origin_info.descrips.clear();
        let mut results = vec![(0, "".to_string()); ctl.origin_info.n()];
        for i in 0..ctl.origin_info.n() {
            results[i].0 = i;
        }
        results.par_iter_mut().for_each(|res| {
            let i = res.0;
            let mut d = ctl.origin_info.dataset_id[i].clone();
            let mut dir = ctl.origin_info.dataset_path[i].clone();
            if dir.ends_with("/outs") {
                dir = dir.rev_before("/outs").to_string();
            }
            let invo = format!("{}/_invocation", dir);
            if path_exists(&invo) {
                let f = open_for_read![invo];
                for line in f.lines() {
                    let s = line.unwrap();
                    // Leave sample_desc alone for internal architecture!
                    if s.contains("sample_desc ") {
                        d = s.between("\"", "\"").to_string();
                    }
                }
            }
            res.1 = d;
        });
        for i in 0..ctl.origin_info.dataset_path.len() {
            ctl.origin_info.descrips.push(results[i].1.clone());
        }
        if ctl.gen_opt.descrip {
            println!("");
            for i in 0..ctl.origin_info.n() {
                if i > 0 {
                    println!("");
                }
                println!(
                    "dataset {} ==> origin {} ==> donor {} ==> dataset descrip = {}",
                    ctl.origin_info.dataset_id[i],
                    // origin_id and donor_id don't make sense if bc specified in META
                    ctl.origin_info.origin_id[i],
                    ctl.origin_info.donor_id[i],
                    ctl.origin_info.descrips[i]
                );
                println!("vdj path = {}", ctl.origin_info.dataset_path[i]);
                if !ctl.origin_info.gex_path.is_empty() {
                    println!("gex path = {}", ctl.origin_info.gex_path[i]);
                }
            }
        }
    }
    ctl.perf_stats(&tall, "in proc_args_tail");
}
