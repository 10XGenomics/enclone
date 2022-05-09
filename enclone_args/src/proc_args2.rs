// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::EncloneControl;
use io_utils::{open_userfile_for_read, path_exists};
use rayon::prelude::*;
use std::fs::{remove_file, File};
use std::{io::BufRead, time::Instant};
use string_utils::TextUtils;
use vector_utils::next_diff;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test a file for writeability by writing and then deleting it.

pub fn test_writeable(val: &str, evil_eye: bool) -> Result<(), String> {
    if evil_eye {
        println!("creating file {} to test writability", val);
    }
    let f = File::create(&val);
    if f.is_err() {
        let mut msgx = format!(
            "\nYou've specified an output file\n{}\nthat cannot be written.\n",
            val
        );
        if val.contains('/') {
            let dir = val.rev_before("/");
            let msg;
            if path_exists(dir) {
                msg = "exists";
            } else {
                msg = "does not exist";
            }
            msgx += &mut format!("Note that the path {} {}.\n", dir, msg);
        }
        return Err(msgx);
    }
    if evil_eye {
        println!("removing file {}", val);
    }
    remove_file(&val).unwrap_or_else(|_| panic!("could not remove file {}", val));
    if evil_eye {
        println!("removal of file {} complete", val);
    }
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Simple arguments.  We test for e.g. PLAIN or PLAIN=, the latter to allow for the case
// where the argument has been set by an environment variable.

pub fn is_simple_arg(arg: &str, x: &str) -> Result<bool, String> {
    if arg == x || arg == &format!("{}=", x) {
        return Ok(true);
    } else if arg.starts_with(&format!("{}=", x)) {
        return Err(format!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}\".\n",
            arg, x
        ));
    }
    Ok(false)
}

// Usize arguments.  We require that these are nonnegative integers.

pub fn is_usize_arg(arg: &str, x: &str) -> Result<bool, String> {
    if arg == x {
        return Err(format!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}=n\", where n >= 0 is an integer.\n",
            arg, x
        ));
    } else if arg.starts_with(&format!("{}=", x)) {
        let val = arg.after(&format!("{}=", x)).parse::<usize>();
        if val.is_ok() {
            return Ok(true);
        } else {
            return Err(format!(
                "\nYour command line includes \"{}\", which is not a valid argument.\n\
                 Perhaps you meant \"{}=n\", where n >= 0 is an integer.\n",
                arg, x
            ));
        }
    }
    Ok(false)
}

// Usize arguments.  We require that these are nonnegative integers.

pub fn is_i32_arg(arg: &str, x: &str) -> Result<bool, String> {
    if arg == x {
        return Err(format!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}=n\", where n >= 0 is an integer.\n",
            arg, x
        ));
    } else if arg.starts_with(&format!("{}=", x)) {
        let val = arg.after(&format!("{}=", x)).parse::<i32>();
        if val.is_ok() {
            return Ok(true);
        } else {
            return Err(format!(
                "\nYour command line includes \"{}\", which is not a valid argument.\n\
                 Perhaps you meant \"{}=n\", where n is an integer.\n",
                arg, x
            ));
        }
    }
    Ok(false)
}

pub fn is_f64_arg(arg: &str, x: &str) -> Result<bool, String> {
    if arg == x {
        return Err(format!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}=n\", where n is a floating point number.\n",
            arg, x
        ));
    } else if arg.starts_with(&format!("{}=", x)) {
        let val = arg.after(&format!("{}=", x)).parse::<f64>();
        if val.is_ok() {
            return Ok(true);
        } else {
            return Err(format!(
                "\nYour command line includes \"{}\", which is not a valid argument.\n\
                 Perhaps you meant \"{}=n\", where n is a floating point number.\n",
                arg, x
            ));
        }
    }
    Ok(false)
}

pub fn is_string_arg(arg: &str, x: &str) -> Result<bool, String> {
    if arg == x {
        return Err(format!(
            "\nYour command line includes \"{}\", which is not a valid argument.\n\
             Perhaps you meant \"{}=s\" for some string s.\n",
            arg, x
        ));
    } else if arg.starts_with(&format!("{}=", x)) {
        return Ok(true);
    }
    Ok(false)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_args_tail(ctl: &mut EncloneControl, args: &Vec<String>) -> Result<(), String> {
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
    if ctl.gen_opt.mouse && !ctl.gen_opt.refname.is_empty() {
        return Err(
            "\nIf you specify REF, please do not also specify MOUSE.  It is enough to\n\
             set REF to a mouse reference sequence.\n"
                .to_string(),
        );
    }

    // Remove "datasets" from lvars if there is only one dataset and LVARS not specified.

    if !lvars_specified && ctl.origin_info.dataset_path.len() == 1 {
        ctl.clono_print_opt.lvars.remove(0);
    }

    // Print command line arguments and dataset summary.

    if !ctl.silent {
        println!();
        for i in 0..args.len() {
            let mut x = args[i].clone();
            if i == 0 && x.contains('/') {
                x = x.rev_after("/").to_string();
            }
            if i > 0 {
                print!(" ");
            }
            print!("{}", x);
        }
        println!();
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
            return Err(format!("\nInput dataset path {} is duplicated.\n", dp[i]));
        }
        i = j;
    }
    if !ctl.silent {
        println!();
    }

    // Get origin descriptions.  Flaky and particularly flaky when internal origin args are paths,
    // since it will look in outs for the file.

    if ctl.gen_opt.internal_run || ctl.gen_opt.descrip || ctl.visual_mode || ctl.gen_opt.vis_dump {
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
            let mut invo = format!("{}/_invocation", dir);
            if !path_exists(&invo) {
                invo = format!("{}/../../_invocation", dir);
            }
            if !path_exists(&invo) {
                invo = format!("{}/../../../_invocation", dir);
            }
            if path_exists(&invo) {
                let f = open_userfile_for_read(&invo);
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
            println!();
            for i in 0..ctl.origin_info.n() {
                if i > 0 {
                    println!();
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
    Ok(())
}
