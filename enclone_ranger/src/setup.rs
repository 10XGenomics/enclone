// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use crate::USING_PAGER;
use enclone_args::proc_args::proc_args;
use enclone_args::proc_args2::is_simple_arg;
use enclone_core::defs::EncloneControl;
use std::sync::atomic::Ordering::SeqCst;
use string_utils::TextUtils;
use vector_utils::erase_if;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Process some arguments.  The order is delicate.

pub fn critical_args(args: &Vec<String>, ctl: &mut EncloneControl) -> Result<Vec<String>, String> {
    // Form the combined set of command-line arguments and "command-line" arguments
    // implied by environment variables.

    let args = args.clone();

    // Check for CELLRANGER.

    for i in 1..args.len() {
        if is_simple_arg(&args[i], "CELLRANGER")? {
            ctl.gen_opt.cellranger = true;
        }
    }

    for i in 1..args.len() {
        if args[i] == *"FORCE_EXTERNAL" {
            ctl.gen_opt.internal_run = false;
        }
    }

    // Determine PRE.

    for i in 1..args.len() {
        if args[i].starts_with("PRE=") {
            let pre = args[i].after("PRE=").split(',').collect::<Vec<&str>>();
            ctl.gen_opt.pre.clear();
            for x in pre.iter() {
                ctl.gen_opt.pre.push(x.to_string());
            }
        }
    }
    Ok(args)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn setup(
    mut ctl: &mut EncloneControl,
    args: &Vec<String>,
    argsx: &mut Vec<String>,
    args_orig: &Vec<String>,
) -> Result<(), String> {
    let using_pager = false;
    // Provide help if requested.

    {
        let mut args = args.clone();
        let mut to_delete = vec![false; args.len()];
        for i in 1..args.len() {
            if args[i] == "NOPAGER" {
                ctl.gen_opt.nopager = true;
                to_delete[i] = true;
            } else if args[i] == "FORCE_EXTERNAL" {
                to_delete[i] = true;
            } else if args[i].starts_with("MAX_CORES=") {
                to_delete[i] = true;
            } else if args[i].starts_with("PRE=") {
                to_delete[i] = true;
            }
        }

        // Proceed.

        if ctl.gen_opt.html && ctl.gen_opt.svg {
            return Err("\nBoth HTML and SVG cannot be used at the same time.\n".to_string());
        }
        erase_if(&mut args, &to_delete);
        *argsx = args.clone();
        let mut argsx = Vec::<String>::new();
        for i in 0..args_orig.len() {
            if args_orig[i] != "HTML"
                && args_orig[i] != "NOPAGER"
                && args_orig[i] != "FORCE_EXTERNAL"
                && args_orig[i] != "NO_KILL"
                && !args_orig[i].starts_with("PRE=")
                && !args_orig[i].starts_with("MAX_CORES=")
            {
                argsx.push(args_orig[i].clone());
            }
        }
    }

    // Pretest for some options.

    ctl.pretty = true;
    ctl.gen_opt.h5 = true;
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "PLAIN")? {
            ctl.pretty = false;
        }
        if is_simple_arg(&args[i], "NH5")? {
            ctl.gen_opt.h5 = false;
        }
    }

    // Turn on pretty trace.

    USING_PAGER.store(using_pager, SeqCst);

    // Process args (and set defaults for them).

    proc_args(&mut ctl, args)?;
    Ok(())
}
