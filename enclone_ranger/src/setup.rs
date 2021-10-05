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

pub fn critical_args(args: &Vec<String>, ctl: &mut EncloneControl) -> Result<Vec<String>, String> {
    let args = args.clone();
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
) -> Result<(), String> {
    let using_pager = false;
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
        erase_if(&mut args, &to_delete);
        *argsx = args.clone();
    }

    // Pretest for some options.

    ctl.pretty = true;
    ctl.gen_opt.h5 = true;
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "PLAIN")? {
            ctl.pretty = false;
        }
    }
    USING_PAGER.store(using_pager, SeqCst);

    // Process args (and set defaults for them).

    proc_args(&mut ctl, args)?;
    Ok(())
}
