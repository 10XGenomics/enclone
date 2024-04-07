// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use enclone::misc1::setup_pager;
use enclone_args::proc_args::proc_args;
use enclone_args::proc_args2::is_simple_arg;
use enclone_build::set_panic_handler;
use enclone_core::defs::EncloneControl;
use enclone_core::{require_readable_file, tilde_expand_me};
use enclone_help::help1::help1;
use enclone_help::help2::help2;
use enclone_help::help3::help3;
use enclone_help::help4::help4;
use enclone_help::help5::help5;
use enclone_help::help_utils::HelpDesk;
use io_utils::{open_for_read, path_exists};
use itertools::Itertools;
use std::env;
use std::io::BufRead;
use std::sync::atomic::Ordering::SeqCst;

use string_utils::TextUtils;
use vector_utils::erase_if;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Process some arguments.  The order is delicate.

pub fn critical_args(args: &Vec<String>, ctl: &mut EncloneControl) -> Result<Vec<String>, String> {
    // Form the combined set of command-line arguments and "command-line" arguments
    // implied by environment variables.

    let mut args = args.clone();
    for (key, value) in env::vars() {
        if key.starts_with("ENCLONE_") && !key.starts_with("ENCLONE_VIS_") {
            args.push(format!("{}={}", key.after("ENCLONE_"), value));
        }
    }

    // Check for EVIL_EYE.

    for i in 1..args.len() {
        if args[i] == "EVIL_EYE" {
            ctl.gen_opt.evil_eye = true;
            if ctl.gen_opt.evil_eye {
                println!("the evil eye is on");
            }
        }
    }

    // Determine PRE.
    if !ctl.cr_opt.cellranger {
        let home = dirs::home_dir().unwrap().to_str().unwrap().to_string();
        ctl.gen_opt.pre = vec![
            format!("{}/enclone/datasets_me", home),
            format!("{}/enclone/datasets", home),
            format!("{}/enclone/datasets2", home),
        ];
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
    for i in (1..args.len()).rev() {
        if args[i].starts_with("PREPOST=") {
            let prepost = args[i].after("PREPOST=");
            let mut pre_plus = Vec::<String>::new();
            for p in ctl.gen_opt.pre.iter() {
                pre_plus.push(format!("{p}/{prepost}"));
            }
            ctl.gen_opt.pre.append(&mut pre_plus);
            break;
        }
    }

    // Process SOURCE.

    let mut args2 = args.clone();
    for i in 1..args.len() {
        if args[i].starts_with("SOURCE=") {
            let f = args[i].after("SOURCE=");
            let mut f2 = f.to_string();
            tilde_expand_me(&mut f2);
            let mut f2s = vec![f2.clone()];
            for pre in ctl.gen_opt.pre.iter() {
                f2s.push(format!("{}/{}", pre, f2));
            }
            let mut found = false;
            for f2 in f2s.iter() {
                if path_exists(f2) {
                    found = true;
                    require_readable_file(f2, "SOURCE")?;
                    let f = open_for_read![&f2];
                    for line in f.lines() {
                        let s = line.unwrap();
                        if !s.starts_with('#') {
                            let fields = s.split(' ').collect::<Vec<&str>>();
                            for j in 0..fields.len() {
                                if !fields[j].is_empty() {
                                    args2.insert(1, fields[j].to_string());
                                }
                            }
                        }
                    }
                }
            }
            if !found {
                return Err(format!(
                    "\nUnable to find SOURCE file {}.\n\
                    This was using PRE={}.\n",
                    f,
                    ctl.gen_opt.pre.iter().format(","),
                ));
            }
        }
    }
    args = args2;
    Ok(args)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn setup(
    ctl: &mut EncloneControl,
    args: &Vec<String>,
    argsx: &mut Vec<String>,
    args_orig: &Vec<String>,
) -> Result<(), String> {
    let mut using_pager = false;

    {
        for i in 2..args.len() {
            if args[i] == "help" {
                return Err("\nThe help argument, if used, must be the first argument \
                    to enclone.\n"
                    .to_string());
            }
        }
        let mut args = args.clone();
        let mut to_delete = vec![false; args.len()];
        let mut plain = false;
        let mut long_help = false;
        for i in 1..args.len() {
            if args[i] == "EVIL_EYE" {
                to_delete[i] = true;
            } else if args[i] == "HTML" {
                ctl.gen_opt.html = true;
                ctl.gen_opt.html_title = "enclone output".to_string();
                to_delete[i] = true;
            } else if args[i].starts_with("HTML=") {
                ctl.gen_opt.html = true;
                let mut title = args[i].after("HTML=").to_string();
                if title.starts_with('\"') && title.ends_with('\"') {
                    title = title.between("\"", "\"").to_string();
                }
                ctl.gen_opt.html_title = title;
                to_delete[i] = true;
            } else if args[i] == "SVG" {
                ctl.gen_opt.svg = true;
                to_delete[i] = true;
            } else if args[i] == "STABLE_DOC" {
                ctl.gen_opt.stable_doc = true;
                to_delete[i] = true;
            } else if args[i] == "LONG_HELP" {
                long_help = true;
                to_delete[i] = true;
            } else if args[i].starts_with("MAX_CORES=")
                || args[i].starts_with("PRE=")
                || args[i].starts_with("PREPOST=")
            {
                to_delete[i] = true;
            } else if args[i] == "PLAIN" {
                to_delete[i] = true;
                plain = true;
            } else if args[i] == "SPLIT" {
                ctl.gen_opt.split = true;
            }
        }

        // Proceed.

        if ctl.gen_opt.html && ctl.gen_opt.svg {
            return Err("\nBoth HTML and SVG cannot be used at the same time.\n".to_string());
        }
        erase_if(&mut args, &to_delete);
        *argsx = args.clone();
        if !ctl.cr_opt.nopager && (args.len() == 1 || args.contains(&"help".to_string())) {
            using_pager = true;
            setup_pager(true);
        }
        let mut help_all = false;
        if args.len() >= 3 && args[1] == "help" && args[2] == "all" {
            help_all = true;
        }
        let mut argsx = Vec::<String>::new();
        for i in 0..args_orig.len() {
            if args_orig[i] != "HTML"
                && args_orig[i] != "STABLE_DOC"
                && args_orig[i] != "NOPAGER"
                && args_orig[i] != "LONG_HELP"
                && !args_orig[i].starts_with("PRE=")
                && !args_orig[i].starts_with("PREPOST=")
                && !args_orig[i].starts_with("MAX_CORES=")
            {
                argsx.push(args_orig[i].clone());
            }
        }
        let mut h = HelpDesk::new(plain, help_all, long_help, ctl.gen_opt.html);
        help1(&argsx, &mut h)?;
        help2(&argsx, ctl, &mut h)?;
        help3(&argsx, &mut h)?;
        help4(&argsx, &mut h)?;
        help5(&argsx, ctl, &mut h)?;
        if argsx.len() == 1 || (argsx.len() > 1 && argsx[1] == "help") {
            return Ok(());
        }
    }

    // Pretest for some options.

    ctl.pretty = true;
    let mut nopretty = false;
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "PLAIN")? {
            ctl.pretty = false;
        }
        if is_simple_arg(&args[i], "NOPRETTY")? {
            nopretty = true;
        }
    }

    if !nopretty && !ctl.cr_opt.cellranger {
        set_panic_handler(args_orig);
        let mut nopager = false;
        for i in 1..args_orig.len() {
            if args_orig[i] == "NOPAGER" {
                nopager = true;
            }
        }
        if !nopager {
            using_pager = true;
            setup_pager(!nopager);
        }
    }

    // Process args (and set defaults for them).

    proc_args(ctl, args)?;
    if ctl.gen_opt.split {
        return Ok(());
    }
    Ok(())
}
