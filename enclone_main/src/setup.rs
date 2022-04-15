// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use crate::USING_PAGER;
use enclone::misc1::setup_pager;
use enclone_args::proc_args::proc_args;
use enclone_args::proc_args2::is_simple_arg;
use enclone_core::defs::{get_config, EncloneControl};
use enclone_core::prepare_for_apocalypse::prepare_for_apocalypse;
use enclone_core::testlist::TEST_FILES_VERSION;
use enclone_core::{require_readable_file, tilde_expand_me, REMOTE_HOST};
use enclone_help::help1::help1;
use enclone_help::help2::help2;
use enclone_help::help3::help3;
use enclone_help::help4::help4;
use enclone_help::help5::help5;
use enclone_help::help_utils::{HelpDesk, HELP_ALL, PLAIN};
use io_utils::{open_for_read, path_exists};
use itertools::Itertools;
use pretty_trace::{new_thread_message, PrettyTrace};
use std::env;
use std::io::BufRead;
use std::sync::atomic::Ordering::SeqCst;
use std::time::Instant;
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

    // Check for CONFIG and EVIL_EYE and CELLRANGER.

    for i in 1..args.len() {
        if args[i].starts_with("CONFIG=") {
            ctl.gen_opt.config_file = args[i].after("CONFIG=").to_string();
        } else if args[i] == "EVIL_EYE" {
            ctl.gen_opt.evil_eye = true;
            if ctl.gen_opt.evil_eye {
                println!("the evil eye is on");
            }
        } else if is_simple_arg(&args[i], "CELLRANGER")? {
            ctl.gen_opt.cellranger = true;
        }
    }

    // Test for internal run and get configuration.

    if ctl.gen_opt.evil_eye {
        println!("testing for internal run");
    }
    for (key, value) in env::vars() {
        if ctl.gen_opt.evil_eye {
            print!("see env {} = {}", key, value);
            if key.starts_with("ENCLONE_") {
                print!(" ***** SETTING ENCLONE VAR ***** ");
            }
            println!("");
        }
        if key.ends_with("USER") && value.ends_with("10xgenomics.com") {
            if ctl.gen_opt.evil_eye {
                println!("getting config");
            }
            if get_config(&ctl.gen_opt.config_file, &mut ctl.gen_opt.config) {
                ctl.gen_opt.internal_run = true;
            }
            if ctl.gen_opt.evil_eye {
                println!("got config");
            }
            break;
        }
    }
    for i in 1..args.len() {
        if args[i] == *"FORCE_EXTERNAL" {
            ctl.gen_opt.internal_run = false;
        }
    }

    // Get configuration info.

    if ctl.gen_opt.internal_run {
        if ctl.gen_opt.evil_eye {
            println!("detected internal run");
        }
        let earth_path = format!(
            "{}/current{}",
            ctl.gen_opt.config["earth"], TEST_FILES_VERSION
        );
        ctl.gen_opt.internal_data_dir = earth_path;
        let cloud_path = format!(
            "{}/current{}",
            ctl.gen_opt.config["cloud"], TEST_FILES_VERSION
        );
        if path_exists(&cloud_path) {
            ctl.gen_opt.internal_data_dir = cloud_path;
        }
    }
    if ctl.gen_opt.config_file.contains(':') {
        let remote_host = ctl.gen_opt.config_file.before(":").to_string();
        REMOTE_HOST.lock().unwrap().clear();
        REMOTE_HOST.lock().unwrap().push(remote_host);
    }

    // Determine PRE.

    if ctl.gen_opt.internal_run {
        ctl.gen_opt.pre = vec![
            ctl.gen_opt.internal_data_dir.clone(),
            "enclone/test/inputs".to_string(),
            "enclone_exec".to_string(),
        ];
    } else if !ctl.gen_opt.cellranger {
        let home = dirs::home_dir().unwrap().to_str().unwrap().to_string();
        ctl.gen_opt.pre = vec![
            format!("{}/enclone/datasets_me", home),
            format!("{}/enclone/datasets", home),
            format!("{}/enclone/datasets2", home),
        ];
    }
    if ctl.gen_opt.config.contains_key("prep") {
        let prep = ctl.gen_opt.config["prep"].split(',').collect::<Vec<&str>>();
        for pre in prep.clone() {
            ctl.gen_opt.pre.push(pre.to_string());
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
    mut ctl: &mut EncloneControl,
    args: &Vec<String>,
    argsx: &mut Vec<String>,
    args_orig: &Vec<String>,
) -> Result<(), String> {
    let t = Instant::now();
    let mut using_pager = false;
    // Provide help if requested.

    let mut bug_reports = "enclone@10xgenomics.com".to_string();
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
        let mut nopager = false;
        let mut plain = false;
        let mut long_help = false;
        if ctl.gen_opt.profile {
            nopager = true;
            ctl.gen_opt.profile = true;
        }
        for i in 1..args.len() {
            if args[i] == "NOPAGER" || args[i] == "EVIL_EYE" || args[i] == "TOY_COM" {
                nopager = true;
                ctl.gen_opt.nopager = true;
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
            } else if args[i] == "FORCE_EXTERNAL" {
                to_delete[i] = true;
            } else if args[i] == "LONG_HELP" {
                long_help = true;
                to_delete[i] = true;
            } else if args[i].starts_with("MAX_CORES=") {
                to_delete[i] = true;
            } else if args[i].starts_with("PRE=") {
                to_delete[i] = true;
            } else if args[i].starts_with("PREPOST=") {
                to_delete[i] = true;
            } else if args[i] == "PLAIN" {
                to_delete[i] = true;
                plain = true;
                unsafe {
                    PLAIN = true;
                }
            } else if args[i] == "SPLIT" {
                ctl.gen_opt.split = true;
            } else if args[i] == "BUG_REPORTS" {
                bug_reports = "".to_string();
            } else if args[i].starts_with("BUG_REPORTS=") {
                bug_reports = args[i].after("BUG_REPORTS=").to_string();
            } else if args[i] == "NO_KILL" {
                to_delete[i] = true;
            }
        }
        for (key, value) in env::vars() {
            if key == "ENCLONE_BUG_REPORTS" {
                bug_reports = value.to_string();
            }
        }
        if ctl.gen_opt.internal_run && ctl.gen_opt.config.contains_key("bug_reports") {
            bug_reports = ctl.gen_opt.config["bug_reports"].clone();
        }

        // Proceed.

        if ctl.gen_opt.html && ctl.gen_opt.svg {
            return Err("\nBoth HTML and SVG cannot be used at the same time.\n".to_string());
        }
        erase_if(&mut args, &to_delete);
        *argsx = args.clone();
        if args.len() == 1 || args.contains(&"help".to_string()) {
            PrettyTrace::new().on();
            if !nopager && !ctl.gen_opt.profile && !ctl.gen_opt.toy_com {
                using_pager = true;
                setup_pager(true);
            }
        }
        let mut help_all = false;
        if args.len() >= 3 && args[1] == "help" && args[2] == "all" {
            unsafe {
                HELP_ALL = true;
            }
            help_all = true;
        }
        let mut argsx = Vec::<String>::new();
        for i in 0..args_orig.len() {
            if args_orig[i] != "HTML"
                && args_orig[i] != "STABLE_DOC"
                && args_orig[i] != "NOPAGER"
                && args_orig[i] != "FORCE_EXTERNAL"
                && args_orig[i] != "NO_KILL"
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
    ctl.gen_opt.h5 = true;
    let mut visual = false;
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "PLAIN")? {
            ctl.pretty = false;
        }
        if is_simple_arg(&args[i], "NOPRETTY")? {
            nopretty = true;
        }
        if is_simple_arg(&args[i], "COMP")? || args[i] == "COMPE" {
            ctl.perf_opt.comp = true;
        }
        if is_simple_arg(&args[i], "COMP2")? {
            ctl.perf_opt.comp = true;
            ctl.perf_opt.comp2 = true;
        }
        if is_simple_arg(&args[i], "NH5")? {
            ctl.gen_opt.h5 = false;
        }
        if is_simple_arg(&args[i], "VISUAL")? {
            visual = true;
        }
    }

    // Turn on pretty trace.

    if ctl.gen_opt.evil_eye {
        println!("about to turn on pretty trace");
    }
    if !nopretty && !ctl.gen_opt.cellranger {
        let mut ctrlc = false;
        for i in 1..args.len() {
            if is_simple_arg(&args[i], "CTRLC")? {
                if visual {
                    return Err(format!("Sorry CTRLC can't be used in visual mode."));
                }
                ctrlc = true;
            }
        }
        let thread_message = new_thread_message();
        if ctrlc {
            PrettyTrace::new().message(thread_message).ctrlc().on();
        } else {
            prepare_for_apocalypse(
                args_orig,
                (ctl.gen_opt.internal_run || REMOTE_HOST.lock().unwrap().len() > 0)
                    && !bug_reports.is_empty(),
                &bug_reports,
            );
            let mut nopager = false;
            for i in 1..args_orig.len() {
                if args_orig[i] == "NOPAGER" || args_orig[i] == "TOY_COM" {
                    nopager = true;
                }
            }
            if !nopager && !ctl.gen_opt.profile {
                using_pager = true;
                setup_pager(!nopager && !ctl.gen_opt.profile);
            }
        }
    }
    USING_PAGER.store(using_pager, SeqCst);
    ctl.perf_stats(&t, "in first part of setup");

    // Process args (and set defaults for them).

    proc_args(&mut ctl, args)?;
    if ctl.gen_opt.split {
        return Ok(());
    }
    Ok(())
}
