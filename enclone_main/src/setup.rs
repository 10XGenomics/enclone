// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use crate::USING_PAGER;
use enclone::misc1::*;
use enclone_args::proc_args::*;
use enclone_args::proc_args2::*;
use enclone_core::defs::*;
use enclone_core::prepare_for_apocalypse::*;
use enclone_core::testlist::TEST_FILES_VERSION;
use enclone_core::*;
use enclone_help::help1::*;
use enclone_help::help2::*;
use enclone_help::help3::*;
use enclone_help::help4::*;
use enclone_help::help5::*;
use enclone_help::help_utils::*;
use io_utils::*;
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::atomic::Ordering::SeqCst;
use std::time::Instant;
use string_utils::*;
use tilde_expand::tilde_expand;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Process SOURCE args.

pub fn process_source(args: &Vec<String>) -> Result<Vec<String>, String> {
    let mut args2 = vec![args[0].clone()];
    for i in 1..args.len() {
        if args[i].starts_with("SOURCE=") {
            let f = args[i].after("SOURCE=");
            let f2 = stringme(&tilde_expand(&f.as_bytes()));
            if !path_exists(&f2) {
                return Err(format!("\nCan't find {}.\n", f));
            }
            let f = open_for_read![&f];
            for line in f.lines() {
                let s = line.unwrap();
                if !s.starts_with('#') {
                    let fields = s.split(' ').collect::<Vec<&str>>();
                    for j in 0..fields.len() {
                        if fields[j].len() > 0 {
                            args2.push(fields[j].to_string());
                        }
                    }
                }
            }
        } else {
            args2.push(args[i].clone());
        }
    }
    Ok(args2)
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
                return Err(format!(
                    "\nThe help argument, if used, must be the first argument \
                    to enclone.\n"
                ));
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
            } else if args[i] == "EVIL_EYE" {
                ctl.evil_eye = true;
            } else if args[i] == "HTML" {
                ctl.gen_opt.html = true;
                ctl.gen_opt.html_title = "enclone output".to_string();
                to_delete[i] = true;
            } else if args[i].starts_with("HTML=") {
                ctl.gen_opt.html = true;
                let mut title = args[i].after("HTML=").to_string();
                if title.starts_with("\"") && title.ends_with("\"") {
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
            } else if args[i].starts_with("CONFIG=") {
                ctl.gen_opt.config_file = args[i].after("CONFIG=").to_string();
            } else if args[i] == "NO_KILL" {
                to_delete[i] = true;
            }
        }
        for (key, value) in env::vars() {
            if key == "ENCLONE_CONFIG" {
                ctl.gen_opt.config_file = value.to_string();
            }
            if key == "ENCLONE_BUG_REPORTS" {
                bug_reports = value.to_string();
            }
        }

        // Test for internal run.

        if ctl.evil_eye {
            println!("testing for internal run");
        }
        for (key, value) in env::vars() {
            if key.contains("USER") && value.ends_with("10xgenomics.com") {
                if ctl.evil_eye {
                    println!("getting config");
                }
                if get_config(&ctl.gen_opt.config_file, &mut ctl.gen_opt.config) {
                    ctl.gen_opt.internal_run = true;
                }
                if ctl.evil_eye {
                    println!("got config");
                }
                break;
            }
        }
        for i in 1..args.len() {
            if args[i] == "FORCE_EXTERNAL".to_string() {
                ctl.gen_opt.internal_run = false;
            }
        }
        if ctl.gen_opt.internal_run {
            if ctl.evil_eye {
                println!("detected internal run");
            }
            if ctl.gen_opt.config.contains_key("bug_reports") {
                bug_reports = ctl.gen_opt.config["bug_reports"].clone();
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
            ctl.gen_opt.pre = vec![
                ctl.gen_opt.internal_data_dir.clone(),
                format!("enclone/test/inputs"),
                format!("enclone_exec"),
            ];
        } else if !ctl.gen_opt.cellranger {
            let home = dirs::home_dir().unwrap().to_str().unwrap().to_string();
            ctl.gen_opt.pre = vec![
                format!("{}/enclone/datasets", home),
                format!("{}/enclone/datasets2", home),
            ];
        }
        if ctl.gen_opt.config_file.contains(":") {
            let remote_host = ctl.gen_opt.config_file.before(":").to_string();
            REMOTE_HOST.lock().unwrap().clear();
            REMOTE_HOST.lock().unwrap().push(remote_host);
        }

        // Proceed.

        if ctl.gen_opt.html && ctl.gen_opt.svg {
            return Err(format!(
                "\nBoth HTML and SVG cannot be used at the same time.\n"
            ));
        }
        erase_if(&mut args, &to_delete);
        *argsx = args.clone();
        if args.len() == 1 || args.contains(&"help".to_string()) {
            PrettyTrace::new().on();
            if !nopager && !ctl.gen_opt.profile && !ctl.gen_opt.toy_com {
                using_pager = true;
                eprintln!("calling pager");
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
        let mut h = HelpDesk::new(plain, help_all, long_help, ctl.gen_opt.html);
        help1(&args, &mut h)?;
        help2(&args, &ctl, &mut h)?;
        help3(&args, &mut h)?;
        help4(&args, &mut h)?;
        help5(&args, &ctl, &mut h)?;
        if args.len() == 1 || (args.len() > 1 && args[1] == "help") {
            return Ok(());
        }
    }

    // Pretest for some options.

    ctl.pretty = true;
    let mut nopretty = false;
    ctl.gen_opt.h5 = true;
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "PLAIN")? {
            ctl.pretty = false;
        }
        if is_simple_arg(&args[i], "NOPRETTY")? {
            nopretty = true;
        }
        if is_simple_arg(&args[i], "COMP")? || args[i] == "COMPE" {
            ctl.comp = true;
        }
        if is_simple_arg(&args[i], "COMP2")? {
            ctl.comp = true;
            ctl.comp2 = true;
        }
        if is_simple_arg(&args[i], "CELLRANGER")? {
            ctl.gen_opt.cellranger = true;
        }
        if is_simple_arg(&args[i], "NH5")? {
            ctl.gen_opt.h5 = false;
        }
    }

    // Turn on pretty trace.

    if ctl.evil_eye {
        println!("about to turn on pretty trace");
    }
    if !nopretty && !ctl.gen_opt.cellranger {
        let mut ctrlc = false;
        for i in 1..args.len() {
            if is_simple_arg(&args[i], "CTRLC")? {
                ctrlc = true;
            }
        }
        let thread_message = new_thread_message();
        if ctrlc {
            PrettyTrace::new().message(&thread_message).ctrlc().on();
        } else {
            prepare_for_apocalypse(
                &args_orig,
                (ctl.gen_opt.internal_run || REMOTE_HOST.lock().unwrap().len() > 0)
                    && bug_reports.len() == 0,
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

    proc_args(&mut ctl, &args)?;
    if ctl.gen_opt.split {
        return Ok(());
    }
    Ok(())
}
