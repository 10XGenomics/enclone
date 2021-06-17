// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use crate::USING_PAGER;
use chrono::{TimeZone, Utc};
use enclone::misc1::*;
use enclone_args::proc_args::*;
use enclone_args::proc_args2::*;
use enclone_core::defs::*;
use enclone_core::testlist::TEST_FILES_VERSION;
use enclone_core::*;
use enclone_help::help1::*;
use enclone_help::help2::*;
use enclone_help::help3::*;
use enclone_help::help4::*;
use enclone_help::help5::*;
use enclone_help::help_utils::*;
use io_utils::*;
use itertools::Itertools;
use lazy_static::lazy_static;
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::process::{Command, Stdio};
use std::sync::atomic::Ordering::SeqCst;
use std::sync::Mutex;
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

lazy_static! {
    pub static ref REMOTE_HOST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
}

pub fn setup(
    mut ctl: &mut EncloneControl,
    args: &Vec<String>,
    argsx: &mut Vec<String>,
    args_orig: &Vec<String>,
) -> Result<(), String> {
    let t = Instant::now();
    let mut using_pager = false;
    // Provide help if requested.

    let mut no_bug_reports = false;
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
            } else if args[i] == "NO_BUG_REPORTS" {
                no_bug_reports = true;
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
            if key == "ENCLONE_NO_BUG_REPORTS" {
                no_bug_reports = true;
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
            let now = Utc::now().naive_utc().timestamp();
            let build_date = version_string().after(":").between(": ", " :").to_string();
            let build_datetime = format!("{} 00:00:00", build_date);
            let then = Utc
                .datetime_from_str(&build_datetime, "%Y-%m-%d %H:%M:%S")
                .unwrap()
                .timestamp();
            let days_since_build = (now - then) / (60 * 60 * 24);
            let mut elapsed_message = String::new();
            if days_since_build > 30 {
                elapsed_message = format!(
                    "Your build is {} days old.  You might want to check \
                    to see if there is a newer build now.\n\n",
                    days_since_build
                );
            }

            if (!ctl.gen_opt.internal_run && REMOTE_HOST.lock().unwrap().len() == 0)
                || no_bug_reports
            {
                let exit_message = format!(
                    "Something has gone badly wrong.  You have probably encountered an internal \
                    error in enclone.\n\n\
                    Please email us at enclone@10xgenomics.com, including the traceback shown\n\
                    above and also the following version information:\n\
                    {} : {}.\n\n\
                    Your command was:\n\n{}\n\n\
                    {}\
                    🌸 Thank you so much for finding a bug and have a nice day! 🌸",
                    env!("CARGO_PKG_VERSION"),
                    version_string(),
                    args_orig.iter().format(" "),
                    elapsed_message,
                );
                PrettyTrace::new().exit_message(&exit_message).on();
            } else {
                // Set up to email bug report on panic.  This is only for internal users!

                let exit_message = format!(
                    "Something has gone badly wrong.  You have probably encountered an internal \
                    error in enclone.\n\n\
                    Here is the version information:\n\
                    {} : {}.\n\n\
                    Your command was:\n\n{}\n\n\
                    {}\
                    Thank you for being a happy internal enclone user.  All of this information \
                    is being\nemailed to enclone@10xgenomics.com, for the developers to \
                    contemplate.\n\n\
                    🌸 Thank you so much for finding a bug and have a nice day! 🌸",
                    env!("CARGO_PKG_VERSION"),
                    version_string(),
                    args_orig.iter().format(" "),
                    elapsed_message,
                );
                fn exit_function(msg: &str) {
                    let msg = format!("{}\n.\n", msg);
                    if !version_string().contains("macos") {
                        let process = Command::new("mail")
                            .arg("-s")
                            .arg("internal automated bug report")
                            .arg("enclone@10xgenomics.com")
                            .stdin(Stdio::piped())
                            .stdout(Stdio::piped())
                            .spawn();
                        let process = process.unwrap();
                        process.stdin.unwrap().write_all(msg.as_bytes()).unwrap();
                        let mut _s = String::new();
                        process.stdout.unwrap().read_to_string(&mut _s).unwrap();
                    } else if REMOTE_HOST.lock().unwrap().len() > 0 {
                        let remote_host = &REMOTE_HOST.lock().unwrap()[0];
                        let process = Command::new("ssh")
                            .arg(&remote_host)
                            .arg("mail")
                            .arg("-s")
                            .arg("\"internal ug report\"")
                            .arg("enclone@10xgenomics.com")
                            .stdin(Stdio::piped())
                            .stdout(Stdio::piped())
                            .spawn();
                        let process = process.unwrap();
                        process.stdin.unwrap().write_all(msg.as_bytes()).unwrap();
                        let mut _s = String::new();
                        process.stdout.unwrap().read_to_string(&mut _s).unwrap();
                    }
                }
                PrettyTrace::new()
                    .exit_message(&exit_message)
                    .run_this(exit_function)
                    .on();
            }
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
