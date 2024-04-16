// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use enclone_args::proc_args::proc_args;
use enclone_args::proc_args2::is_simple_arg;
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

#[cfg(not(target_os = "windows"))]
use pager::Pager;

use std::env;
use std::io::BufRead;

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
        ctl.cr_opt.pre = vec![
            format!("{}/enclone/datasets_me", home),
            format!("{}/enclone/datasets", home),
            format!("{}/enclone/datasets2", home),
        ];
    }
    for i in (1..args.len()).rev() {
        if args[i].starts_with("PREPOST=") {
            let prepost = args[i].after("PREPOST=");
            let mut pre_plus = Vec::<String>::new();
            for p in ctl.cr_opt.pre.iter() {
                pre_plus.push(format!("{p}/{prepost}"));
            }
            ctl.cr_opt.pre.append(&mut pre_plus);
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
            for pre in ctl.cr_opt.pre.iter() {
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
                    ctl.cr_opt.pre.iter().format(","),
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
    // Provide help if requested.
    let mut nopager = false;
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
            if args[i] == "NOPAGER" || args[i] == "EVIL_EYE" {
                nopager = true;
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
        if !nopager && (args.len() == 1 || args.contains(&"help".to_string())) {
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
    for i in 1..args.len() {
        if is_simple_arg(&args[i], "PLAIN")? {
            ctl.pretty = false;
        }
    }

    if !nopager {
        setup_pager(true);
    }

    // Process args (and set defaults for them).

    proc_args(ctl, args)?;
    if ctl.gen_opt.split {
        return Ok(());
    }
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This section contains a function that supports paging.  It does not work under Windows, and
// we describe here all the *known* problems with getting enclone to work under Windows.
// 1. It does not compile for us.  When we tried, there was a problem with libhdf-5.
// 2. Paging is turned off, because the pager crate doesn't compile under Windows, and porting
//    it to Windows appears nontrivial.
// 3. ANSI escape characters are not handled correctly, at least by default.
// In addition, we have some concerns about what it would mean to properly test enclone on Windows,
// given that some users might have older OS installs, and support for ANSI escape characters
// appears to have been changed in 2018. This is not made easier by the Windows Subsystem for
// Linux.

#[cfg(not(target_os = "windows"))]
pub fn setup_pager(pager: bool) {
    // If the output is going to a terminal, set up paging so that output is in effect piped to
    // "less -R -F -X -K".
    //
    // ∙ The option -R is used to render ANSI escape characters correctly.  We do not use
    //   -r instead because if you navigate backwards in less -r, stuff gets screwed up,
    //   which is consistent with the scary stuff in the man page for less at -r.  However -R will
    //   not display all unicode characters correctly, so those have to be picked carefully,
    //   by empirically testing that e.g. "echo ◼ | less -R -F -X" renders correctly.
    //
    // ∙ The -F option makes less exit immediately if all the output can be seen in one screen.
    //
    // ∙ The -X option is needed because we found that in full screen mode on OSX Catalina, output
    //   was sent to the alternate screen, and hence it appeared that one got no output at all
    //   from enclone.  This is really bad, so do not turn off this option!

    if pager {
        Pager::with_pager("less -R -F -X -K").setup();
    }
}

#[cfg(target_os = "windows")]
pub fn setup_pager(_pager: bool) {}
