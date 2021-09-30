// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Set up response on panic.  If the email argument is sent, then an email is sent to the
// argument bug_reports.  The caller only sets these for internal users.

use crate::version_string;
use crate::{BUG_REPORT_ADDRESS, REMOTE_HOST};
use chrono::{TimeZone, Utc};
use itertools::Itertools;
use pretty_trace::*;
use std::env;
use std::io::{Read, Write};
use std::process::{Command, Stdio};
use string_utils::*;

pub fn prepare_for_apocalypse(args: &Vec<String>, email: bool, bug_reports: &str) {
    if email {
        assert!(!bug_reports.is_empty());
    }
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
    if !email {
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
            args.iter().format(" "),
            elapsed_message,
        );
        PrettyTrace::new().exit_message(&exit_message).on();
    } else {
        // Set up to email bug report on panic.  This is only for internal users!

        let mut contemplate = "".to_string();
        if bug_reports == "enclone@10xgenomics.com" {
            contemplate = ", for the developers to contemplate".to_string();
        }
        let exit_message = format!(
            "Something has gone badly wrong.  You have probably encountered an internal \
            error in enclone.\n\n\
            Here is the version information:\n\
            {} : {}.\n\n\
            Your command was:\n\n{}\n\n\
            {}\
            Thank you for being a happy internal enclone user.  All of this information \
            is being\nemailed to {}{}.\n\n\
            🌸 Thank you so much for finding a bug and have a nice day! 🌸",
            env!("CARGO_PKG_VERSION"),
            version_string(),
            args.iter().format(" "),
            elapsed_message,
            bug_reports,
            contemplate,
        );
        BUG_REPORT_ADDRESS
            .lock()
            .unwrap()
            .push(bug_reports.to_string());
        fn exit_function(msg: &str) {
            let msg = format!("{}\n.\n", msg);
            let bug_report_address = &BUG_REPORT_ADDRESS.lock().unwrap()[0];
            if !version_string().contains("macos") {
                let process = Command::new("mail")
                    .arg("-s")
                    .arg("internal automated bug report")
                    .arg(&bug_report_address)
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
                    .arg("\"internal automated bug report\"")
                    .arg(&bug_report_address)
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
}
