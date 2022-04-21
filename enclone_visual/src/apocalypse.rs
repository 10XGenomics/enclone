// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::*;
use chrono::{TimeZone, Utc};
use enclone_core::version_string;
use enclone_core::{BUG_REPORT_ADDRESS, REMOTE_HOST};
use itertools::Itertools;
use pretty_trace::*;
use std::env;
use std::io::{Read, Write};
use std::process::Stdio;
use std::sync::atomic::Ordering::SeqCst;

pub fn prepare_for_apocalypse_visual() {
    let email = INTERNAL.load(SeqCst);
    let ctrlc = CTRLC.load(SeqCst);
    let bug_reports = &BUG_REPORTS.lock().unwrap()[0];
    let args: Vec<String> = std::env::args().collect();
    if email {
        assert!(bug_reports.len() > 0);
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
        if !ctrlc {
            PrettyTrace::new().exit_message(&exit_message).on();
        } else {
            PrettyTrace::new().exit_message(&exit_message).ctrlc().on();
        }
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
            let mut msg = msg.to_string();

            // Get messages and shrink.

            let mut messages = compressed_message_history();
            let mut messages2 = Vec::<String>::new();
            for i in 0..messages.len() {
                if i == messages.len() - 1
                    || !messages[i].starts_with("ArchiveName(")
                    || !messages[i + 1].starts_with("ArchiveName(")
                {
                    messages2.push(messages[i].clone());
                }
            }
            messages = messages2;
            let mut messages2 = Vec::<String>::new();
            let mut i = 0;
            while i < messages.len() {
                if i < messages.len() - 1
                    && messages[i] == "SubmitButtonPressed(Ok(()))"
                    && messages[i + 1] == "ComputationDone(Ok(()))"
                {
                    messages2.push("Submit button pressed".to_string());
                    i += 1;
                } else {
                    messages2.push(messages[i].clone());
                }
                i += 1;
            }
            messages = messages2;
            let mut messages2 = Vec::<String>::new();
            let mut i = 0;
            while i < messages.len() {
                if i < messages.len() - 1
                    && messages[i] == "Save"
                    && messages[i + 1] == "CompleteSave(Ok(()))"
                {
                    messages2.push("Save button pressed".to_string());
                    i += 1;
                } else {
                    messages2.push(messages[i].clone());
                }
                i += 1;
            }
            messages = messages2;
            let mut messages2 = Vec::<String>::new();
            let mut i = 0;
            while i < messages.len() {
                if i < messages.len() - 1
                    && messages[i] == "ArchiveRefresh"
                    && messages[i + 1] == "ArchiveRefreshComplete(Ok(()))"
                {
                    messages2.push("Refresh button pressed".to_string());
                    i += 1;
                } else {
                    messages2.push(messages[i].clone());
                }
                i += 1;
            }
            messages = messages2;
            let mut messages2 = Vec::<String>::new();
            let mut i = 0;
            while i < messages.len() {
                if i < messages.len() - 1
                    && messages[i] == "DoShare(true)"
                    && messages[i + 1] == "CompleteDoShare(Ok(()))"
                {
                    messages2.push("sharing".to_string());
                    i += 1;
                } else {
                    messages2.push(messages[i].clone());
                }
                i += 1;
            }
            messages = messages2;
            for i in 0..messages.len() {
                if messages[i] == "ArchiveOpen(Ok(()))" {
                    messages[i] = "Archive button pressed".to_string();
                } else if messages[i] == "ArchiveClose" {
                    messages[i] = "Dismiss button pressed".to_string();
                } else if messages[i] == "DelButtonPressed(Ok(()))" {
                    messages[i] = "Del button pressed".to_string();
                }
            }

            // Proceed.

            println!("enclone visual message history:\n");
            msg += "enclone visual message history:\n\n";
            for i in 0..messages.len() {
                println!("[{}] {}", i + 1, messages[i]);
                msg += &mut format!("[{}] {}\n", i + 1, messages[i]);
            }
            println!("");
            let msg = format!("{}\n.\n", msg);
            let bug_report_address = &BUG_REPORT_ADDRESS.lock().unwrap()[0];
            if !version_string().contains("macos") {
                let process = std::process::Command::new("mail")
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
                let process = std::process::Command::new("ssh")
                    .arg(&remote_host)
                    .arg("mail")
                    .arg("-s")
                    .arg("\"internal automated bug report\"")
                    .arg(&bug_report_address)
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn();
                let mut process = process.unwrap();
                process
                    .stdin
                    .as_ref()
                    .unwrap()
                    .write_all(msg.as_bytes())
                    .unwrap();
                if !process.wait().unwrap().success() {
                    eprintln!("\nAttempt to send email failed.\n");
                }
            }
        }
        if !ctrlc {
            PrettyTrace::new()
                .exit_message(&exit_message)
                .run_this(exit_function)
                .on();
        } else {
            PrettyTrace::new()
                .exit_message(&exit_message)
                .run_this(exit_function)
                .ctrlc()
                .on();
        }
    }
}
