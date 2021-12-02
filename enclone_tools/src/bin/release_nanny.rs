// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Runs continuously, waiting for a new release to appear, and when that happens, updates the
// remote linux binary, and exits.  This is launched by start_release.

use enclone_core::defs::get_config;
use io_utils::*;
use pretty_trace::PrettyTrace;
use std::collections::HashMap;
use std::env;
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::process::{Command, Stdio};
use std::thread;
use std::time::Duration;
use string_utils::{strme, TextUtils};

fn mail(address: &str, title: &str) {
    let process = std::process::Command::new("mail")
        .arg("-s")
        .arg(&title)
        .arg(&address)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn();
    let process = process.unwrap();
    let msg = String::new();
    process.stdin.unwrap().write_all(msg.as_bytes()).unwrap();
}

fn main() {
    PrettyTrace::new().on();

    // Don't run on a Mac.

    if cfg!(any(target_os = "macos", target_os = "ios")) {
        eprintln!("\nrelease_nanny should not be run from a Mac\n");
        std::process::exit(1);
    }

    // Determine email address.

    let mut user = None;
    for (key, value) in env::vars() {
        if key == "USER" {
            user = Some(value.clone());
        }
    }
    if user.is_none() {
        eprintln!("\nrelease nanny unable to determine user name\n");
        std::process::exit(1);
    }
    let address = format!("{}@10xgenomics.com", user.unwrap());

    // Loop until enclone executable updates.

    loop {
        // Sleep for a minute.

        thread::sleep(Duration::from_millis(60 * 1000));

        // Get version number of latest enclone.

        let mut github_version = None;
        let o = Command::new("curl")
            .arg("-sI")
            .arg(
                "https://github.com/10XGenomics/enclone/\
                releases/latest/download/enclone_linux",
            )
            .output()
            .expect("failed to execute curl");
        let out = strme(&o.stdout);
        for line in out.lines() {
            if line.starts_with("Location:") {
                github_version = Some(line.between("/download/", "/enclone_linux").to_string());
            }
        }
        assert!(github_version.is_some());
        let github_version = github_version.unwrap();

        // Get version number of enclone on remote server.

        let mut config = HashMap::<String, String>::new();
        let mut config_file = String::new();
        for (key, value) in env::vars() {
            if key == "ENCLONE_CONFIG" {
                config_file = value.to_string();
            }
        }
        let mut remote_version_file = String::new();
        let mut remote_version = None;
        if get_config(&config_file, &mut config) {
            let bin = &config["enclone_linux_bin"];
            remote_version_file = format!("{}/version", bin);
            let version = std::fs::read_to_string(&remote_version_file).unwrap();
            remote_version = Some(version.before("\n").to_string());
        }
        assert!(remote_version.is_some());
        let remote_version = remote_version.unwrap();

        // Test for change.

        if github_version != remote_version {
            let msg = format!(
                "release nanny start update from {} to {}",
                github_version, remote_version,
            );
            mail(&address, &msg);
            let bin = &config["enclone_linux_bin"];
            let mut f = open_for_write_new![&remote_version_file];
            fwrite!(f, "{}\n", github_version);
            let current = format!("{}/enclone", bin);
            let last = format!("{}/enclone_last", bin);
            if path_exists(&last) {
                std::fs::remove_file(&last).unwrap();
            }
            if path_exists(&current) {
                std::fs::rename(&current, &last).unwrap();
            }
            let o = Command::new("curl")
                .arg("-s")
                .arg("-L")
                .arg(
                    "https://github.com/10XGenomics/enclone/\
                    releases/latest/download/enclone_linux",
                )
                .arg("--output")
                .arg(&current)
                .output()
                .expect("failed to execute curl");
            if o.status.code() != Some(0) {
                eprintln!(
                    "Update failed with the following error message:\n{}",
                    strme(&o.stderr)
                );
                std::process::exit(1);
            }
            let perms = std::fs::Permissions::from_mode(0o775);
            std::fs::set_permissions(&current, perms).unwrap();
            let msg = format!(
                "release nanny end update from {} to {}",
                github_version, remote_version,
            );
            mail(&address, &msg);

            break;
        }
    }
}
