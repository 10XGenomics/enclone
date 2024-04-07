// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Note: if enclone is run from the command line, and fails, it will still return exit status
// zero.  As far as we know, in all other cases where it is not run from the command line, it
// returns exit status zero.

use enclone_main::main_enclone::main_enclone;

use io_utils::*;

use std::env;

use std::process::Command;
use string_utils::*;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    let mut update = false;
    for i in 1..args.len() {
        if args[i] == "UPDATE" {
            update = true;
        }
    }

    // Update mode.

    if update {
        if args.len() != 2 {
            eprintln!(
                "\nYou've specified UPDATE, but in update mode, we expect only a single \
                argument.\nIf you'd like to update, please type just enclone UPDATE.\n"
            );
            std::process::exit(1);
        }
        let mut home = String::new();
        for (key, value) in env::vars() {
            if key == "HOME" {
                home = value.clone();
            }
        }
        if home.is_empty() {
            eprintln!("Weird, unable to determine your home directory.\n");
            std::process::exit(1);
        }
        let datasets = format!("{}/enclone/datasets", home);
        if !path_exists(&datasets) || !std::path::Path::new(&datasets).is_dir() {
            eprintln!(
                "\nSomething odd has happened.  There should be a directory ~/enclone and \
                inside that, a directory datasets.\n"
            );
            std::process::exit(1);
        }
        let list = dir_list(&datasets);
        let size;
        if list.len() <= 3 {
            size = "small";
        } else if list.len() <= 40 {
            size = "medium";
        } else {
            size = "large";
        }
        println!("updating enclone using size = {}", size);
        let o = Command::new("bash")
            .arg("-c")
            .arg(&format!(
                "curl -sSf -L bit.ly/enclone_install | bash -s {}",
                size
            ))
            .output()
            .expect("failed to execute curl");
        print!("{}{}", strme(&o.stdout), strme(&o.stderr));
        std::process::exit(0);
    }

    // Standard run of enclone.

    if args.len() < 2 || args[1] != "SERVER" {
        // Test for error.

        if let Err(err) = main_enclone(args) {
            eprintln!("{err}");
            std::process::exit(1);
        }

        // Done.

        std::process::exit(0);
    }

    Ok(())
}
