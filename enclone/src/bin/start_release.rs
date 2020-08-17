// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Start a release, assuming that you want to bump the z in x.y.z.
//
// See enclone/release_instructions.

use io_utils::*;
use itertools::Itertools;
use pretty_trace::*;
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();

    // Step 1.  Bump version x.y.z to x.y.z+1 in every Cargo.toml.

    let all = read_dir(".").unwrap();
    let mut versions = Vec::<String>::new();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        let toml = format!("{}/Cargo.toml", f);
        if path_exists(&toml) {
            let g = open_for_read![&toml];
            let mut found_version = false;
            for line in g.lines() {
                let s = line.unwrap();
                if s.starts_with("version = \"") {
                    versions.push(s.between("\"", "\"").to_string());
                    found_version = true;
                }
            }
            if !found_version {
                eprintln!("\nFailed to find version in {}.\n", toml);
                std::process::exit(1);
            }
        }
    }
    unique_sort(&mut versions);
    if versions.len() > 1 {
        eprintln!(
            "\nFound multiple versions: {}",
            versions.iter().format(", ")
        );
        std::process::exit(1);
    }
    let old_version = versions[0].clone();
    let version = format!(
        "{}.{}",
        versions[0].rev_before("."),
        versions[0].rev_after(".").force_usize() + 1
    );
    let all = read_dir(".").unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        let toml = format!("{}/Cargo.toml", f);
        if path_exists(&toml) {
            let mut newg = Vec::<String>::new();
            {
                let g = open_for_read![&toml];
                for line in g.lines() {
                    let s = line.unwrap();
                    if s.starts_with("version = ") {
                        newg.push(format!("version = \"{}\"", version));
                    } else {
                        newg.push(s.clone());
                    }
                }
            }
            let mut g = open_for_write_new![&toml];
            fwrite!(g, "{}\n", newg.iter().format("\n"));
        }
    }

    // 2. Edit README.md to reflect the upcoming version.

    {
        let readme = include_str!["../../../README.md"].replace(&old_version, &version);
        let mut g = open_for_write_new!["README.md"];
        fwriteln!(g, "{}", readme);
    }

    // 3. Commit and push changes.

    let new = Command::new("git")
        .arg("commit")
        .arg("-a")
        .output()
        .expect(&format!("failed to execute git commit"));
    if new.status.code() != Some(0) {
        eprintln!("\ngit commit failed\n");
        std::process::exit(1);
    }
    let new = Command::new("git")
        .arg("push")
        .output()
        .expect(&format!("failed to execute git push"));
    if new.status.code() != Some(0) {
        eprintln!("\ngit push failed\n");
        std::process::exit(1);
    }

    // 4. Tag the commit.

    let new = Command::new("git")
        .arg("tag")
        .arg(&version)
        .output()
        .expect(&format!("failed to execute git tag"));
    if new.status.code() != Some(0) {
        eprintln!("\ngit tag failed\n");
        std::process::exit(1);
    }

    // 5. Trigger the release.

    let new = Command::new("git")
        .arg("push")
        .arg("origin")
        .arg("--tags")
        .output()
        .expect(&format!("failed to trigger release"));
    if new.status.code() != Some(0) {
        eprintln!("\nattempt to trigger release failed\n");
        std::process::exit(1);
    }

    // 6. Done.

    println!("\nAll done, looks like it worked!\n");
    println!("Please read enclone/release_instructions.\n");
}
