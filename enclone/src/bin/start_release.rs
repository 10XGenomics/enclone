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

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 0. Verify that we're in the top level directory.

    if !path_exists("Cargo.lock") {
        eprintln!("\nLooks like you're not in the top level directory.\n");
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 1. Verify that we're on master.

    println!("\nchecking git status");
    let new = Command::new("git")
        .arg("status")
        .output()
        .expect(&format!("failed to execute git status"));
    if new.status.code() != Some(0) {
        eprintln!("\ngit status failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }
    let s = strme(&new.stdout);
    if !s.contains("Your branch is up to date with 'origin/master'.") {
        eprintln!(
            "\nExpected to see message:\n\
            Your branch is up to date with 'origin/master'.\n"
        );
        std::process::exit(1);
    }
    if !s.contains("nothing to commit, working tree clean") {
        eprintln!(
            "\nExpected to see message:\n\
            nothing to commit, working tree clean\n"
        );
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 2. Pull.

    println!("running git pull");
    let new = Command::new("git")
        .arg("pull")
        .output()
        .expect(&format!("failed to execute git pull"));
    if new.status.code() != Some(0) {
        eprintln!("\ngit pull failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 3. Test "cargo b".

    println!("running cargo b");
    let new = Command::new("cargo")
        .arg("b")
        .output()
        .expect(&format!("failed to execute cargo b"));
    if new.status.code() != Some(0) {
        eprintln!("\ncargo b failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 4. Test "cargo t".

    println!("running cargo t");
    let new = Command::new("cargo")
        .arg("b")
        .output()
        .expect(&format!("failed to execute cargo t"));
    if new.status.code() != Some(0) {
        eprintln!("\ncargo t failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 5.  Bump version x.y.z to x.y.z+1 in every Cargo.toml.

    println!("bumping version");
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

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 6. Edit README.md to reflect the upcoming version.

    {
        let readme = include_str!["../../../README.md"];
        let new_readme = readme.replace(&old_version, &version);
        if new_readme == readme {
            eprintln!("\nFailed to update version in README.md.");
            eprintln!(
                "Could not change version from {} to {}.",
                old_version, version
            );
            eprintln!("Please do this now and continue manually with the release.");
            eprintln!("And figure out how this happened.\n");
            std::process::exit(1);
        }
        let mut g = open_for_write_new!["README.md"];
        fwrite!(g, "{}", new_readme);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 7. Run "cargo b" again to update Cargo.lock.

    println!("running cargo b");
    let new = Command::new("cargo")
        .arg("b")
        .output()
        .expect(&format!("failed to execute cargo b"));
    if new.status.code() != Some(0) {
        eprintln!("\ncargo b (2) failed \n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 8. Commit and push changes.

    println!("committing changes");
    let new = Command::new("git")
        .arg("commit")
        .arg("-a")
        .arg("-m")
        .arg("bump version")
        .output()
        .expect(&format!("failed to execute git commit"));
    if new.status.code() != Some(0) {
        eprintln!("\ngit commit failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }
    println!("pushing changes");
    let new = Command::new("git")
        .arg("push")
        .output()
        .expect(&format!("failed to execute git push"));
    if new.status.code() != Some(0) {
        eprintln!("\ngit push failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 9. Tag the commit.

    println!("tagging commit");
    let new = Command::new("git")
        .arg("tag")
        .arg(&format!("v{}", version))
        .output()
        .expect(&format!("failed to execute git tag"));
    if new.status.code() != Some(0) {
        eprintln!("\ngit tag failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 10. Trigger the release.

    println!("triggering release");
    let new = Command::new("git")
        .arg("push")
        .arg("origin")
        .arg("--tags")
        .output()
        .expect(&format!("failed to trigger release"));
    if new.status.code() != Some(0) {
        eprintln!("\nattempt to trigger release failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Step 11. Done.

    println!("\nAll done, looks like it worked!\n");
    println!("GitHub should now be making a release.\n");
    println!("Please read enclone/release_instructions.\n");
}
