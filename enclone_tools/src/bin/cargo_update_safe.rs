// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// For each dependent crate, cargo update it and run cargo b.  If that succeeds,
// commit the change.  Otherwise, undo.
//
// Limitations:
// * Only looks at crates in master.toml having simple specification crate = "version".
// * Assumes that master.toml is complete.
// * Updates crates even if they have changed very recently.
// * Should check that results are unchanged.

use io_utils::*;
use perf_stats::*;
use pretty_trace::*;
use std::collections::HashMap;
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use std::time::Instant;
use string_utils::*;

fn print_dot(dots: &mut usize) {
    if *dots > 0 && *dots % 10 == 0 {
        print!(" ");
    }
    print!(".");
    std::io::stdout().flush().unwrap();
    *dots += 1;
    if *dots == 90 {
        println!();
        std::io::stdout().flush().unwrap();
        *dots = 0;
    }
}

fn reset() {
    let o = Command::new("git")
        .arg("reset")
        .arg("--hard")
        .output()
        .unwrap_or_else(|_| panic!("{}", "\n\nfailed to execute git reset".to_string()));
    if o.status.code() != Some(0) {
        eprintln!("\ngit reset failed\n");
        std::process::exit(1);
    }
}

fn main() {
    let t = Instant::now();
    PrettyTrace::new().on();
    println!();

    // Test for synced to master.  This code is essentially identical to code in an enclone test.

    let mut version = HashMap::<String, String>::new();
    let f = open_for_read!["master.toml"];
    for line in f.lines() {
        let s = line.unwrap();
        if !s.starts_with('#') && s.contains('=') {
            version.insert(s.before(" = ").to_string(), s.after(" = ").to_string());
        }
    }
    let all = read_dir(".").unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        let toml = format!("{}/Cargo.toml", f);
        if path_exists(&toml) {
            let g = open_for_read![&toml];
            for line in g.lines() {
                let s = line.unwrap();
                if s.contains(" =") {
                    let cratex = s.before(" =").to_string();
                    if version.contains_key(&cratex) {
                        let t = format!("{} = {}", cratex, version[&cratex]);
                        if t != s {
                            eprintln!("\nFound change in {}.\nold: {}\nnew: {}", toml, s, t);
                            eprintln!("You probably need to run sync_to_master\n");
                            std::process::exit(1);
                        }
                    }
                }
            }
        }
    }

    // Loop until nothing changes.

    loop {
        let mut changed = false;

        // Get crate versions from master.toml.

        let mut master = Vec::<(String, String)>::new();
        let f = open_for_read!["master.toml"];
        for line in f.lines() {
            let s = line.unwrap();
            if !s.starts_with('#') && s.contains('=') {
                let cratex = s.before(" = ").to_string();
                let version = s.after(" = ").to_string();
                if version.starts_with('"') && version.ends_with('"') {
                    let version = version.after("\"").rev_before("\"").to_string();
                    let fields = version.split('.').collect::<Vec<&str>>();
                    let mut ok = true;
                    for i in 0..fields.len() {
                        if fields[i].parse::<usize>().is_err() {
                            ok = false;
                        }
                    }
                    if ok {
                        master.push((cratex, version));
                    }
                }
            }
        }

        // Look for possible raises to crate versions in master.toml.

        println!("examining {} directly dependent crates", master.len());
        for x in master.iter() {
            let (cratex, old) = (&x.0, &x.1);
            let o = Command::new("cargo")
                .arg("search")
                .arg("--limit")
                .arg("1")
                .arg(&cratex)
                .output()
                .unwrap_or_else(|_| panic!("{}", "\n\nfailed to execute cargo search".to_string()));
            if o.status.code() != Some(0) {
                eprintln!("\ncargo search failed\n");
                std::process::exit(1);
            }
            let out = strme(&o.stdout);
            let mut new = String::new();
            for line in out.lines() {
                let c = line.before(" =");
                assert_eq!(c, cratex);
                new = line.between("\"", "\"").to_string();
                break;
            }
            let old_dots = old.matches('.').count();
            let mut new_dots = new.matches('.').count();
            if new_dots == old_dots + 1 {
                new = new.rev_before(".").to_string();
                new_dots -= 1;
            } else if new_dots == old_dots + 2 {
                new = new.rev_before(".").rev_before(".").to_string();
                new_dots -= 2;
            }
            if old_dots == new_dots {
                if old_dots >= 1 && old.before(".") != new.before(".") {
                    new = new.before(".").to_string();
                } else if old_dots == 2 && old.between(".", ".") != new.between(".", ".") {
                    new = new.rev_before(".").to_string();
                }
            }
            if *old != new {
                println!("trying to update {} from {} to {}", cratex, old, new);
                let mut new_lines = Vec::<String>::new();
                let f = open_for_read!["master.toml"];
                for line in f.lines() {
                    let s = line.unwrap();
                    let mut saved = false;
                    if !s.starts_with('#') && s.contains('=') {
                        let cratey = s.before(" = ").to_string();
                        if cratey == *cratex {
                            new_lines.push(format!("{} = \"{}\"", cratex, new));
                            saved = true;
                        }
                    }
                    if !saved {
                        new_lines.push(s);
                    }
                }
                {
                    let mut f = open_for_write_new!["master.toml"];
                    for line in new_lines.iter() {
                        fwriteln!(f, "{}", line);
                    }
                }
                let o = Command::new("sync_to_master").output().unwrap_or_else(|_| {
                    panic!("{}", "\n\nfailed to execute sync_to_master".to_string())
                });
                if o.status.code() != Some(0) {
                    eprintln!("\nsync_to_master failed");
                    std::process::exit(1);
                }
                let o = Command::new("cargo").arg("b").output().unwrap_or_else(|_| {
                    panic!("{}", "\n\nfailed to execute cargo b 0".to_string())
                });
                if o.status.code() != Some(0) {
                    reset();
                    continue;
                }
                println!("succeeded, committing change to {}", cratex);
                changed = true;
                let new = Command::new("git")
                    .arg("commit")
                    .arg("-a")
                    .arg("-m")
                    .arg(&format!("raise version of crate {}", cratex))
                    .output()
                    .unwrap_or_else(|_| {
                        panic!("{}", "\n\nfailed to execute git commit".to_string())
                    });
                if new.status.code() != Some(0) {
                    println!("\n\ngit commit failed, something is wrong");
                    std::process::exit(1);
                }
            }
        }

        // Get complete list of crate:version instances from Cargo.lock.

        let f = open_for_read!["Cargo.lock"];
        let mut crates = Vec::<String>::new();
        let mut versions = Vec::<String>::new();
        let mut lines = Vec::<String>::new();
        for line in f.lines() {
            let s = line.unwrap();
            lines.push(s);
        }
        let mut i = 0;
        while i < lines.len() - 1 {
            if lines[i].starts_with("name = \"") {
                let cratex = lines[i].between("name = \"", "\"");
                if !lines[i + 1].starts_with("version = \"") {
                    eprintln!("\nProblem with Cargo.lock entry for crate {}.\n", cratex);
                    std::process::exit(1);
                }
                let version = lines[i + 1].between("\"", "\"");
                crates.push(cratex.to_string());
                versions.push(version.to_string());
                i += 1;
            }
            i += 1;
        }
        println!("\nupdating all {} crates", crates.len());

        // Attempt to update each crate.

        let mut dots = 0;
        println!();
        for i in 0..crates.len() {
            let cratex = &crates[i];
            let version = &versions[i];

            // Update crate.

            let new = Command::new("cargo")
                .arg("update")
                .arg("-p")
                .arg(&format!("{}:{}", cratex, version))
                .output()
                .unwrap_or_else(|_| panic!("{}", "\n\nfailed to execute cargo update".to_string()));
            if new.status.code() != Some(0) {
                print_dot(&mut dots);
                continue;
            }

            // See what changed in Cargo.lock.  In order to declare a change, we require that the
            // crate in question changed.  Other things can change, and that's flaky, and if it
            // happens, we don't do anything.

            let new = Command::new("git")
                .arg("diff")
                .output()
                .unwrap_or_else(|_| panic!("{}", "\n\nfailed to execute git diff".to_string()));
            if new.status.code() != Some(0) {
                println!("\n\ngit diff failed, something is wrong");
                std::process::exit(1);
            }
            let updated = strme(&new.stdout).contains(&format!("+ \"{} ", cratex));
            if !updated {
                let new = Command::new("git")
                    .arg("checkout")
                    .arg("Cargo.lock")
                    .output()
                    .unwrap_or_else(|_| {
                        panic!("{}", "\n\nfailed to execute git checkout".to_string())
                    });
                if new.status.code() != Some(0) {
                    println!("\n\ngit checkout failed, something is wrong");
                    println!("\ntry checking for disk quota exceeded\n");
                    std::process::exit(1);
                }
                print_dot(&mut dots);
                continue;
            }

            // Now try to compile.

            let new = Command::new("cargo")
                .arg("b")
                .output()
                .unwrap_or_else(|_| panic!("{}", "\n\nfailed to execute cargo b".to_string()));
            if new.status.code() != Some(0) {
                let new = Command::new("git")
                    .arg("checkout")
                    .arg("Cargo.lock")
                    .output()
                    .unwrap_or_else(|_| {
                        panic!("{}", "\n\nfailed to execute git checkout".to_string())
                    });
                if new.status.code() != Some(0) {
                    println!("\n\ngit checkout failed, something is wrong");
                    println!("err =\n{}", strme(&new.stderr));
                    std::process::exit(1);
                }
                print_dot(&mut dots);
                continue;
            }

            // Finally, commit the change.

            if dots > 0 {
                println!();
                dots = 0;
            }
            println!("committing change to {}", cratex);
            changed = true;
            let new = Command::new("git")
                .arg("commit")
                .arg("-a")
                .arg("-m")
                .arg(&format!("update crate {}", cratex))
                .output()
                .unwrap_or_else(|_| panic!("{}", "\n\nfailed to execute git commit".to_string()));
            if new.status.code() != Some(0) {
                println!("\n\ngit commit failed, something is wrong");
                std::process::exit(1);
            }
        }
        if dots > 0 {
            println!();
        }
        println!();
        if !changed {
            break;
        }
    }
    println!("done, used {:.1} minutes\n", elapsed(&t) / 60.0);
}
