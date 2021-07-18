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
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use std::time::Instant;
use string_utils::*;
use vector_utils::*;

fn print_dot(dots: &mut usize) {
    if *dots > 0 && *dots % 10 == 0 {
        print!(" ");
    }
    print!(".");
    std::io::stdout().flush().unwrap();
    *dots += 1;
    if *dots == 90 {
        println!("");
        std::io::stdout().flush().unwrap();
        *dots = 0;
    }
}

fn reset() {
    let o = Command::new("git")
        .arg("reset")
        .arg("--hard")
        .output()
        .expect(&format!("\n\nfailed to execute git reset"));
    if o.status.code() != Some(0) {
        eprintln!("\ngit reset failed\n");
        std::process::exit(1);
    }
}

fn main() {
    let t = Instant::now();
    PrettyTrace::new().on();

    // Get crate versions from master.toml.

    let mut master = Vec::<(String, String)>::new();
    let f = open_for_read!["master.toml"];
    for line in f.lines() {
        let s = line.unwrap();
        if !s.starts_with('#') && s.contains("=") {
            let cratex = s.before(" = ").to_string();
            let version = s.after(" = ").to_string();
            if version.starts_with('"') && version.ends_with('"') {
                let version = version.after("\"").rev_before("\"").to_string();
                let fields = version.split('.').collect::<Vec<&str>>();
                let mut ok = true;
                for i in 0..fields.len() {
                    if !fields[i].parse::<usize>().is_ok() {
                        ok = false;
                    }
                }
                if ok {
                    master.push((cratex, version));
                }
            }
        }
    }

    // Loop until nothing changes.

    loop {
        let mut changed = false;

        // Look for possible raises to crate versions in master.toml.

        println!("\nexamining {} directly dependent crates", master.len());
        for x in master.iter() {
            let (cratex, old) = (&x.0, &x.1);
            let o = Command::new("cargo")
                .arg("search")
                .arg("--limit")
                .arg("1")
                .arg(&cratex)
                .output()
                .expect(&format!("\n\nfailed to execute cargo search"));
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
                    if !s.starts_with('#') && s.contains("=") {
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
                let o = Command::new("sync_to_master")
                    .output()
                    .expect(&format!("\n\nfailed to execute sync_to_master"));
                if o.status.code() != Some(0) {
                    eprintln!("\nsync_to_master failed");
                    std::process::exit(1);
                }
                let o = Command::new("cargo")
                    .arg("update")
                    .arg("-p")
                    .arg(&cratex)
                    .output()
                    .expect(&format!("\n\nfailed to execute cargo update 0"));
                if o.status.code() != Some(0) {
                    reset();
                    continue;
                }
                let o = Command::new("cargo")
                    .arg("b")
                    .output()
                    .expect(&format!("\n\nfailed to execute cargo b 0"));
                if o.status.code() != Some(0) {
                    reset();
                    continue;
                }
                {
                    let mut f = open_for_write_new!["master.toml"];
                    for line in new_lines.iter() {
                        fwriteln!(f, "{}", line);
                    }
                }
                println!("succeeded, committing change to {}", cratex);
                changed = true;
                let new = Command::new("git")
                    .arg("commit")
                    .arg("-a")
                    .arg("-m")
                    .arg(&format!("raise version of crate {}", cratex))
                    .output()
                    .expect(&format!("\n\nfailed to execute git commit"));
                if new.status.code() != Some(0) {
                    println!("\n\ngit commit failed, something is wrong");
                    std::process::exit(1);
                }
            }
        }

        // Get complete list of crates from Cargo.lock.

        let f = open_for_read!["Cargo.lock"];
        let mut crates = Vec::<String>::new();
        for line in f.lines() {
            let s = line.unwrap();
            if s.starts_with("name = \"") {
                let cratex = s.between("name = \"", "\"");
                crates.push(cratex.to_string());
            }
        }
        unique_sort(&mut crates);
        println!("\nupdating all {} crates", crates.len());

        // Attempt to update each crate.

        let mut dots = 0;
        println!("");
        for cratex in crates.iter() {
            let cratex = &*cratex;

            // Update crate.

            let new = Command::new("cargo")
                .arg("update")
                .arg("-p")
                .arg(&cratex)
                .output()
                .expect(&format!("\n\nfailed to execute cargo update"));
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
                .expect(&format!("\n\nfailed to execute git diff"));
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
                    .expect(&format!("\n\nfailed to execute git checkout"));
                if new.status.code() != Some(0) {
                    println!("\n\ngit checkout failed, something is wrong");
                    std::process::exit(1);
                }
                print_dot(&mut dots);
                continue;
            }

            // Now try to compile.

            let new = Command::new("cargo")
                .arg("b")
                .output()
                .expect(&format!("\n\nfailed to execute cargo b"));
            if new.status.code() != Some(0) {
                let new = Command::new("git")
                    .arg("checkout")
                    .arg("Cargo.lock")
                    .output()
                    .expect(&format!("\n\nfailed to execute git checkout"));
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
                println!("");
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
                .expect(&format!("\n\nfailed to execute git commit"));
            if new.status.code() != Some(0) {
                println!("\n\ngit commit failed, something is wrong");
                std::process::exit(1);
            }
        }
        if dots > 0 {
            println!("");
        }
        println!("");
        if !changed {
            break;
        }
    }
    println!("done, used {:.1} minutes\n", elapsed(&t) / 60.0);
}
