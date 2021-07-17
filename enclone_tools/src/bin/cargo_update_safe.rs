// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// For each dependent crate, cargo update it and run cargo b.  If that succeeds,
// commit the change.  Otherwise, undo.

use io_utils::*;
use perf_stats::*;
use pretty_trace::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
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

fn main() {
    let t = Instant::now();
    PrettyTrace::new().on();
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

    // Loop until nothing changes.

    loop {
        let mut dots = 0;
        println!("");
        let mut changed = false;
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
