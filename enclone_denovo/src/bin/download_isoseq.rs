// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Download the isoseq data in the file ../isoseq_list, that are not already downloaded.
//
// This process is extremely slow.  It may emit messages like this:
//
// 2020-11-12T14:31:01 fastq-dump.2.9.2 sys: timeout exhausted while reading file within network
// system module - mbedtls_ssl_read returned -76 ( NET - Reading information from the
// socket failed )
//
// which appear to be harmless.
//
// It has also failed on occasion.

use pretty_trace::*;
use std::fs::read_dir;
use std::process::Command;
use string_utils::*;
use vector_utils::*;

fn main() {
    PrettyTrace::new().on();
    let dir = "/mnt/assembly/isoseq";
    let all = read_dir(&dir).unwrap();
    let mut owned = Vec::<String>::new();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        if f.contains('.') {
            owned.push(f.rev_after("/").before(".").to_string());
        }
    }
    owned.sort();
    let f = include_str!("../isoseq_list");
    let _ = std::env::set_current_dir(&dir);
    for line in f.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let line = line.replace(" ", "");
        let acc = line.after(":").split(',').collect::<Vec<&str>>();
        for i in 0..acc.len() {
            if !bin_member(&owned, &acc[i].to_string()) {
                println!("downloading {}", acc[i]);
                let o = Command::new(
                    "/mnt/customer1/swops/03_SOFTWARE/\
                    SRA/sratoolkit.2.9.2-centos_linux64/bin/fastq-dump",
                )
                .arg("--fasta")
                .arg(&acc[i])
                .output()
                .expect("failed to execute fastq-dump");
                if o.status.code().unwrap() != 0 {
                    eprintln!("\nFAILED\n");
                    eprintln!("stderr:\n{}\n", strme(&o.stderr));
                    std::process::exit(1);
                }
            }
        }
    }
}
