// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
//
// Download one isoseq dataset.

use pretty_trace::PrettyTrace;
use std::env;
use std::process::Command;
use string_utils::strme;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let acc = &args[1];
    let o = Command::new(
        "/mnt/customer1/swops/03_SOFTWARE/\
        SRA/sratoolkit.2.9.2-centos_linux64/bin/fastq-dump",
    )
    .arg("--fasta")
    .arg(&acc)
    .output()
    .expect("failed to execute fastq-dump");
    if o.status.code().unwrap() != 0 {
        eprintln!("\nFAILED\n");
        eprintln!("stderr:\n{}\n", strme(&o.stderr));
        std::process::exit(1);
    }
}
