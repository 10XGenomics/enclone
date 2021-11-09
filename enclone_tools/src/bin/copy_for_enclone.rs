// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// copy_for_enclone SOURCE=... TARGET=...

use enclone_tools::copy_for_enclone::copy_for_enclone;
use pretty_trace::PrettyTrace;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = std::env::args().collect();
    let mut source = String::new();
    let mut target = String::new();
    for i in 1..args.len() {
        if args[i].starts_with("SOURCE=") {
            source = args[i].after("SOURCE=").to_string();
        } else if args[i].starts_with("TARGET=") {
            target = args[i].after("TARGET=").to_string();
        } else {
            eprint!("\nUnknown arg.\n");
            std::process::exit(1);
        }
    }
    assert!(source.len() > 0);
    assert!(target.len() > 0);
    copy_for_enclone(&source, &target);
}
