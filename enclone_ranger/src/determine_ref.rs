// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Start of code to determine the reference sequence that is to be used.

use enclone_core::defs::EncloneControl;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn determine_ref(ctl: &mut EncloneControl, refx: &mut String) -> Result<(), String> {
    let fx = File::open(&ctl.gen_opt.refname);
    let f = BufReader::new(fx.unwrap());
    for line in f.lines() {
        let s = line.unwrap();
        *refx += &s;
        *refx += "\n";
    }
    Ok(())
}
