// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use std::env;

use enclone_main::main_enclone::*;

fn main() {
    let args: Vec<String> = env::args().collect();
    main_enclone(&args);
}
