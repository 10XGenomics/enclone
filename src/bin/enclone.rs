// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

extern crate enclone;
use enclone::*;
use std::env;

use main_enclone::*;

fn main() {
    let args: Vec<String> = env::args().collect();
    main_enclone( &args );
}
