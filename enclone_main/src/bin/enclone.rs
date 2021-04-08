// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use std::env;

use enclone_main::main_enclone::*;

#[tokio::main]
async fn main() {
    let args: Vec<String> = env::args().collect();
    main_enclone(&args);
}
