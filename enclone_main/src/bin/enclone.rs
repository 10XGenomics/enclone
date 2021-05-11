// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use std::env;

use enclone_main::main_enclone::*;

#[tokio::main]
async fn main() {
    let mut args: Vec<String> = env::args().collect();
    let res = main_enclone(&mut args).await;
    if res.is_err() {
        eprintln!("{}", res.unwrap_err());
        std::process::exit(1);
    }
    std::process::exit(0);
}
