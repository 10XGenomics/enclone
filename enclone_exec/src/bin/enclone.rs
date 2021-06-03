// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_main::main_enclone::main_enclone;
use enclone_visual::enclone_client::enclone_client;
use enclone_visual::enclone_server::enclone_server;
use pretty_trace::*;
use std::env;
use std::time::Instant;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let t = Instant::now();
    PrettyTrace::new().on();
    let mut args: Vec<String> = env::args().collect();

    // Client run of enclone.

    for i in 0..args.len() {
        if args[i] == "VIS" || args[i].starts_with("VIS=") {
            enclone_client(&t).await?;
        }
    }

    // Standard run of enclone.

    if args.len() < 2 || args[1] != "SERVER" {
        let res = main_enclone(&mut args).await;
        if res.is_err() {
            eprintln!("{}", res.unwrap_err());
            std::process::exit(1);
        }
        std::process::exit(0);
    }

    // Server run of enclone.

    enclone_server().await?;
    Ok(())
}
