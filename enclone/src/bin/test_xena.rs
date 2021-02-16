// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Code to test responsiveness of xena.  Uses hardcoded ids.  It was observed that in some
// cases an individual https call would take about a minute.  This happened sporadically

// Cargo.toml:
// attohttpc = { version = "0.12", default-features = false, features = ["compress", "tls-rustls"] }
// perf_stats = "0.1.2"
// rayon = "1.3"

// use attohttpc::*;
use perf_stats::*;
use rayon::prelude::*;
use std::time::{Duration, Instant};

fn main() {
    let mut ids = Vec::<usize>::new();
    for i in 1096913..=1096926 {
        ids.push(i);
    }

    ids.par_iter_mut().for_each(|q| {
        let t = Instant::now();
        let url = format!("https://xena.fuzzplex.com/api/analyses/{}", q);
        const TIMEOUT: u64 = 120; // timeout in seconds
        let req = attohttpc::get(url).read_timeout(Duration::new(TIMEOUT, 0));
        let response = req.send();
        if response.is_err() {
            eprintln!("\n1 failed to execute access to xena {}\n", q);
            std::process::exit(1);
        }
        let response = response.unwrap();
        if !response.is_success() {
            eprintln!("\n2 failed to execute access to xena {}\n", q);
            std::process::exit(1);
        }
        println!("xena call for {} took {:.2} seconds", q, elapsed(&t));
    });

    /*
        use std::process::Command;
        ids.par_iter_mut().for_each(|q| {
            let t = Instant::now();
            let url = format!("https://xena.fuzzplex.com/api/analyses/{}", q);
            let o = Command::new("curl")
                .arg(url.clone())
                .output()
                .expect("failed to execute xena http");
            let m = String::from_utf8(o.stdout).unwrap();
            if o.status.code() != Some(0) {
                let merr = String::from_utf8(o.stderr).unwrap();
                eprintln!(
                    "\nSomething went wrong. the URL \
                    \nhttp://xena.fuzzplex.com/api/analyses/{}\n\
                    failed with the following stderr:\n{}\n\
                    and the following stdout:\n{}\n",
                    q, merr, m
                );
                std::process::exit(1);
            }
            println!("xena call for {} took {:.2} seconds", q, elapsed(&t));
        });
    */
}
