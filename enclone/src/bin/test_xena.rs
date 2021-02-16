// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Code to test responsiveness of xena.  Uses hardcoded ids.  It was observed that in some
// cases an individual https call would take about a minute.  This happened sporadically

// Cargo.toml:
// attohttpc = { version = "0.12", default-features = false, features = ["compress", "tls-rustls"] }
// perf_stats = "0.1.2"
// rayon = "1.3"

use perf_stats::*;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

fn main() {
    let mut ids = Vec::<usize>::new();
    for i in 1096800..=1097000 {
        ids.push(i);
    }
    let tall = Instant::now();
    let spinlock: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));
    ids.par_iter_mut().for_each(|q| {
        let t = Instant::now();
        while spinlock.load(Ordering::SeqCst) != 0 {}
        let url = format!("https://xena.fuzzplex.com/api/analyses/{}", q);
        const TIMEOUT: u64 = 120; // timeout in seconds
        spinlock.store(1, Ordering::SeqCst);
        let req = attohttpc::get(url).read_timeout(Duration::new(TIMEOUT, 0));
        let response = req.send();
        if response.is_err() || !response.unwrap().is_success() {
            eprintln!("\nfailed to execute access to xena {}\n", q);
            std::process::exit(1);
        }
        spinlock.store(0, Ordering::SeqCst);
        println!("xena call for {} took {:.2} seconds", q, elapsed(&t));
    });
    eprintln!("used {:.2} seconds total", elapsed(&tall));
}
