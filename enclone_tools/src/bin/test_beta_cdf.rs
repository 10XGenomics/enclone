// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Play with beta cdf.
//
// usage test_beta_cdf signal noise

use pretty_trace::PrettyTrace;
use statrs::distribution::ContinuousCDF;
use std::env;
use string_utils::*;

fn beta_cdf(x: f64, a: f64, b: f64) -> f64 {
    let n = statrs::distribution::Beta::new(a, b).unwrap();
    n.cdf(x)
}

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let s = args[1].force_f64();
    let n = args[2].force_f64();
    let p = 1.0 - beta_cdf(0.9, s + 2.0, n + 2.0);
    println!("p = {:.1}%", 100.0 * p);
}
