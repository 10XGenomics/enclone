// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use io_utils::*;
use itertools::Itertools;
use pretty_trace::*;

// Given x ≠ 0, find r and s such that x = r * 10^s, and 1.0 <= |r| < 10.

fn normalize_f32(x: f32, r: &mut f32, s: &mut isize) {
    assert!(x != 0.0);
    let y = x.abs();
    *s = y.log10().floor() as isize;
    *r = x * 10.0_f32.powi(-(*s as i32));
}

// Return x^p more accurately then x.powi(p).

/*
fn powix(x: f64, p: i32) -> f64 {
    let mut y = 1.0_f64;
    let mut p = p;
    if p >= 0 {
        while p > 0 {
            y *= x;
            p -= 1;
        }
    } else {
        while p < 0 {
            y /= x;
            p += 1;
        }
    }
    y
}
*/

/*
fn round(x: f64, p: i32) {
    let y = x * 10.0_f64.powi(-p);
    let y = y.round() as isize;
    if p < 0 {
        let mut tenp = 1;
        for j in 0..-p {
            tenp += 10;
        }
        let a = p / tenp;
        let b = p % tenp;
*/


// Return n * 10^p, more exactly than n as f64 * 10.0_f64.powi(p).

// fn exact_rep(n: usize, p: u32) -> f64 {


fn main() {
    PrettyTrace::new().on();

    // Always choose tick marks that lie within the range.

    let mut examples = Vec::<(f32, f32)>::new();

    examples.push((-10.961007, -5.8206754));    // -10, -9, -8, -7, -6
    examples.push((0.05962117, 1.02));          // 0.2, 0.4, 0.6, 0.8, 1.0
    examples.push((1.2301,     1.68));          // 1.3, 1.4, 1.5, 1.6 (not in control set)
    examples.push((0.0,        462.06));        // 100, 200, 300, 400
    examples.push((0.0,        8.16));          // 2, 4, 6, 8
    examples.push((0.98,       4.08));          // 1, 2, 3, 4
    examples.push((0.0,        0.0016109196));  // 0.0005, 0.0010, 0.0015

    let max_ticks = 5;

    println!("");
    for ex in 0..examples.len() {

        println!("==========================================================================\n");


        let low = examples[ex].0;
        let high = examples[ex].1;

        let mut low = low;
        if low == 0.0 {
            low += 0.00001;
        }

        // Find ri and si such that:
        // low  = r1 x 10^s1 where -10 < r1 < +10
        // high = r2 x 10^s2 where -10 < r2 < +10.
    
        let (mut r1, mut s1) = (0.0, 0);
        let (mut r2, mut s2) = (0.0, 0);
        normalize_f32(low, &mut r1, &mut s1);
        normalize_f32(high, &mut r2, &mut s2);
    
        // Find the first decimal position p where r and s disagree.
    
        let mut p;
        if s1 > s2 {
            p = s1 as i32;
        } else if s2 > s1 {
            p = s2 as i32;
        } else {
            p = s1 as i32;
            let (mut x1, mut x2) = (r1, r2);
            while x1.floor() == x2.floor() {
                p -= 1;
                x1 *= 10.0;
                x2 *= 10.0;
            }
        }
    
        printme!(low, high);
        // printme!(r1, s1);
        // printme!(r2, s2);
        // printme!(p);
    

        let mut best = 0;
        for q in [p, p - 1].iter() {
            let p = *q;
            // println!("using p = {}", p);
    
            // p -= 1; // OPTIONAL
            // printme!(p);
    
            let s1 = s1 as i32;
            let s2 = s2 as i32;
            let p = p as i32;
    
            let r1 = r1 * 10.0_f32.powi(s1 - p);
            let r2 = r2 * 10.0_f32.powi(s2 - p);
        
            // println!("{} * 10^{}", r1, p);
            // println!("{} * 10^{}", r2, p);
        
            let n1 = r1.ceil() as i32;
            let n2 = r2.floor() as i32;
            
            let mut ns = Vec::<i32>::new();
            for n in n1..=n2 {
                ns.push(n);
            }
            if ns.len() > 1 {
                if ns.len() > best && ns.len() <= max_ticks {
                    best = ns.len();
                    println!("\n{} = [{}] x 10^{}", ns.len(), ns.iter().format(", "), p);

                    let mut nsx = Vec::<f64>::new();
                    for i in 0..ns.len() {
                        nsx.push(ns[i] as f64 * 10.0_f64.powi(p));
                    }
                    println!("\n{} = [{}]", nsx.len(), nsx.iter().format(", "));

                } else {
                    let mut ns2 = Vec::<i32>::new();
                    for n in ns.iter() {
                        let n = *n;
                        if n % 2 == 0 {
                            ns2.push(n);
                        }
                    }
                    if ns2.len() > best && ns2.len() <= max_ticks {
                        let ns = ns2;
                        best = ns.len();
                        println!("\n{} = [{}] x 10^{}", ns.len(), ns.iter().format(", "), p);
                    }
                    let mut ns5 = Vec::<i32>::new();
                    for n in ns.iter() {
                        let n = *n;
                        if n % 5 == 0 {
                            ns5.push(n);
                        }
                    }
                    if ns5.len() > best && ns5.len() <= max_ticks {
                        let ns = ns5;
                        best = ns.len();
                        println!("\n{} = [{}] x 10^{}", ns.len(), ns.iter().format(", "), p);
                    }
                }
            }
            if ns.len() >= max_ticks {
                break;
            }
        
        }
        println!("");
    }
    println!("==========================================================================\n");
}
