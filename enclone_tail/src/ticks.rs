// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Find tick marks, given a range on an axis, and a maximum number of ticks.
// The ticks lie in the range, and if the left end of the range is zero, zero is not included.
//
// This has been designed to approximately match the behavior of plotters 0.3.0 on a
// set of test cases.
//
// Searching the internet for tick marks algorithm yields some other methods.

use string_utils::*;

// Given x ≠ 0, find r and s such that x = r * 10^s, and 1.0 <= |r| < 10.

fn normalize_f32(x: f32, r: &mut f32, s: &mut isize) {
    assert!(x != 0.0);
    let y = x.abs();
    *s = y.log10().floor() as isize;
    *r = x * 10.0_f32.powi(-(*s as i32));
}

pub fn ticks(low: f32, high: f32, max_ticks: usize, verbose: bool) -> Vec<String> {
    assert!(low <= high);
    if verbose {
        println!("\nlow = {}, high = {}", low, high);
    }
    if low == high {
        return vec![format!("{}", low)];
    }
    let mut low = low;
    if low == 0.0 {
        low = std::f32::MIN_POSITIVE * 1000.0; // tests fail without some multiplier
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

    // Find best ticks.  The result of this is a an integer vector ns and an integer p.
    // The ticks are ns x 10^p.

    let mut best = 0;
    let mut best_ns = Vec::<i32>::new();
    let mut best_p = 0_i32;
    for q in [p, p - 1, p - 2].iter() {
        let (s1, s2) = (s1 as i32, s2 as i32);
        let p = *q as i32;
        let r1 = r1 * 10.0_f32.powi(s1 - p);
        let r2 = r2 * 10.0_f32.powi(s2 - p);
        let n1 = r1.ceil() as i32;
        let n2 = r2.floor() as i32;
        let mut ns = Vec::<i32>::new();
        for n in n1..=n2 {
            ns.push(n);
        }
        if ns.len() > 1 {
            if ns.len() > best && ns.len() <= max_ticks {
                best = ns.len();
                best_ns = ns.clone();
                best_p = p;
            } else {
                for div in [2, 5].iter() {
                    let mut nsx = Vec::<i32>::new();
                    for n in ns.iter() {
                        let n = *n;
                        if n % div == 0 {
                            nsx.push(n);
                        }
                    }
                    if nsx.len() > best && nsx.len() <= max_ticks {
                        best = nsx.len();
                        best_ns = nsx.clone();
                        best_p = p
                    }
                }
            }
        }
        if ns.len() >= max_ticks {
            break;
        }
    }
    let ns = best_ns;
    let p = best_p;

    // Format the ticks.

    let mut ticks = Vec::<String>::new();
    if verbose {
        use itertools::Itertools;
        println!("ticks = [{}] x 10^{}", ns.iter().format(", "), p);
    }
    for x in ns.iter() {
        let mut x = *x;
        let mut neg = false;
        if x < 0 {
            neg = true;
            x = -x;
        }
        let mut tick = format!("{}", x);
        if p >= 0 {
            for _ in 0..p {
                tick += "0";
            }
        } else {
            let q = -p as usize;
            if q < tick.len() {
                let t = tick.as_bytes();
                let r = t.len() - q;
                tick = format!("{}.{}", strme(&t[0..r]), strme(&t[r..t.len()]));
            } else {
                let mut zeros = String::new();
                for _ in 0..q - tick.len() {
                    zeros += "0";
                }
                tick = format!("0.{}{}", zeros, tick);
            }
        }
        if neg {
            tick = format!("-{}", tick);
        }
        ticks.push(tick);
    }
    ticks
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_trace::*;
    #[test]
    fn test_ticks() {
        PrettyTrace::new().on();
        let max_ticks = 5;
        let mut examples = Vec::<(f32, f32, Vec<&str>)>::new();

        // Examples designed to match plotters behavior.

        examples.push((-10.961007, -5.8206754, vec!["-10", "-9", "-8", "-7", "-6"]));
        examples.push((0.05962117, 1.02, vec!["0.2", "0.4", "0.6", "0.8", "1.0"]));
        examples.push((0.0, 462.06, vec!["100", "200", "300", "400"]));
        examples.push((0.0, 8.16, vec!["2", "4", "6", "8"]));
        examples.push((0.98, 4.08, vec!["1", "2", "3", "4"]));
        examples.push((0.0, 0.0016109196, vec!["0.0005", "0.0010", "0.0015"]));
        examples.push((96.236, 101.08199, vec!["97", "98", "99", "100", "101"]));

        // Other examples.

        examples.push((1.2301, 1.68, vec!["1.3", "1.4", "1.5", "1.6"]));
        examples.push((2.1266, 2.2134001, vec!["2.14", "2.16", "2.18", "2.20"]));
        examples.push((-2.0, 1.0, vec!["-2", "-1", "0", "1"]));
        examples.push((-0.01, 0.1, vec!["0.00", "0.05", "0.10"]));
        examples.push((-0.03, 0.01, vec!["-0.03", "-0.02", "-0.01", "0.00", "0.01"]));
        examples.push((-0.2, 0.1, vec!["-0.2", "-0.1", "0.0", "0.1"]));

        // Test examples.

        for x in examples.iter() {
            let mut y = Vec::<String>::new();
            for m in x.2.iter() {
                y.push(m.to_string());
            }
            if ticks(x.0, x.1, max_ticks, false) != y {
                assert_eq!(ticks(x.0, x.1, max_ticks, true), y);
            }
        }
    }
}
