// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

pub fn median(x: &[usize]) -> f64 {
    let h = x.len() / 2;
    if x.len() % 2 == 1 {
        x[h] as f64
    } else {
        (x[h - 1] + x[h]) as f64 / 2.0
    }
}

pub fn rounded_median(x: &[usize]) -> usize {
    let h = x.len() / 2;
    if x.len() % 2 == 1 {
        x[h]
    } else {
        let s = x[h - 1] + x[h];
        if s % 2 == 0 {
            s / 2
        } else {
            s / 2 + 1
        }
    }
}

pub fn median_f64(x: &[f64]) -> f64 {
    let h = x.len() / 2;
    if x.len() % 2 == 1 {
        x[h]
    } else {
        (x[h - 1] + x[h]) / 2.0
    }
}
