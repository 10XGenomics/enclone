// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::*;

pub fn gene_scan_test(
    ctl: &EncloneControl,
    stats: &Vec<(String, Vec<f64>)>,
    stats_orig: &Vec<(String, Vec<f64>)>,
    nexacts: usize,
    n: usize,
    in_test: &mut Vec<bool>,
    in_control: &mut Vec<bool>,
) {
    // See if we're in the test and control sets for gene scan (non-exact case).

    if ctl.gen_opt.gene_scan_test.is_some() && !ctl.gen_opt.gene_scan_exact {
        let x = ctl.gen_opt.gene_scan_test.clone().unwrap();
        let mut means = Vec::<f64>::new();
        for i in 0..x.n() {
            let mut vals = Vec::<f64>::new();
            for j in 0..stats.len() {
                if stats[j].0 == x.var[i] {
                    vals.append(&mut stats[j].1.clone());
                    break;
                }
            }
            let mut mean = 0.0;
            for j in 0..vals.len() {
                mean += vals[j];
            }
            mean /= n as f64;
            means.push(mean);
        }
        in_test.push(x.satisfied(&means));
        let x = ctl.gen_opt.gene_scan_control.clone().unwrap();
        let mut means = Vec::<f64>::new();
        for i in 0..x.n() {
            let mut vals = Vec::<f64>::new();
            for j in 0..stats.len() {
                if stats[j].0 == x.var[i] {
                    vals.append(&mut stats[j].1.clone());
                    break;
                }
            }
            let mut mean = 0.0;
            for j in 0..vals.len() {
                mean += vals[j];
            }
            mean /= n as f64;
            means.push(mean);
        }
        in_control.push(x.satisfied(&means));
    }

    // See if we're in the test and control sets for gene scan (exact case).

    if ctl.gen_opt.gene_scan_test.is_some() && ctl.gen_opt.gene_scan_exact {
        let x = ctl.gen_opt.gene_scan_test.clone().unwrap();
        for k in 0..nexacts {
            let mut means = Vec::<f64>::new();
            for i in 0..x.n() {
                let mut vals = Vec::<f64>::new();
                let mut count = 0;
                for j in 0..stats_orig.len() {
                    if stats_orig[j].0 == x.var[i] {
                        if count == k {
                            vals.append(&mut stats_orig[j].1.clone());
                            break;
                        } else {
                            count += 1;
                        }
                    }
                }
                let mut mean = 0.0;
                for j in 0..vals.len() {
                    mean += vals[j];
                }
                mean /= vals.len() as f64;
                means.push(mean);
            }
            in_test.push(x.satisfied(&means));
            let x = ctl.gen_opt.gene_scan_control.clone().unwrap();
            let mut means = Vec::<f64>::new();
            for i in 0..x.n() {
                let mut vals = Vec::<f64>::new();
                let mut count = 0;
                for j in 0..stats_orig.len() {
                    if stats_orig[j].0 == x.var[i] {
                        if count == k {
                            vals.append(&mut stats_orig[j].1.clone());
                            break;
                        } else {
                            count += 1;
                        }
                    }
                }
                let mut mean = 0.0;
                for j in 0..vals.len() {
                    mean += vals[j];
                }
                mean /= n as f64;
                means.push(mean);
            }
            in_control.push(x.satisfied(&means));
        }
    }
}
