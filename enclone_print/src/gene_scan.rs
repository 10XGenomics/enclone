// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::EncloneControl;
use string_utils::TextUtils;

pub fn gene_scan_test(
    ctl: &EncloneControl,
    stats: &Vec<(String, Vec<String>)>,
    stats_orig: &Vec<(String, Vec<String>)>,
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
                    for k in 0..stats[j].1.len() {
                        if stats[j].1[k].parse::<f64>().is_ok() {
                            vals.push(stats[j].1[k].force_f64());
                        }
                    }
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
                    for k in 0..stats[j].1.len() {
                        if stats[j].1[k].parse::<f64>().is_ok() {
                            vals.push(stats[j].1[k].force_f64());
                        }
                    }
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
                            for l in 0..stats_orig[j].1.len() {
                                if stats_orig[j].1[l].parse::<f64>().is_ok() {
                                    vals.push(stats_orig[j].1[l].force_f64());
                                }
                            }
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
                            for l in 0..stats_orig[j].1.len() {
                                if stats_orig[j].1[l].parse::<f64>().is_ok() {
                                    vals.push(stats_orig[j].1[l].force_f64());
                                }
                            }
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
