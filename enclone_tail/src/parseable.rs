// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Generate parseable output.

use enclone_core::defs::{justification, EncloneControl, ExactClonotype, POUT_SEP};
use io_utils::{fwrite, fwriteln};
use itertools::Itertools;
use std::collections::HashMap;
use std::io::Write;
use string_utils::strme;
use tables::print_tabular;

pub fn generate_parseable_output(
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    group_ncells: usize,
    i: usize,
    j: usize,
    oo: usize,
    pcols: &Vec<String>,
    pcols_show: &Vec<String>,
    pout: &mut Box<dyn std::io::Write>,
    glog: &mut Vec<u8>,
) {
    if !ctl.parseable_opt.pout.is_empty() {
        let mut rows = Vec::<Vec<String>>::new();
        for m in 0..out_datas[oo].len() {
            out_datas[oo][m].insert("group_id".to_string(), format!("{}", i + 1));
            out_datas[oo][m].insert("group_ncells".to_string(), format!("{}", group_ncells));
            out_datas[oo][m].insert("clonotype_id".to_string(), format!("{}", j + 1));
        }
        if !ctl.parseable_opt.pno_header {
            if ctl.parseable_opt.pout == *"stdout" && (!ctl.gen_opt.noprint || (i == 0 && j == 0)) {
                fwriteln!(glog, "{}", pcols_show.iter().format(","));
            }
            if ctl.parseable_opt.pout == *"stdouth" {
                rows.push(pcols_show.clone());
            }
        }
        let x = &out_datas[oo];
        for (u, y) in x.iter().enumerate() {
            if !ctl.parseable_opt.pbarcode {
                if ctl.parseable_opt.pout != *"stdouth" {
                    for (i, c) in pcols.iter().enumerate() {
                        if i > 0 {
                            if ctl.parseable_opt.pout != *"stdout" {
                                fwrite!(pout, ",");
                            } else {
                                fwrite!(glog, ",");
                            }
                        }
                        if y.contains_key(c) {
                            let val = &y[c];
                            let val = val.replace(POUT_SEP, ",");
                            if !val.contains(',') {
                                if ctl.parseable_opt.pout != *"stdout" {
                                    fwrite!(pout, "{}", val);
                                } else {
                                    fwrite!(glog, "{}", val);
                                }
                            } else if ctl.parseable_opt.pout != *"stdout" {
                                fwrite!(pout, "\"{}\"", val);
                            } else {
                                fwrite!(glog, "\"{}\"", val);
                            }
                        } else if ctl.parseable_opt.pout != *"stdout" {
                            fwrite!(pout, "");
                        } else {
                            fwrite!(glog, "");
                        }
                    }
                    if ctl.parseable_opt.pout != *"stdout" {
                        fwriteln!(pout, "");
                    } else {
                        fwriteln!(glog, "");
                    }
                } else {
                    let mut row = Vec::<String>::new();
                    for c in pcols.iter() {
                        if y.contains_key(c) {
                            let val = &y[c];
                            let val = val.replace(POUT_SEP, ",");
                            row.push(val.clone());
                        } else {
                            row.push("".to_string());
                        }
                    }
                    rows.push(row);
                }
            } else {
                let ex = &exact_clonotypes[exacts[oo][u]];
                let n = ex.ncells();
                if ctl.parseable_opt.pout != *"stdouth" {
                    // Traverse the cells in the exact subclonotype.
                    for m in 0..n {
                        // Traverse the parseable fields to be displayed.
                        for (i, c) in pcols.iter().enumerate() {
                            if i > 0 {
                                if ctl.parseable_opt.pout != *"stdout" {
                                    fwrite!(pout, ",");
                                } else {
                                    fwrite!(glog, ",");
                                }
                            }
                            // Test for whether the out_data contain the field.
                            if y.contains_key(c) {
                                let mut id = 0;
                                let vals = y[c].split(POUT_SEP).collect::<Vec<&str>>();
                                if vals.len() > 1 {
                                    id = m;
                                }
                                if id >= vals.len() {
                                    panic!(
                                        "id >= vals.len() where id = {} and vals.len() \
                                        = {},\nparseable variable = {}, barcodes include \
                                        {}, n = {}, y[c] = {}",
                                        id,
                                        vals.len(),
                                        c,
                                        ex.clones[0][0].barcode,
                                        n,
                                        y[c],
                                    );
                                }
                                let val = vals[id];
                                if !val.contains(',') {
                                    if ctl.parseable_opt.pout != *"stdout" {
                                        fwrite!(pout, "{}", val);
                                    } else {
                                        fwrite!(glog, "{}", val);
                                    }
                                } else if ctl.parseable_opt.pout != *"stdout" {
                                    fwrite!(pout, "\"{}\"", val);
                                } else {
                                    fwrite!(glog, "\"{}\"", val);
                                }
                            } else if ctl.parseable_opt.pout != *"stdout" {
                                fwrite!(pout, "");
                            } else {
                                fwrite!(glog, "");
                            }
                        }
                        if ctl.parseable_opt.pout != *"stdout" {
                            fwriteln!(pout, "");
                        } else {
                            fwriteln!(glog, "");
                        }
                    }
                } else {
                    for m in 0..n {
                        let mut row = Vec::<String>::new();
                        for c in pcols.iter() {
                            if y.contains_key(c) {
                                let mut id = 0;
                                let vals = y[c].split(POUT_SEP).collect::<Vec<&str>>();
                                if vals.len() > 1 {
                                    id = m;
                                }
                                let val = vals[id];
                                row.push(val.to_string());
                            } else {
                                row.push("".to_string());
                            }
                        }
                        rows.push(row);
                    }
                }
            }
        }
        if ctl.parseable_opt.pout == *"stdouth" {
            if ctl.gen_opt.noprint {
                for (k, r) in rows.iter().enumerate() {
                    if k > 0 || (i == 0 && j == 0) {
                        let s = format!("{}\n", r.iter().format("\t"));
                        glog.append(&mut s.as_bytes().to_vec());
                    }
                }
            } else {
                let mut log = Vec::<u8>::new();
                let mut justify = Vec::<u8>::new();
                for x in rows[0].iter() {
                    justify.push(justification(x));
                }
                print_tabular(&mut log, &rows, 2, Some(justify));
                if ctl.gen_opt.noprint && (i > 0 || j > 0) {
                    let mut x = String::new();
                    for (k, line) in strme(&log).lines().enumerate() {
                        if k > 0 {
                            x += &mut format!("{}\n", line);
                        }
                    }
                    log = x.as_bytes().to_vec();
                }
                fwrite!(glog, "{}", strme(&log));
            }
        }
    }
}
