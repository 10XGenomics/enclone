// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Group and print clonotypes.  For now, limited grouping functionality.

use crate::group::group_and_print_clonotypes;

use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype, GexInfo};

use enclone_core::enclone_structs::{BarcodeFates, GexReaders, JoinInfo};
use enclone_proto::types::DonorReferenceItem;
use io_utils::fwrite;
use ndarray::s;
use rayon::prelude::*;

use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;
use string_utils::{strme, TextUtils};
use tables::print_tabular;
use vdj_ann::refx::RefData;
use vector_utils::bin_position;

use hdf5::Reader;

pub fn tail_code(
    tall: &Instant,
    refdata: &RefData,
    pics: &Vec<String>,
    group_pics: &mut Vec<String>,
    last_widths: &mut Vec<u32>,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    join_info: &Vec<JoinInfo>,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    fate: &[BarcodeFates],
    tests: &Vec<usize>,
    controls: &Vec<usize>,
    gex_readers: &[Option<GexReaders<'_>>],
    dref: &Vec<DonorReferenceItem>,
    groups: &Vec<Vec<(i32, String)>>,
    opt_d_val: &Vec<(usize, Vec<Vec<Vec<usize>>>)>,
    svgs: &mut Vec<String>,
    summary: &mut String,
) -> Result<(), String> {
    // Print clonotypes.

    group_and_print_clonotypes(
        tall,
        refdata,
        pics,
        group_pics,
        last_widths,
        exacts,
        rsi,
        exact_clonotypes,
        ctl,
        out_datas,
        join_info,
        gex_info,
        vdj_cells,
        fate,
        dref,
        groups,
        opt_d_val,
        svgs,
        summary,
    )?;

    // Do gene scan.

    if ctl.gen_opt.gene_scan.is_some() {
        println!("\nFEATURE SCAN\n");
        let mut test_cells = 0;
        if !ctl.gen_opt.gene_scan_exact {
            for i in tests.iter() {
                for u in exacts[*i].iter() {
                    test_cells += exact_clonotypes[*u].ncells();
                }
            }
            println!(
                "{} clonotypes containing {} cells in test set",
                tests.len(),
                test_cells
            );
        } else {
            for u in tests.iter() {
                test_cells += exact_clonotypes[*u].ncells();
            }
            println!(
                "{} exact subclonotypes containing {} cells in test set",
                tests.len(),
                test_cells
            );
        }
        let mut control_cells = 0;
        if !ctl.gen_opt.gene_scan_exact {
            for i in controls.iter() {
                for u in exacts[*i].iter() {
                    control_cells += exact_clonotypes[*u].ncells();
                }
            }
            println!(
                "{} clonotypes containing {} cells in control set\n",
                controls.len(),
                control_cells
            );
        } else {
            for u in controls.iter() {
                control_cells += exact_clonotypes[*u].ncells();
            }
            println!(
                "{} exact subclonotypes containing {} cells in control set\n",
                controls.len(),
                control_cells
            );
        }
        if tests.is_empty() {
            if !ctl.gen_opt.gene_scan_exact {
                return Err("Gene scan failed, no test clonotypes.\n".to_string());
            } else {
                return Err("Gene scan failed, no test exact subclonotypes.\n".to_string());
            }
        }
        if controls.is_empty() {
            if !ctl.gen_opt.gene_scan_exact {
                return Err("Gene scan failed, no control clonotypes.\n".to_string());
            } else {
                return Err("Gene scan failed, no control exact subclonotypes.\n".to_string());
            }
        }
        println!("enriched features\n");
        let mut results = Vec::<(usize, Vec<u8>, f64, f64, f64)>::new();
        let nf = gex_info.gex_features[0].len();
        for fid in 0..nf {
            results.push((fid, Vec::<u8>::new(), 0.0, 0.0, 0.0));
        }
        results.par_iter_mut().for_each(|res| {
            let fid = res.0;
            // NOT SURE THIS IS BACKWARD COMPATIBLE!
            let gene = gex_info.gex_features[0][fid]
                .after("\t")
                .after("\t")
                .contains("Gene");
            let mut test_values = Vec::<f64>::new();
            let mut control_values = Vec::<f64>::new();
            for pass in 1..=2 {
                let tc;
                let vals;
                if pass == 1 {
                    tc = &tests;
                    vals = &mut test_values;
                } else {
                    tc = &controls;
                    vals = &mut control_values;
                }
                if !ctl.gen_opt.gene_scan_exact {
                    for j in 0..tc.len() {
                        for m in 0..exacts[tc[j]].len() {
                            let ex = &exact_clonotypes[exacts[tc[j]][m]];
                            for l in 0..ex.clones.len() {
                                let li = ex.clones[l][0].dataset_index;
                                let bc = ex.clones[l][0].barcode.clone();
                                let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                if p >= 0 {
                                    let mut raw_count = 0 as f64;

                                    let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                                    // p+1 OK?
                                    let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize;
                                    let d: Vec<u32>;
                                    let ind: Vec<u32>;
                                    let (d, ind) = gex_readers[li]
                                        .as_ref()
                                        .unwrap()
                                        .get_range(z1..z2)
                                        .unwrap();
                                    for j in 0..d.len() {
                                        if ind[j] == fid as u32 {
                                            raw_count = d[j] as f64;
                                            break;
                                        }
                                    }

                                    let mult: f64;
                                    if gene {
                                        mult = gex_info.gex_mults[li];
                                    } else {
                                        mult = gex_info.fb_mults[li];
                                    }
                                    if !ctl.gen_opt.full_counts {
                                        vals.push(raw_count * mult);
                                    } else {
                                        vals.push(raw_count);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    for j in 0..tc.len() {
                        let ex = &exact_clonotypes[tc[j]];
                        for l in 0..ex.clones.len() {
                            let li = ex.clones[l][0].dataset_index;
                            let bc = ex.clones[l][0].barcode.clone();
                            let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                            if p >= 0 {
                                let mut raw_count = 0 as f64;
                                let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                                // p+1 OK?
                                let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize;
                                let d: Vec<u32>;
                                let ind: Vec<u32>;
                                let (d, ind) =
                                    gex_readers[li].as_ref().unwrap().get_range(z1..z2).unwrap();
                                for j in 0..d.len() {
                                    if ind[j] == fid as u32 {
                                        raw_count = d[j] as f64;
                                        break;
                                    }
                                }

                                let mult: f64;
                                if gene {
                                    mult = gex_info.gex_mults[li];
                                } else {
                                    mult = gex_info.fb_mults[li];
                                }
                                if !ctl.gen_opt.full_counts {
                                    vals.push(raw_count * mult);
                                } else {
                                    vals.push(raw_count);
                                }
                            }
                        }
                    }
                }
            }
            let mut test_mean = 0.0;
            for i in 0..test_values.len() {
                test_mean += test_values[i];
            }
            test_mean /= test_values.len() as f64;
            let mut control_mean = 0.0;
            for i in 0..control_values.len() {
                control_mean += control_values[i];
            }
            control_mean /= control_values.len() as f64;
            let mut vals = Vec::<f64>::new();
            let threshold = ctl.gen_opt.gene_scan.as_ref().unwrap().threshold.clone();
            for i in 0..threshold.var.len() {
                if threshold.var[i] == *"t" {
                    vals.push(test_mean);
                } else {
                    vals.push(control_mean);
                }
            }
            if threshold.satisfied(&vals) {
                fwrite!(res.1, "{}", gex_info.gex_features[0][fid]);
                res.2 = test_mean;
                res.3 = control_mean;
                res.4 = test_mean / control_mean;
            }
        });
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "id".to_string(),
            "name".to_string(),
            "library_type".to_string(),
            "test".to_string(),
            "control".to_string(),
            "enrichment".to_string(),
        ];
        rows.push(row);
        for fid in 0..nf {
            if !results[fid].1.is_empty() {
                let stuff = strme(&results[fid].1);
                let fields = stuff.split('\t').collect::<Vec<&str>>();
                let mut row = Vec::<String>::new();
                row.push(fields[0].to_string());
                row.push(fields[1].to_string());
                row.push(fields[2].to_string());
                row.push(format!("{:.2}", results[fid].2));
                row.push(format!("{:.2}", results[fid].3));
                row.push(format!("{:.2}", results[fid].4));
                rows.push(row);
            }
        }
        let mut log = Vec::<u8>::new();
        print_tabular(&mut log, &rows, 2, Some(b"lllrrr".to_vec()));
        print!("{}", strme(&log));
    }

    // Report time.

    Ok(())
}
