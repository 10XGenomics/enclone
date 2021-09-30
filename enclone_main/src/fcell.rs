// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Filter using constraints imposed by FCELL.

use enclone_core::defs::*;
use enclone_print::proc_lvar1::*;
use evalexpr::*;
use io_utils::*;
use ndarray::s;
use rayon::prelude::*;
use std::env;
use std::thread;
use std::time;
use std::time::Instant;
use string_utils::*;
use vector_utils::*;

pub fn filter_by_fcell(
    ctl: &EncloneControl,
    orbits: &mut Vec<Vec<i32>>,
    info: &Vec<CloneInfo>,
    exact_clonotypes: &mut Vec<ExactClonotype>,
    gex_info: &GexInfo,
) -> Result<(), String> {
    if !ctl.clono_filt_opt_def.fcell.is_empty() {
        // Load the GEX and FB data.  This is quite horrible: the code and computation are
        // duplicated verbatim in stop.rs.

        let tdi = Instant::now();
        let mut d_readers = Vec::<Option<hdf5::Reader>>::new();
        let mut ind_readers = Vec::<Option<hdf5::Reader>>::new();
        for li in 0..ctl.origin_info.n() {
            if !ctl.origin_info.gex_path[li].is_empty() && !gex_info.gex_matrices[li].initialized()
            {
                let x = gex_info.h5_data[li].as_ref();
                if x.is_none() {
                    // THIS FAILS SPORADICALLY, OBSERVED MULTIPLE TIMES,
                    // CAUSING PUSH TO D_READERS BELOW TO FAIL.
                    eprintln!("\nWeird, gex_info.h5_data[li].as_ref() is None.");
                    eprintln!("Path = {}.", ctl.origin_info.gex_path[li]);
                    let current = env::current_dir().unwrap();
                    println!(
                        "The current working directory is {}",
                        current.canonicalize().unwrap().display()
                    );
                    if path_exists(&ctl.origin_info.gex_path[li]) {
                        eprintln!(
                            "The directory that is supposed to contain \
                            raw_feature_bc_matrix.h5 exists."
                        );
                        let list = dir_list(&ctl.origin_info.gex_path[li]);
                        eprintln!(
                            "This directory is {} and its contents are:",
                            ctl.origin_info.gex_path[li]
                        );
                        for i in 0..list.len() {
                            eprintln!("{}.  {}", i + 1, list[i]);
                        }
                        let h5_path =
                            format!("{}/raw_feature_bc_matrix.h5", ctl.origin_info.gex_path[li]);
                        eprintln!("H5 path = {}.", h5_path);
                        if !path_exists(&h5_path) {
                            let mut msg = format!("H5 path {} does not exist.\n", h5_path);
                            msg += "Retrying a few times to see if it appears.\n";
                            for _ in 0..5 {
                                msg += "Sleeping for 0.1 seconds.";
                                thread::sleep(time::Duration::from_millis(100));
                                if !path_exists(&h5_path) {
                                    msg += "Now h5 path does not exist.\n";
                                } else {
                                    msg += "Now h5 path exists.\n";
                                    break;
                                }
                            }
                            msg += "Aborting.\n";
                            return Err(msg);
                        } else {
                            println!("h5 path exists.");
                        }
                    } else {
                        println!("Path exists.");
                    }
                    println!();
                }
                d_readers.push(Some(x.unwrap().as_reader()));
                ind_readers.push(Some(gex_info.h5_indices[li].as_ref().unwrap().as_reader()));
            } else {
                d_readers.push(None);
                ind_readers.push(None);
            }
        }
        let mut h5_data = Vec::<(usize, Vec<u32>, Vec<u32>)>::new();
        for li in 0..ctl.origin_info.n() {
            h5_data.push((li, Vec::new(), Vec::new()));
        }
        h5_data.par_iter_mut().for_each(|res| {
            let li = res.0;
            if !ctl.origin_info.gex_path[li].is_empty()
                && !gex_info.gex_matrices[li].initialized()
                && ctl.gen_opt.h5_pre
            {
                res.1 = d_readers[li].as_ref().unwrap().read_raw().unwrap();
                res.2 = ind_readers[li].as_ref().unwrap().read_raw().unwrap();
            }
        });
        ctl.perf_stats(&tdi, "setting up readers, zero");

        // Proceed.

        let mut orbits2 = Vec::<Vec<i32>>::new();
        for i in 0..orbits.len() {
            let mut o = orbits[i].clone();
            let mut to_deletex = vec![false; o.len()];
            for j in 0..o.len() {
                let x: &CloneInfo = &info[o[j] as usize];
                let ex = &mut exact_clonotypes[x.clonotype_index];
                let mut to_delete = vec![false; ex.ncells()];
                let mut d_all = vec![Vec::<u32>::new(); ex.clones.len()];
                let mut ind_all = vec![Vec::<u32>::new(); ex.clones.len()];
                for l in 0..ex.clones.len() {
                    let li = ex.clones[l][0].dataset_index;
                    let bc = ex.clones[l][0].barcode.clone();
                    if !gex_info.gex_barcodes.is_empty() {
                        let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                        if p >= 0 && !gex_info.gex_matrices[li].initialized() {
                            let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                            let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize; // p+1 OK?
                            if ctl.gen_opt.h5_pre {
                                d_all[l] = h5_data[li].1[z1..z2].to_vec();
                                ind_all[l] = h5_data[li].2[z1..z2].to_vec();
                            } else {
                                d_all[l] = d_readers[li]
                                    .as_ref()
                                    .unwrap()
                                    .read_slice(s![z1..z2])
                                    .unwrap()
                                    .to_vec();
                                ind_all[l] = ind_readers[li]
                                    .as_ref()
                                    .unwrap()
                                    .read_slice(s![z1..z2])
                                    .unwrap()
                                    .to_vec();
                            }
                        }
                    }
                }
                for l in 0..ex.ncells() {
                    let li = ex.clones[l][0].dataset_index;
                    let bc = &ex.clones[l][0].barcode;
                    let mut keep = true;
                    for x in ctl.clono_filt_opt_def.fcell.iter() {
                        let alt = &ctl.origin_info.alt_bc_fields[li];
                        let vars = x.iter_variable_identifiers().collect::<Vec<&str>>();
                        let mut vals = Vec::<String>::new();
                        for m in 0..vars.len() {
                            let var = vars[m].to_string();
                            let mut val = String::new();
                            let mut found = false;
                            'uloop: for u in 0..alt.len() {
                                if alt[u].0 == var && alt[u].1.contains_key(&bc.clone()) {
                                    val = alt[u].1[&bc.clone()].clone();
                                    found = true;
                                    break 'uloop;
                                }
                            }
                            if !found && gex_info.feature_id[li].contains_key(&var) {
                                let p = bin_position(&gex_info.gex_barcodes[li], bc);
                                if p >= 0 {
                                    let fid = gex_info.feature_id[li][&var];
                                    let raw_count = get_gex_matrix_entry(
                                        ctl, gex_info, fid, &d_all, &ind_all, li, l, p as usize,
                                        &var,
                                    );
                                    val = format!("{:.2}", raw_count);
                                }
                            }
                            vals.push(val);
                        }
                        let mut c = HashMapContext::new();
                        for m in 0..vars.len() {
                            if vals[m].parse::<i64>().is_ok() {
                                c.set_value(
                                    vars[m].into(),
                                    evalexpr::Value::from(vals[m].force_i64()),
                                )
                                .unwrap();
                            } else if vals[m].parse::<f64>().is_ok() {
                                c.set_value(
                                    vars[m].into(),
                                    evalexpr::Value::from(vals[m].force_f64()),
                                )
                                .unwrap();
                            } else {
                                c.set_value(vars[m].into(), vals[m].clone().into()).unwrap();
                            }
                        }
                        let res = x.eval_with_context(&c);
                        let ok = res == Ok(evalexpr::Value::from(true));
                        if !ok {
                            keep = false;
                        }
                    }
                    if !keep {
                        to_delete[l] = true;
                    }
                }
                erase_if(&mut ex.clones, &to_delete);
                if ex.ncells() == 0 {
                    to_deletex[j] = true;
                }
            }
            erase_if(&mut o, &to_deletex);
            if !o.is_empty() {
                orbits2.push(o.clone());
            }
        }
        *orbits = orbits2;
    }
    Ok(())
}
