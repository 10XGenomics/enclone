// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Filter using constraints imposed by FCELL.

use enclone_core::defs::*;
use enclone_print::proc_lvar1::*;
use evalexpr::*;
use ndarray::s;
use string_utils::*;
use vector_utils::*;

pub fn filter_by_fcell(
    ctl: &EncloneControl,
    orbits: &mut Vec<Vec<i32>>,
    info: &Vec<CloneInfo>,
    exact_clonotypes: &mut Vec<ExactClonotype>,
    gex_info: &GexInfo,
    d_readers: &Vec<Option<hdf5::Reader>>,
    ind_readers: &Vec<Option<hdf5::Reader>>,
    h5_data: &Vec<(usize, Vec<u32>, Vec<u32>)>,
) {
    if !ctl.clono_filt_opt_def.fcell.is_empty() {
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
                        if p >= 0 { 
                            if !gex_info.gex_matrices[li].initialized() {
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
                                if alt[u].0 == var {
                                    if alt[u].1.contains_key(&bc.clone()) {
                                        val = alt[u].1[&bc.clone()].clone();
                                        found = true;
                                        break 'uloop;
                                    }
                                }
                            }
                            if !found {
                                let mut ux = Vec::<usize>::new();
                                if ctl.clono_print_opt.regex_match[li].contains_key(&var) {
                                    ux = ctl.clono_print_opt.regex_match[li][&var].clone();
                                }
                                if ux.len() > 0 {
                                    let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                    if p >= 0 {
                                        let mut raw_count = 0.0;
                                        for fid in ux.iter() {
                                            let raw_counti = get_gex_matrix_entry(
                                                &ctl, &gex_info, *fid, &d_all, &ind_all, li, l, 
                                                p as usize, &var,
                                            );
                                            raw_count += raw_counti;
                                        }
                                        val = format!("{:.2}", raw_count);
                                    }
                                } else {
                                    if gex_info.feature_id[li].contains_key(&var) {
                                        let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                                        if p >= 0 {
                                            let fid = gex_info.feature_id[li][&var];
                                            let raw_count = get_gex_matrix_entry(
                                                &ctl, &gex_info, fid, &d_all, &ind_all, li, l, 
                                                p as usize, &var,
                                            );
                                            val = format!("{:.2}", raw_count);
                                        }
                                    }
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
}
