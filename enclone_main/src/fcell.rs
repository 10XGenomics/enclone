// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Filter using constraints imposed by FCELL.

use enclone_core::defs::*;
use evalexpr::*;
use string_utils::*;
use vector_utils::*;

pub fn filter_by_fcell(
    ctl: &EncloneControl,
    orbits: &mut Vec<Vec<i32>>,
    info: &Vec<CloneInfo>,
    exact_clonotypes: &mut Vec<ExactClonotype>,
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
                for k in 0..ex.ncells() {
                    let li = ex.clones[k][0].dataset_index;
                    let bc = &ex.clones[k][0].barcode;
                    let mut keep = true;
                    for x in ctl.clono_filt_opt_def.fcell.iter() {
                        let alt = &ctl.origin_info.alt_bc_fields[li];
                        let vars = x.iter_variable_identifiers().collect::<Vec<&str>>();
                        let mut vals = Vec::<String>::new();
                        for m in 0..vars.len() {
                            let mut val = String::new();
                            'uloop: for u in 0..alt.len() {
                                if alt[u].0 == vars[m] {
                                    if alt[u].1.contains_key(&bc.clone()) {
                                        val = alt[u].1[&bc.clone()].clone();
                                        break 'uloop;
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
                        to_delete[k] = true;
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
