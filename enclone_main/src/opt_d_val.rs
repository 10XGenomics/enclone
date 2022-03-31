// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use enclone_core::opt_d::opt_d;
use enclone_proto::types::DonorReferenceItem;
use rayon::prelude::*;
use std::time::Instant;
use vdj_ann::refx::RefData;

// Assign a D segment to each "left" column in a clonotype (if we need this information).
// The assignments are to exact subclonotypes, and might differ across a clonotype, even
// though the true values have to be the same.  This is also true for V and J segments,
// although they are less likely to vary.

pub fn make_opt_d_val(
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    refdata: &RefData,
    drefs: &Vec<DonorReferenceItem>,
    opt_d_val: &mut Vec<(usize, Vec<Vec<Vec<usize>>>)>,
) {
    let t = Instant::now();
    let mut need_opt_d_val =
        ctl.clono_group_opt.vdj_refname || ctl.clono_group_opt.vdj_heavy_refname;
    for x in ctl.gen_opt.gvars.iter() {
        if x.starts_with("d_inconsistent_") {
            need_opt_d_val = true;
        }
    }
    if need_opt_d_val {
        for i in 0..exacts.len() {
            opt_d_val.push((i, Vec::new()));
        }
        opt_d_val.par_iter_mut().for_each(|res| {
            let i = res.0;
            res.1 = vec![Vec::<Vec<usize>>::new(); rsi[i].mat.len()];
            for col in 0..rsi[i].mat.len() {
                let mut dvotes = Vec::<Vec<usize>>::new();
                for u in 0..exacts[i].len() {
                    let ex = &exact_clonotypes[exacts[i][u]];
                    let m = rsi[i].mat[col][u];
                    if m.is_some() {
                        let m = m.unwrap();
                        if ex.share[m].left {
                            let mut scores = Vec::<f64>::new();
                            let mut ds = Vec::<Vec<usize>>::new();
                            let mid = rsi[i].mat[col][u].unwrap();
                            opt_d(
                                ex.share[mid].v_ref_id,
                                ex.share[mid].j_ref_id,
                                &ex.share[mid].seq_del,
                                &ex.share[mid].annv,
                                &ex.share[mid].cdr3_aa,
                                refdata,
                                drefs,
                                &mut scores,
                                &mut ds,
                                ctl.gen_opt.jscore_match,
                                ctl.gen_opt.jscore_mismatch,
                                ctl.gen_opt.jscore_gap_open,
                                ctl.gen_opt.jscore_gap_extend,
                                ctl.gen_opt.jscore_bits_multiplier,
                                rsi[i].vpids[col],
                            );
                            let mut opt = Vec::new();
                            if !ds.is_empty() {
                                opt = ds[0].clone();
                            }
                            dvotes.push(opt);
                        }
                    } else {
                        dvotes.push(Vec::new());
                    }
                }
                res.1[col] = dvotes;
            }
        });
    }
    ctl.perf_stats(&t, "computing opt_d");
}
