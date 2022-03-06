// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

use crate::align_to_vdj_ref::align_to_vdj_ref;
use crate::defs::{EncloneControl, ExactClonotype};
use crate::opt_d::{jflank, opt_d};
use bio_edit::alignment::AlignmentOperation::{Del, Ins, Subst};
use enclone_proto::types::DonorReferenceItem;
use rayon::prelude::*;
use vdj_ann::refx::RefData;

// This is largely copied from align_n.

pub fn heavy_complexity(
    refdata: &RefData,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    dref: &Vec<DonorReferenceItem>,
) -> Vec<usize> {
    let mut results = Vec::<(usize, usize)>::new();
    for i in 0..exact_clonotypes.len() {
        results.push((i, 0));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let ex = &exact_clonotypes[i];
        for r in 0..ex.share.len() {
            if ex.share[r].left && ex.share.len() == 2 && !ex.share[1 - r].left {
                let mut seq = ex.share[r].seq_del.clone();
                let mut vref = refdata.refs[ex.share[r].v_ref_id].to_ascii_vec();
                if ex.share[r].v_ref_id_donor.is_some() {
                    vref = dref[ex.share[r].v_ref_id_donor_alt_id.unwrap()]
                        .nt_sequence
                        .clone();
                }
                let mut vstart = ex.share[r].cdr3_start - 2;

                // Compensate for indel.  Code here and next work imperfectly and
                // there would be value in investigating the error cases.

                if ex.share[r].ins.len() > 0 {
                    vstart -= ex.share[r].ins[0].1.len();
                } else if ex.share[r].seq.len() < ex.share[r].seq_del.len() {
                    vstart += ex.share[r].seq_del.len() - ex.share[r].seq.len();
                }

                // Prevent crash (working around bug).

                if vstart > vref.len() {
                    vstart = vref.len();
                }

                // Keep going.

                vref = vref[vstart..vref.len()].to_vec();
                let mut concat = vref.clone();
                let mut drefx = Vec::<u8>::new();
                let mut d2ref = Vec::<u8>::new();
                let mut drefname = String::new();
                let mut scores = Vec::<f64>::new();
                let mut ds = Vec::<Vec<usize>>::new();
                opt_d(
                    ex,
                    r,
                    refdata,
                    dref,
                    &mut scores,
                    &mut ds,
                    ctl,
                    ex.share[r].v_ref_id_donor_alt_id,
                );
                let mut opt = Vec::new();
                if !ds.is_empty() {
                    opt = ds[0].clone();
                }
                for j in 0..opt.len() {
                    let d = opt[j];
                    if j == 0 {
                        drefx = refdata.refs[d].to_ascii_vec();
                    } else {
                        d2ref = refdata.refs[d].to_ascii_vec();
                    }
                    if j > 0 {
                        drefname += ":";
                    }
                    drefname += &mut refdata.name[d].clone();
                }
                concat.append(&mut drefx.clone());
                concat.append(&mut d2ref.clone());
                let mut jref = refdata.refs[ex.share[r].j_ref_id].to_ascii_vec();
                let jend = jflank(&seq, &jref);
                let mut seq_start = vstart as isize;
                // probably not exactly right
                if ex.share[r].annv.len() > 1 {
                    let q1 = ex.share[r].annv[0].0 + ex.share[r].annv[0].1;
                    let q2 = ex.share[r].annv[1].0;
                    seq_start += q2 as isize - q1 as isize;
                }
                let mut seq_end = seq.len() - (jref.len() - jend);
                // very flaky bug workaround
                // asserted on BCR=180030 CDR3=CARERDLIWFGPW JALIGN1
                if seq_start as usize > seq_end {
                    seq_start = vstart as isize;
                }
                if seq_end <= seq_start as usize {
                    seq_end = seq.len(); // bug fix for problem found by customer,
                                         // couldn't reproduce internally
                }
                seq = seq[seq_start as usize..seq_end].to_vec();
                jref = jref[0..jend].to_vec();
                concat.append(&mut jref.clone());
                let (ops, _score) =
                    align_to_vdj_ref(&seq, &vref, &drefx, &d2ref, &jref, &drefname, true, ctl);
                for i in 0..ops.len() {
                    if ops[i] == Subst {
                        res.1 += 1;
                    } else if ops[i] == Ins {
                        res.1 += 1;
                    } else if ops[i] == Del && (i == 0 || ops[i - 1] != Del) {
                        res.1 += 1;
                    }
                }
            }
        }
    });
    let mut comp = Vec::<usize>::new();
    for i in 0..results.len() {
        comp.push(results[i].1);
    }
    comp
}
