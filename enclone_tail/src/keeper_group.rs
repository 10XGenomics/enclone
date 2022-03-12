// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use enclone_proto::types::DonorReferenceItem;
use vdj_ann::refx::RefData;
use vector_utils::unique_sort;
use vector_utils::VecUtils;

pub fn keeper_group(
    o: &Vec<i32>, // clonotypes in group
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    rsi: &Vec<ColInfo>,
    dref: &Vec<DonorReferenceItem>,
) -> bool {
    if o.len() < ctl.clono_group_opt.min_group {
        return false;
    }
    if ctl.clono_group_opt.min_group_donors > 1 {
        let mut donors = Vec::<usize>::new();
        for j in 0..o.len() {
            let x = o[j] as usize;
            let s = &exacts[x];
            for k in 0..s.len() {
                let ex = &exact_clonotypes[s[k]];
                for l in 0..ex.clones.len() {
                    let d = ex.clones[l][0].donor_index;
                    if d.is_some() {
                        donors.push(d.unwrap());
                    }
                }
            }
        }
        unique_sort(&mut donors);
        if donors.len() < ctl.clono_group_opt.min_group_donors {
            return false;
        }
    }
    if ctl.clono_group_opt.cdr3h_len_var {
        let mut lens = Vec::<usize>::new();
        for j in 0..o.len() {
            let x = o[j] as usize;
            let s = &exacts[x];
            for k in 0..s.len() {
                let ex = &exact_clonotypes[s[k]];
                for l in 0..ex.share.len() {
                    if ex.share[l].left {
                        lens.push(ex.share[l].cdr3_aa.len());
                    }
                }
            }
        }
        unique_sort(&mut lens);
        if lens.solo() {
            return false;
        }
    }
    if ctl.clono_group_opt.cdr3.len() > 0 {
        let mut found = false;
        for j in 0..o.len() {
            let x = o[j] as usize;
            let s = &exacts[x];
            for k in 0..s.len() {
                let ex = &exact_clonotypes[s[k]];
                for l in 0..ex.share.len() {
                    if ex.share[l].cdr3_aa == ctl.clono_group_opt.cdr3 {
                        found = true;
                    }
                }
            }
        }
        if !found {
            return false;
        }
    }
    if ctl.clono_group_opt.naive || ctl.clono_group_opt.no_naive {
        let mut found = false;
        // for each clonotype
        for j in 0..o.len() {
            let x = o[j] as usize;
            let rsi = &rsi[x];
            let s = &exacts[x];
            // for each exact subclonotype
            for u in 0..s.len() {
                let mut diffs = 0;
                let ex = &exact_clonotypes[s[u]];
                for m in 0..rsi.mat.len() {
                    if rsi.mat[m][u].is_some() {
                        let r = rsi.mat[m][u].unwrap();
                        let seq = &ex.share[r].seq_del_amino;
                        let mut vref = refdata.refs[rsi.vids[m]].to_ascii_vec();
                        if rsi.vpids[m].is_some() {
                            vref = dref[rsi.vpids[m].unwrap()].nt_sequence.clone();
                        }
                        let jref = refdata.refs[rsi.jids[m]].to_ascii_vec();
                        let z = seq.len();
                        for p in 0..z {
                            let b = seq[p];
                            if p < vref.len() - ctl.heur.ref_v_trim && b != vref[p] {
                                diffs += 1;
                            }
                            if p >= z - (jref.len() - ctl.heur.ref_j_trim)
                                && b != jref[jref.len() - (z - p)]
                            {
                                diffs += 1;
                            }
                        }
                    }
                }
                if diffs == 0 {
                    found = true;
                }
            }
        }
        if !found && ctl.clono_group_opt.naive {
            return false;
        }
        if found && ctl.clono_group_opt.no_naive {
            return false;
        }
    }
    if ctl.clono_group_opt.donor.len() > 0 {
        for d in ctl.clono_group_opt.donor.iter() {
            let mut found = false;
            for j in 0..o.len() {
                let x = o[j] as usize;
                let s = &exacts[x];
                for k in 0..s.len() {
                    let ex = &exact_clonotypes[s[k]];
                    for l in 0..ex.clones.len() {
                        let dx = &ctl.origin_info.donor_id[ex.clones[l][0].dataset_index];
                        if dx == d {
                            found = true;
                        }
                    }
                }
            }
            if !found {
                return false;
            }
        }
    }
    true
}
