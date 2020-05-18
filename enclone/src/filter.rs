// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Test a clonotype to see if it passes the filters.

use vdj_ann::*;

use self::refx::*;
use enclone_core::defs::*;
use std::cmp::*;
use string_utils::*;
use vector_utils::*;

pub fn survives_filter(
    exacts: &Vec<usize>,
    rsi: &ColInfo,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    refdata: &RefData,
) -> bool {
    let mut mults = Vec::<usize>::new();
    for i in 0..exacts.len() {
        mults.push(exact_clonotypes[exacts[i]].clones.len());
    }
    let n: usize = mults.iter().sum();
    if n == 0 {
        return false;
    }
    if n < ctl.clono_filt_opt.ncells_low {
        return false;
    }
    if ctl.clono_filt_opt.marked {
        let mut marked = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            for i in 0..ex.clones.len() {
                if ex.clones[i][0].marked {
                    marked = true;
                }
            }
        }
        if !marked {
            return false;
        }
    }
    let cols = rsi.vids.len();
    if ctl.clono_filt_opt.barcode.len() > 0 {
        let mut ok = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            for i in 0..ex.clones.len() {
                for j in 0..ctl.clono_filt_opt.barcode.len() {
                    if ex.clones[i][0].barcode == ctl.clono_filt_opt.barcode[j] {
                        ok = true;
                    }
                }
            }
        }
        if !ok {
            return false;
        }
    }
    if ctl.clono_filt_opt.del {
        let mut ok = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            for m in 0..ex.share.len() {
                if ex.share[m].seq_del.contains(&b'-') {
                    ok = true;
                }
            }
        }
        if !ok {
            return false;
        }
    }
    if ctl.clono_filt_opt.vdup {
        let mut dup = false;
        let mut x = rsi.vids.clone();
        x.sort();
        let mut i = 0;
        while i < x.len() {
            let j = next_diff(&x, i);
            if j - i > 1 {
                dup = true;
            }
            i = j;
        }
        if !dup {
            return false;
        }
    }
    if ctl.clono_filt_opt.cdiff {
        let mut cdiff = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            for m in 0..ex.share.len() {
                let cstart = ex.share[m].j_stop;
                let clen = ex.share[m].full_seq.len() - cstart;
                let cid = ex.share[m].c_ref_id;
                if cid.is_some() {
                    let r = &refdata.refs[cid.unwrap()];
                    for i in 0..min(clen, r.len()) {
                        let tb = ex.share[m].full_seq[cstart + i];
                        let rb = r.to_ascii_vec()[i];
                        if tb != rb {
                            cdiff = true;
                        }
                    }
                }
            }
        }
        if !cdiff {
            return false;
        }
    }
    if ctl.clono_filt_opt.have_onesie {
        let mut have = false;
        for i in 0..exacts.len() {
            if exact_clonotypes[exacts[i]].share.len() == 1 {
                have = true;
            }
        }
        if !have {
            return false;
        }
    }
    if !ctl.clono_filt_opt.vj.is_empty() {
        let mut have_vj = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            for j in 0..ex.share.len() {
                if ex.share[j].seq == ctl.clono_filt_opt.vj {
                    have_vj = true;
                }
            }
        }
        if !have_vj {
            return false;
        }
    }
    if n > ctl.clono_filt_opt.ncells_high {
        return false;
    }
    if exacts.len() < ctl.clono_filt_opt.min_exacts {
        return false;
    }
    if !ctl.clono_filt_opt.seg.is_empty() {
        let mut hit = false;
        for j in 0..ctl.clono_filt_opt.seg.len() {
            for cx in 0..cols {
                if refdata.name[rsi.vids[cx]] == ctl.clono_filt_opt.seg[j] {
                    hit = true;
                }
                let did = rsi.dids[cx];
                if did.is_some() {
                    let did = did.unwrap();
                    if refdata.name[did] == ctl.clono_filt_opt.seg[j] {
                        hit = true;
                    }
                }
                if refdata.name[rsi.jids[cx]] == ctl.clono_filt_opt.seg[j] {
                    hit = true;
                }
                if rsi.cids[cx].is_some() {
                    if refdata.name[rsi.cids[cx].unwrap()] == ctl.clono_filt_opt.seg[j] {
                        hit = true;
                    }
                }
            }
        }
        if !hit {
            return false;
        }
    }
    if !ctl.clono_filt_opt.segn.is_empty() {
        let mut hit = false;
        for j in 0..ctl.clono_filt_opt.segn.len() {
            for cx in 0..cols {
                if refdata.id[rsi.vids[cx]] == ctl.clono_filt_opt.segn[j].force_i32() {
                    hit = true;
                }
                let did = rsi.dids[cx];
                if did.is_some() {
                    let did = did.unwrap();
                    if refdata.id[did] == ctl.clono_filt_opt.segn[j].force_i32() {
                        hit = true;
                    }
                }
                if refdata.id[rsi.jids[cx]] == ctl.clono_filt_opt.segn[j].force_i32() {
                    hit = true;
                }
                if rsi.cids[cx].is_some() {
                    if refdata.id[rsi.cids[cx].unwrap()] == ctl.clono_filt_opt.segn[j].force_i32() {
                        hit = true;
                    }
                }
            }
        }
        if !hit {
            return false;
        }
    }
    if mults.iter().sum::<usize>() < ctl.clono_filt_opt.ncells_low {
        return false;
    }
    let mut numi = 0;
    for i in 0..exacts.len() {
        let ex = &exact_clonotypes[exacts[i]];
        numi = max(numi, ex.max_umi_count());
    }
    if numi < ctl.clono_filt_opt.min_umi {
        return false;
    }
    let mut lis = Vec::<usize>::new();
    for s in exacts.iter() {
        let mut z = exact_clonotypes[*s].dataset_indices();
        lis.append(&mut z);
    }
    unique_sort(&mut lis);
    if lis.len() < ctl.clono_filt_opt.min_datasets {
        return false;
    }
    if lis.len() > ctl.clono_filt_opt.max_datasets {
        return false;
    }
    if cols < ctl.clono_filt_opt.min_chains || cols > ctl.clono_filt_opt.max_chains {
        return false;
    }
    if ctl.clono_filt_opt.cdr3.is_some() {
        let mut ok = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            for j in 0..ex.share.len() {
                if ctl
                    .clono_filt_opt
                    .cdr3
                    .as_ref()
                    .unwrap()
                    .is_match(&ex.share[j].cdr3_aa)
                {
                    ok = true;
                }
            }
        }
        if !ok {
            return false;
        }
    }
    let mut donors = Vec::<usize>::new();
    for u in 0..exacts.len() {
        let ex = &exact_clonotypes[exacts[u]];
        for m in 0..ex.clones.len() {
            if ex.clones[m][0].donor_index.is_some() {
                donors.push(ex.clones[m][0].donor_index.unwrap());
            }
        }
    }
    unique_sort(&mut donors);
    if ctl.clono_filt_opt.fail_only && donors.len() <= 1 {
        return false;
    }
    return true;
}
