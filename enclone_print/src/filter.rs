// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Test a clonotype to see if it passes the filters.
// See also enclone_core src for a list of these filters and
// the related struct.

use vdj_ann::refx;

use self::refx::RefData;
use edit_distance::edit_distance;
use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype, GexInfo};
use enclone_core::opt_d::opt_d;
use enclone_proto::types::DonorReferenceItem;
use std::cmp::{max, min};
use string_utils::TextUtils;
use vector_utils::{make_freq, next_diff, unique_sort};

pub fn survives_filter(
    exacts: &Vec<usize>,
    rsi: &ColInfo,
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    refdata: &RefData,
    gex_info: &GexInfo,
    dref: &Vec<DonorReferenceItem>,
) -> bool {
    let mut mults = Vec::<usize>::new();
    for i in 0..exacts.len() {
        mults.push(exact_clonotypes[exacts[i]].clones.len());
    }
    let n: usize = mults.iter().sum();
    if n == 0 {
        return false;
    }

    // Clonotypes with at least n cells

    if n < ctl.clono_filt_opt.ncells_low {
        return false;
    }

    // Clonotypes having iNKT or MAIT evidence

    if ctl.clono_filt_opt.inkt {
        let mut evidence = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            if ex.share[0].inkt_alpha_chain_gene_match
                || ex.share[0].inkt_alpha_chain_junction_match
                || ex.share[0].inkt_beta_chain_gene_match
                || ex.share[0].inkt_beta_chain_junction_match
            {
                evidence = true;
            }
        }
        if !evidence {
            return false;
        }
    }
    if ctl.clono_filt_opt.mait {
        let mut evidence = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            if ex.share[0].mait_alpha_chain_gene_match
                || ex.share[0].mait_alpha_chain_junction_match
                || ex.share[0].mait_beta_chain_gene_match
                || ex.share[0].mait_beta_chain_junction_match
            {
                evidence = true;
            }
        }
        if !evidence {
            return false;
        }
    }

    // Clonotypes marked by heuristics

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

    // Marked clonotypes which are also B cells by annotation

    if ctl.clono_filt_opt_def.marked_b {
        let mut marked_b = false;
        for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            for i in 0..ex.clones.len() {
                if ex.clones[i][0].marked {
                    let li = ex.clones[i][0].dataset_index;
                    let bc = &ex.clones[i][0].barcode;
                    if gex_info.cell_type[li].contains_key(&bc.clone())
                        && gex_info.cell_type[li][&bc.clone()].starts_with('B')
                    {
                        marked_b = true;
                    }
                }
            }
        }
        if !marked_b {
            return false;
        }
    }
    let cols = rsi.vids.len();

    // Barcode required

    if !ctl.clono_filt_opt.barcode.is_empty() {
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

    // Clonotypes with deletions

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

    // Clonotypes with same V gene in 2 chains

    if ctl.clono_filt_opt.vdup {
        let mut dup = false;
        let mut x = rsi.vids.clone();
        x.sort_unstable();
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

    // Clonotypes with constant region differences

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

    // Clonotypes with onesie exact subclonotypes

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

    // Clonotypes with full length V..J

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

    // Clonotypes with no more than n cells

    if n > ctl.clono_filt_opt.ncells_high {
        return false;
    }

    // Clonotypes with at least n chains or at most n chains

    if exacts.len() < ctl.clono_filt_opt.min_exacts {
        return false;
    }
    if exacts.len() > ctl.clono_filt_opt.max_exacts {
        return false;
    }

    // Clonotypes with given V gene name

    for i in 0..ctl.clono_filt_opt.seg.len() {
        let mut hit = false;
        for j in 0..ctl.clono_filt_opt.seg[i].len() {
            for cx in 0..cols {
                if refdata.name[rsi.vids[cx]] == ctl.clono_filt_opt.seg[i][j] {
                    hit = true;
                }
                let did = rsi.dids[cx];
                if did.is_some() {
                    let did = did.unwrap();
                    if refdata.name[did] == ctl.clono_filt_opt.seg[i][j] {
                        hit = true;
                    }
                }
                if refdata.name[rsi.jids[cx]] == ctl.clono_filt_opt.seg[i][j] {
                    hit = true;
                }
                if rsi.cids[cx].is_some()
                    && refdata.name[rsi.cids[cx].unwrap()] == ctl.clono_filt_opt.seg[i][j]
                {
                    hit = true;
                }
            }
        }
        if !hit {
            return false;
        }
    }
    for i in 0..ctl.clono_filt_opt.nseg.len() {
        let mut hit = false;
        for j in 0..ctl.clono_filt_opt.nseg[i].len() {
            for cx in 0..cols {
                if refdata.name[rsi.vids[cx]] == ctl.clono_filt_opt.nseg[i][j] {
                    hit = true;
                }
                let did = rsi.dids[cx];
                if did.is_some() {
                    let did = did.unwrap();
                    if refdata.name[did] == ctl.clono_filt_opt.nseg[i][j] {
                        hit = true;
                    }
                }
                if refdata.name[rsi.jids[cx]] == ctl.clono_filt_opt.nseg[i][j] {
                    hit = true;
                }
                if rsi.cids[cx].is_some()
                    && refdata.name[rsi.cids[cx].unwrap()] == ctl.clono_filt_opt.nseg[i][j]
                {
                    hit = true;
                }
            }
        }
        if hit {
            return false;
        }
    }

    // Clonotypes with given V gene number/allele

    for i in 0..ctl.clono_filt_opt.segn.len() {
        let mut hit = false;
        for j in 0..ctl.clono_filt_opt.segn[i].len() {
            for cx in 0..cols {
                if refdata.id[rsi.vids[cx]] == ctl.clono_filt_opt.segn[i][j].force_i32() {
                    hit = true;
                }
                let did = rsi.dids[cx];
                if did.is_some() {
                    let did = did.unwrap();
                    if refdata.id[did] == ctl.clono_filt_opt.segn[i][j].force_i32() {
                        hit = true;
                    }
                }
                if refdata.id[rsi.jids[cx]] == ctl.clono_filt_opt.segn[i][j].force_i32() {
                    hit = true;
                }
                if rsi.cids[cx].is_some()
                    && refdata.id[rsi.cids[cx].unwrap()]
                        == ctl.clono_filt_opt.segn[i][j].force_i32()
                {
                    hit = true;
                }
            }
        }
        if !hit {
            return false;
        }
    }
    for i in 0..ctl.clono_filt_opt.nsegn.len() {
        let mut hit = false;
        for j in 0..ctl.clono_filt_opt.nsegn[i].len() {
            for cx in 0..cols {
                if refdata.id[rsi.vids[cx]] == ctl.clono_filt_opt.nsegn[i][j].force_i32() {
                    hit = true;
                }
                let did = rsi.dids[cx];
                if did.is_some() {
                    let did = did.unwrap();
                    if refdata.id[did] == ctl.clono_filt_opt.nsegn[i][j].force_i32() {
                        hit = true;
                    }
                }
                if refdata.id[rsi.jids[cx]] == ctl.clono_filt_opt.nsegn[i][j].force_i32() {
                    hit = true;
                }
                if rsi.cids[cx].is_some()
                    && refdata.id[rsi.cids[cx].unwrap()]
                        == ctl.clono_filt_opt.nsegn[i][j].force_i32()
                {
                    hit = true;
                }
            }
        }
        if hit {
            return false;
        }
    }

    // Clonotypes with at least n cells

    if mults.iter().sum::<usize>() < ctl.clono_filt_opt.ncells_low {
        return false;
    }
    let mut numi = 0;
    for i in 0..exacts.len() {
        let ex = &exact_clonotypes[exacts[i]];
        numi = max(numi, ex.max_umi_count());
    }

    // Clonotypes with at least n UMIs for contig

    if numi < ctl.clono_filt_opt.min_umi {
        return false;
    }
    let mut lis = Vec::<usize>::new();
    for s in exacts.iter() {
        let mut z = exact_clonotypes[*s].dataset_indices();
        lis.append(&mut z);
    }
    unique_sort(&mut lis);

    // Clonotypes found in at least n datasets

    if lis.len() < ctl.clono_filt_opt.min_datasets {
        return false;
    }

    // Clonotypes found in at least n origins

    let mut origins = Vec::<String>::new();
    for id in lis.iter() {
        origins.push(ctl.origin_info.origin_id[*id].clone());
    }
    unique_sort(&mut origins);
    if origins.len() < ctl.clono_filt_opt.min_origins {
        return false;
    }

    // Clonotypes found in at least n donors

    let mut donors = Vec::<String>::new();
    for id in lis.iter() {
        donors.push(ctl.origin_info.donor_id[*id].clone());
    }
    unique_sort(&mut donors);
    if donors.len() < ctl.clono_filt_opt.min_donors {
        return false;
    }

    // Clonotypes in no more than n datasets

    if lis.len() > ctl.clono_filt_opt.max_datasets {
        return false;
    }

    // Implement MIN_DATASET_RATIO.

    if ctl.clono_filt_opt.min_dataset_ratio > 0 {
        let mut datasets = Vec::<usize>::new();
        for i in 0..exacts.len() {
            let ex = &exact_clonotypes[exacts[i]];
            for j in 0..ex.ncells() {
                datasets.push(ex.clones[j][0].dataset_index);
            }
        }
        datasets.sort_unstable();
        let mut freq = Vec::<(u32, usize)>::new();
        make_freq(&datasets, &mut freq);
        if freq.len() == 1
            || freq[0].0 < ctl.clono_filt_opt.min_dataset_ratio as u32 * max(1, freq[1].0)
        {
            return false;
        }
    }

    // Clonotypes with no more and no less than min and max chains

    if cols < ctl.clono_filt_opt.min_chains || cols > ctl.clono_filt_opt.max_chains {
        return false;
    }

    // Clonotypes with given junction AA sequence

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

    // Clonotypes having given CDR3 (given by Levenshtein distance pattern).

    if !ctl.clono_filt_opt.cdr3_lev.is_empty() {
        let fields = ctl
            .clono_filt_opt
            .cdr3_lev
            .split('|')
            .collect::<Vec<&str>>();
        let mut cdr3 = Vec::<String>::new();
        let mut dist = Vec::<usize>::new();
        for i in 0..fields.len() {
            cdr3.push(fields[i].before("~").to_string());
            dist.push(fields[i].after("~").force_usize());
        }
        let mut ok = false;
        'exact_loop: for s in exacts.iter() {
            let ex = &exact_clonotypes[*s];
            for j in 0..ex.share.len() {
                for k in 0..cdr3.len() {
                    if edit_distance(&ex.share[j].cdr3_aa, &cdr3[k]) <= dist[k] {
                        ok = true;
                        break 'exact_loop;
                    }
                }
            }
        }
        if !ok {
            return false;
        }
    }

    // Donors.

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

    // Inconsistent D genes.

    if ctl.clono_filt_opt.d_inconsistent {
        let mut inconsistent = false;
        for col in 0..rsi.mat.len() {
            let mut dvotes = Vec::<Vec<usize>>::new();
            for u in 0..exacts.len() {
                let ex = &exact_clonotypes[exacts[u]];
                let m = rsi.mat[col][u];
                if m.is_some() {
                    let m = m.unwrap();
                    if ex.share[m].left {
                        let mut scores = Vec::<f64>::new();
                        let mut ds = Vec::<Vec<usize>>::new();
                        opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
                        let mut opt = Vec::new();
                        if !ds.is_empty() {
                            opt = ds[0].clone();
                        }
                        dvotes.push(opt);
                    }
                }
            }
            unique_sort(&mut dvotes);
            if dvotes.len() >= 2 {
                inconsistent = true;
            }
        }
        if !inconsistent {
            return false;
        }
    }

    // None D genes.

    if ctl.clono_filt_opt.d_none {
        let mut none = false;
        for col in 0..rsi.mat.len() {
            for u in 0..exacts.len() {
                let ex = &exact_clonotypes[exacts[u]];
                let m = rsi.mat[col][u];
                if m.is_some() {
                    let m = m.unwrap();
                    if ex.share[m].left {
                        let mut scores = Vec::<f64>::new();
                        let mut ds = Vec::<Vec<usize>>::new();
                        opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
                        let mut opt = Vec::new();
                        if !ds.is_empty() {
                            opt = ds[0].clone();
                        }
                        if opt.is_empty() {
                            none = true;
                        }
                    }
                }
            }
        }
        if !none {
            return false;
        }
    }

    // Second D genes.

    if ctl.clono_filt_opt.d_second {
        let mut second = false;
        for col in 0..rsi.mat.len() {
            for u in 0..exacts.len() {
                let ex = &exact_clonotypes[exacts[u]];
                let m = rsi.mat[col][u];
                if m.is_some() {
                    let m = m.unwrap();
                    if ex.share[m].left {
                        let mut scores = Vec::<f64>::new();
                        let mut ds = Vec::<Vec<usize>>::new();
                        opt_d(ex, col, u, rsi, refdata, dref, &mut scores, &mut ds, ctl);
                        let mut opt = Vec::new();
                        if !ds.is_empty() {
                            opt = ds[0].clone();
                        }
                        if opt.len() == 2 {
                            second = true;
                        }
                    }
                }
            }
        }
        if !second {
            return false;
        }
    }

    // Done.

    true
}
