// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Miscellaneous functions.

use crate::innate::mark_innate;
use crate::misc3::study_consensus;
use amino::aa_seq;
use debruijn::dna_string::DnaString;
use enclone_core::defs::{EncloneControl, ExactClonotype, TigData, TigData0, TigData1};
use io_utils::{fwriteln, open_for_write_new};
use rayon::prelude::*;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;
use string_utils::strme;
use vdj_ann::refx::RefData;
use vector_utils::{
    erase_if, next_diff, next_diff12_4, next_diff1_2, next_diff1_3, reverse_sort, unique_sort,
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Filter out putative gel bead contamination.  We look for cases where inside a
// given exact subclonotype, the same first or last half of a barcode is reused, and one
// instance has at least 10-fold higher UMI count.  If the fraction of the "bad"
// clones is at least 20%, delete them.

pub fn filter_gelbead_contamination(
    ctl: &EncloneControl,
    mut clones: &mut Vec<Vec<TigData0>>,
    fate: &mut Vec<(usize, String, String)>,
) {
    const GB_UMI_MULT: usize = 10;
    const GB_MIN_FRAC: f64 = 0.2;
    let mut bch = vec![Vec::<(usize, String, usize, usize)>::new(); 2];
    for l in 0..clones.len() {
        let li = clones[l][0].dataset_index;
        let bc = &clones[l][0].barcode;
        let mut numi = 0;
        for j in 0..clones[l].len() {
            numi += clones[l][j].umi_count;
        }
        bch[0].push((li, bc[0..8].to_string(), numi, l));
        bch[1].push((li, bc[8..16].to_string(), numi, l));
    }
    let mut bad = vec![false; clones.len()];
    for l in 0..2 {
        bch[l].sort();
        let mut m = 0;
        while m < bch[l].len() {
            let n = next_diff12_4(&bch[l], m as i32) as usize;
            let mut count = 0;
            for u1 in m..n {
                for u2 in m..n {
                    if bch[l][u1].2 >= GB_UMI_MULT * bch[l][u2].2 {
                        count += 1;
                    }
                }
            }
            if count as f64 / clones.len() as f64 >= GB_MIN_FRAC {
                for u1 in m..n {
                    for u2 in m..n {
                        if bch[l][u1].2 >= GB_UMI_MULT * bch[l][u2].2 {
                            bad[bch[l][u2].3] = true;
                        }
                    }
                }
            }
            m = n;
        }
    }
    for i in 0..bad.len() {
        if bad[i] {
            fate.push((
                clones[i][0].dataset_index,
                clones[i][0].barcode.clone(),
                "failed WHITEF filter".to_string(),
            ));
        }
    }
    if !ctl.gen_opt.nwhitef {
        erase_if(&mut clones, &bad);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn create_exact_subclonotype_core(
    // inputs:
    tig_bc: &Vec<Vec<TigData>>,
    r: usize,
    s: usize,
    to_delete: &Vec<bool>,
    // outputs:
    share: &mut Vec<TigData1>,
    clones: &mut Vec<Vec<TigData0>>,
) {
    for m in 0..tig_bc[r].len() {
        // Form 5'-UTR consensus sequence.

        let mut utr = Vec::<u8>::new();
        let mut pos = 0;
        let mut last_calls = 0;
        loop {
            let mut calls = Vec::<(u8, u8)>::new(); // (base,qual)
            for t in r..s {
                if !to_delete[t - r] && tig_bc[t][m].v_start > pos {
                    let p = tig_bc[t][m].v_start - pos - 1;
                    calls.push((tig_bc[t][m].full_seq[p], tig_bc[t][m].full_quals[p]));
                }
            }
            if calls.is_empty() || calls.len() < last_calls / 10 {
                break;
            }
            last_calls = calls.len();
            calls.sort_unstable();
            let mut callsx = Vec::<(usize, u8)>::new(); // (qual,base)
            let mut i = 0;
            while i < calls.len() {
                let j = next_diff1_2(&calls, i as i32) as usize;
                let mut q = 0;
                for k in i..j {
                    q += calls[k].1 as usize;
                }
                callsx.push((q, calls[i].0));
                i = j;
            }
            reverse_sort(&mut callsx);
            utr.push(callsx[0].1);
            pos += 1;
        }
        utr.reverse();

        // Form constant region consensus sequence.

        let mut constx = Vec::<u8>::new();
        let mut pos = 0;
        let mut last_calls = 0;
        loop {
            let mut calls = Vec::<(u8, u8)>::new(); // (base,qual)
            for t in r..s {
                if !to_delete[t - r] && tig_bc[t][m].j_stop + pos < tig_bc[t][m].full_seq.len() {
                    let p = tig_bc[t][m].j_stop + pos;
                    calls.push((tig_bc[t][m].full_seq[p], tig_bc[t][m].full_quals[p]));
                }
            }
            if calls.is_empty() || calls.len() < last_calls / 10 {
                break;
            }
            last_calls = calls.len();
            calls.sort_unstable();
            let mut callsx = Vec::<(usize, u8)>::new(); // (qual,base)
            let mut i = 0;
            while i < calls.len() {
                let j = next_diff1_2(&calls, i as i32) as usize;
                let mut q = 0;
                for k in i..j {
                    q += calls[k].1 as usize;
                }
                callsx.push((q, calls[i].0));
                i = j;
            }
            reverse_sort(&mut callsx);
            constx.push(callsx[0].1);
            pos += 1;
        }

        // Form full sequence.

        let mut full = utr.clone();
        let mut z = tig_bc[r][m].seq().to_vec();
        full.append(&mut z);
        full.append(&mut constx);

        // Note that here we are taking the first entry (r), sort of assuming
        // that all the entries are the same, which in principle they should be.

        let aa = aa_seq(tig_bc[r][m].seq(), 0);
        let mut d_start = None;
        if tig_bc[r][m].d_start.is_some() {
            d_start = Some(tig_bc[r][m].d_start.unwrap() + utr.len() - tig_bc[r][m].v_start);
        }
        share.push(TigData1 {
            cdr3_dna: tig_bc[r][m].cdr3_dna.clone(),
            seq: tig_bc[r][m].seq().to_vec(),
            seq_del: tig_bc[r][m].seq().to_vec(), // may get changed later
            seq_del_amino: tig_bc[r][m].seq().to_vec(), // may get changed later
            ins: Vec::new(),                      // may get changed later
            aa_mod_indel: aa,                     // may get changed later
            full_seq: full,
            v_start: utr.len(),
            v_stop: tig_bc[r][m].v_stop + utr.len() - tig_bc[r][m].v_start,
            v_stop_ref: tig_bc[r][m].v_stop_ref,
            d_start,
            j_start: tig_bc[r][m].j_start + utr.len() - tig_bc[r][m].v_start,
            j_start_ref: tig_bc[r][m].j_start_ref,
            j_stop: tig_bc[r][m].j_stop + utr.len() - tig_bc[r][m].v_start,
            u_ref_id: tig_bc[r][m].u_ref_id,
            v_ref_id: tig_bc[r][m].v_ref_id,
            v_ref_id_donor: None,
            v_ref_id_donor_alt_id: None,
            v_ref_id_donor_donor: None,
            d_ref_id: tig_bc[r][m].d_ref_id,
            j_ref_id: tig_bc[r][m].j_ref_id,
            c_ref_id: tig_bc[r][m].c_ref_id,
            fr1_start: tig_bc[r][m].fr1_start,
            fr2_start: tig_bc[r][m].fr2_start,
            fr3_start: tig_bc[r][m].fr3_start,
            cdr1_start: tig_bc[r][m].cdr1_start,
            cdr2_start: tig_bc[r][m].cdr2_start,
            cdr3_aa: tig_bc[r][m].cdr3_aa.clone(),
            cdr3_start: tig_bc[r][m].cdr3_start,
            left: tig_bc[r][m].left,
            chain_type: tig_bc[r][m].chain_type.clone(),
            annv: tig_bc[r][m].annv.clone(),
            // these get set when making CloneInfo objects:
            vs: DnaString::new(),
            vs_notesx: String::new(),
            js: DnaString::new(),
            // iNKT and MAIT annotations
            inkt_alpha_chain_gene_match: false,
            inkt_alpha_chain_junction_match: false,
            inkt_beta_chain_gene_match: false,
            inkt_beta_chain_junction_match: false,
            mait_alpha_chain_gene_match: false,
            mait_alpha_chain_junction_match: false,
            mait_beta_chain_gene_match: false,
            mait_beta_chain_junction_match: false,
        });
    }
    for t in r..s {
        if !to_delete[t - r] {
            let mut x = Vec::<TigData0>::new();
            for m in 0..tig_bc[t].len() {
                x.push(TigData0 {
                    quals: tig_bc[t][m].quals.clone(),
                    v_start: tig_bc[t][m].v_start,
                    j_stop: tig_bc[t][m].j_stop,
                    c_start: tig_bc[t][m].c_start,
                    full_seq: tig_bc[t][m].full_seq.clone(),
                    barcode: tig_bc[t][m].barcode.clone(),
                    tigname: tig_bc[t][m].tigname.clone(),
                    dataset_index: tig_bc[t][m].dataset_index,
                    origin_index: tig_bc[t][m].origin_index,
                    donor_index: tig_bc[t][m].donor_index,
                    tag_index: tig_bc[t][m].tag_index,
                    umi_count: tig_bc[t][m].umi_count,
                    read_count: tig_bc[t][m].read_count,
                    validated_umis: tig_bc[t][m].validated_umis.clone(),
                    non_validated_umis: tig_bc[t][m].non_validated_umis.clone(),
                    invalidated_umis: tig_bc[t][m].invalidated_umis.clone(),
                    frac_reads_used: tig_bc[t][m].frac_reads_used,
                    marked: false,
                });
            }
            clones.push(x);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Find exact subclonotypes.

pub fn find_exact_subclonotypes(
    ctl: &EncloneControl,
    tig_bc: &Vec<Vec<TigData>>,
    refdata: &RefData,
    fate: &mut Vec<HashMap<String, String>>,
) -> Vec<ExactClonotype> {
    let mut exact_clonotypes = Vec::<ExactClonotype>::new();
    let mut r = 0;
    let mut groups = Vec::<(usize, usize)>::new();
    let t = Instant::now();
    while r < tig_bc.len() {
        let mut s = r + 1;
        while s < tig_bc.len() {
            let mut ok = true;
            if tig_bc[s].len() != tig_bc[r].len() {
                break;
            }
            for m in 0..tig_bc[r].len() {
                let (cid1, cid2) = (tig_bc[r][m].c_ref_id, tig_bc[s][m].c_ref_id);
                if tig_bc[s][m].cdr3_dna != tig_bc[r][m].cdr3_dna
                    || tig_bc[s][m].seq() != tig_bc[r][m].seq()

                    // Working around a bug here.  See above for explanation.

                    // || cid1 != cid2 {

                    || ( cid1.is_none() && cid2.is_some() ) || ( cid1.is_some() && cid2.is_none() )
                    || ( cid1.is_some() && cid2.is_some()
                        && refdata.name[cid1.unwrap()] != refdata.name[cid2.unwrap()] )

                    || ( cid1.is_some() && cid2.is_some()
                        && tig_bc[r][m].c_start.unwrap() + tig_bc[s][m].j_stop
                        < tig_bc[s][m].c_start.unwrap() + tig_bc[r][m].j_stop )

                    // Check for different donors if MIX_DONORS specified on command line.
                    // Note funky redundancy in checking each chain

                    || ( !ctl.clono_filt_opt_def.donor
                        && tig_bc[r][m].donor_index != tig_bc[s][m].donor_index )
                {
                    ok = false;
                    break;
                }
            }
            if !ok {
                break;
            }
            s += 1;
        }
        groups.push((r, s));
        r = s;
    }
    ctl.perf_stats(&t, "finding exact subclonotypes one");
    let t = Instant::now();
    let mut results = Vec::<(usize, Vec<ExactClonotype>, Vec<(usize, String, String)>)>::new();
    for i in 0..groups.len() {
        results.push((i, Vec::new(), Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let r = groups[i].0;
        let s = groups[i].1;

        // Print reused barcodes.

        if ctl.gen_opt.reuse {
            for t1 in r..s {
                for t2 in t1 + 1..s {
                    if tig_bc[t1][0].barcode == tig_bc[t2][0].barcode {
                        println!("see reuse of barcode {}", tig_bc[t1][0].barcode);
                        print!(
                            "{}: numis =",
                            ctl.origin_info.dataset_id[tig_bc[t1][0].dataset_index]
                        );
                        for m in 0..tig_bc[t1].len() {
                            print!(" {}", tig_bc[t1][m].umi_count);
                        }
                        println!();
                        print!(
                            "{}: numis =",
                            ctl.origin_info.dataset_id[tig_bc[t2][0].dataset_index]
                        );
                        for m in 0..tig_bc[t2].len() {
                            print!(" {}", tig_bc[t2][m].umi_count);
                        }
                        println!("\n");
                    }
                }
            }
        }

        // Delete reused barcodes.  In principle we could instead choose the instance having
        // higher UMI counts.  Also we might test for more evidence of concurrence, to avoid
        // the case where a barcode was accidentally reused.

        let mut to_delete = vec![false; s - r];
        let mut bc = Vec::<(String, usize)>::new();
        for t in r..s {
            bc.push((tig_bc[t][0].barcode.clone(), t));
        }
        bc.sort();
        let mut i = 0;
        while i < bc.len() {
            let j = next_diff1_2(&bc, i as i32) as usize;
            if j - i >= 2 {
                for k in i..j {
                    let t = bc[k].1;
                    if ctl.clono_filt_opt_def.bc_dup {
                        to_delete[t - r] = true;
                    }
                    res.2.push((
                        tig_bc[t][0].dataset_index,
                        tig_bc[t][0].barcode.clone(),
                        "failed BC_DUP filter".to_string(),
                    ));
                }
            }
            i = j;
        }

        // Create the exact subclonotype.

        let mut share = Vec::<TigData1>::new();
        let mut clones = Vec::<Vec<TigData0>>::new();
        create_exact_subclonotype_core(tig_bc, r, s, &to_delete, &mut share, &mut clones);

        // Explore consensus.

        let mut _count = 0;
        study_consensus(
            &mut _count,
            ctl,
            &share,
            &clones,
            &exact_clonotypes,
            refdata,
        );

        // Filter out putative gel bead contamination.

        filter_gelbead_contamination(ctl, &mut clones, &mut res.2);

        // Save exact subclonotype.

        if !clones.is_empty() {
            res.1.push(ExactClonotype { share, clones });
        }
    });
    ctl.perf_stats(&t, "finding exact subclonotypes two");
    let t = Instant::now();
    let mut max_exact = 0;
    for i in 0..results.len() {
        if !results[i].1.is_empty() {
            max_exact = max(max_exact, results[i].1[0].ncells());
            exact_clonotypes.append(&mut results[i].1);
        }
        for j in 0..results[i].2.len() {
            fate[results[i].2[j].0].insert(results[i].2[j].1.clone(), results[i].2[j].2.clone());
        }
    }
    if ctl.gen_opt.utr_con || ctl.gen_opt.con_con {
        println!();
        return Vec::new();
    }
    if !ctl.silent {
        println!(
            "found {} exact subclonotypes from {} productive pairs",
            exact_clonotypes.len(),
            tig_bc.len()
        );
        println!("max exact subclonotype size = {}", max_exact);
    }

    // Edit if NWEAK_ONESIES not specified.

    /*
    if ctl.clono_filt_opt.weak_onesies {
        let mut total_cells = 0;
        for i in 0..exact_clonotypes.len() {
            total_cells += exact_clonotypes[i].ncells();
        }
        let mut exacts2 = Vec::<ExactClonotype>::new();
        for i in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[i];
            if ex.share.len() == 1 && ex.ncells() > 1 && ex.ncells() * 1000 < total_cells {
                for j in 0..ex.clones.len() {
                    exacts2.push(ExactClonotype {
                        share: ex.share.clone(),
                        clones: vec![ex.clones[j].clone()],
                    });
                }
            } else {
                exacts2.push(exact_clonotypes[i].clone());
            }
        }
        exact_clonotypes = exacts2;
    }
    */

    // Fill in iNKT and MAIT annotations.

    mark_innate(refdata, &mut exact_clonotypes);

    // Do other stuff.

    if !ctl.gen_opt.fasta.is_empty() {
        let mut f = open_for_write_new![&ctl.gen_opt.fasta];
        for i in 0..exact_clonotypes.len() {
            let x = &exact_clonotypes[i];
            for j in 0..x.share.len() {
                fwriteln!(
                    f,
                    ">exact_clonotype{}.chain{}.VJ\n{}",
                    i,
                    j + 1,
                    strme(&x.share[j].seq)
                );
            }
        }
    }
    if ctl.gen_opt.exact.is_some() {
        let ex = &exact_clonotypes[ctl.gen_opt.exact.unwrap()];
        println!("\nEXACT CLONOTYPE {}", ctl.gen_opt.exact.unwrap());
        for i in 0..ex.share.len() {
            let vid = ex.share[i].v_ref_id;
            let jid = ex.share[i].j_ref_id;
            println!(
                "chain {} = {} + {} = {}",
                i + 1,
                refdata.name[vid],
                refdata.name[jid],
                ex.share[i].cdr3_aa
            );
        }
        for i in 0..ex.clones.len() {
            let x = &ex.clones[i][0];
            println!(
                "clone {} = {}.{}",
                i + 1,
                ctl.origin_info.dataset_id[x.dataset_index],
                x.barcode
            );
        }
        println!();
    }
    ctl.perf_stats(&t, "finding exact subclonotypes three");
    exact_clonotypes
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Search for SHM indels.  Exploratory.

pub fn search_for_shm_indels(ctl: &EncloneControl, tig_bc: &Vec<Vec<TigData>>) {
    if ctl.gen_opt.indels {
        println!("CDR3s associated with possible SHM indels");
        let mut cs = Vec::<((String, usize), usize, String)>::new();
        for i in 0..tig_bc.len() {
            for j in 0..tig_bc[i].len() {
                let x = &tig_bc[i][j];
                cs.push((
                    (x.cdr3_dna.clone(), x.v_ref_id),
                    x.seq().len(),
                    x.cdr3_aa.clone(),
                ));
            }
        }
        unique_sort(&mut cs);
        let mut i = 0;
        while i < cs.len() {
            let j = next_diff1_3(&cs, i as i32) as usize;
            if j - i > 1 {
                println!("{}", cs[i].2);
            }
            i = j;
        }
        println!();
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Look for barcode reuse.  The primary purpose of this is to detect instances where two
// datasets were obtained from the same cDNA (from the same GEM well).

pub fn check_for_barcode_reuse(
    ctl: &EncloneControl,
    tig_bc: &Vec<Vec<TigData>>,
) -> Result<(), String> {
    if !ctl.gen_opt.accept_reuse {
        const MIN_REUSE_FRAC_TO_SHOW: f64 = 0.25;
        let mut all = Vec::<(String, usize, usize)>::new();
        let mut total = vec![0; ctl.origin_info.dataset_id.len()];
        for i in 0..tig_bc.len() {
            all.push((tig_bc[i][0].barcode.clone(), tig_bc[i][0].dataset_index, i));
            total[tig_bc[i][0].dataset_index] += 1;
        }
        all.par_sort();
        let mut reuse = Vec::<(usize, usize)>::new();
        let mut i = 0;
        while i < all.len() {
            let j = next_diff1_3(&all, i as i32) as usize;
            for k1 in i..j {
                for k2 in k1 + 1..j {
                    // We require identity on one cdr3_aa.  That seems to be about the right amount
                    // of concurrence that should be required.  If two datasets arise from the same
                    // barcode in the same cDNA (from the same GEM well), there can still be
                    // assembly differences.

                    let mut ok = false;
                    let (u1, u2) = (all[k1].2, all[k2].2);
                    for v1 in 0..tig_bc[u1].len() {
                        for v2 in 0..tig_bc[u2].len() {
                            if tig_bc[u1][v1].cdr3_aa == tig_bc[u2][v2].cdr3_aa {
                                ok = true;
                            }
                        }
                    }
                    if ok {
                        reuse.push((all[k1].1, all[k2].1));
                    }
                }
            }
            i = j;
        }
        reuse.sort_unstable();
        let mut i = 0;
        let mut found = false;
        let mut msg = String::new();
        while i < reuse.len() {
            let j = next_diff(&reuse, i);
            let n = j - i;
            let (l1, l2) = (reuse[i].0, reuse[i].1);
            let (n1, n2) = (total[l1], total[l2]);
            let frac = n as f64 / min(n1, n2) as f64;
            if frac >= MIN_REUSE_FRAC_TO_SHOW {
                if !found {
                    found = true;
                    msg += &mut "\nSignificant barcode reuse detected.  If at least 25% of the barcodes \
                        in one dataset\nare present in another dataset, is is likely that two datasets \
                        arising from the\nsame library were included as input to enclone.  Since this \
                        would normally occur\nonly by accident, enclone exits.  \
                        If you wish to override this behavior,\nplease rerun with the argument \
                        ACCEPT_REUSE.\n\nHere are the instances of reuse that were observed:\n\n".to_string();
                }
                msg += &mut format!(
                    "{}, {} ==> {} of {}, {} barcodes ({:.1}%)\n",
                    ctl.origin_info.dataset_id[l1],
                    ctl.origin_info.dataset_id[l2],
                    n,
                    n1,
                    n2,
                    100.0 * frac
                );
            }
            i = j;
        }
        if found {
            return Err(msg);
        }
    }
    Ok(())
}
