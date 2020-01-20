// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Miscellaneous functions.

use crate::defs::*;
use crate::misc3::*;
use debruijn::dna_string::*;
use io_utils::*;
use perf_stats::*;
use rayon::prelude::*;
use std::cmp::{max, min};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Filter out putative gel bead contamination.  We look for cases where inside a
// given exact subclonotype, the same first or last half barcode is reused, and one
// instance has at least 10-fold higher UMI count.  If the fraction of the "bad"
// clones is at least 20%, delete them.

pub fn filter_gelbead_contamination(ctl: &EncloneControl, mut clones: &mut Vec<Vec<TigData0>>) {
    if !ctl.gen_opt.nwhitef {
        const GB_UMI_MULT: usize = 10;
        const GB_MIN_FRAC: f64 = 0.2;
        let mut bch = vec![Vec::<(usize, String, usize, usize)>::new(); 2];
        for l in 0..clones.len() {
            let li = clones[l][0].lena_index;
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
                if !to_delete[t - r] && tig_bc[t][m].v_start >= pos + 1 {
                    let p = tig_bc[t][m].v_start - pos - 1;
                    calls.push((tig_bc[t][m].full_seq[p], tig_bc[t][m].full_quals[p]));
                }
            }
            if calls.is_empty() || calls.len() < last_calls / 10 {
                break;
            }
            last_calls = calls.len();
            calls.sort();
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
            calls.sort();
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
        let mut z = tig_bc[r][m].seq.clone();
        full.append(&mut z);
        full.append(&mut constx);

        // Note that here we are taking the first entry (r), sort of assuming
        // that all the entries are the same, which in principle they should be.

        share.push(TigData1 {
            cdr3_dna: tig_bc[r][m].cdr3_dna.clone(),
            seq: tig_bc[r][m].seq.clone(),
            seq_del: tig_bc[r][m].seq.clone(), // may get changed later
            seq_del_amino: tig_bc[r][m].seq.clone(), // may get changed later
            full_seq: full,
            v_start: utr.len(),
            v_stop: tig_bc[r][m].v_stop + utr.len() - tig_bc[r][m].v_start,
            v_stop_ref: tig_bc[r][m].v_stop_ref,
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
            cdr3_aa: tig_bc[r][m].cdr3_aa.clone(),
            cdr3_start: tig_bc[r][m].cdr3_start,
            left: tig_bc[r][m].left,
            chain_type: tig_bc[r][m].chain_type.clone(),
            annv: tig_bc[r][m].annv.clone(),
            // these get set when making CloneInfo objects:
            vs: DnaString::new(),
            vs_notesx: String::new(),
            js: DnaString::new(),
        });
    }
    for t in r..s {
        if !to_delete[t - r] {
            let mut x = Vec::<TigData0>::new();
            for m in 0..tig_bc[t].len() {
                x.push(TigData0 {
                    quals: tig_bc[t][m].quals.clone(),
                    v_start: tig_bc[t][m].v_start.clone(),
                    j_stop: tig_bc[t][m].j_stop.clone(),
                    full_seq: tig_bc[t][m].full_seq.clone(),
                    barcode: tig_bc[t][m].barcode.clone(),
                    tigname: tig_bc[t][m].tigname.clone(),
                    lena_index: tig_bc[t][m].lena_index,
                    umi_count: tig_bc[t][m].umi_count,
                    read_count: tig_bc[t][m].read_count,
                });
            }
            clones.push(x);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Find exact subclonotypes.  We use a strict definition of exact subclonotypes, not
// allowing for errors.

pub fn find_exact_subclonotypes(
    ctl: &EncloneControl,
    tig_bc: &Vec<Vec<TigData>>,
    refdata: &RefData,
) -> Vec<ExactClonotype> {
    let texact = Instant::now();
    let mut exact_clonotypes = Vec::<ExactClonotype>::new();
    let mut r = 0;
    let mut max_exact = 0;
    let mut _count = 0;
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
                    || tig_bc[s][m].seq != tig_bc[r][m].seq

                    // Working around a bug here.  See above for explanation.

                    // || cid1 != cid2 {

                    || ( cid1.is_none() && cid2.is_some() ) || ( cid1.is_some() && cid2.is_none() )
                    || ( cid1.is_some() && cid2.is_some()
                        && refdata.name[cid1.unwrap()] != refdata.name[cid2.unwrap()] )
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

        // Print reused barcodes.

        if ctl.gen_opt.reuse {
            for t1 in r..s {
                for t2 in t1 + 1..s {
                    if tig_bc[t1][0].barcode == tig_bc[t2][0].barcode {
                        println!("see reuse of barcode {}", tig_bc[t1][0].barcode);
                        print!(
                            "{}: numis =",
                            ctl.sample_info.dataset_id[tig_bc[t1][0].lena_index]
                        );
                        for m in 0..tig_bc[t1].len() {
                            print!(" {}", tig_bc[t1][m].umi_count);
                        }
                        println!("");
                        print!(
                            "{}: numis =",
                            ctl.sample_info.dataset_id[tig_bc[t2][0].lena_index]
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
        // higher UMI counts.

        let mut to_delete = vec![false; s - r];
        for t1 in r..s {
            for t2 in t1 + 1..s {
                if tig_bc[t1][0].barcode == tig_bc[t2][0].barcode {
                    to_delete[t1 - r] = true;
                    to_delete[t2 - r] = true;
                }
            }
        }

        // Create the exact subclonotype.

        let mut share = Vec::<TigData1>::new();
        let mut clones = Vec::<Vec<TigData0>>::new();
        create_exact_subclonotype_core(&tig_bc, r, s, &to_delete, &mut share, &mut clones);

        // Explore consensus.

        study_consensus(
            &mut _count,
            &ctl,
            &share,
            &clones,
            &exact_clonotypes,
            &refdata,
        );

        // Filter out putative gel bead contamination.

        filter_gelbead_contamination(&ctl, &mut clones);

        // Save exact subclonotype.

        if clones.len() > 0 {
            exact_clonotypes.push(ExactClonotype {
                share: share,
                clones: clones,
            });
            max_exact = max(max_exact, s - r);
        }
        r = s;
    }
    if ctl.gen_opt.utr_con || ctl.gen_opt.con_con {
        println!("");
        std::process::exit(0);
    }
    if !ctl.silent {
        println!(
            "found {} exact subclonotypes from {} productive pairs",
            exact_clonotypes.len(),
            tig_bc.len()
        );
        println!("max exact subclonotype size = {}", max_exact);
    }
    if ctl.comp {
        println!(
            "used {:.1} seconds finding exact subclonotypes, peak mem = {:.2} GB",
            elapsed(&texact),
            peak_mem_usage_gb()
        );
    }
    if ctl.gen_opt.fasta.len() > 0 {
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
                ctl.sample_info.dataset_id[x.lena_index],
                x.barcode
            );
        }
        println!("");
    }
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
                    x.seq.len(),
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
        println!("");
        std::process::exit(0);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Look for barcode reuse.  The primary purpose of this is to detect instances where two
// datasets were obtained from the same cDNA (from the same GEM well).

pub fn check_for_barcode_reuse(ctl: &EncloneControl, tig_bc: &Vec<Vec<TigData>>) {
    if ctl.gen_opt.no_reuse {
        const MIN_REUSE_FRAC_TO_SHOW: f64 = 0.1;
        let mut all = Vec::<(String, usize, usize)>::new();
        let mut total = vec![0; ctl.sample_info.dataset_id.len()];
        for i in 0..tig_bc.len() {
            all.push((tig_bc[i][0].barcode.clone(), tig_bc[i][0].lena_index, i));
            total[tig_bc[i][0].lena_index] += 1;
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
        reuse.sort();
        let mut found = false;
        let mut i = 0;
        while i < reuse.len() {
            let j = next_diff(&reuse, i);
            let n = j - i;
            let (l1, l2) = (reuse[i].0, reuse[i].1);
            let (n1, n2) = (total[l1], total[l2]);
            let frac = n as f64 / min(n1, n2) as f64;
            if frac >= MIN_REUSE_FRAC_TO_SHOW {
                if !found {
                    eprintln!("\nBARCODE REUSE DETECTED:\n");
                    found = true;
                }
                eprintln!(
                    "{}, {} ==> {} of {}, {} ({:.1}%)",
                    ctl.sample_info.dataset_id[l1],
                    ctl.sample_info.dataset_id[l2],
                    n,
                    n1,
                    n2,
                    100.0 * frac
                );
            }
            i = j;
        }
        if found {
            eprintln!("");
            std::process::exit(1);
        }
    }
}
