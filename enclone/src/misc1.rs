// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Miscellaneous functions.

use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype, TigData};
use equiv::EquivRel;
use itertools::Itertools;
#[cfg(not(target_os = "windows"))]
use pager::Pager;
use perf_stats::elapsed;
use std::collections::HashMap;
use std::time::Instant;
use string_utils::stringme;
use vector_utils::{
    bin_member, bin_position, erase_if, next_diff, next_diff1_3, unique_sort, VecUtils,
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This section contains a function that supports paging.  It does not work under Windows, and
// we describe here all the *known* problems with getting enclone to work under Windows.
// 1. It does not compile for us.  When we tried, there was a problem with libhdf-5.
// 2. Paging is turned off, because the pager crate doesn't compile under Windows, and porting
//    it to Windows appears nontrivial.
// 3. ANSI escape characters are not handled correctly, at least by default.
// In addition, we have some concerns about what it would mean to properly test enclone on Windows,
// given that some users might have older OS installs, and support for ANSI escape characters
// appears to have been changed in 2018. This is not made easier by the Windows Subsystem for
// Linux.

#[cfg(not(target_os = "windows"))]
pub fn setup_pager(pager: bool) {
    // If the output is going to a terminal, set up paging so that output is in effect piped to
    // "less -R -F -X -K".
    //
    // ∙ The option -R is used to render ANSI escape characters correctly.  We do not use
    //   -r instead because if you navigate backwards in less -r, stuff gets screwed up,
    //   which is consistent with the scary stuff in the man page for less at -r.  However -R will
    //   not display all unicode characters correctly, so those have to be picked carefully,
    //   by empirically testing that e.g. "echo ◼ | less -R -F -X" renders correctly.
    //
    // ∙ The -F option makes less exit immediately if all the output can be seen in one screen.
    //
    // ∙ The -X option is needed because we found that in full screen mode on OSX Catalina, output
    //   was sent to the alternate screen, and hence it appeared that one got no output at all
    //   from enclone.  This is really bad, so do not turn off this option!

    if pager {
        Pager::with_pager("less -R -F -X -K").setup();
    }
}

#[cfg(target_os = "windows")]
pub fn setup_pager(pager: bool) {}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Lookup for heavy chain reuse (special purpose experimental option).
// This is interesting but not likely to yield interesting examples of heavy chain reuse
// because biologically it doesn't make sense that one would have both H-L1 and H-L2 expanded.
// However, what the code does suggest is ways that heavy chain concordance might be used
// to improve the join algorithm.

pub fn lookup_heavy_chain_reuse(
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    info: &Vec<CloneInfo>,
    eq: &EquivRel,
) {
    if ctl.gen_opt.heavy_chain_reuse {
        let t = Instant::now();
        println!("\nheavy chain reuse by cdr3_aa:\n");
        let mut cdr3 = Vec::<(String, String, usize)>::new();
        for i in 0..info.len() {
            let ex = &exact_clonotypes[info[i].clonotype_id];

            // Assume exact subclonotype has structure H+L.
            // OMG do not assume that heavy comes before light.  It's not true!!!!!!!!!!!!!!!!!!!!!

            if ex.share.len() != 2 {
                continue;
            }
            let mut left = false;
            let mut right = false;
            let mut lid = 0;
            for j in 0..ex.share.len() {
                if ex.share[j].left {
                    left = true;
                    lid = j;
                } else {
                    right = true;
                }
            }
            if !left || !right {
                continue;
            }

            // Assume there are at least two cells in the exact subclonotype.

            if ex.clones.len() < 2 {
                continue;
            }

            // Currently set up to work on CDR3_AA but could switch to CDR3_DNA and just use
            // CDR3_AA for tracking.

            cdr3.push((
                ex.share[lid].cdr3_aa.clone(),
                ex.share[lid].cdr3_aa.clone(),
                eq.class_id(i as i32) as usize,
            ));
        }

        // For each of two positions, ignore the amino acid at that position by putting an
        // asterisk there.  Could make this number of positions configurable.

        let mut dio = Vec::<Vec<String>>::new();
        println!("cdr3 has size {}", cdr3.len());
        for z1 in 0..25 {
            for z2 in z1 + 1..25 {
                let mut xcdr3 = cdr3.clone();
                for i in 0..xcdr3.len() {
                    if xcdr3[i].0.len() > z2 {
                        let mut t = xcdr3[i].0.as_bytes().to_vec();
                        t[z1] = b'*';
                        t[z2] = b'*';
                        xcdr3[i].0 = stringme(&t);
                    }
                }

                // Look for heavy chain CDR3_AAs that are identical except for one position, and
                // lie in different orbits.

                xcdr3.sort();
                let mut i = 0;
                while i < xcdr3.len() {
                    let j = next_diff1_3(&xcdr3, i as i32) as usize;
                    let mut ids = Vec::<usize>::new();
                    for k in i..j {
                        ids.push(xcdr3[k].2);
                    }
                    unique_sort(&mut ids);
                    if !ids.solo() {
                        let mut x = Vec::<String>::new();
                        for k in i..j {
                            x.push(xcdr3[k].1.clone());
                        }
                        unique_sort(&mut x);
                        dio.push(x);
                    }
                    i = j;
                }
            }
        }
        unique_sort(&mut dio);
        for i in 0..dio.len() {
            println!("{} = {}", i + 1, dio[i].iter().format(", "));
        }
        println!(
            "\nused {:.2} seconds in heavy chain reuse calculation\n",
            elapsed(&t)
        );
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// If a V..J segment appears in exactly one dataset, with frequency n, let x be the total
// number of productive pairs for that dataset, and let y be the total number of productive
// pairs for all datasets from the same origin.  If (x/y)^n <= 10^-6, i.e. the probability
// that assuming even distribution, all instances of that V..J ended up in that one dataset,
// delete all the productive pairs for that V..J segment that do not have at least 100
// supporting UMIs.  (Note no attempt to do Bonferroni correction.)
//
// For the case of two datasets for one origin, with equal numbers of productive pairs in
// each, this corresponds roughly to the case n = 20.
//
// Note that we could modify this to allow *some* occurrences in other datasets.
//
// There are only certain ways that these misdistribution events could happen:
//
// 1. A cell (and particularly a plasma cell or plasmablast) bursts after drawing cells to
//    make libraries, leaving behind cell fragments that seed separate GEMs
//    (probably most likely).
// 2. Multiple gel beads end up in one GEM.
// 3. Something involving like cells sticking together and subsequently separating.
// 4. Physical contamination of libraries.
// 5. Informatic mixup of libraries.
// 6. Nothing other than a low probability event (unlikely).
//
// Note that in case 1, we have evidence that a plasma cell or plasmablast existed in the
// original cells that were drawn (perhaps breaking up in the process of drawing), and was
// subsequently distintegrated.

pub fn cross_filter(
    ctl: &EncloneControl,
    tig_bc: &mut Vec<Vec<TigData>>,
    fate: &mut Vec<HashMap<String, String>>,
) {
    // Get the list of dataset origins.  Here we allow the same origin name to have been used
    // for more than one donor, as we haven't explicitly prohibited that.

    let mut origins = Vec::<(String, String)>::new();
    for i in 0..ctl.origin_info.n() {
        origins.push((
            ctl.origin_info.donor_id[i].clone(),
            ctl.origin_info.origin_id[i].clone(),
        ));
    }
    unique_sort(&mut origins);
    let mut to_origin = vec![0; ctl.origin_info.n()];
    for i in 0..ctl.origin_info.n() {
        to_origin[i] = bin_position(
            &origins,
            &(
                ctl.origin_info.donor_id[i].clone(),
                ctl.origin_info.origin_id[i].clone(),
            ),
        ) as usize;
    }

    // For each dataset index, and each origin, compute the total number of productive pairs.

    let mut n_dataset_index = vec![0; ctl.origin_info.n()];
    let mut n_origin = vec![0; origins.len()];
    for i in 0..tig_bc.len() {
        for j in 0..tig_bc[i].len() {
            let x = &tig_bc[i][j];
            n_dataset_index[x.dataset_index] += 1;
            n_origin[to_origin[x.dataset_index]] += 1;
        }
    }

    // Find all the V..J segments, and for each, the number of times it appears in each
    // dataset ID.
    //
    // Note that there is no point running this unless we have at least two dataset IDs, and in
    // fact unless there is an origin with at least two dataset IDs.  Better: just gather data
    // for the origin for which there are at least two dataset IDs.  Also no point if NCROSS.

    let mut vjx = Vec::<(Vec<u8>, usize, usize)>::new(); // (V..J, dataset index, count)
    {
        for i in 0..tig_bc.len() {
            for j in 0..tig_bc[i].len() {
                let x = &tig_bc[i][j];
                vjx.push((x.seq().to_vec(), x.dataset_index, 1));
            }
        }
        vjx.sort();
        let mut to_delete = vec![false; vjx.len()];
        let mut i = 0;
        while i < vjx.len() {
            let j = next_diff(&vjx, i); // actually only need to check first two fields
            vjx[i].2 = j - i;
            for k in i + 1..j {
                to_delete[k] = true;
            }
            i = j;
        }
        erase_if(&mut vjx, &to_delete);
    }

    // Now do the cross filter.

    let mut blacklist = Vec::<Vec<u8>>::new();
    let mut i = 0;
    while i < vjx.len() {
        let j = next_diff1_3(&vjx, i as i32) as usize;
        if j - i == 1 {
            let dataset_index = vjx[i].1;
            let n = vjx[i].2;
            let x = n_dataset_index[dataset_index];
            let y = n_origin[to_origin[dataset_index]];
            if y > 0 {
                let p = (x as f64 / y as f64).powi(n as i32);
                if p <= 1.0e-6 {
                    blacklist.push(vjx[i].0.clone());
                }
            }
        }
        i = j;
    }
    blacklist.sort();
    let mut to_delete = vec![false; tig_bc.len()];
    const UMIS_SAVE: usize = 100;
    for i in 0..tig_bc.len() {
        for j in 0..tig_bc[i].len() {
            if tig_bc[i][j].umi_count < UMIS_SAVE
                && bin_member(&blacklist, &tig_bc[i][j].seq().to_vec())
            {
                fate[tig_bc[i][0].dataset_index].insert(
                    tig_bc[i][0].barcode.clone(),
                    "failed CROSS filter".to_string(),
                );
                if !ctl.clono_filt_opt_def.ncross {
                    to_delete[i] = true;
                }
                break;
            }
        }
    }
    erase_if(tig_bc, &to_delete);
}
