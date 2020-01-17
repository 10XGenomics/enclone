// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Miscellaneous functions.

use ansi_escape::*;
use crate::defs::*;
use equiv::*;
use itertools::*;
use libc;
use pager::Pager;
use perf_stats::*;
use std::time::Instant;
use string_utils::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn setup_pager() {

    // If output is going to a terminal, emit the ANSI escape character that disables the
    // alternate screen buffer.  We did this because we found that in full screen mode on
    // mac OSX Catalina, enclone would appear to produce no output, because the output was
    // being sent to the alternate screen.

    if unsafe { libc::isatty(libc::STDOUT_FILENO) } != 0 {
        let mut log = Vec::<u8>::new();
        emit_disable_alternate_screen_buffer_escape(&mut log);
        print!( "{}", strme(&log) );
    }

    // If the output is going to a terminal, set up paging so that output is in effect piped to
    // "less -r -F".  The option -r is used to render ANSI escape characters correctly and also
    // to properly display special unicode characters, including at least the red dot we show
    // when a clonotype contains cells from two donors.  We do not use the -R option because it
    // incorrectly handles this.  The -F option makes less exit immediately if all the output can
    // be seen in one screen.

    Pager::with_pager("less -r -F").setup();
}

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
        std::process::exit(0);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// If a V..J segment appears in exactly one lena id, with frequency n, let x be the total
// number of productive pairs for that lena, and let y be the total number of productive
// pairs for all lena ids from the same sample.  If (x/y)^n <= 10^-6, i.e. the probability
// that assuming even distribution, all instances of that V..J ended up in that one lena,
// delete all the productive pairs for that V..J segment that do not have at least 100
// supporting UMIs.  (Note no attempt to do bonferroni correction.)
//
// For the case of two lena ids for one sample, with equal numbers of productive pairs in
// each, this corresponds roughly to the case n = 20.
//
// Note that we could modify this to allow *some* occurrences in other lena ids.
//
// There are only certain ways that these misdistribution events could happen:
//
// 1. A cell (and particularly a plasma cell or plasmablast) bursts after drawing cells to
//    make libraries, leaving behind cell fragments that seed separate GEMs
// (probably most likely).
// 2. Multiple gel beads end up in one GEM.
// 3. Something involving like cells sticking together and subsequently separating.
// 4. Physical contamination of libraries.
// 5. Informatic mixup of libraries.
// 6. Nothing other than a low probability event (unlikely).
//
// Note that in case 1, we have evidence that a plasma cell or permblast existed in the
// original cells that were drawn (perhaps breaking up in the process of drawing), and was
// subsequently distintegrated.

pub fn cross_filter(ctl: &EncloneControl, mut tig_bc: &mut Vec<Vec<TigData>>) {
    if !ctl.clono_filt_opt.ncross {
        // Get the list of samples.  Here we allow the same sample name to have been used for
        // more than one donor, as we haven't explicitly prohibited that.

        let mut samples = Vec::<(String, String)>::new();
        for i in 0..ctl.sample_info.n() {
            samples.push((
                ctl.sample_info.donor_id[i].clone(),
                ctl.sample_info.sample_id[i].clone(),
            ));
        }
        unique_sort(&mut samples);
        let mut to_sample = vec![0; ctl.sample_info.n()];
        for i in 0..ctl.sample_info.n() {
            to_sample[i] = bin_position(
                &samples,
                &(
                    ctl.sample_info.donor_id[i].clone(),
                    ctl.sample_info.sample_id[i].clone(),
                ),
            ) as usize;
        }

        // For each lena index, and each sample, compute the total number of productive pairs.

        let mut n_lena_index = vec![0; ctl.sample_info.n()];
        let mut n_sample = vec![0; samples.len()];
        for i in 0..tig_bc.len() {
            for j in 0..tig_bc[i].len() {
                let x = &tig_bc[i][j];
                n_lena_index[x.lena_index] += 1;
                n_sample[to_sample[x.lena_index]] += 1;
            }
        }

        // Find all the V..J segments, and for each, the number of times it appears in each lena id.
        //
        // Note that there is no point running this unless we have at least two lena ids, and in
        // fact unless there is a sample with at least two lena ids.  Better: just gather data for
        // the sample for which there are at least two lena ids.  Also no point if NCROSS.

        let mut vjx = Vec::<(Vec<u8>, usize, usize)>::new(); // (V..J, lena index, count)
        {
            for i in 0..tig_bc.len() {
                for j in 0..tig_bc[i].len() {
                    let x = &tig_bc[i][j];
                    vjx.push((x.seq.clone(), x.lena_index, 1));
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
                let lena_index = vjx[i].1;
                let n = vjx[i].2;
                let x = n_lena_index[lena_index];
                let y = n_sample[to_sample[lena_index]];
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
                if tig_bc[i][j].umi_count < UMIS_SAVE && bin_member(&blacklist, &tig_bc[i][j].seq) {
                    to_delete[i] = true;
                }
            }
        }
        erase_if(&mut tig_bc, &to_delete);
    }
}
