// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Print statistics.

use crate::alluvial_fb::*;
use crate::fate::print_fate;
use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype, GexInfo};
use enclone_core::median::median;
use io_utils::{fwrite, fwriteln};
use perf_stats::elapsed;
#[cfg(not(target_os = "macos"))]
use perf_stats::peak_mem_usage_gb;
use stats_utils::percent_ratio;
use std::cmp::max;
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;
use string_utils::{add_commas, TextUtils};
use tables::print_tabular_vbox;
use vector_utils::*;

pub fn print_stats(
    tall: &Instant,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    groups: &Vec<Vec<(i32, String)>>,
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    fate: &Vec<HashMap<String, String>>,
    logx: &mut Vec<u8>,
    nclono2: &mut usize,
    two_chain: &mut usize,
    three_chain: &mut usize,
    four_chain: &mut usize,
    opt_d_val: &Vec<(usize, Vec<Vec<Vec<usize>>>)>,
) {
    // Compute some umi stats.

    let nclono = exacts.len();
    let mut umish = Vec::<usize>::new();
    let mut umisl = Vec::<usize>::new();
    let mut umis = Vec::<usize>::new();
    for i in 0..nclono {
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            for k in 0..ex.clones.len() {
                let mut nu = 0;
                for l in 0..ex.share.len() {
                    if ex.share[l].left {
                        umish.push(ex.clones[k][l].umi_count);
                    } else {
                        umisl.push(ex.clones[k][l].umi_count);
                    }
                    nu += ex.clones[k][l].umi_count;
                }
                if ex.share.len() == 2 {
                    umis.push(nu);
                }
            }
        }
    }
    umish.sort_unstable();
    umisl.sort_unstable();
    umis.sort_unstable();
    let (mut middleh, mut denomh) = (0, 0);
    for j in umish.len() / 3..(2 * umish.len()) / 3 {
        middleh += umish[j];
        denomh += 1;
    }
    let mut middle_mean_umish = 0.0;
    if denomh > 0 {
        middle_mean_umish = (middleh as f64) / (denomh as f64);
    }
    let (mut middlel, mut denoml) = (0, 0);
    for j in umisl.len() / 3..(2 * umisl.len()) / 3 {
        middlel += umisl[j];
        denoml += 1;
    }
    let mut middle_mean_umisl = 0.0;
    if denoml > 0 {
        middle_mean_umisl = (middlel as f64) / (denoml as f64);
    }
    let (mut middle, mut denom) = (0, 0);
    for j in umis.len() / 3..(2 * umis.len()) / 3 {
        middle += umis[j];
        denom += 1;
    }
    let mut middle_mean_umis = 0.0;
    if denom > 0 {
        middle_mean_umis = (middle as f64) / (denom as f64);
    }

    // Compute n1 and n2 and n23 and n4.

    let mut n1 = 0;
    let mut n2 = 0;
    let mut n23 = 0;
    let mut n4 = 0;
    for i in 0..nclono {
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            if ex.nchains() == 1 {
                n1 += ex.ncells();
            }
            if ex.nchains() == 2 {
                n2 += ex.ncells();
            }
            if ex.nchains() == 2 || ex.nchains() == 3 {
                n23 += ex.ncells();
            }
            if ex.nchains() == 4 {
                n4 += ex.ncells();
            }
        }
    }

    // Print summary stats.

    let mut ncells = 0;
    *nclono2 = 0;
    *two_chain = 0;
    *three_chain = 0;
    *four_chain = 0;
    let mut ncc = Vec::<(usize, usize)>::new();
    let mut sd = Vec::<(Option<usize>, Option<usize>)>::new();
    let mut merges = 0;
    let mut merges2 = 0;
    let (mut numis, mut nreads) = (0, 0);
    let mut nreads_adjusted = 0.0;
    let mut numis2 = 0;
    let mut ncells2 = 0;
    let mut cells_by_donor = vec![0; ctl.origin_info.donor_list.len()];
    let mut mixes = 0;
    for i in 0..nclono {
        let mut cells_by_donor_this = vec![0; ctl.origin_info.donor_list.len()];
        if rsi[i].mat.len() == 2 {
            *two_chain += 1;
        } else if rsi[i].mat.len() == 3 {
            *three_chain += 1;
        } else if rsi[i].mat.len() == 4 {
            *four_chain += 1;
        }
        let mut n = 0;
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            if ex.share.len() == 2 {
                ncells2 += ex.ncells();
            }
            n += ex.ncells();
            for k in 0..ex.clones.len() {
                let x = &ex.clones[k][0];
                sd.push((x.origin_index, x.donor_index));
                if x.donor_index.is_some() {
                    cells_by_donor[x.donor_index.unwrap()] += 1;
                    cells_by_donor_this[x.donor_index.unwrap()] += 1;
                }
                for m in 0..ex.clones[k].len() {
                    numis += ex.clones[k][m].umi_count;
                    if ex.share.len() == 2 {
                        numis2 += ex.clones[k][m].umi_count;
                    }
                    let n = ex.clones[k][m].read_count;
                    nreads += n;
                    let mut x = n as f64;
                    if ex.clones[k][0].frac_reads_used.is_some() {
                        x /= ex.clones[k][0].frac_reads_used.unwrap() as f64 / 1_000_000.0;
                    }
                    nreads_adjusted += x;
                }
            }
        }
        if ctl.origin_info.donor_list.len() > 1 && ctl.clono_filt_opt_def.donor {
            for j1 in 0..exacts[i].len() {
                let ex1 = &exact_clonotypes[exacts[i][j1]];
                for j2 in j1..exacts[i].len() {
                    let ex2 = &exact_clonotypes[exacts[i][j2]];
                    for k1 in 0..ex1.clones.len() {
                        let x1 = &ex1.clones[k1][0];
                        for k2 in 0..ex2.clones.len() {
                            if (j1, k1) < (j2, k2) {
                                let x2 = &ex2.clones[k2][0];
                                if x1.donor_index.is_some() && x2.donor_index.is_some() {
                                    if x1.donor_index.unwrap() != x2.donor_index.unwrap() {
                                        mixes += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if n >= 2 {
            *nclono2 += 1;
        }
        for j in 0..cells_by_donor_this.len() {
            let n = cells_by_donor_this[j];
            if n > 1 {
                merges += n - 1;
                merges2 += (n * (n - 1)) / 2;
            }
        }
        ncells += n;
        ncc.push((rsi[i].mat.len(), n));
    }
    sd.sort();
    let mut sdx = Vec::<(Option<usize>, Option<usize>, usize)>::new();
    let mut i = 0;
    while i < sd.len() {
        let j = next_diff(&sd, i);
        sdx.push((sd[i].0, sd[i].1, j - i));
        i = j;
    }
    if ctl.gen_opt.summary {
        fwriteln!(logx, "\nSUMMARY STATISTICS");
        fwriteln!(logx, "1. overall");
        fwriteln!(logx, "   • number of datasets = {}", ctl.origin_info.n());
        fwriteln!(logx, "   • number of donors = {}", ctl.origin_info.donors);
        let mut vcells = 0;
        for i in 0..vdj_cells.len() {
            vcells += vdj_cells[i].len();
        }
        fwriteln!(logx, "   • original number of cells = {}", vcells);

        // Print mean reads per cell if known.

        let mut known = true;
        for i in 0..ctl.origin_info.n() {
            if ctl.origin_info.cells_cellranger[i].is_none() {
                known = false;
            } else if ctl.origin_info.mean_read_pairs_per_cell_cellranger[i].is_none() {
                known = false;
            }
        }
        let (mut cells, mut read_pairs) = (0, 0);
        if known {
            for i in 0..ctl.origin_info.n() {
                let c = ctl.origin_info.cells_cellranger[i].unwrap();
                let rpc = ctl.origin_info.mean_read_pairs_per_cell_cellranger[i].unwrap();
                cells += c;
                read_pairs += cells * rpc;
            }
            let rpc = ((read_pairs as f64) / (cells as f64)).round();
            fwriteln!(logx, "   • read pairs per cell = {}", rpc);
        }

        // Compute marking stats.

        let (mut nmarked, mut nmarked_good, mut ndubious) = (0, 0, 0);
        let (mut nfake, mut nfake_marked, mut ngood, mut ngood_marked) = (0, 0, 0, 0);
        if ctl.gen_opt.mark_stats || ctl.gen_opt.mark_stats2 {
            for i in 0..nclono {
                let mut datasets = Vec::<usize>::new();
                let mut ncells = 0;
                for j in 0..exacts[i].len() {
                    let ex = &exact_clonotypes[exacts[i][j]];
                    ncells += ex.ncells();
                    for l in 0..ex.ncells() {
                        datasets.push(ex.clones[l][0].dataset_index);
                    }
                }
                datasets.sort_unstable();
                let mut freq = Vec::<(u32, usize)>::new();
                make_freq(&datasets, &mut freq);
                let mut fake = false;
                let mut di = -1;
                if freq.len() == 1 || freq[0].0 >= 10 * freq[1].0 {
                    di = freq[0].1 as isize;
                }
                if ncells >= 10 && di >= 0 {
                    nfake += ncells - 1;
                    fake = true;
                }
                if ncells >= 2 {
                    for j in 0..exacts[i].len() {
                        let ex = &exact_clonotypes[exacts[i][j]];

                        // Determine if cell is called a B cell.

                        for l in 0..ex.ncells() {
                            let mut b = false;
                            let li = ex.clones[l][0].dataset_index;
                            let bc = &ex.clones[l][0].barcode;
                            if gex_info.cell_type[li].contains_key(&bc.clone())
                                && gex_info.cell_type[li][&bc.clone()].starts_with('B')
                            {
                                b = true;
                            }

                            // Record accordingly.

                            if ex.clones[l][0].dataset_index as isize == di || !b {
                                ndubious += 1;
                            }
                            if !fake
                                && ncells >= 10
                                && b
                                && (ex.share.len() == 2 || ex.share.len() == 3)
                            {
                                ngood += 1;
                                if ex.clones[l][0].marked {
                                    ngood_marked += 1;
                                }
                            }
                        }
                    }
                }
                if fake {
                    ngood += 1;
                }
                let mut fake_marks = 0;
                for j in 0..exacts[i].len() {
                    let ex = &exact_clonotypes[exacts[i][j]];
                    for l in 0..ex.ncells() {
                        if ex.clones[l][0].marked {
                            nmarked += 1;
                            if fake {
                                nfake_marked += 1;
                                fake_marks += 1;
                            }
                            let chains_ok = ex.nchains() >= 2 && ex.nchains() <= 3;
                            let mut b = false;
                            let li = ex.clones[l][0].dataset_index;
                            let bc = &ex.clones[l][0].barcode;
                            if gex_info.cell_type[li].contains_key(&bc.clone())
                                && gex_info.cell_type[li][&bc.clone()].starts_with('B')
                            {
                                b = true;
                            }
                            if chains_ok && freq.len() >= 2 && b {
                                nmarked_good += 1;
                            }
                        }
                    }
                }
                if fake_marks == ncells {
                    nfake_marked -= 1;
                }
            }
        }

        // Print computational performance stats.

        if !ctl.gen_opt.summary_clean {
            fwriteln!(
                logx,
                "   • total elapsed time = {:.1} seconds",
                elapsed(tall)
            );
            #[cfg(not(target_os = "macos"))]
            fwriteln!(logx, "   • peak memory = {:.1} GB", peak_mem_usage_gb());
        }

        // Print barcode fate.

        print_fate(&ctl, &fate, logx);

        // Print other stats.

        fwriteln!(logx, "3. for the selected clonotypes");

        // Print summary table for chains / clonotypes / cells.

        ncc.sort_unstable();
        let mut i = 0;
        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "chains".to_string(),
            "clonotypes with this".to_string(),
            "cells in these".to_string(),
            "%".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "".to_string(),
            "number of chains".to_string(),
            "clonotypes".to_string(),
            "".to_string(),
        ];
        rows.push(row);
        let row = vec!["\\hline".to_string(); 4];
        rows.push(row);
        while i < ncc.len() {
            let j = next_diff1_2(&ncc, i as i32) as usize;
            let nchains_this = ncc[i].0;
            let nclono_this = j - i;
            let mut ncells_this = 0;
            for k in i..j {
                ncells_this += ncc[k].1;
            }
            let row = vec![
                format!("{}", nchains_this),
                format!("{}", nclono_this),
                format!("{}", ncells_this),
                format!("{:.1}", percent_ratio(ncells_this, ncells)),
            ];
            rows.push(row);
            i = j;
        }
        let row = vec![
            "total".to_string(),
            format!("{}", nclono),
            format!("{}", ncells),
            "100.0".to_string(),
        ];
        rows.push(row);
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"l|r|r|r".to_vec(), false, false);
        log = log.replace("\n", "\n   ");
        fwrite!(logx, "   {}", log);

        // Print other cell/clonotype stats.

        fwriteln!(
            logx,
            "• number of clonotypes having at least two cells = {}",
            nclono2
        );
        fwriteln!(
            logx,
            "   • number of intradonor cell-cell merges = {}",
            add_commas(merges)
        );
        fwriteln!(
            logx,
            "   • number of intradonor cell-cell merges (quadratic) = {}",
            add_commas(merges2)
        );
        if cells_by_donor.len() > 1 && ctl.clono_filt_opt_def.donor {
            let mut cross = 0;
            let mut intra = 0;
            for i1 in 0..cells_by_donor.len() {
                if cells_by_donor[i1] > 1 {
                    intra += cells_by_donor[i1] * (cells_by_donor[i1] - 1) / 2;
                }
                for i2 in i1 + 1..cells_by_donor.len() {
                    cross += cells_by_donor[i1] * cells_by_donor[i2];
                }
            }
            fwriteln!(
                logx,
                "   • number of intradonor comparisons = {}",
                add_commas(intra)
            );
            fwriteln!(
                logx,
                "   • number of cross-donor comparisons = {}",
                add_commas(cross)
            );
            fwriteln!(
                logx,
                "   • number of cross-donor comparisons that mix donors = {}",
                add_commas(mixes)
            );
            let rate = (mixes as f64) * 1_000_000_000.0 / (cross as f64);
            fwriteln!(
                logx,
                "   • rate of cross donor mixing = {:.2} x 10^-9",
                rate
            );
            let bogus = (intra as f64) * (mixes as f64) / (cross as f64);
            let bogus = bogus.round() as usize;
            fwriteln!(
                logx,
                "   • estimated number of false intradonor merges = {}",
                add_commas(bogus)
            );

            let adjusted = if bogus <= merges2 { merges2 - bogus } else { 0 };
            fwriteln!(
                logx,
                "   • adjusted cell-cell merges (quadratic) = {}",
                add_commas(adjusted)
            );
        }
        fwriteln!(logx, "   • number of cells having 1 chain = {}", n1);
        fwriteln!(logx, "   • number of cells having 2 or 3 chains = {}", n23);
        let mut doublet_rate = 0.0;
        if n2 > 0 || n4 > 0 {
            doublet_rate = n4 as f64 / (n2 + n4) as f64;
        }
        let celltype;
        if ctl.gen_opt.bcr {
            celltype = "B";
        } else {
            celltype = "T";
        }
        fwrite!(
            logx,
            "   • estimated {}-{} doublet rate = {:.1}% = {}/{}",
            celltype,
            celltype,
            100.0 * doublet_rate,
            n4,
            n2 + n4
        );
        fwriteln!(logx, " = cells with 4 chains / cells with 2 or 4 chains");

        // Print UMI stats.

        let hchain;
        let lchain;
        if ctl.gen_opt.bcr {
            hchain = "heavy chain";
            lchain = "light chain";
        } else {
            hchain = "TRB";
            lchain = "TRA";
        }
        fwriteln!(
            logx,
            "   • mean over middle third of contig UMI counts ({}) = {:.2}",
            hchain,
            middle_mean_umish,
        );
        fwriteln!(
            logx,
            "   • mean over middle third of contig UMI counts ({}) = {:.2}",
            lchain,
            middle_mean_umisl,
        );
        fwriteln!(
            logx,
            "   • mean over middle third of cell UMI counts for cells having two chains = {:.2}",
            middle_mean_umis,
        );
        fwriteln!(
            logx,
            "   • mean UMIs per cell = {:.2}",
            numis as f64 / ncells as f64,
        );
        fwriteln!(
            logx,
            "   • mean UMIs per cell having two chains = {:.2}",
            numis2 as f64 / ncells2 as f64,
        );
        fwriteln!(
            logx,
            "   • for reads contributing to UMIs in reported chains, mean reads per UMI = {:.2}",
            nreads as f64 / numis as f64,
        );
        if known && ctl.gen_opt.internal_run {
            if ctl.gen_opt.no_uncap_sim {
                nreads_adjusted = nreads as f64;
            }
            fwriteln!(
                logx,
                "   • read utilization = {:.1}%\n     (please see notes in the file UNDOCUMENTED)",
                100.0 * nreads_adjusted / read_pairs as f64
            );
        }

        // Print validated UMI stats.

        let mut missing_valid = false;
        let mut left_valids = Vec::<usize>::new();
        let mut right_valids = Vec::<usize>::new();
        let mut lefts = 0;
        let mut rights = 0;
        let mut have_both = 0;
        let mut have_both_denom = 0;
        for i in 0..nclono {
            for j in 0..exacts[i].len() {
                let ex = &exact_clonotypes[exacts[i][j]];
                for k in 0..ex.clones.len() {
                    if ex.clones[k][0].validated_umis.is_some() && ex.clones[k].len() == 2 {
                        have_both_denom += 1;
                        if !ex.clones[k][0].validated_umis.as_ref().unwrap().is_empty()
                            && !ex.clones[k][1].validated_umis.as_ref().unwrap().is_empty()
                        {
                            have_both += 1;
                        }
                    }
                    for m in 0..ex.clones[k].len() {
                        if ex.clones[k][m].validated_umis.is_none() {
                            missing_valid = true;
                        } else if ex.share[m].left {
                            left_valids
                                .push(ex.clones[k][m].validated_umis.as_ref().unwrap().len());
                            lefts += ex.clones[k][m].umi_count;
                        } else {
                            right_valids
                                .push(ex.clones[k][m].validated_umis.as_ref().unwrap().len());
                            rights += ex.clones[k][m].umi_count;
                        }
                    }
                }
            }
        }
        if !missing_valid && !left_valids.is_empty() && !right_valids.is_empty() {
            left_valids.sort_unstable();
            right_valids.sort_unstable();
            let (mut left_sum, mut right_sum) = (0, 0);
            for i in 0..left_valids.len() {
                left_sum += left_valids[i];
            }
            for i in 0..right_valids.len() {
                right_sum += right_valids[i];
            }
            fwriteln!(
                logx,
                "   • median validated UMIs: {} = {} ({:.1}%), {} = {} ({:.1}%)",
                hchain,
                median(&left_valids),
                100.0 * left_sum as f64 / lefts as f64,
                lchain,
                median(&right_valids),
                100.0 * right_sum as f64 / rights as f64,
            );
            fwriteln!(
                logx,
                "   • fraction of two-chain cells having validated UMIs for both chains: {:.1}%",
                100.0 * have_both as f64 / have_both_denom as f64
            );
        }

        // Print marking stats.

        if ctl.gen_opt.mark_stats {
            fwriteln!(logx, "   --------------------------------");
            fwriteln!(logx, "   • number of dubious cells = {}", ndubious);
            fwriteln!(logx, "   • number of marked cells = {}", nmarked);
            fwriteln!(logx, "   • number of good marked cells = {}", nmarked_good);
        }
        if ctl.gen_opt.mark_stats2 {
            fwriteln!(logx, "   --------------------------------");
            fwriteln!(
                logx,
                "   • number of fake expanded clonotype cells = {}",
                nfake
            );
            fwriteln!(
                logx,
                "   • number of these that are marked = {}",
                nfake_marked
            );
            fwriteln!(logx, "   • residual = {}", nfake - nfake_marked);
            fwriteln!(
                logx,
                "   • number of good expanded clonotype cells = {}",
                ngood
            );
            fwriteln!(
                logx,
                "   • number of these that are marked = {}",
                ngood_marked
            );
        }

        // Print origin (sample)/donor table, but only if there is more than one.

        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "origin".to_string(),
            "donor".to_string(),
            "cells".to_string(),
        ];
        rows.push(row);
        let row = vec!["\\hline".to_string(); 3];
        rows.push(row);
        for i in 0..sdx.len() {
            let mut row = Vec::<String>::new();
            if sdx[i].0.is_some() {
                row.push(ctl.origin_info.origin_list[sdx[i].0.unwrap()].to_string());
            } else {
                row.push("?".to_string());
            }
            if sdx[i].1.is_some() {
                row.push(ctl.origin_info.donor_list[sdx[i].1.unwrap()].to_string());
            } else {
                row.push("?".to_string());
            }
            row.push(format!("{}", sdx[i].2));
            rows.push(row);
        }
        if rows.len() > 3 {
            let mut log = String::new();
            print_tabular_vbox(&mut log, &rows, 2, &b"llr".to_vec(), false, false);
            log = log.replace("\n", "\n   ");
            fwrite!(logx, "   {}", log);
        }

        // Print gene expression cluster analysis, but only if there is one library and the
        // cell types file was provided (which would be for PD cellranger runs with GEX data).

        if ctl.origin_info.n() == 1 && gex_info.cell_type_specified[0] {
            let cell_type = &gex_info.cell_type[0];
            let cluster = &gex_info.cluster[0];
            if !cluster.is_empty() {
                let mut nclust = 0;
                for x in cluster.iter() {
                    nclust = max(*x.1, nclust);
                }
                let mut ncells = vec![0; nclust];
                let mut ncells_type = vec![0; nclust];
                let mut ncells_vdj = vec![0; nclust];
                let mut ncells_vdj_max = vec![0; nclust];
                for x in cluster.iter() {
                    ncells[x.1 - 1] += 1;
                    let bc = &x.0;
                    let typex = &cell_type[*bc];
                    if ctl.gen_opt.bcr && typex.starts_with('B') {
                        ncells_type[x.1 - 1] += 1;
                    } else if ctl.gen_opt.tcr && typex.starts_with('T') {
                        ncells_type[x.1 - 1] += 1;
                    }
                }
                for i in 0..nclono {
                    let mut cs = Vec::<usize>::new();
                    for j in 0..exacts[i].len() {
                        let ex = &exact_clonotypes[exacts[i][j]];
                        for k in 0..ex.ncells() {
                            let bc = &ex.clones[k][0].barcode;
                            if cluster.contains_key(&bc.clone()) {
                                cs.push(cluster[&bc.clone()] - 1);
                            }
                        }
                    }
                    cs.sort_unstable();
                    let mut freq = Vec::<(u32, usize)>::new();
                    make_freq(&cs, &mut freq);
                    for x in freq.iter() {
                        ncells_vdj[x.1] += x.0;
                        ncells_vdj_max[x.1] = max(x.0, ncells_vdj_max[x.1]);
                    }
                }
                let mut rows = Vec::<Vec<String>>::new();
                let bt;
                if ctl.gen_opt.bcr {
                    bt = "B";
                } else {
                    bt = "T";
                }
                rows.push(vec![
                    "gex cluster".to_string(),
                    "cells".to_string(),
                    format!("%{}", bt),
                    "%in clonotypes".to_string(),
                    "max % in one clonotype".to_string(),
                ]);
                rows.push(vec!["\\hline".to_string(); 5]);
                for i in 0..nclust {
                    let n = ncells[i];
                    let mut row = Vec::<String>::new();
                    row.push(format!("{}", i + 1));
                    row.push(format!("{}", n));
                    row.push(format!("{:.1}", percent_ratio(ncells_type[i], n)));
                    row.push(format!("{:.1}", percent_ratio(ncells_vdj[i] as usize, n)));
                    row.push(format!(
                        "{:.1}",
                        percent_ratio(ncells_vdj_max[i] as usize, n)
                    ));
                    rows.push(row);
                }
                let mut log = String::new();
                print_tabular_vbox(&mut log, &rows, 2, &b"r|r|r|r|r".to_vec(), false, false);
                log = log.replace("\n", "\n   ");
                fwrite!(logx, "   {}", log);
            }
        }

        // Print dataset-level variable values.

        if !ctl.gen_opt.dvars.is_empty() {
            fwriteln!(logx, "\nDATASET-LEVEL METRICS");
            let mut row = vec!["dataset".to_string()];
            for j in 0..ctl.gen_opt.dvars.len() {
                let var = ctl.gen_opt.dvars[j].clone();
                let mut display_var = var.clone();
                if var.contains(':') {
                    display_var = var.before(":").to_string();
                }
                row.push(display_var);
            }
            let mut rows = vec![row];
            for i in 0..ctl.origin_info.n() {
                let mut row = Vec::<String>::new();
                let dataset_name = &ctl.origin_info.dataset_id[i];
                row.push(dataset_name.clone());
                for j in 0..ctl.gen_opt.dvars.len() {
                    let mut var = ctl.gen_opt.dvars[j].clone();
                    if var.contains(':') {
                        var = var.after(":").to_string();
                    }
                    let mut value = String::new();
                    if gex_info.json_metrics[i].contains_key(&var.to_string()) {
                        value = format!("{:.2}", gex_info.json_metrics[i][&var.to_string()]);
                    }
                    if value.is_empty() {
                        let mut feature = String::new();
                        let mut typex = String::new();
                        let mut fail = false;
                        if var.ends_with("_cellular_r") {
                            feature = var.before("_cellular_r").to_string();
                            typex = "r".to_string();
                        } else if var.ends_with("_cellular_u") {
                            feature = var.before("_cellular_u").to_string();
                            typex = "u".to_string();
                        } else {
                            fail = true;
                        }
                        if fail {
                            value = "undefined".to_string();
                        } else if typex == "r" {
                            if !gex_info.feature_metrics[i]
                                .contains_key(&(feature.clone(), "num_reads".to_string()))
                                || !gex_info.feature_metrics[i]
                                    .contains_key(&(feature.clone(), "num_reads_cells".to_string()))
                            {
                                value = "undefined".to_string();
                            } else {
                                let num = gex_info.feature_metrics[i]
                                    [&(feature.clone(), "num_reads_cells".to_string())]
                                    .force_usize();
                                let den = gex_info.feature_metrics[i]
                                    [&(feature.clone(), "num_reads".to_string())]
                                    .force_usize();
                                if den == 0 {
                                    value = "0/0".to_string();
                                } else {
                                    value = format!("{:.1}", 100.0 * num as f64 / den as f64);
                                }
                            }
                        } else if !gex_info.feature_metrics[i]
                            .contains_key(&(feature.clone(), "num_umis".to_string()))
                        {
                            value = "undefined".to_string();
                        } else if !gex_info.feature_metrics[i]
                            .contains_key(&(feature.clone(), "num_umis_cells".to_string()))
                        {
                            value = "undefined".to_string();
                        } else {
                            let num = gex_info.feature_metrics[i]
                                [&(feature.clone(), "num_umis_cells".to_string())]
                                .force_usize();
                            let den = gex_info.feature_metrics[i]
                                [&(feature.clone(), "num_umis".to_string())]
                                .force_usize();
                            if den == 0 {
                                value = "0/0".to_string();
                            } else {
                                value = format!("{:.1}", 100.0 * num as f64 / den as f64);
                            }
                        }
                    }
                    row.push(value);
                }
                rows.push(vec!["\\hline".to_string(); row.len()]);
                rows.push(row);
            }
            let mut just = vec![b'l'];
            for _ in 0..ctl.gen_opt.dvars.len() {
                just.push(b'|');
                just.push(b'r');
            }
            let mut log = String::new();
            print_tabular_vbox(&mut log, &rows, 2, &just, false, false);
            fwrite!(logx, "{}", log);
        }

        // Print global variable values.

        if !ctl.gen_opt.gvars.is_empty() {
            fwriteln!(logx, "\nGLOBAL VARIABLES\n");
            let (mut total, mut bads) = (0, 0);
            let mut need_inc = false;
            for x in ctl.gen_opt.gvars.iter() {
                if x.starts_with("d_inconsistent_") {
                    need_inc = true;
                }
            }
            if need_inc {
                for i in 0..exacts.len() {
                    for col in 0..rsi[i].mat.len() {
                        for u1 in 0..exacts[i].len() {
                            let m1 = rsi[i].mat[col][u1];
                            if m1.is_some() {
                                let m1 = m1.unwrap();
                                let ex1 = &exact_clonotypes[exacts[i][u1]];
                                if ex1.share[m1].left {
                                    for u2 in (u1 + 1)..exacts[i].len() {
                                        let m2 = rsi[i].mat[col][u2];
                                        if m2.is_some() {
                                            total += 1;
                                            if opt_d_val[i].1[col][u1] != opt_d_val[i].1[col][u2] {
                                                bads += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for var in ctl.gen_opt.gvars.iter() {
                let mut val = String::new();
                if *var == "d_inconsistent_%" {
                    val = format!("{:.2}", 100.0 * bads as f64 / total as f64);
                } else if *var == "d_inconsistent_n" {
                    val = format!("{}", total);
                }
                fwriteln!(logx, "{} = {}", var, val);
            }
        }

        // Print group stats, in the symmetric grouping case.

        if ctl.clono_group_opt.style == "symmetric" {
            fwriteln!(logx, "\nsymmetric grouping statistics");
            let mut rows = Vec::<Vec<String>>::new();
            let row = vec!["group size".to_string(), "number of clonotypes".to_string()];
            rows.push(row);
            rows.push(vec!["\\hline".to_string(); 2]);
            let mut lens = Vec::<usize>::new();
            for i in 0..groups.len() {
                lens.push(groups[i].len());
            }
            lens.sort();
            let mut i = 0;
            while i < lens.len() {
                let j = next_diff(&lens, i);
                let row = vec![format!("{}", lens[i]), format!("{}", j - i)];
                rows.push(row);
                i = j;
            }
            let mut log = String::new();
            print_tabular_vbox(&mut log, &rows, 2, &b"r|r".to_vec(), false, false);
            logx.append(&mut log.as_bytes().to_vec());
        }
    }

    // Print summary csv stats.

    if ctl.gen_opt.summary_csv {
        println!("\nmiddle_mean_umis_heavy,middle_mean_umis_light,n_twothreesie");
        println!("{:.2},{:.2},{}", middle_mean_umish, middle_mean_umisl, n23);
    }

    // Make alluvial tables for feature barcode data.  We determine cellular using vdj_cells,
    // which is not the only way of doing it.

    if ctl.gen_opt.summary {
        description_table(ctl, logx);
        alluvial_fb_reads(ctl, gex_info, vdj_cells, logx);
        alluvial_fb(ctl, gex_info, vdj_cells, logx);
    }
}
