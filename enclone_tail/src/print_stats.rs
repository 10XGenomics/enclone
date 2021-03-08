// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Print statistics.

use enclone_core::defs::*;
use io_utils::*;
use perf_stats::*;
use stats_utils::*;
use std::cmp::max;
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;
use string_utils::*;
use tables::*;
use vector_utils::*;

fn median(x: &[usize]) -> f64 {
    let h = x.len() / 2;
    if x.len() % 2 == 1 {
        x[h] as f64
    } else {
        (x[h - 1] + x[h]) as f64 / 2.0
    }
}

pub fn print_stats(
    tall: &Instant,
    exacts: &Vec<Vec<usize>>,
    rsi: &Vec<ColInfo>,
    exact_clonotypes: &Vec<ExactClonotype>,
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    fate: &Vec<HashMap<String, String>>,
    logx: &mut Vec<u8>,
    nclono2: &mut usize,
    two_chain: &mut usize,
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
    umish.sort();
    umisl.sort();
    umis.sort();
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
    let mut ncc = Vec::<(usize, usize)>::new();
    let mut sd = Vec::<(Option<usize>, Option<usize>)>::new();
    let mut merges = 0;
    let (mut numis, mut nreads) = (0, 0);
    let mut nreads_adjusted = 0.0;
    let mut numis2 = 0;
    let mut ncells2 = 0;
    for i in 0..nclono {
        if rsi[i].mat.len() == 2 {
            *two_chain += 1;
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
        if n >= 2 {
            *nclono2 += 1;
        }
        if n >= 1 {
            // not sure how n = 0 can happen but it does, maybe should trap this
            merges += n - 1;
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
                datasets.sort();
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
                            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                                if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                                    b = true;
                                }
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
                            if gex_info.cell_type[li].contains_key(&bc.clone()) {
                                if gex_info.cell_type[li][&bc.clone()].starts_with('B') {
                                    b = true;
                                }
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
                elapsed(&tall)
            );
            #[cfg(not(target_os = "macos"))]
            fwriteln!(logx, "   • peak memory = {:.1} GB", peak_mem_usage_gb());
        }

        // Print barcode fate.

        fwriteln!(logx, "2. barcode fate");
        let mut fates = Vec::<String>::new();
        for i in 0..fate.len() {
            for f in fate[i].iter() {
                if f.1.contains(" GEX ") && ctl.clono_filt_opt.ngex {
                    continue;
                }
                if f.1.contains(" CROSS ") && ctl.clono_filt_opt.ncross {
                    continue;
                }
                if f.1.contains(" UMI ") && !ctl.clono_filt_opt.umi_filt {
                    continue;
                }
                if f.1.contains(" UMI_RATIO ") && !ctl.clono_filt_opt.umi_ratio_filt {
                    continue;
                }
                if f.1.contains(" GRAPH_FILTER ") && ctl.gen_opt.ngraph_filter {
                    continue;
                }
                if f.1.contains(" QUAL") && !ctl.clono_filt_opt.qual_filter {
                    continue;
                }
                if f.1.contains(" WEAK_CHAINS ") && !ctl.clono_filt_opt.weak_chains {
                    continue;
                }
                if f.1.contains(" FOURSIE_KILL ") && !ctl.clono_filt_opt.weak_foursies {
                    continue;
                }
                if f.1.contains(" WHITEF ") && ctl.gen_opt.nwhitef {
                    continue;
                }
                if f.1.contains(" BC_DUP ") && !ctl.clono_filt_opt.bc_dup {
                    continue;
                }
                if f.1.contains(" IMPROPER ") && ctl.merge_all_impropers {
                    continue;
                }
                fates.push(f.1.clone());
            }
        }
        fates.sort();
        let mut freq = Vec::<(u32, String)>::new();
        make_freq(&fates, &mut freq);
        let mut rows = Vec::<Vec<String>>::new();
        rows.push(vec!["barcodes".to_string(), "why deleted".to_string()]);
        rows.push(vec!["\\hline".to_string(); 2]);
        for i in 0..freq.len() {
            rows.push(vec![format!("{}", freq[i].0), freq[i].1.clone()]);
        }
        rows.push(vec![format!("{}", fates.len()), "total".to_string()]);
        let mut log = String::new();
        print_tabular_vbox(&mut log, &rows, 2, &b"r|l".to_vec(), false, false);
        log.truncate(log.len() - 1);
        log = log.replace("\n", "\n   ");
        fwrite!(logx, "   {}\n", log);

        // Print other stats.

        fwriteln!(logx, "3. for the selected clonotypes");

        // Print summary table for chains / clonotypes / cells.

        ncc.sort();
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
            "   • number of cell-cell merges = {}",
            add_commas(merges)
        );
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
                        if ex.clones[k][0].validated_umis.as_ref().unwrap().len() > 0
                            && ex.clones[k][1].validated_umis.as_ref().unwrap().len() > 0
                        {
                            have_both += 1;
                        }
                    }
                    for m in 0..ex.clones[k].len() {
                        if ex.clones[k][m].validated_umis.is_none() {
                            missing_valid = true;
                        } else {
                            if ex.share[m].left {
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
        }
        if !missing_valid && left_valids.len() > 0 && right_valids.len() > 0 {
            left_valids.sort();
            right_valids.sort();
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
                row.push(format!(
                    "{}",
                    ctl.origin_info.origin_list[sdx[i].0.unwrap()]
                ));
            } else {
                row.push("?".to_string());
            }
            if sdx[i].1.is_some() {
                row.push(format!("{}", ctl.origin_info.donor_list[sdx[i].1.unwrap()]));
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
            if cluster.len() > 0 {
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
                    let typex = &cell_type[bc.clone()];
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
                    cs.sort();
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
    }

    // Print summary csv stats.

    if ctl.gen_opt.summary_csv {
        println!("\nmiddle_mean_umis_heavy,middle_mean_umis_light,n_twothreesie");
        println!("{:.2},{:.2},{}", middle_mean_umish, middle_mean_umisl, n23);
    }
}
