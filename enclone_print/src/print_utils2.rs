// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function row_fill,
// plus a small helper function get_gex_matrix_entry.

use crate::print_utils1::*;
use crate::proc_cvar1::*;
use crate::proc_cvar2::*;
use crate::proc_lvar::*;
use amino::*;
use enclone_core::allowed_vars::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use itertools::*;
use ndarray::s;
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// The following code creates a row in the enclone output table for a clonotype.  Simultaneously
// it generates a row of parseable output.  And it does some other things that are not described
// here.
//
// TODO: Awful interface, should work to improve.

pub fn row_fill(
    pass: usize,
    u: usize,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    mults: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    gex_info: &GexInfo,
    refdata: &RefData,
    varmat: &Vec<Vec<Vec<u8>>>,
    fp: &Vec<Vec<usize>>,
    vars_amino: &Vec<Vec<usize>>,
    show_aa: &Vec<Vec<usize>>,
    field_types: &Vec<Vec<u8>>,
    bads: &mut Vec<bool>,
    gex_low: &mut usize,
    row: &mut Vec<String>,                       // row of human-readable output
    out_data: &mut Vec<HashMap<String, String>>, // row of parseable output
    cx: &mut Vec<Vec<String>>,
    d_all: &mut Vec<Vec<u32>>,
    ind_all: &mut Vec<Vec<u32>>,
    rsi: &ColInfo,
    dref: &Vec<DonorReferenceItem>,
    groups: &HashMap<usize, Vec<usize>>,
    d_readers: &Vec<Option<hdf5::Reader>>,
    ind_readers: &Vec<Option<hdf5::Reader>>,
    h5_data: &Vec<(usize, Vec<u32>, Vec<u32>)>,
    stats: &mut Vec<(String, Vec<f64>)>,
    vdj_cells: &Vec<Vec<String>>,
    n_vdj_gex: &Vec<usize>,
    lvarsc: &Vec<String>,
    nd_fields: &Vec<String>,
    peer_groups: &Vec<Vec<(usize, u8, u32)>>,
    extra_parseables: &Vec<String>,
) {
    // Redefine some things to reduce dependencies.

    let mat = &rsi.mat;
    let cvars = &ctl.clono_print_opt.cvars;
    let lvars = lvarsc.clone();
    let clonotype_id = exacts[u];
    let ex = &exact_clonotypes[clonotype_id];
    let mut pcols_sort = ctl.parseable_opt.pcols_sort.clone();
    for i in 0..pcols_sort.len() {
        pcols_sort[i] = pcols_sort[i].replace("_Σ", "_sum");
        pcols_sort[i] = pcols_sort[i].replace("_μ", "_mean");
    }
    pcols_sort.sort();
    let mut extra_args = ctl.gen_opt.tree.clone();
    if ctl.plot_opt.plot_xy_filename.len() > 0 {
        extra_args.push(ctl.plot_opt.plot_xy_xvar.clone());
        extra_args.push(ctl.plot_opt.plot_xy_yvar.clone());
    }
    unique_sort(&mut extra_args);
    macro_rules! speakc {
        ($u:expr, $col:expr, $var:expr, $val:expr) => {
            if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && $col + 1 <= ctl.parseable_opt.pchains
            {
                let mut v = $var.clone();
                v = v.replace("_Σ", "_sum");
                v = v.replace("_μ", "_mean");

                // Strip escape character sequences from val.  Can happen in notes, maybe
                // other places.

                let mut val_clean = String::new();
                let mut chars = Vec::<char>::new();
                let valx = format!("{}", $val);
                for c in valx.chars() {
                    chars.push(c);
                }
                let mut escaped = false;
                for l in 0..chars.len() {
                    if chars[l] == '' {
                        escaped = true;
                    }
                    if escaped {
                        if chars[l] == 'm' {
                            escaped = false;
                        }
                        continue;
                    }
                    val_clean.push(chars[l]);
                }

                // Proceed.

                let varc = format!("{}{}", v, $col + 1);
                if pcols_sort.is_empty() || bin_member(&pcols_sort, &varc) {
                    out_data[$u].insert(varc, val_clean);
                }
            }
        };
    }
    let cols = varmat[0].len();
    if ctl.gen_opt.row_fill_verbose {
        eprintln!("");
    }

    // Compute dataset indices, gex, gex_min, gex_max, gex_mean, gex_sum,
    // n_gex_cell, n_gex, entropy.

    let mut dataset_indices = Vec::<usize>::new();
    for l in 0..ex.clones.len() {
        dataset_indices.push(ex.clones[l][0].dataset_index);
    }
    unique_sort(&mut dataset_indices);
    let mut lenas = Vec::<String>::new();
    for l in dataset_indices.iter() {
        lenas.push(ctl.origin_info.dataset_id[*l].clone());
    }
    row.push("".to_string()); // row number (#), filled in below
    let mut counts = Vec::<usize>::new();
    let mut fcounts = Vec::<f64>::new();
    let mut n_gex = 0;
    let mut n_gexs = Vec::<usize>::new();
    let mut total_counts = Vec::<usize>::new();
    // It might be possible to speed this up a lot by pulling part of the "let d" and
    // "let ind" constructs out of the loop.
    if lvars.contains(&"entropy".to_string()) {
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = ex.clones[l][0].barcode.clone();
            if !gex_info.gex_barcodes.is_empty() {
                let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                if p >= 0 {
                    let mut raw_count = 0;
                    if gex_info.gex_matrices[li].initialized() {
                        let row = gex_info.gex_matrices[li].row(p as usize);
                        for j in 0..row.len() {
                            let f = row[j].0;
                            let n = row[j].1;
                            if gex_info.is_gex[li][f] {
                                raw_count += n;
                            }
                        }
                    } else {
                        let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                        let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize; // is p+1 OK??
                        let d: Vec<u32>;
                        let ind: Vec<u32>;
                        if ctl.gen_opt.h5_pre {
                            d = h5_data[li].1[z1..z2].to_vec();
                            ind = h5_data[li].2[z1..z2].to_vec();
                        } else {
                            d = d_readers[li]
                                .as_ref()
                                .unwrap()
                                .read_slice(s![z1..z2])
                                .unwrap()
                                .to_vec();
                            ind = ind_readers[li]
                                .as_ref()
                                .unwrap()
                                .read_slice(s![z1..z2])
                                .unwrap()
                                .to_vec();
                        }
                        for j in 0..d.len() {
                            if gex_info.is_gex[li][ind[j] as usize] {
                                raw_count += d[j] as usize;
                            }
                        }
                        d_all[l] = d;
                        ind_all[l] = ind;
                    }
                    total_counts.push(raw_count);
                }
            }
        }
    }
    let mut entropies = Vec::<f64>::new();
    for l in 0..ex.clones.len() {
        let li = ex.clones[l][0].dataset_index;
        let bc = ex.clones[l][0].barcode.clone();
        if !gex_info.gex_barcodes.is_empty() {
            if bin_member(&gex_info.gex_cell_barcodes[li], &bc) {
                n_gex += 1;
                n_gexs.push(1);
            } else {
                n_gexs.push(0);
            }
            let mut count = 0;
            let mut fcount = 0.0;
            let mut entropy = 0.0;
            let p = bin_position(&gex_info.gex_barcodes[li], &bc);
            if p >= 0 {
                let mut raw_count = 0;
                if gex_info.gex_matrices[li].initialized() {
                    let row = gex_info.gex_matrices[li].row(p as usize);
                    for j in 0..row.len() {
                        let f = row[j].0;
                        let n = row[j].1;
                        if gex_info.is_gex[li][f] {
                            if lvars.contains(&"entropy".to_string()) {
                                let q = n as f64 / total_counts[l] as f64;
                                entropy -= q * q.log2();
                            }
                            raw_count += n;
                        }
                    }
                } else {
                    let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                    let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize; // is p+1 OK??
                    let d: Vec<u32>;
                    let ind: Vec<u32>;
                    if ctl.gen_opt.h5_pre {
                        d = h5_data[li].1[z1..z2].to_vec();
                        ind = h5_data[li].2[z1..z2].to_vec();
                    } else {
                        d = d_readers[li]
                            .as_ref()
                            .unwrap()
                            .read_slice(s![z1..z2])
                            .unwrap()
                            .to_vec();
                        ind = ind_readers[li]
                            .as_ref()
                            .unwrap()
                            .read_slice(s![z1..z2])
                            .unwrap()
                            .to_vec();
                    }
                    for j in 0..d.len() {
                        if gex_info.is_gex[li][ind[j] as usize] {
                            let n = d[j] as usize;
                            if lvars.contains(&"entropy".to_string()) {
                                let q = n as f64 / total_counts[l] as f64;
                                entropy -= q * q.log2();
                            }
                            raw_count += n;
                        }
                    }
                    d_all[l] = d;
                    ind_all[l] = ind;
                }
                if !ctl.gen_opt.full_counts {
                    count = (raw_count as f64 * gex_info.gex_mults[li]).round() as usize;
                    fcount = raw_count as f64 * gex_info.gex_mults[li];
                } else {
                    count = (raw_count as f64).round() as usize;
                    fcount = raw_count as f64;
                }
            }
            counts.push(count);
            fcounts.push(fcount);
            entropies.push(entropy);
        }
    }
    let count_unsorted = counts.clone();
    counts.sort();
    for n in counts.iter() {
        if *n < 100 {
            *gex_low += 1;
        }
    }
    let (mut gex_median, mut gex_min, mut gex_max, mut gex_mean, mut gex_sum) = (0, 0, 0, 0.0, 0.0);
    if counts.len() > 0 {
        gex_median = rounded_median(&counts);
        gex_min = counts[0];
        gex_max = counts[counts.len() - 1];
        gex_sum = fcounts.iter().sum::<f64>();
        gex_mean = gex_sum / fcounts.len() as f64;
    }
    let entropies_unsorted = entropies.clone();
    entropies.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut entropy = 0.0;
    if entropies.len() > 0 {
        entropy = median_f64(&entropies);
    }

    // Output lead variable columns.
    // WARNING!  If you add lead variables, you may need to add them to the function
    // LinearCondition::require_valid_variables.

    let mut all_lvars = lvars.clone();
    if ctl.parseable_opt.pout.len() == 0 {
    } else if ctl.parseable_opt.pcols.is_empty() {
        for i in 0..LVARS_ALLOWED.len() {
            if !lvars.contains(&LVARS_ALLOWED[i].to_string()) {
                all_lvars.push(LVARS_ALLOWED[i].to_string());
            }
        }
    } else {
        for i in 0..ctl.parseable_opt.pcols.len() {
            if !lvars.contains(&ctl.parseable_opt.pcols[i].to_string()) {
                all_lvars.push(ctl.parseable_opt.pcols[i].to_string());
            }
        }
    }
    for x in extra_args.iter() {
        if !lvars.contains(&x) {
            all_lvars.push(x.clone());
        }
    }
    let mut alt_bcs = Vec::<String>::new();
    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.clone());
        }
    }
    unique_sort(&mut alt_bcs);
    for i in 0..all_lvars.len() {
        let x = &all_lvars[i];
        proc_lvar(
            i,
            &x,
            pass,
            u,
            &ctl,
            &exacts,
            &mults,
            &exact_clonotypes,
            &gex_info,
            &refdata,
            &varmat,
            &fp,
            row,
            out_data,
            d_all,
            ind_all,
            &rsi,
            &dref,
            &groups,
            stats,
            &vdj_cells,
            &n_vdj_gex,
            &nd_fields,
            &lvars,
            &lenas,
            &alt_bcs,
            n_gex,
            &n_gexs,
            gex_min,
            gex_max,
            gex_mean,
            gex_sum,
            gex_median,
            &count_unsorted,
            entropy,
            &entropies_unsorted,
            &fcounts,
        );
    }

    // Sanity check.  It's here because if it fails and that failure was not detected, something
    // exceptionally cryptic would happen downstream.

    if row.len() != lvars.len() + 1 {
        let mut msg = format!(
            "Oops, row.len() != lvars.len() + 1, as in fact we have\n\
            row.len() = {} and lvars.len() = {}, and in more detail,\n\
            row = {}\n\
            and lvars = {}.\nThis happened on a clonotype that included the barcode {}.",
            row.len(),
            lvars.len(),
            row.iter().format(","),
            lvars.iter().format(","),
            ex.clones[0][0].barcode,
        );
        if !ctl.gen_opt.row_fill_verbose {
            msg += &format!(
                "\n\nYou may find it helpful to add the options\n\
                BARCODE={} ROW_FILL_VERBOSE\n\
                to the command line.  Depending on other arguments, you might also need to \
                add MAX_CORES=1.",
                ex.clones[0][0].barcode
            );
        }
        panic!(msg);
    }

    // Get the relevant barcodes.

    let mut bli = Vec::<(String, usize, usize)>::new();
    for l in 0..ex.clones.len() {
        bli.push((
            ex.clones[l][0].barcode.clone(),
            ex.clones[l][0].dataset_index,
            l,
        ));
    }
    bli.sort();

    // Traverse the chains.

    for col in 0..cols {
        let mid = mat[col][u];
        if mid.is_none() {
            continue;
        }
        let mid = mid.unwrap();
        let ex = &exact_clonotypes[clonotype_id];
        let seq_amino = rsi.seqss_amino[col][u].clone();

        // Get UMI and read stats.

        let mut numis = Vec::<usize>::new();
        let mut nreads = Vec::<usize>::new();
        for j in 0..ex.clones.len() {
            numis.push(ex.clones[j][mid].umi_count);
            nreads.push(ex.clones[j][mid].read_count);
        }
        numis.sort();
        let median_numis = rounded_median(&numis);
        let utot: usize = numis.iter().sum();
        let u_mean = (utot as f64 / numis.len() as f64).round() as usize;
        let u_min = *numis.iter().min().unwrap();
        let u_max = *numis.iter().max().unwrap();
        nreads.sort();
        let rtot: usize = nreads.iter().sum();
        let r_mean = (rtot as f64 / nreads.len() as f64).round() as usize;
        let r_min = *nreads.iter().min().unwrap();
        let r_max = *nreads.iter().max().unwrap();
        let median_nreads = rounded_median(&nreads);

        // Speak quality score column entries.

        if ctl.parseable_opt.pout.len() > 0 && col + 1 <= ctl.parseable_opt.pchains {
            for i in 0..pcols_sort.len() {
                if pcols_sort[i].starts_with('q')
                    && pcols_sort[i].ends_with(&format!("_{}", col + 1))
                {
                    let n = pcols_sort[i].after("q").rev_before("_").force_usize();
                    if n < ex.share[mid].seq.len() {
                        let mut quals = Vec::<u8>::new();
                        for j in 0..ex.clones.len() {
                            quals.push(ex.clones[j][mid].quals[n]);
                        }
                        let q = format!("{}", quals.iter().format(","));
                        out_data[u].insert(pcols_sort[i].clone(), q);
                    }
                }
            }
        }

        // Speak some other column entries.

        let xm = &ex.share[mid];
        speakc!(u, col, "vj_aa".to_string(), stringme(&aa_seq(&xm.seq, 0)));
        speakc!(
            u,
            col,
            "vj_aa_nl".to_string(),
            stringme(&aa_seq(&xm.seq, xm.fr1_start))
        );
        speakc!(u, col, "vj_seq".to_string(), stringme(&xm.seq));
        speakc!(
            u,
            col,
            "vj_seq_nl".to_string(),
            stringme(&xm.seq[xm.fr1_start..])
        );
        speakc!(u, col, "seq".to_string(), stringme(&xm.full_seq));
        speakc!(u, col, "v_start".to_string(), xm.v_start);
        let (mut d_start, mut d_frame) = (String::new(), String::new());
        if xm.d_start.is_some() {
            d_start = format!("{}", xm.d_start.unwrap());
            d_frame = format!("{}", (xm.d_start.unwrap() - xm.v_start) % 3);
        }
        speakc!(u, col, "d_start".to_string(), d_start);
        speakc!(u, col, "d_frame".to_string(), d_frame);
        let cid = xm.c_ref_id;
        if cid.is_some() {
            let cid = cid.unwrap();
            speakc!(u, col, "const_id".to_string(), refdata.id[cid]);
        }
        let uid = ex.share[mid].u_ref_id;
        if uid.is_some() {
            let uid = uid.unwrap();
            speakc!(u, col, "utr_id".to_string(), refdata.id[uid]);
            speakc!(u, col, "utr_name".to_string(), refdata.name[uid]);
        }
        speakc!(u, col, "cdr3_start".to_string(), xm.cdr3_start);
        let mut vv = Vec::<usize>::new();
        for x in vars_amino[col].iter() {
            vv.push(*x / 3);
        }
        unique_sort(&mut vv);
        let mut varaa = Vec::<u8>::new();
        for p in vv.iter() {
            // what does it mean if this fails?
            if 3 * p + 3 <= seq_amino.len() {
                if seq_amino[3 * p..3 * p + 3].to_vec() == b"---".to_vec() {
                    varaa.push(b'-');
                } else {
                    varaa.push(codon_to_aa(&seq_amino[3 * p..3 * p + 3]));
                }
            }
        }
        speakc!(u, col, "var_aa".to_string(), strme(&varaa));

        // Define all_vars.

        let rsi_vars = &rsi.cvars[col];
        let mut all_vars = rsi_vars.clone();
        for j in 0..CVARS_ALLOWED.len() {
            let var = &CVARS_ALLOWED[j];
            if !rsi_vars.contains(&var.to_string()) {
                all_vars.push(var.to_string());
            }
        }
        for j in 0..CVARS_ALLOWED_PCELL.len() {
            let var = &CVARS_ALLOWED_PCELL[j];
            if !rsi_vars.contains(&var.to_string()) {
                all_vars.push(var.to_string());
            }
        }
        all_vars.append(&mut extra_parseables.clone());
        for x in extra_args.iter() {
            if !rsi_vars.contains(&x) {
                all_vars.push(x.clone());
            }
        }

        // Create column entry.

        let mut somelist = vec![false; all_vars.len()];
        for j in 0..all_vars.len() {
            let var = &all_vars[j];
            let varc = format!("{}{}", var, col + 1);
            if j < rsi.cvars[col].len() && cvars.contains(&var) {
                somelist[j] = true;
            } else if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && col + 1 <= ctl.parseable_opt.pchains
                && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc))
            {
                somelist[j] = true;
            }
        }
        for j in 0..all_vars.len() {
            // Decide if there is nothing to compute.  This is almost certainly not optimal.

            let mut needed = false;
            let var = &all_vars[j];
            let varc = format!("{}{}", var, col + 1);
            if j < rsi.cvars[col].len() && cvars.contains(&var) {
                needed = true;
            } else if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && col + 1 <= ctl.parseable_opt.pchains
                && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc))
            {
                needed = true;
            } else if *var == "amino".to_string() {
                needed = true;
            } else if *var == "u_cell".to_string() || *var == "r_cell".to_string() {
                needed = true;
            } else if *var == "white".to_string() || ctl.clono_filt_opt.whitef {
                needed = true;
            }
            if extra_args.contains(&varc) {
                needed = true;
            }
            if !needed {
                continue;
            }
            let col_var = j < rsi_vars.len();
            if !col_var && ctl.parseable_opt.pout.len() == 0 && extra_args.is_empty() {
                continue;
            }

            // Compute.

            if !proc_cvar1(
                &var,
                j,
                col,
                mid,
                pass,
                u,
                &ex,
                &ctl,
                &exacts,
                &exact_clonotypes,
                &refdata,
                &varmat,
                out_data,
                &rsi,
                &dref,
                &peer_groups,
                &show_aa,
                &field_types,
                col_var,
                &pcols_sort,
                bads,
                cx,
                u_min,
                u_max,
                u_mean,
                median_numis,
                utot,
                median_nreads,
                r_min,
                r_max,
                r_mean,
                rtot,
            ) {
                let _ = proc_cvar2(
                    &var,
                    j,
                    col,
                    mid,
                    pass,
                    u,
                    &ex,
                    &ctl,
                    &exacts,
                    &exact_clonotypes,
                    &refdata,
                    &varmat,
                    out_data,
                    &rsi,
                    &dref,
                    &peer_groups,
                    &show_aa,
                    &field_types,
                    col_var,
                    &pcols_sort,
                    bads,
                    cx,
                    u_min,
                    u_max,
                    u_mean,
                    median_numis,
                    utot,
                    median_nreads,
                    r_min,
                    r_max,
                    r_mean,
                    rtot,
                );
            }
        }
    }
}
