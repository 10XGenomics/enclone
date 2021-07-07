// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function row_fill,
// plus a small helper function get_gex_matrix_entry.

use crate::proc_cvar1::*;
use crate::proc_cvar2::*;
use crate::proc_lvar1::*;
use crate::proc_lvar2::*;
use amino::*;
use enclone_core::allowed_vars::*;
use enclone_core::defs::*;
use enclone_core::median::*;
use enclone_proto::types::*;
use itertools::*;
use ndarray::s;
use std::collections::{HashMap, HashSet};
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
    stats: &mut Vec<(String, Vec<String>)>,
    vdj_cells: &Vec<Vec<String>>,
    n_vdj_gex: &Vec<usize>,
    lvarsc: &Vec<String>,
    lvarsh: &HashSet<String>,
    nd_fields: &Vec<String>,
    peer_groups: &Vec<Vec<(usize, u8, u32)>>,
    extra_args: &Vec<String>,
    all_vars: &Vec<String>,
    need_gex: bool,
    fate: &Vec<HashMap<String, String>>,
) -> Result<(), String> {
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
    macro_rules! speakc {
        ($u:expr, $col:expr, $var:expr, $val:expr) => {
            if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && (ctl.parseable_opt.pchains == "max"
                    || $col + 1 <= ctl.parseable_opt.pchains.force_usize())
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
    let mut count_unsorted = Vec::<usize>::new();
    let mut fcounts = Vec::<f64>::new();
    let mut n_gex = 0;
    let mut n_gexs = Vec::<usize>::new();
    let mut total_counts = Vec::<usize>::new();
    // It might be possible to speed this up a lot by pulling part of the "let d" and
    // "let ind" constructs out of the loop.
    let have_entropy = lvarsh.contains(&"entropy".to_string());
    if have_entropy {
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
    let (mut gex_median, mut gex_min, mut gex_max, mut gex_mean, mut gex_sum) = (0, 0, 0, 0.0, 0.0);
    if have_entropy || need_gex {
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
                                if have_entropy {
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
                                if lvarsh.contains(&"entropy".to_string()) {
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
        count_unsorted = counts.clone();
        counts.sort();
        for n in counts.iter() {
            if *n < 100 {
                *gex_low += 1;
            }
        }
        if counts.len() > 0 {
            gex_median = rounded_median(&counts);
            gex_min = counts[0];
            gex_max = counts[counts.len() - 1];
            gex_sum = fcounts.iter().sum::<f64>();
            gex_mean = gex_sum / fcounts.len() as f64;
        }
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
            if !lvarsh.contains(&LVARS_ALLOWED[i].to_string()) {
                all_lvars.push(LVARS_ALLOWED[i].to_string());
            }
        }
    } else {
        for i in 0..ctl.parseable_opt.pcols.len() {
            if !lvarsh.contains(&ctl.parseable_opt.pcols[i].to_string()) {
                all_lvars.push(ctl.parseable_opt.pcols[i].to_string());
            }
        }
    }
    for x in extra_args.iter() {
        if !lvarsh.contains(&*x) {
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
        if !proc_lvar1(
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
            &extra_args,
            &fate,
        ) {
            let _ = proc_lvar2(
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
                &extra_args,
                &fate,
            );
        }
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
        panic!("{}", msg);
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
        // Set up chain variable macro.  This is the mechanism for generating
        // both human-readable and parseable output for chain variables.

        macro_rules! cvar_stats1 {
            ($i: expr, $var:expr, $val:expr) => {
                if $i < rsi.cvars[col].len() && cvars.contains(&$var) {
                    cx[col][$i] = $val.clone();
                }
                speakc!(u, col, $var, $val);
                let varc = format!("{}{}", $var, col + 1);
                stats.push((varc, vec![$val.to_string(); ex.ncells()]));
            };
        }

        // Process variables that need to be computed even if the chain entry is empty.
        // NO: WHY?  WHY WOULD WE WANT TO PRINT THESE?  BEHAVIOR CHANGED.  DON'T KNOW WHY
        // WE EVER DID THIS>

        let rsi_vars = &ctl.clono_print_opt.cvars;
        let have_notes = rsi.cvars[col].contains(&"notes".to_string());
        let mut notes_pos = 0;
        let mut notes_in = false;
        for j in 0..rsi_vars.len() {
            if all_vars[j] == "notes" {
                notes_pos = j;
                notes_in = true;
            }
        }
        // these lines moved to prevent printing if chain is absent
        let mid = mat[col][u];
        if mid.is_none() {
            continue;
        }
        for j in 0..all_vars.len() {
            let mut jj = j;
            if !have_notes && notes_in && j >= notes_pos {
                jj -= 1;
            }

            // Decide if there is nothing to compute.  This is almost certainly not optimal.
            // Also largely duplicated below.

            let mut needed = false;
            let var = &all_vars[j];
            let varc = format!("{}{}", var, col + 1);
            if jj < rsi.cvars[col].len() && cvars.contains(&var) {
                needed = true;
            } else if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && (ctl.parseable_opt.pchains == "max"
                    || col + 1 <= ctl.parseable_opt.pchains.force_usize())
                && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc))
            {
                needed = true;
            }
            if extra_args.contains(&varc) {
                needed = true;
            }
            if !needed {
                continue;
            }
            let col_var = jj < rsi_vars.len();
            if !col_var && ctl.parseable_opt.pout.len() == 0 && extra_args.is_empty() {
                continue;
            }

            // Process some variables.

            if *var == "v_name" {
                cvar_stats1![j, var, refdata.name[rsi.vids[col]]];
            } else if *var == "v_id" {
                cvar_stats1![j, var, format!("{}", refdata.id[rsi.vids[col]])];
            } else if *var == "d_name" {
                let mut dname = String::new();
                if rsi.dids[col].is_some() {
                    dname = refdata.name[rsi.dids[col].unwrap()].clone();
                }
                cvar_stats1![j, var, dname];
            } else if *var == "d_id" {
                let mut did = String::new();
                if rsi.dids[col].is_some() {
                    did = format!("{}", refdata.id[rsi.dids[col].unwrap()]);
                }
                cvar_stats1![j, var, did];
            } else if *var == "j_name" {
                cvar_stats1![j, var, refdata.name[rsi.jids[col]]];
            } else if *var == "j_id" {
                cvar_stats1![j, var, format!("{}", refdata.id[rsi.jids[col]])];
            } else if *var == "v_name" {
                cvar_stats1![j, var, refdata.name[rsi.vids[col]]];
            } else if *var == "v_id" {
                cvar_stats1![j, var, format!("{}", refdata.id[rsi.vids[col]])];
            } else if *var == "d_name" {
                let mut dname = String::new();
                if rsi.dids[col].is_some() {
                    dname = refdata.name[rsi.dids[col].unwrap()].clone();
                }
                cvar_stats1![j, var, dname];
            } else if *var == "d_id" {
                let mut did = String::new();
                if rsi.dids[col].is_some() {
                    did = format!("{}", refdata.id[rsi.dids[col].unwrap()]);
                }
                cvar_stats1![j, var, did];
            } else if *var == "j_name" {
                cvar_stats1![j, var, refdata.name[rsi.jids[col]]];
            } else if *var == "j_id" {
                cvar_stats1![j, var, format!("{}", refdata.id[rsi.jids[col]])];
            }
        }

        // Keep going.

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

        // Create column entry.

        for j in 0..all_vars.len() {
            let mut jj = j;
            if !have_notes && notes_in && j >= notes_pos {
                jj -= 1;
            }
            if all_vars[j] == "notes" && !have_notes {
                continue;
            }

            // Decide if there is nothing to compute.  This is almost certainly not optimal.

            let mut needed = false;
            let var = &all_vars[j];
            if !ex.share[mid].left
                && (*var == "d1_name"
                    || *var == "d2_name"
                    || *var == "d_delta"
                    || *var == "d_Δ"
                    || *var == "d1_score"
                    || *var == "d2_score")
            {
                continue;
            }
            let varc = format!("{}{}", var, col + 1);
            if jj < rsi.cvars[col].len() && cvars.contains(&var) {
                needed = true;
            } else if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && (ctl.parseable_opt.pchains == "max"
                    || col + 1 <= ctl.parseable_opt.pchains.force_usize())
                && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc))
            {
                needed = true;
            } else if *var == "amino" {
                needed = true;
            } else if *var == "u_cell" || *var == "r_cell" {
                needed = true;
            } else if *var == "white" || ctl.clono_filt_opt_def.whitef {
                needed = true;
            }
            if extra_args.contains(&varc) {
                needed = true;
            }
            if !needed {
                continue;
            }
            let col_var = jj < rsi_vars.len();
            if !col_var && ctl.parseable_opt.pout.len() == 0 && extra_args.is_empty() {
                continue;
            }

            // Compute.

            if !proc_cvar1(
                &var,
                jj,
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
                &extra_args,
                stats,
            )? {
                let _ = proc_cvar2(
                    &var,
                    jj,
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
                    &extra_args,
                    stats,
                );
            }
        }
    }
    Ok(())
}
