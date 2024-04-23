// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the single function row_fill,
// plus a small helper function get_gex_matrix_entry.

use crate::print_clonotypes::print_utils1::color_codon;

use super::proc_cvar_auto::proc_cvar_auto;
use super::proc_lvar2::proc_lvar2;
use super::proc_lvar_auto::proc_lvar_auto;
use amino::{codon_to_aa, nucleotide_to_aminoacid_sequence};
use enclone_core::allowed_vars::LVARS_ALLOWED;
use enclone_core::defs::{AlleleData, ColInfo, EncloneControl, ExactClonotype, GexInfo, POUT_SEP};
use enclone_core::enclone_structs::{BarcodeFates, GexReaders};
use enclone_core::median::median_f64;
use enclone_proto::types::DonorReferenceItem;
use enclone_vars::decode_arith;
use expr_tools::{define_evalexpr_context, vars_of_node};
use itertools::Itertools;
use stats_utils::percent_ratio;
use std::collections::{HashMap, HashSet};
use std::fmt::Write;
use string_utils::{stringme, strme, TextUtils};
use vdj_ann::refx::RefData;
use vector_utils::next_diff12_4;
use vector_utils::{bin_member, bin_position, unique_sort};

/// The following code creates a row in the enclone output table for a clonotype.  Simultaneously
/// it generates a row of parseable output.  And it does some other things that are not described
/// here.
///
/// TODO: Awful interface, should work to improve.
pub fn row_fill(
    pass: usize,
    u: usize,
    ctl: &EncloneControl,
    exacts: &[usize],
    mults: &[usize],
    exact_clonotypes: &[ExactClonotype],
    gex_info: &GexInfo,
    refdata: &RefData,
    varmat: &[Vec<Vec<u8>>],
    fp: &[Vec<usize>],
    vars_amino: &[Vec<usize>],
    show_aa: &[Vec<usize>],
    ref_diff_pos: &[Vec<Vec<usize>>],
    field_types: &[Vec<u8>],
    bads: &mut [bool],
    row: &mut Vec<String>,                    // row of human-readable output
    out_data: &mut [HashMap<String, String>], // row of parseable output
    cx: &mut [Vec<String>],
    d_all: &mut [Vec<u32>],
    ind_all: &mut [Vec<u32>],
    rsi: &ColInfo,
    dref: &[DonorReferenceItem],
    groups: &HashMap<usize, Vec<usize>>,
    gex_readers: &[Option<GexReaders<'_>>],
    stats: &mut Vec<(String, Vec<String>)>,
    stats_pass1: &[Vec<(String, Vec<String>)>],
    vdj_cells: &[Vec<String>],
    n_vdj_gex: &[usize],
    lvarsc: &[String],
    lvarsh: &HashSet<String>,
    nd_fields: &[String],
    peer_groups: &[Vec<(usize, u8, u32)>],
    extra_args: &[String],
    all_vars: &[String],
    need_gex: bool,
    fate: &[BarcodeFates],
    cdr3_con: &[Vec<u8>],
    allele_data: &AlleleData,
) -> Result<(), String> {
    // Redefine some things to reduce dependencies.

    let mat = &rsi.mat;
    let cvars = &ctl.clono_print_opt.cvars;
    let lvars = lvarsc;
    let clonotype_id = exacts[u];
    let ex = &exact_clonotypes[clonotype_id];
    let mut pcols_sort = ctl.parseable_opt.pcols_sort.clone();
    for pcol in &mut pcols_sort {
        *pcol = pcol.replace("_Σ", "_sum").replace("_μ", "_mean");
    }
    pcols_sort.sort();
    macro_rules! speakc {
        ($u:expr, $col:expr, $var:expr, $val:expr) => {
            if pass == 2
                && ctl.parseable_opt.pout.len() > 0
                && (ctl.parseable_opt.pchains == "max"
                    || $col < ctl.parseable_opt.pchains.force_usize())
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
        eprintln!();
    }

    // Compute dataset indices, gex, gex_mean, gex_sum,
    // n_gex_cell, n_gex.

    let mut dataset_indices = Vec::<usize>::new();
    for l in 0..ex.clones.len() {
        dataset_indices.push(ex.clones[l][0].dataset_index);
    }
    unique_sort(&mut dataset_indices);
    let mut lenas = Vec::<String>::new();
    for l in &dataset_indices {
        lenas.push(ctl.origin_info.dataset_id[*l].clone());
    }
    row.push(String::new()); // row number (#), filled in below
    let mut counts = Vec::<usize>::new();
    let mut gex_counts_unsorted = Vec::<usize>::new();
    let mut gex_fcounts_unsorted = Vec::<f64>::new();
    let mut n_gexs = Vec::<usize>::new();

    // It may not make any sense at all for this code to be here.

    if ctl.clono_filt_opt_def.whitef {
        let mut bch = vec![Vec::<(usize, String, usize, usize)>::new(); 2];
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = &ex.clones[l][0].barcode;
            let mut numi = 0;
            for j in 0..ex.clones[l].len() {
                numi += ex.clones[l][j].umi_count;
            }
            bch[0].push((li, bc[0..8].to_string(), numi, l));
            bch[1].push((li, bc[8..16].to_string(), numi, l));
        }
        let mut junk = 0;
        let mut bad = vec![false; ex.clones.len()];
        for l in 0..2 {
            bch[l].sort();
            let mut m = 0;
            while m < bch[l].len() {
                let n = next_diff12_4(&bch[l], m);
                for u1 in m..n {
                    for u2 in m..n {
                        if bch[l][u1].2 >= 10 * bch[l][u2].2 {
                            bad[bch[l][u2].3] = true;
                        }
                    }
                }
                m = n;
            }
        }
        for b in bad {
            if b {
                junk += 1;
            }
        }
        let junk_rate = percent_ratio(junk, ex.clones.len());
        // WRONG!  THIS IS SUPPOSED TO BE EXECUTED ON PASS 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if junk_rate == 0.0 {
            bads[u] = true;
        }
    }

    // It might be possible to speed this up a lot by pulling part of the "let d" and
    // "let ind" constructs out of the loop.

    let (mut gex_mean, mut gex_sum) = (0.0, 0.0);
    if need_gex {
        for l in 0..ex.clones.len() {
            let li = ex.clones[l][0].dataset_index;
            let bc = ex.clones[l][0].barcode.clone();
            if !gex_info.gex_barcodes.is_empty() {
                if bin_member(&gex_info.gex_cell_barcodes[li], &bc) {
                    n_gexs.push(1);
                } else {
                    n_gexs.push(0);
                }
                let mut count = 0;
                let mut fcount = 0.0;
                let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                if p >= 0 {
                    let mut raw_count = 0;

                    let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                    let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize; // is p+1 OK??
                    let (d, ind) = gex_readers[li].as_ref().unwrap().get_range(z1..z2).unwrap();
                    for j in 0..d.len() {
                        if gex_info.is_gex[li][ind[j] as usize] {
                            let n = d[j] as usize;
                            raw_count += n;
                        }
                    }
                    d_all[l] = d;
                    ind_all[l] = ind;

                    if !ctl.gen_opt.full_counts {
                        count = (raw_count as f64 * gex_info.gex_mults[li]).round() as usize;
                        fcount = raw_count as f64 * gex_info.gex_mults[li];
                    } else {
                        count = (raw_count as f64).round() as usize;
                        fcount = raw_count as f64;
                    }
                }
                counts.push(count);
                gex_fcounts_unsorted.push(fcount);
            }
        }
        gex_counts_unsorted = counts.clone();
        counts.sort_unstable();
        if !counts.is_empty() {
            gex_sum = gex_fcounts_unsorted.iter().sum::<f64>();
            gex_mean = gex_sum / gex_fcounts_unsorted.len() as f64;
        }
    }

    // Output lead variable columns.
    // WARNING!  If you add lead variables, you may need to add them to the function
    // LinearCondition::require_valid_variables.

    let mut all_lvars = lvars.iter().map(String::as_str).collect::<Vec<_>>();
    if ctl.parseable_opt.pout.is_empty() {
    } else if ctl.parseable_opt.pcols.is_empty() {
        for var in LVARS_ALLOWED {
            if !lvarsh.contains(&var.to_string()) {
                all_lvars.push(var);
            }
        }
    } else {
        for i in 0..ctl.parseable_opt.pcols.len() {
            if !lvarsh.contains(&ctl.parseable_opt.pcols[i].to_string()) {
                all_lvars.push(ctl.parseable_opt.pcols[i].as_str());
            }
        }
    }
    for x in extra_args {
        if !lvarsh.contains(x) {
            all_lvars.push(x.as_str());
        }
    }
    let mut alt_bcs = Vec::<&str>::new();
    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.as_str());
        }
    }
    unique_sort(&mut alt_bcs);

    macro_rules! speak {
        ($u:expr, $var:expr, $val:expr) => {
            if pass == 2 && (ctl.parseable_opt.pout.len() > 0 || extra_args.len() > 0) {
                let mut v = $var.to_string();
                v = v.replace("_Σ", "_sum");
                v = v.replace("_μ", "_mean");
                if ctl.parseable_opt.pcols.is_empty()
                    || bin_member(&ctl.parseable_opt.pcols_sortx, &v)
                    || bin_member(&extra_args, &v)
                {
                    out_data[$u].insert(v, $val);
                }
            }
        };
    }

    'lvar_loop: for (i, &x) in all_lvars.iter().enumerate() {
        // Process VAR_DEF variables.

        for var_def in &ctl.gen_opt.var_def {
            if x == var_def.0 && i < lvars.len() {
                if pass == 2 {
                    let comp = &var_def.2;
                    let vars = vars_of_node(comp); // computing this here might be inefficient
                    let mut out_vals = Vec::<String>::new();
                    for k in 0..ex.clones.len() {
                        let mut in_vals = Vec::<String>::new();
                        for var in &vars {
                            let var = decode_arith(var);
                            let mut found = false;
                            for stat in stats.iter() {
                                if stat.0 == var {
                                    in_vals.push(stat.1[k].clone());
                                    found = true;
                                    break;
                                }
                            }
                            if !found {
                                for stat in &stats_pass1[u] {
                                    if stat.0 == var {
                                        in_vals.push(stat.1[k].clone());
                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if !found {
                                in_vals.push(String::new());
                            }
                        }
                        let c = define_evalexpr_context(&vars, &in_vals);
                        let res = comp.eval_with_context(&c);
                        // if res.is_err() {
                        //     eprintln!("\nInternal error, failed to compute {}.\n", x);
                        //     std::process::exit(1);
                        // }
                        let val = res.unwrap();
                        let val = val.as_number();
                        out_vals.push(match val {
                            Err(_) => String::new(),
                            Ok(val) => format!("{val:.1}"),
                        });
                    }
                    let mut median = String::new();
                    let mut out_valsf = Vec::<f64>::new();
                    let mut all_float = true;
                    for y in &out_vals {
                        if y.parse::<f64>().is_err() {
                            all_float = false;
                        } else {
                            out_valsf.push(y.force_f64());
                        }
                    }
                    if all_float {
                        out_valsf.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        median = format!("{:.1}", median_f64(&out_valsf));
                    }
                    if i < lvars.len() {
                        row.push(median.clone());
                    }
                    if pass == 2 {
                        if ctl.parseable_opt.pbarcode {
                            speak!(
                                u,
                                x.to_string(),
                                format!("{}", out_vals.iter().format(POUT_SEP))
                            );
                        } else {
                            speak!(u, x.to_string(), median.clone());
                        }
                    }
                    stats.push((x.to_string(), out_vals.clone()));
                } else if i < lvars.len() {
                    row.push(String::new());
                }
                continue 'lvar_loop;
            }
        }

        // Process other lvars.

        if !proc_lvar_auto(
            i,
            pass,
            x,
            exacts,
            exact_clonotypes,
            u,
            rsi,
            refdata,
            ctl,
            extra_args,
            out_data,
            stats,
            lvars,
            row,
            fate,
            dref,
            varmat,
            fp,
            n_vdj_gex,
            vdj_cells,
            gex_info,
            groups,
            mults,
            nd_fields,
            &gex_counts_unsorted,
            &gex_fcounts_unsorted,
            &n_gexs,
            gex_readers,
            &alt_bcs,
        )? {
            let _ = proc_lvar2(
                i,
                x,
                pass,
                u,
                ctl,
                exacts,
                exact_clonotypes,
                gex_info,
                row,
                out_data,
                d_all,
                ind_all,
                stats,
                lvars,
                &alt_bcs,
                gex_mean,
                gex_sum,
                &gex_fcounts_unsorted,
                extra_args,
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
            write!(
                msg,
                "\n\nYou may find it helpful to add the options\n\
                BARCODE={} ROW_FILL_VERBOSE\n\
                to the command line.  Depending on other arguments, you might also need to \
                add MAX_CORES=1.",
                ex.clones[0][0].barcode
            )
            .unwrap();
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
        // Process variables that need to be computed even if the chain entry is empty.
        // NO: WHY?  WHY WOULD WE WANT TO PRINT THESE?  BEHAVIOR CHANGED.  DON'T KNOW WHY
        // WE EVER DID THIS>

        let rsi_vars = &ctl.clono_print_opt.cvars;
        let have_notes = rsi.cvars[col].contains(&"notes".to_string());
        let mut notes_pos = 0;
        let mut notes_in = false;
        for (j, var) in all_vars.iter().take(rsi_vars.len()).enumerate() {
            if var == "notes" {
                notes_pos = j;
                notes_in = true;
            }
        }
        // these lines moved to prevent printing if chain is absent
        let mid = mat[col][u];
        if mid.is_none() {
            continue;
        }
        for (j, var) in all_vars.iter().enumerate() {
            let mut jj = j;
            if !have_notes && notes_in && j >= notes_pos {
                jj -= 1;
            }

            // Decide if there is nothing to compute.  This is almost certainly not optimal.
            // Also largely duplicated below.

            let mut needed = false;
            let varc = format!("{var}{}", col + 1);
            if jj < rsi.cvars[col].len() && cvars.contains(&var.to_string())
                || pass == 2
                    && !ctl.parseable_opt.pout.is_empty()
                    && (ctl.parseable_opt.pchains == "max"
                        || col < ctl.parseable_opt.pchains.force_usize())
                    && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc))
                || extra_args.contains(&varc)
            {
                needed = true;
            }
            if !needed {
                continue;
            }
            let col_var = jj < rsi_vars.len();
            if !col_var && ctl.parseable_opt.pout.is_empty() && extra_args.is_empty() {
                continue;
            }
        }

        // Keep going.

        let mid = mid.unwrap();
        let ex = &exact_clonotypes[clonotype_id];
        let seq_amino = rsi.seqss_amino[col][u].clone();

        // Speak some other column entries.

        let xm = &ex.share[mid];
        speakc!(
            u,
            col,
            "vj_aa".to_string(),
            stringme(&nucleotide_to_aminoacid_sequence(&xm.seq, 0))
        );
        speakc!(u, col, "vj_seq".to_string(), stringme(&xm.seq));
        let mut dna = Vec::<u8>::new();
        for p in xm.fr1_start..xm.seq_del_amino.len() {
            for j in 0..xm.ins.len() {
                if xm.ins[j].0 == p {
                    let mut z = xm.ins[j].1.clone();
                    dna.append(&mut z);
                }
            }
            if xm.seq_del_amino[p] != b'-' {
                dna.push(xm.seq_del_amino[p]);
            }
        }
        speakc!(
            u,
            col,
            "vj_aa_nl".to_string(),
            stringme(&nucleotide_to_aminoacid_sequence(&dna, 0))
        );
        speakc!(u, col, "vj_seq_nl".to_string(), stringme(&dna));
        speakc!(u, col, "seq".to_string(), stringme(&xm.full_seq));
        let mut vv = Vec::<usize>::new();
        for x in &vars_amino[col] {
            vv.push(*x / 3);
        }
        unique_sort(&mut vv);
        let mut varaa = Vec::<u8>::new();
        for p in &vv {
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

        for (j, var) in all_vars.iter().enumerate() {
            let mut jj = j;
            if !have_notes && notes_in && j >= notes_pos {
                jj -= 1;
            }
            if var == "notes" && !have_notes {
                continue;
            }

            // Decide if there is nothing to compute.  This is almost certainly not optimal.

            let mut needed = false;
            if !ex.share[mid].left
                && (var == "d1_name"
                    || var == "d2_name"
                    || var == "d_delta"
                    || var == "d_Δ"
                    || var == "d1_score"
                    || var == "d2_score")
            {
                continue;
            }
            let varc = format!("{var}{}", col + 1);
            if jj < rsi.cvars[col].len() && cvars.contains(&var.to_string())
                || pass == 2
                    && !ctl.parseable_opt.pout.is_empty()
                    && (ctl.parseable_opt.pchains == "max"
                        || col < ctl.parseable_opt.pchains.force_usize())
                    && (pcols_sort.is_empty() || bin_member(&pcols_sort, &varc))
                || var == "amino"
                || var == "u_cell"
                || var == "r_cell"
                || var == "white"
                || ctl.clono_filt_opt_def.whitef
                || extra_args.contains(&varc)
            {
                needed = true;
            }
            if !needed {
                continue;
            }
            let col_var = jj < rsi_vars.len();
            if !col_var && ctl.parseable_opt.pout.is_empty() && extra_args.is_empty() {
                continue;
            }

            // Compute.

            if !proc_cvar_auto(
                jj,
                pass,
                var,
                ex,
                exacts,
                exact_clonotypes,
                mid,
                col,
                u,
                rsi,
                refdata,
                dref,
                ctl,
                extra_args,
                &pcols_sort,
                cx,
                varmat,
                out_data,
                stats,
                allele_data,
            )? && var == "amino"
                && col_var
            {
                let mut last_color = "black".to_string();
                for k in 0..show_aa[col].len() {
                    let p = show_aa[col][k];
                    if k > 0
                        && field_types[col][k] != field_types[col][k - 1]
                        && !ctl.gen_opt.nospaces
                    {
                        cx[col][jj] += " ";
                    }
                    if 3 * p + 3 <= seq_amino.len()
                        && seq_amino[3 * p..3 * p + 3].to_vec() == b"---".to_vec()
                    {
                        cx[col][jj] += "-";
                    } else if 3 * p + 3 > seq_amino.len()
                        || seq_amino[3 * p..3 * p + 3].contains(&b'-')
                    {
                        cx[col][jj] += "*";
                    } else {
                        let x = &peer_groups[rsi.vids[col]];
                        let last = k == show_aa[col].len() - 1;
                        let log = color_codon(
                            ctl,
                            &seq_amino,
                            ref_diff_pos,
                            x,
                            col,
                            mid,
                            p,
                            u,
                            &mut last_color,
                            last,
                            cdr3_con,
                            exacts,
                            exact_clonotypes,
                        );
                        cx[col][jj] += strme(&log);
                    }
                }
            }
        }
    }
    Ok(())
}
