// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Make alluvial tables for feature barcode data.  We determine cellular using vdj_cells,
// which is not the only way of doing it.

use enclone_core::defs::{EncloneControl, GexInfo};
use enclone_core::stringulate::*;
use io_utils::fwrite;
use itertools::Itertools;
use stats_utils::percent_ratio;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::io::Write;
use string_utils::*;
use tables::print_tabular_vbox;
use vector_utils::{bin_member, unique_sort};

pub fn alluvial_fb(
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    vdj_cells: &Vec<Vec<String>>,
    logx: &mut Vec<u8>,
) {
    let mut fs = Vec::<FeatureBarcodeAlluvialTable>::new();
    let mut have_some = false;
    for li in 0..ctl.origin_info.n() {
        let m = &gex_info.fb_top_matrices[li];
        if m.initialized() {
            have_some = true;
            let mut keep = 4;
            let mut specials = Vec::<String>::new();
            if !ctl.gen_opt.fb_show.is_empty() {
                let fb_show = ctl.gen_opt.fb_show.split(',').collect::<Vec<&str>>();
                for x in fb_show.iter() {
                    if x.parse::<usize>().is_ok() {
                        keep = x.force_usize();
                    } else {
                        specials.push(x.to_string());
                    }
                }
            }
            let mut cells = Vec::<String>::new();
            for i in 0..vdj_cells[li].len() {
                cells.push(vdj_cells[li][i].before("-").to_string());
            }
            let brn = &gex_info.fb_brn[li];
            let total = gex_info.fb_total_umis[li] as usize;
            let mut seq_to_id = HashMap::<String, String>::new();
            {
                let fref = &gex_info.feature_refs[li];
                let (mut id_pos, mut seq_pos, mut type_pos) = (0, 0, 0);
                for (i, line) in fref.lines().enumerate() {
                    let fields = parse_csv(line);
                    if i == 0 {
                        for j in 0..fields.len() {
                            if fields[j] == "id" {
                                id_pos = j;
                            } else if fields[j] == "sequence" {
                                seq_pos = j;
                            } else if fields[j] == "feature_type" {
                                type_pos = j;
                            }
                        }
                    } else if fields[type_pos] == "Antibody Capture" {
                        seq_to_id.insert(fields[seq_pos].to_string(), fields[id_pos].to_string());
                    }
                }
            }
            let (mut cellular_ref, mut cellular_nref) = (0, 0);
            let (mut ncellular_ref, mut ncellular_nref) = (0, 0);
            {
                for i in 0..brn.len() {
                    if bin_member(&cells, &brn[i].0) {
                        cellular_ref += brn[i].1 as usize;
                        cellular_nref += brn[i].2 as usize;
                    } else {
                        ncellular_ref += brn[i].1 as usize;
                        ncellular_nref += brn[i].2 as usize;
                    }
                }
            }
            let mut top_ref = Vec::<(usize, String)>::new();
            let mut top_nref = Vec::<(usize, String)>::new();
            for i in 0..m.ncols() {
                let bc = m.col_label(i);
                if seq_to_id.contains_key(&bc) {
                    if top_ref.len() < keep {
                        top_ref.push((i, m.col_label(i)));
                    }
                } else if top_nref.len() < keep {
                    top_nref.push((i, m.col_label(i)));
                }
            }
            for i in 0..specials.len() {
                let mut found = false;
                for j in 0..m.ncols() {
                    let bc = m.col_label(j);
                    if specials[i] == bc {
                        if seq_to_id.contains_key(&bc) {
                            top_ref.push((j, m.col_label(j)));
                        } else {
                            top_nref.push((j, m.col_label(j)));
                        }
                        found = true;
                    }
                }
                if !found {
                    if seq_to_id.contains_key(&specials[i]) {
                        top_ref.push((1000000, specials[i].clone()));
                    } else {
                        top_nref.push((1000000, specials[i].clone()));
                    }
                }
            }
            unique_sort(&mut top_ref);
            unique_sort(&mut top_nref);
            let xr = max(top_ref.len(), 1);
            let xnr = max(top_nref.len(), 1);
            let nrows = 4 * (xr + xnr) - 1;
            let ncols = 4;
            let mut rows = vec![vec![String::new(); ncols]; nrows];
            let mut csv_rows = vec![vec![String::new(); 6]; nrows];
            rows[2 * (xr + xnr) - 1][0] = "100.0".to_string();
            for j in 1..ncols {
                rows[2 * (xr + xnr) - 1][j] = "\\hline".to_string();
            }
            let mut count = 0;
            for pass in 0..2 {
                for i in 0..xr {
                    count += 1;
                    rows[count][3] = "\\hline".to_string();
                    if i == xr - 1 {
                        rows[count][2] = "\\hline".to_string();
                    }
                    count += 1
                }
                for i in 0..xnr {
                    count += 1;
                    if pass == 0 || i < xnr - 1 {
                        rows[count][3] = "\\hline".to_string();
                    }
                    count += 1;
                }
            }
            fn pr(x: usize, y: usize) -> String {
                format!("{:>4.1}", percent_ratio(x, y))
            }
            fn pr0(x: usize, y: usize) -> String {
                format!("{:.1}", percent_ratio(x, y))
            }
            for pass in 0..2 {
                for i in 0..top_ref.len() {
                    let c = top_ref[i].0;
                    let seq = &top_ref[i].1;
                    let (mut cell, mut ncell) = (0, 0);
                    let label = format!("{} = {}", seq, seq_to_id[seq]);
                    if c < 1000000 {
                        for j in 0..m.nrows() {
                            if bin_member(&vdj_cells[li], &m.row_label(j)) {
                                cell += m.value(j, c as usize);
                            } else {
                                ncell += m.value(j, c as usize);
                            }
                        }
                    }
                    if pass == 0 {
                        let r = 2 * i;
                        rows[r][3] = format!("{} {}", pr(cell, total), label);
                        csv_rows[r][0] = ctl.origin_info.dataset_id[li].clone();
                        csv_rows[r][1] = "cellular".to_string();
                        csv_rows[r][2] = "reference".to_string();
                        csv_rows[r][3] = format!("{}", pr0(cell, total));
                        csv_rows[r][4] = seq.clone();
                        csv_rows[r][5] = seq_to_id[seq].clone();
                    } else {
                        let r = 2 * (xr + xnr) + 2 * i;
                        rows[r][3] = format!("{} {}", pr(ncell, total), label);
                        csv_rows[r][0] = ctl.origin_info.dataset_id[li].clone();
                        csv_rows[r][1] = "noncellular".to_string();
                        csv_rows[r][2] = "reference".to_string();
                        csv_rows[r][3] = format!("{}", pr0(ncell, total));
                        csv_rows[r][4] = seq.clone();
                        csv_rows[r][5] = seq_to_id[seq].clone();
                    }
                }
                for i in 0..top_nref.len() {
                    let c = top_nref[i].0;
                    let seq = &top_nref[i].1;
                    let (mut cell, mut ncell) = (0, 0);
                    if c < 1000000 {
                        for j in 0..m.nrows() {
                            if bin_member(&vdj_cells[li], &m.row_label(j)) {
                                cell += m.value(j, c as usize);
                            } else {
                                ncell += m.value(j, c as usize);
                            }
                        }
                    }
                    if pass == 0 {
                        let r = 2 * (i + xr);
                        rows[r][3] = format!("{} {}", pr(cell, total), seq);
                        csv_rows[r][0] = ctl.origin_info.dataset_id[li].clone();
                        csv_rows[r][1] = "cellular".to_string();
                        csv_rows[r][2] = "nonreference".to_string();
                        csv_rows[r][3] = format!("{}", pr0(cell, total));
                        csv_rows[r][4] = seq.clone();
                    } else {
                        let r = 2 * (xr + xnr) + 2 * (i + xr);
                        csv_rows[r][0] = ctl.origin_info.dataset_id[li].clone();
                        csv_rows[r][1] = "noncellular".to_string();
                        csv_rows[r][2] = "nonreference".to_string();
                        csv_rows[r][3] = format!("{}", pr0(ncell, total));
                        csv_rows[r][4] = seq.clone();
                        rows[r][3] = format!("{} {}", pr(ncell, total), seq);
                    }
                }
            }
            rows[xr - 1][2] = format!("{} reference", pr(cellular_ref, total));
            rows[2 * (xr + xnr) + xr - 1][2] = format!("{} reference", pr(ncellular_ref, total));
            rows[2 * xr + xnr - 1][2] = format!("{} nonreference", pr(cellular_nref, total));
            rows[2 * xr + 2 * (xr + xnr) + xnr - 1][2] =
                format!("{} nonreference", pr(ncellular_nref, total));
            rows[xr + xnr - 1][1] = format!("{} cellular", pr(cellular_ref + cellular_nref, total));
            rows[2 * (xr + xnr) + xr + xnr - 1][1] =
                format!("{} noncellular", pr(ncellular_ref + ncellular_nref, total));
            let mut display_text = String::new();
            print_tabular_vbox(
                &mut display_text,
                &rows,
                0,
                &b"l|l|l|l".to_vec(),
                false,
                false,
            );
            if !ctl.visual_mode {
                fwrite!(
                    logx,
                    "\nfeature barcode UMI distribution for {}\n{}",
                    ctl.origin_info.dataset_id[li],
                    display_text
                );
            } else {
                let mut spreadsheet_text = String::new();
                for (i, r) in csv_rows.iter().enumerate() {
                    if i % 2 == 0 {
                        spreadsheet_text += &mut format!("{}\n", r.iter().format("\t "));
                    }
                }
                let f = FeatureBarcodeAlluvialTable {
                    id: ctl.origin_info.dataset_id[li].clone(),
                    display_text: display_text.clone(),
                    spreadsheet_text: spreadsheet_text.clone(),
                };
                fs.push(f);
            }
        }
    }
    if ctl.visual_mode && have_some {
        let tables = FeatureBarcodeAlluvialTableSet { s: fs };
        logx.append(&mut tables.to_string().as_bytes().to_vec());
    }
}

pub fn common_gumis(ctl: &EncloneControl, gex_info: &GexInfo, logx: &mut Vec<u8>) {
    let mut fs = Vec::<FeatureBarcodeCommonGumisTable>::new();
    let mut have_some = false;
    for li in 0..ctl.origin_info.n() {
        let m = &gex_info.fb_top_matrices[li];
        if m.initialized() {
            have_some = true;
            let freqs = &gex_info.fb_common_gumis[li].0;
            let umis = &gex_info.fb_common_gumis[li].1;
            let mut rows = Vec::<Vec<String>>::new();
            let mut csv_rows = Vec::<Vec<String>>::new();
            let row = vec!["percent".to_string(), "UMI".to_string()];
            rows.push(row);
            rows.push(vec!["\\hline".to_string(); 2]);
            for i in 0..min(20, freqs.len()) {
                let row = vec![format!("{:.1}", 100.0 * freqs[i]), stringme(&umis[i])];
                rows.push(row);
                csv_rows.push(vec![
                    ctl.origin_info.dataset_id[li].clone(),
                    format!("{:.1}", 100.0 * freqs[i]),
                    stringme(&umis[i]),
                ]);
            }
            let mut display_text = String::new();
            print_tabular_vbox(&mut display_text, &rows, 0, &b"r|l".to_vec(), false, false);
            if !ctl.visual_mode {
                fwrite!(
                    logx,
                    "\nfrequent UMIs having feature barcode GGGGGGGGGGGGGGG for {}\n{}",
                    ctl.origin_info.dataset_id[li],
                    display_text
                );
            } else {
                let mut spreadsheet_text = String::new();
                for r in csv_rows.iter() {
                    spreadsheet_text += &mut format!("{}\n", r.iter().format("\t "));
                }
                let f = FeatureBarcodeCommonGumisTable {
                    id: ctl.origin_info.dataset_id[li].clone(),
                    display_text: display_text.clone(),
                    spreadsheet_text: spreadsheet_text.clone(),
                };
                fs.push(f);
            }
        }
    }
    if ctl.visual_mode && have_some {
        let tables = FeatureBarcodeCommonGumisTableSet { s: fs };
        logx.append(&mut tables.to_string().as_bytes().to_vec());
    }
}
