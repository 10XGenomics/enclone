// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::print_utils1::*;
use crate::print_utils3::*;
use crate::print_utils5::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use string_utils::*;
use vdj_ann::refx::*;

pub fn finish_table(
    n: usize,
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    rsi: &ColInfo,
    vars: &Vec<Vec<usize>>,
    show_aa: &Vec<Vec<usize>>,
    field_types: &Vec<Vec<u8>>,
    lvars: &Vec<String>,
    refdata: &RefData,
    dref: &Vec<DonorReferenceItem>,
    peer_groups: &Vec<Vec<(usize, u8, u32)>>,
    mlog: &Vec<u8>,
    logz: &mut String,
    stats: &Vec<(String, Vec<f64>)>,
    sr: &mut Vec<(Vec<String>, Vec<Vec<String>>, Vec<Vec<u8>>, usize)>,
) {
    // Build table stuff.

    let mut row1 = Vec::<String>::new();
    let mut justify = Vec::<u8>::new();
    let mut rows = Vec::<Vec<String>>::new();
    let mut drows = Vec::<Vec<String>>::new();
    build_table_stuff(
        &ctl,
        &exacts,
        &exact_clonotypes,
        &rsi,
        &vars,
        &show_aa,
        &field_types,
        &mut row1,
        &mut justify,
        &mut drows,
        &mut rows,
        &lvars,
    );

    // Insert universal and donor reference rows.

    insert_reference_rows(
        &ctl,
        &rsi,
        &show_aa,
        &field_types,
        &refdata,
        &dref,
        &row1,
        &mut drows,
        &mut rows,
        &exacts,
        &exact_clonotypes,
        &peer_groups,
    );

    // Insert consensus row.

    insert_consensus_row(
        &ctl,
        &rsi,
        exacts.len(),
        &field_types,
        &show_aa,
        &row1,
        &mut rows,
    );

    // Insert horizontal line.

    let cols = rsi.mat.len();
    if !drows.is_empty() {
        let mut width = 1 + lvars.len();
        for col in 0..cols {
            width += rsi.cvars[col].len();
        }
        rows.push(vec!["\\hline".to_string(); width]);
    }

    // Build the diff row.

    build_diff_row(
        &ctl,
        &rsi,
        &mut rows,
        &mut drows,
        &row1,
        exacts.len(),
        &field_types,
        &show_aa,
    );

    // Finish building table content.

    for j in 0..sr.len() {
        sr[j].0[0] = format!("{}", j + 1); // row number (#)
        rows.push(sr[j].0.clone());
        rows.append(&mut sr[j].1.clone());
    }

    // Add sum and mean rows.

    if ctl.clono_print_opt.sum {
        let mut row = Vec::<String>::new();
        row.push("Σ".to_string());
        for i in 0..lvars.len() {
            let mut x = lvars[i].clone();
            if x.contains(':') {
                x = x.before(":").to_string();
            }
            let mut found = false;
            let mut total = 0.0;
            for j in 0..stats.len() {
                if stats[j].0 == x {
                    found = true;
                    for k in 0..stats[j].1.len() {
                        total += stats[j].1[k];
                    }
                }
            }
            if !found {
                row.push(String::new());
            } else {
                if !lvars[i].ends_with("_%") {
                    row.push(format!("{}", total.round() as usize));
                } else {
                    row.push(format!("{:.2}", total));
                }
            }
        }
        // This is necessary but should not be:
        for cx in 0..cols {
            for _ in 0..rsi.cvars[cx].len() {
                row.push(String::new());
            }
        }
        rows.push(row);
    }
    if ctl.clono_print_opt.mean {
        let mut row = Vec::<String>::new();
        row.push("μ".to_string());
        for i in 0..lvars.len() {
            let mut x = lvars[i].clone();
            if x.contains(':') {
                x = x.before(":").to_string();
            }
            let mut found = false;
            let mut total = 0.0;
            for j in 0..stats.len() {
                if stats[j].0 == x {
                    found = true;
                    for k in 0..stats[j].1.len() {
                        total += stats[j].1[k];
                    }
                }
            }
            let mean = total / n as f64;
            if !found {
                row.push(String::new());
            } else {
                if !lvars[i].ends_with("_%") {
                    row.push(format!("{:.1}", mean));
                } else {
                    row.push(format!("{:.2}", mean));
                }
            }
        }
        // This is necessary but should not be:
        for cx in 0..cols {
            for _ in 0..rsi.cvars[cx].len() {
                row.push(String::new());
            }
        }
        rows.push(row);
    }

    // Make table.

    for i in 0..rows.len() {
        for j in 0..rows[i].len() {
            rows[i][j] = rows[i][j].replace("|TRX", "TRB");
            rows[i][j] = rows[i][j].replace("|TRY", "TRA");
        }
    }
    for cx in 0..cols {
        justify.push(b'|');
        for m in 0..rsi.cvars[cx].len() {
            justify.push(justification(&rsi.cvars[cx][m]));
        }
    }
    make_table(&ctl, &mut rows, &justify, &mlog, logz);
}
