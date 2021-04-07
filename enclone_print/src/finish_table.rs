// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::build_table_stuff::*;
use crate::print_utils1::*;
use crate::print_utils3::*;
use crate::print_utils5::*;
use enclone_core::defs::*;
use enclone_proto::types::*;
use std::collections::HashMap;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

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
    mlog: &mut Vec<u8>,
    logz: &mut String,
    stats: &Vec<(String, Vec<String>)>,
    sr: &mut Vec<(Vec<String>, Vec<Vec<String>>, Vec<Vec<u8>>, usize)>,
    extra_args: &Vec<String>,
    pcols_sort: &Vec<String>,
    out_data: &mut Vec<HashMap<String, String>>,
    rord: &Vec<usize>,
    pass: usize,
) {
    // Fill in exact_subclonotype_id, reorder.

    let nexacts = exacts.len();
    if ctl.parseable_opt.pout.len() > 0 || !extra_args.is_empty() {
        for u in 0..nexacts {
            macro_rules! speak {
                ($u:expr, $var:expr, $val:expr) => {
                    if pass == 2 && (ctl.parseable_opt.pout.len() > 0 || !extra_args.is_empty()) {
                        if pcols_sort.is_empty()
                            || bin_member(&pcols_sort, &$var.to_string())
                            || bin_member(&extra_args, &$var.to_string())
                        {
                            out_data[$u].insert($var.to_string(), $val);
                        }
                    }
                };
            }
            speak![rord[u], "exact_subclonotype_id", format!("{}", u + 1)];
        }
        let mut out_data2 = vec![HashMap::<String, String>::new(); nexacts];
        for v in 0..nexacts {
            out_data2[v] = out_data[rord[v]].clone();
        }
        *out_data = out_data2;
    }

    // Add header text to mlog.

    let mat = &rsi.mat;
    add_header_text(&ctl, &exacts, &exact_clonotypes, &rord, &mat, mlog);

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
                        if stats[j].1[k].parse::<f64>().is_ok() {
                            total += stats[j].1[k].force_f64();
                        }
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
                        if stats[j].1[k].parse::<f64>().is_ok() {
                            total += stats[j].1[k].force_f64();
                        }
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

    // Add phylogeny.

    let nexacts = exacts.len();
    if ctl.toy {
        let mut vrefs = Vec::<Vec<u8>>::new();
        let mut jrefs = Vec::<Vec<u8>>::new();
        for cx in 0..cols {
            let (mut vref, mut jref) = (Vec::<u8>::new(), Vec::<u8>::new());
            for u in 0..nexacts {
                let m = rsi.mat[cx][u];
                if m.is_some() {
                    let m = m.unwrap();
                    jref = exact_clonotypes[exacts[u]].share[m].js.to_ascii_vec();
                }
                let vseq1 = refdata.refs[rsi.vids[cx]].to_ascii_vec();
                if rsi.vpids[cx].is_some() {
                    vref = dref[rsi.vpids[cx].unwrap()].nt_sequence.clone();
                } else {
                    vref = vseq1.clone();
                }
            }
            vrefs.push(vref);
            jrefs.push(jref);
        }
        for u1 in 0..nexacts {
            let ex1 = &exact_clonotypes[exacts[u1]];
            for u2 in u1 + 1..nexacts {
                let ex2 = &exact_clonotypes[exacts[u2]];
                let (mut d1, mut d2) = (0, 0);
                let mut d = 0;
                for cx in 0..cols {
                    let (m1, m2) = (rsi.mat[cx][u1], rsi.mat[cx][u2]);
                    if m1.is_none() || m2.is_none() {
                        continue;
                    }
                    let (m1, m2) = (m1.unwrap(), m2.unwrap());
                    let (s1, s2) = (&ex1.share[m1].seq_del, &ex2.share[m2].seq_del);
                    let n = s1.len();
                    let (vref, jref) = (&vrefs[cx], &jrefs[cx]);
                    for j in 0..vars[cx].len() {
                        let p = vars[cx][j];
                        if s1[p] != s2[p] {
                            if p < vref.len() - ctl.heur.ref_v_trim {
                                if s1[p] == vref[p] {
                                    d1 += 1;
                                } else if s2[p] == vref[p] {
                                    d2 += 1;
                                }
                            } else if p >= n - (jref.len() - ctl.heur.ref_j_trim) {
                                if s1[p] == jref[jref.len() - (n - p)] {
                                    d1 += 1;
                                } else if s2[p] == jref[jref.len() - (n - p)] {
                                    d2 += 1;
                                }
                            } else {
                                d += 1;
                            }
                        }
                    }
                }
                if (d1 == 0) ^ (d2 == 0) {
                    if d1 == 0 {
                        *logz += &format!("{} ==> {}", u1 + 1, u2 + 1);
                    } else {
                        *logz += &format!("{} ==> {}", u2 + 1, u1 + 1);
                    }
                    let s = format!(
                        "; u1 = {}, u2 = {}, d1 = {}, d2 = {}, d = {}\n",
                        u1 + 1,
                        u2 + 1,
                        d1,
                        d2,
                        d
                    );
                    *logz += &s;
                }
            }
        }
    }
}
