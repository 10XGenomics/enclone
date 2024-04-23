// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use ansi_escape::bold;
use enclone_core::defs::{justification, ColInfo, EncloneControl, ExactClonotype};
use itertools::Itertools;
use string_utils::{strme, TextUtils};
use vector_utils::unique_sort;

use super::print_utils1::insert_position_rows;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn build_table_stuff(
    ctl: &EncloneControl,
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    rsi: &ColInfo,
    vars: &[Vec<usize>],
    show_aa: &[Vec<usize>],
    field_types: &[Vec<u8>],
    row1: &mut Vec<String>,
    justify: &mut Vec<u8>,
    drows: &mut Vec<Vec<String>>,
    rows: &mut Vec<Vec<String>>,
    lvars: &[String],
) {
    // Build lead header row and justification to match.

    let cols = rsi.vids.len();
    let nexacts = exacts.len();
    if !ctl.clono_print_opt.bu {
        row1.push("#".to_string());
    } else {
        row1.push("#  barcode".to_string());
    }
    justify.push(b'l');
    for x in lvars {
        let x = if x.contains(':') {
            x.before(":")
        } else {
            x.as_str()
        };
        row1.push(x.to_string());
        justify.push(justification(x));
    }

    // Insert main chain row.  Then insert chain info row if we're using CHAIN_SPLIT.

    let mut row = vec![String::new(); row1.len()];
    for j in 0..cols {
        if rsi.chain_descrip[j].contains(&"IGH".to_string())
            || rsi.chain_descrip[j].contains(&"TRB".to_string())
        {
            row.push(bold(&format!("CHAIN {}", j + 1)));
        } else {
            row.push(format!("CHAIN {}", j + 1));
        }
        for _ in 1..rsi.cvars[j].len() {
            row.push("\\ext".to_string());
        }
    }
    rows.push(row);
    let mut row = vec![String::new(); row1.len()];
    for j in 0..cols {
        if !ctl.gen_opt.fold_headers {
            if rsi.chain_descrip[j].contains(&"IGH".to_string())
                || rsi.chain_descrip[j].contains(&"TRB".to_string())
            {
                row.push(bold(&rsi.chain_descrip[j].to_string()));
            } else {
                row.push(rsi.chain_descrip[j].to_string());
            }
        } else if rsi.chain_descrip[j].contains(&"IGH".to_string())
            || rsi.chain_descrip[j].contains(&"TRB".to_string())
        {
            row.push(bold(rsi.chain_descrip[j].before(" ◆ ")));
        } else {
            row.push(rsi.chain_descrip[j].before(" ◆ ").to_string());
        }
        for _ in 1..rsi.cvars[j].len() {
            row.push("\\ext".to_string());
        }
    }
    rows.push(row);
    if ctl.gen_opt.fold_headers {
        let mut row = vec![String::new(); row1.len()];
        for j in 0..cols {
            let mut next = rsi.chain_descrip[j].after(" ◆ ").to_string();
            if next.contains(" ◆ ") {
                next = next.before(" ◆ ").to_string();
            }
            if rsi.chain_descrip[j].contains(&"IGH".to_string())
                || rsi.chain_descrip[j].contains(&"TRB".to_string())
            {
                row.push(bold(&format!("◆ {next}")));
            } else {
                row.push(format!("◆ {next}"));
            }
            for _ in 1..rsi.cvars[j].len() {
                row.push("\\ext".to_string());
            }
        }
        rows.push(row);
        let mut have_last = false;
        for j in 0..cols {
            if rsi.chain_descrip[j].after(" ◆ ").contains(" ◆ ") {
                have_last = true;
            }
        }
        if have_last {
            let mut row = vec![String::new(); row1.len()];
            for j in 0..cols {
                let mut last = String::new();
                if rsi.chain_descrip[j].after(" ◆ ").contains(" ◆ ") {
                    last = rsi.chain_descrip[j].after(" ◆ ").after(" ◆ ").to_string();
                }
                if last.is_empty() {
                    row.push(String::new());
                } else if rsi.chain_descrip[j].contains(&"IGH".to_string())
                    || rsi.chain_descrip[j].contains(&"TRB".to_string())
                {
                    row.push(bold(&format!("◆ {last}")));
                } else {
                    row.push(format!("◆ {last}"));
                }
                for _ in 1..rsi.cvars[j].len() {
                    row.push("\\ext".to_string());
                }
            }
            rows.push(row);
        }
    }

    // Insert divider row (horizontal line across the chains).

    let mut row = vec![String::new(); lvars.len() + 1];
    let mut ncall = 0;
    for j in 0..cols {
        ncall += rsi.cvars[j].len();
    }
    row.append(&mut vec!["\\hline".to_string(); ncall]);
    rows.push(row);

    // Insert position rows.

    *drows = insert_position_rows(ctl, rsi, show_aa, field_types, vars, row1);
    let mut drows2 = drows.clone();
    rows.append(&mut drows2);

    // Insert main per-chain header row.

    let mut row = vec![String::new(); row1.len()];
    for cx in 0..cols {
        let show = &show_aa[cx];
        for j in 0..rsi.cvars[cx].len() {
            if rsi.cvars[cx][j] != *"amino" {
                if drows.is_empty() {
                    row.push(rsi.cvars[cx][j].to_string());
                } else {
                    row.push(String::new());
                }
                continue;
            }
            for u in 0..nexacts {
                let m = rsi.mat[cx][u];
                if m.is_none() {
                    continue;
                }
                let m = m.unwrap();
                let mut n = show.len();
                for k in 1..show.len() {
                    if field_types[cx][k] != field_types[cx][k - 1] && !ctl.gen_opt.nospaces {
                        n += 1;
                    }
                }
                let mut ch = vec![' '; n];
                let amino = &ctl.clono_print_opt.amino;
                let ex = &exact_clonotypes[exacts[u]];
                let x = &ex.share[m];
                let mut cs1 = 0;
                if rsi.cdr1_starts[cx].is_some() {
                    cs1 = rsi.cdr1_starts[cx].unwrap();
                }
                let mut cs2 = 0;
                if rsi.cdr2_starts[cx].is_some() {
                    cs2 = rsi.cdr2_starts[cx].unwrap();
                }
                let mut fs2 = 0;
                if rsi.fr2_starts[cx].is_some() {
                    fs2 = rsi.fr2_starts[cx].unwrap();
                }
                let mut fs3 = 0;
                if rsi.fr3_starts[cx].is_some() {
                    fs3 = rsi.fr3_starts[cx].unwrap();
                }
                let fields = [
                    (
                        "fwr1",
                        rsi.fr1_starts[cx],
                        cs1,
                        rsi.cdr1_starts[cx].is_some(),
                    ),
                    (
                        "fwr2",
                        fs2,
                        cs2,
                        rsi.fr2_starts[cx].is_some() && rsi.cdr2_starts[cx].is_some(),
                    ),
                    (
                        "fwr3",
                        fs3,
                        rsi.cdr3_starts[cx],
                        rsi.fr3_starts[cx].is_some(),
                    ),
                    (
                        "cdr1",
                        cs1,
                        fs2,
                        rsi.cdr1_starts[cx].is_some() && rsi.fr2_starts[cx].is_some(),
                    ),
                    (
                        "cdr2",
                        cs2,
                        fs3,
                        rsi.cdr2_starts[cx].is_some() && rsi.fr3_starts[cx].is_some(),
                    ),
                    (
                        "cdr3",
                        rsi.cdr3_starts[cx],
                        rsi.cdr3_starts[cx] + x.cdr3_aa.len() * 3,
                        true,
                    ),
                    (
                        "fwr4",
                        rsi.cdr3_starts[cx] + x.cdr3_aa.len() * 3,
                        rsi.seq_del_lens[cx] - 1,
                        true,
                    ),
                ];
                for field in fields {
                    if field.3 && field.1 <= field.2 && amino.contains(&field.0.to_string()) {
                        let cs1 = field.1 / 3;
                        let mut ch_start = 0;
                        for (k, (&s, &t)) in show.iter().zip(field_types[cx].iter()).enumerate() {
                            if k > 0 && t != field_types[cx][k - 1] && !ctl.gen_opt.nospaces {
                                ch_start += 1;
                            }
                            if s == cs1 {
                                break;
                            }
                            ch_start += 1;
                        }
                        let q = (field.2 - field.1) / 3;

                        // Catch an error condition that has happened a few times.

                        if ch_start + q > ch.len() {
                            let mut ds = Vec::<&str>::new();
                            for clone in &ex.clones {
                                let li = clone[m].dataset_index;
                                ds.push(ctl.origin_info.dataset_id[li].as_str());
                            }
                            unique_sort(&mut ds);
                            let fields_msg = format!(
                                "fields[z].0 = {}, fields[z].1 = {}, fields[z].2 = {},",
                                field.0, field.1, field.2,
                            );
                            panic!(
                                "Internal error, out of range in \
                                build_table_stuff, CDR3 = {}, datasets = {},\n\
                                ch_start = {}, q = {}, ch.len() = {},\n\
                                {}\n\
                                show = {},\n\
                                field_types = {}.",
                                x.cdr3_aa,
                                ds.iter().format(","),
                                ch_start,
                                q,
                                ch.len(),
                                fields_msg,
                                show.iter().format(","),
                                field_types[cx].iter().format(","),
                            );
                        }

                        // Do the work.

                        let t = field.0.to_ascii_uppercase();
                        let t = t.as_bytes();
                        // The second form of this gets converted below.
                        let c = if t[0] == b'C' || !ctl.gen_opt.nospaces {
                            "═"
                        } else {
                            "┅"
                        };
                        let mut s = String::new();
                        if q >= 4 {
                            let left = (q - 3) / 2;
                            let right = q - left - 4;
                            s.reserve((left + right) * c.len() + t.len());
                            for _ in 0..left {
                                s += c;
                            }
                            s += strme(t);
                            for _ in 0..right {
                                s += c;
                            }
                        } else if q == 3 {
                            s += strme(&t[0..1]);
                            s += strme(&t[2..4]);
                        } else if q == 2 {
                            s += strme(&t[0..1]);
                            s += strme(&t[3..4]);
                        } else if q == 1 {
                            s += strme(&t[3..4]);
                        }
                        let schars = s.chars().take(q).collect::<Vec<_>>();
                        ch[ch_start..(q + ch_start)].copy_from_slice(&schars);
                    }
                }
                let mut s = String::new();
                for c in ch {
                    s.push(c);
                }
                s = s.trim_end().to_string();

                // Convert ┅ to ═ in different color.

                if ctl.gen_opt.nospaces {
                    let chars = s.chars().collect::<Vec<_>>();
                    s.clear();
                    for (i, &char_i) in chars.iter().enumerate() {
                        if (i == 0 || (i > 0 && chars[i - 1] != '┅')) && char_i == '┅' {
                            s += "[01m[38;5;198m";
                            s.push('═');
                        } else if char_i == '┅' {
                            s.push('═');
                            if i == chars.len() - 1 || (i < chars.len() - 1 && chars[i + 1] != '┅')
                            {
                                s += "[0m";
                            }
                        } else {
                            s.push(char_i);
                        }
                    }
                    s = s.replace("FWR1", "[01m[38;5;198mFWR1[0m");
                    s = s.replace("FWR2", "[01m[38;5;198mFWR2[0m");
                    s = s.replace("FWR3", "[01m[38;5;198mFWR3[0m");
                    s = s.replace("FWR4", "[01m[38;5;198mFWR4[0m");
                }

                // Save.

                row.push(s);
                break;
            }
        }
    }
    rows.push(row);
}
