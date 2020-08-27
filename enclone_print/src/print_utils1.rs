// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use ansi_escape::*;
use enclone_core::defs::*;
use io_utils::*;
use itertools::*;
use std::cmp::max;
use std::collections::HashMap;
use std::io::Write;
use string_utils::*;
use tables::*;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn make_table(
    ctl: &EncloneControl,
    rows: &mut Vec<Vec<String>>,
    justify: &Vec<u8>,
    mlog: &Vec<u8>,
    logz: &mut String,
) {
    // In plain mode, strip escape characters.

    if !ctl.pretty {
        for i in 0..rows.len() {
            for j in 0..rows[i].len() {
                let mut x = Vec::<u8>::new();
                let mut escaped = false;
                let s = rows[i][j].as_bytes();
                for l in 0..s.len() {
                    if s[l] == b'' {
                        escaped = true;
                    }
                    if escaped {
                        if s[l] == b'm' {
                            escaped = false;
                        }
                        continue;
                    }
                    x.push(s[l]);
                }
                rows[i][j] = stringme(&x);
            }
        }
    }

    // Make table.

    let log0 = stringme(&mlog);
    let mut log = String::new();
    if ctl.debug_table_printing {
        for i in 0..rows.len() {
            println!("");
            for j in 0..rows[i].len() {
                println!(
                    "row = {}, col = {}, entry = {}, vis width = {}",
                    i,
                    j,
                    rows[i][j],
                    visible_width(&rows[i][j])
                );
            }
        }
        println!("");
    }
    print_tabular_vbox(
        &mut log,
        &rows,
        2,
        &justify,
        ctl.debug_table_printing,
        false,
    );
    if ctl.debug_table_printing {
        println!("{}", log);
    }
    let mut cs = vec![Vec::<char>::new(); rows.len() + 2];
    let mut row = 0;
    for c in log.chars() {
        if c == '\n' {
            row += 1;
        } else {
            cs[row].push(c);
        }
    }
    log = log0;

    // Process each row.

    for i in 0..cs.len() {
        for j in 0..cs[i].len() {
            log.push(cs[i][j]);
        }
        log.push('\n');
    }

    // Make some character substitutions.

    let mut barcode = false;
    let mut header = false;
    let mut x = Vec::<char>::new();
    for c in log.chars() {
        x.push(c);
    }
    let mut j = 0;
    while j < x.len() {
        // DEFAULT
        /*
        const TEXTCOLOR: usize = 200;
        const BACKGROUND: usize = 229;
        */

        // NOT CRAZY
        /*
        const TEXTCOLOR: usize = 200;
        const BACKGROUND: usize = 225;
        */

        // NEW COLOR SCHEME
        const TEXTCOLOR: usize = 18;
        const BACKGROUND: usize = 255;

        let c = x[j];

        // $ is a placeholder for •, and $ is only in barcodes line if PER_CELL is specified.
        // In plain mode, we just make the substitution, whereas in fancy mode, we change the
        // text and background color for the entire line.
        // *** bullets now off ***
        if c == '$' {
            if ctl.pretty {
                *logz += &format!("[38;5;{}m[48;5;{}m ", TEXTCOLOR, BACKGROUND);
                barcode = true;
            } else {
                logz.push('•');
            }

        // In a barcode line, elide end escapes.  Not exactly sure how these get here.
        } else if barcode && c == '' && x[j + 1] == '[' && x[j + 2] == '0' && x[j + 3] == 'm' {
            j += 3;

        // In a barcode line, hop around │ symbols, which should not be colorized.
        } else if barcode && c == '│' && x[j + 1] != '\n' {
            // *logz += "[0m│";
            *logz += &format!("[0m[48;5;{}m│", BACKGROUND);
            *logz += &format!("[38;5;{}m[48;5;{}m", TEXTCOLOR, BACKGROUND);
        } else if barcode && c == '│' && x[j + 1] == '\n' {
            *logz += "[0m│";
            // *logz += &format!("[0m[48;5;{}m│[0m", BACKGROUND);
            barcode = false;

        // Do similar things for header line, but bold the line instead.
        } else if c == '#' {
            if ctl.pretty {
                *logz += &format!("[01m#");
                header = true;
            } else {
                logz.push('#');
            }

        // In a header line, elide end escapes.  Not exactly sure how these get here.
        // (did not check to see if this does anything)
        } else if header && c == '' && x[j + 1] == '[' && x[j + 2] == '0' && x[j + 3] == 'm' {
            j += 3;

        // In a header line, hop around │ symbols, which should not be colorized.
        } else if header && c == '│' && x[j + 1] != '\n' {
            *logz += "[0m│";
            *logz += &format!("[01m");
        } else if header && c == '│' && x[j + 1] == '\n' {
            *logz += "[0m│";
            header = false;

        // Otherwise just save the character.
        } else {
            logz.push(c);
        }
        j += 1;
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn print_digit(p: usize, i: usize, digits: usize, ds: &mut String) {
    if digits == 1 {
        *ds += &format!("{}", p);
    } else if digits == 2 {
        if i == 0 {
            if p >= 10 {
                *ds += &format!("{}", p / 10);
            } else {
                ds.push(' ');
            }
        } else {
            *ds += &format!("{}", p % 10);
        }
    } else {
        if i == 0 {
            if p >= 100 {
                *ds += &format!("{}", p / 100);
            } else {
                ds.push(' ');
            }
        } else if i == 1 {
            if p >= 10 {
                *ds += &format!("{}", (p % 100) / 10);
            } else {
                ds.push(' ');
            }
        } else {
            *ds += &format!("{}", p % 10);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Determine the number of digits in a nonnegative integer.

pub fn ndigits(n: usize) -> usize {
    let mut d = 1;
    let mut x = n;
    while x >= 10 {
        d += 1;
        x /= 10;
    }
    d
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn make_diff_row(
    ctl: &EncloneControl,
    rsi: &ColInfo,
    cols: usize,
    diff_pos: usize,
    drows: &Vec<Vec<String>>,
    row1: &mut Vec<String>,
    rows: &mut Vec<Vec<String>>,
) {
    let mut xrow = vec!["".to_string(); row1.len()];
    let mut xrow_filled = vec![false; row1.len()];
    for cx in 0..cols {
        for j in 0..rsi.cvars[cx].len() {
            if rsi.cvars[cx][j] != "amino".to_string() {
                xrow.push(rsi.cvars[cx][j].to_string());
                xrow_filled.push(true);
            } else {
                xrow.push("".to_string());
                xrow_filled.push(false);
            }
        }
    }
    let nc = row1.len();
    if !drows.is_empty() {
        let mut ncall = 0;
        for j in 0..cols {
            for z in 0..rsi.cvars[j].len() {
                let mut c = Vec::<Vec<u8>>::new();
                let mut start = 5 + drows.len();
                if drows.len() >= 1 {
                    start += 3;
                }
                if ctl.clono_print_opt.sum {
                    start += 1;
                }
                if ctl.clono_print_opt.mean {
                    start += 1;
                }
                for k in start..rows.len() {
                    if rows[k][0].contains("$") {
                        continue;
                    }
                    if rows[k][ncall + z + nc].len() > 0 {
                        c.push(rows[k][ncall + z + nc].as_bytes().to_vec());
                    }
                }

                // Package characters with ANSI escape codes that come before them.
                // The is a dorky way of identifying codons that are different, by
                // virtue of them being shown as colored amino acids.

                let mut c2 = Vec::<Vec<Vec<u8>>>::new();
                for i in 0..c.len() {
                    c2.push(package_characters_with_escapes(&c[i]));
                }

                // Proceed.

                if (rsi.cvars[j][z] != "amino" && rsi.cvars[j][z] != "var") || c.len() == 0 {
                    row1.push("".to_string());
                    continue;
                }
                let mut dots = Vec::<u8>::new();
                for m in 0..c2[0].len() {
                    let mut digits_or_blanks = true;
                    for l in 0..c2.len() {
                        if c2[l].is_empty() {
                            // needed?
                            continue;
                        }
                        if !(c2[l][m] == b" ".to_vec()
                            || (c2[l][m] >= b"0".to_vec() && c2[l][m] <= b"9".to_vec()))
                        {
                            digits_or_blanks = false;
                        }
                        if c2[l][m].contains(&b' ') {
                            digits_or_blanks = true;
                        }
                    }
                    let mut same = true;
                    for l in 1..c2.len() {
                        if c2[l].is_empty() {
                            // needed?
                            continue;
                        }
                        if c2[l][m] != c2[0][m] {
                            same = false;
                        }
                    }
                    let mut sep = true;
                    for l in 0..c2.len() {
                        if c2[l].is_empty() {
                            // needed?
                            continue;
                        }
                        if c2[l][m] != b"|" {
                            sep = false;
                        }
                    }
                    if sep {
                    } else if digits_or_blanks {
                        dots.push(b' ');
                    } else if same {
                        dots.push(b'.');
                    } else {
                        dots.push(b'x');
                    }
                }
                row1.push(format!("{}", strme(&dots)));
            }
            ncall += rsi.cvars[j].len();
        }
        if !drows.is_empty() {
            for i in 0..row1.len() {
                if xrow_filled[i] {
                    row1[i] = xrow[i].clone();
                }
            }
        }
        rows[diff_pos] = row1.to_vec();
    } else {
        if !drows.is_empty() {
            for i in 0..row1.len() {
                if xrow_filled[i] {
                    row1[i] = xrow[i].clone();
                }
            }
        }
        for i in 0..row1.len() {
            rows[diff_pos - 1][i] = row1[i].clone();
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Start to generate parseable output, warn about multi-donor clonotypes, and other things.

pub fn start_gen(
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    refdata: &RefData,
    rsi: &ColInfo,
    out_data: &mut Vec<HashMap<String, String>>,
    mut mlog: &mut Vec<u8>,
) {
    let pcols_sort = &ctl.parseable_opt.pcols_sort;
    macro_rules! speak {
        ($u:expr, $var:expr, $val:expr) => {
            if ctl.parseable_opt.pout.len() > 0 {
                if pcols_sort.is_empty() || bin_member(&pcols_sort, &$var.to_string()) {
                    out_data[$u].insert($var.to_string(), $val);
                }
            }
        };
    }
    macro_rules! speakc {
        ($u:expr, $col:expr, $var:expr, $val:expr) => {
            if ctl.parseable_opt.pout.len() > 0 && $col + 1 <= ctl.parseable_opt.pchains {
                let varc = format!("{}{}", $var, $col + 1);
                if pcols_sort.is_empty() || bin_member(&pcols_sort, &varc) {
                    out_data[$u].insert(varc, format!("{}", $val));
                }
            }
        };
    }
    let nexacts = exacts.len();
    let mut n = 0;
    for u in 0..nexacts {
        n += exact_clonotypes[exacts[u]].ncells();
    }
    if ctl.parseable_opt.pout.len() > 0 {
        *out_data = vec![HashMap::<String, String>::new(); nexacts];
    }
    let cols = rsi.vids.len();
    let mut ncells = 0;
    for u in 0..exacts.len() {
        ncells += exact_clonotypes[exacts[u]].ncells();
    }
    for u in 0..exacts.len() {
        speak!(u, "nchains", format!("{}", cols));
        speak!(u, "clonotype_ncells", format!("{}", ncells));
        let mut bc = Vec::<String>::new();
        for x in exact_clonotypes[exacts[u]].clones.iter() {
            bc.push(x[0].barcode.clone());
        }
        bc.sort();
        speak!(u, "barcodes", format!("{}", bc.iter().format(",")));
        for d in ctl.origin_info.dataset_list.iter() {
            if d.len() > 0 {
                let mut bc = Vec::<String>::new();
                for i in 0..exact_clonotypes[exacts[u]].clones.len() {
                    let q = &exact_clonotypes[exacts[u]].clones[i];
                    if ctl.origin_info.dataset_id[q[0].dataset_index] == *d {
                        bc.push(q[0].barcode.clone());
                    }
                }
                speak!(
                    u,
                    &format!("{}_barcodes", d),
                    format!("{}", bc.iter().format(","))
                );
            }
        }
        if ctl.parseable_opt.pbarcode {
            let mut bc = Vec::<String>::new();
            for x in exact_clonotypes[exacts[u]].clones.iter() {
                bc.push(x[0].barcode.clone());
            }
            speak!(u, "barcode", format!("{}", bc.iter().format(";")));
            for d in ctl.origin_info.dataset_list.iter() {
                if d.len() > 0 {
                    let mut bc = Vec::<String>::new();
                    for i in 0..exact_clonotypes[exacts[u]].clones.len() {
                        let q = &exact_clonotypes[exacts[u]].clones[i];
                        if ctl.origin_info.dataset_id[q[0].dataset_index] == *d {
                            bc.push(q[0].barcode.clone());
                        } else {
                            bc.push("".to_string());
                        }
                    }
                    speak!(
                        u,
                        &format!("{}_barcode", d),
                        format!("{}", bc.iter().format(";"))
                    );
                }
            }
        }
        for cx in 0..cols {
            let vid = rsi.vids[cx];
            speakc!(u, cx, "v_name", refdata.name[vid]);
            speakc!(u, cx, "v_id", refdata.id[vid]);
            let did = rsi.dids[cx];
            if did.is_some() {
                let did = did.unwrap();
                speakc!(u, cx, "d_name", refdata.name[did]);
                speakc!(u, cx, "d_id", refdata.id[did]);
            }
            let jid = rsi.jids[cx];
            speakc!(u, cx, "j_name", refdata.name[jid]);
            speakc!(u, cx, "j_id", refdata.id[jid]);
        }
    }

    // Start to print the clonotype.

    let mut donors = Vec::<usize>::new();
    for u in 0..exacts.len() {
        let ex = &exact_clonotypes[exacts[u]];
        for m in 0..ex.clones.len() {
            if ex.clones[m][0].donor_index.is_some() {
                let d = ex.clones[m][0].donor_index.unwrap();
                if ctl.origin_info.donor_list[d].len() > 0 {
                    donors.push(d);
                }
            }
        }
    }
    unique_sort(&mut donors);
    fwriteln!(&mut mlog, "CLONOTYPE = {} CELLS", n);
    if donors.len() > 1 && !ctl.gen_opt.nwarn {
        if ctl.pretty {
            // emit_red_escape(&mut mlog);
            // what is below is a brighter red
            mlog.append(&mut b"[38;5;9m".to_vec());
        }
        fwrite!(&mut mlog, "██");
        if ctl.pretty {
            emit_end_escape(&mut mlog);
        }
        fwriteln!(
            &mut mlog,
            " WARNING: This clonotype contains cells from multiple donors."
        );
        let mut donor_names = Vec::<String>::new();
        for i in 0..donors.len() {
            donor_names.push(ctl.origin_info.donor_list[donors[i]].clone());
        }
        fwriteln!(&mut mlog, "donors = {}", donor_names.iter().format(","));
        fwriteln!(&mut mlog, "datasets in which these donors appear:");
        for i in 0..donors.len() {
            let mut datasets = Vec::<String>::new();
            for u in 0..nexacts {
                let ex = &exact_clonotypes[exacts[u]];
                for l in 0..ex.clones.len() {
                    if ex.clones[l][0].donor_index.is_some() {
                        if ex.clones[l][0].donor_index.unwrap() == donors[i] {
                            datasets.push(
                                ctl.origin_info.dataset_id[ex.clones[l][0].dataset_index].clone(),
                            );
                        }
                    }
                }
            }
            unique_sort(&mut datasets);
            // This message is pretty flaky in the case where bc has been specified in META.
            fwriteln!(
                &mut mlog,
                "donor {}: {}",
                i + 1,
                datasets.iter().format(",")
            );
        }
    }

    // Print barcodes.

    if ctl.clono_print_opt.barcodes {
        let mut bc = Vec::<String>::new();
        for u in 0..nexacts {
            let ex = &exact_clonotypes[exacts[u]];
            for l in 0..ex.clones.len() {
                bc.push(ex.clones[l][0].barcode.clone());
            }
        }
        unique_sort(&mut bc);
        fwriteln!(&mut mlog, "• {}", bc.iter().format(","));
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn insert_position_rows(
    rsi: &ColInfo,
    show_aa: &Vec<Vec<usize>>,
    field_types: &Vec<Vec<u8>>,
    vars: &Vec<Vec<usize>>,
    row1: &Vec<String>,
) -> Vec<Vec<String>> {
    let cols = rsi.cdr3_starts.len();
    let mut drows = Vec::<Vec<String>>::new();
    let mut digits = 0;
    for zpass in 1..=2 {
        if zpass == 2 {
            drows = vec![vec![String::new(); row1.len()]; digits];
        }
        for cx in 0..cols {
            for m in 0..rsi.cvars[cx].len() {
                if zpass == 1 {
                    if rsi.cvars[cx][m] == "amino".to_string() {
                        for p in show_aa[cx].iter() {
                            digits = max(digits, ndigits(*p));
                        }
                    } else if rsi.cvars[cx][m] == "var".to_string() {
                        for p in vars[cx].iter() {
                            digits = max(digits, ndigits(*p));
                        }
                    }
                } else {
                    for i in 0..digits {
                        if rsi.cvars[cx][m] == "amino".to_string() {
                            let mut ds = String::new();
                            for (j, p) in show_aa[cx].iter().enumerate() {
                                if j > 0 && field_types[cx][j] != field_types[cx][j - 1] {
                                    ds += " ";
                                }
                                print_digit(*p, i, digits, &mut ds);
                            }
                            drows[i].push(ds);
                        } else if rsi.cvars[cx][m] == "var".to_string() {
                            let mut ds = String::new();
                            for p in vars[cx].iter() {
                                print_digit(*p, i, digits, &mut ds);
                            }
                            drows[i].push(ds);
                        } else {
                            drows[i].push(String::new());
                        }
                    }
                }
            }
        }
    }
    drows
}
