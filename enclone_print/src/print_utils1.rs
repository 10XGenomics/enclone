// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use amino::codon_to_aa;
use ansi_escape::{
    emit_bold_escape, emit_eight_bit_color_escape, emit_end_escape, emit_red_escape,
};
use enclone_core::cell_color::CellColor;
use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype, GexInfo, TigData1, POUT_SEP};
use enclone_core::print_tools::{color_by_property, emit_codon_color_escape};
use enclone_vars::decode_arith;
use expr_tools::vars_of_node;
use io_utils::{fwrite, fwriteln};
use itertools::Itertools;
use std::cmp::max;
use std::collections::HashMap;
use std::io::Write;
use string_utils::{stringme, strme};
use tables::{print_tabular_vbox, visible_width};
use vector_utils::{bin_member, lower_bound1_3, meet_size, unique_sort, upper_bound1_3, VecUtils};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn test_internal_error_seq(seq: &[u8], dna: &[u8], cdr3: &str) -> Result<(), String> {
    let mut found = false;
    for i in 0..seq.len() {
        if seq[i..].starts_with(dna) {
            found = true;
        }
    }
    if !found {
        return Err(format!(
            "\nInternal error, failed to find {}, CDR3 = {}.\n",
            strme(dna),
            cdr3
        ));
    }
    Ok(())
}

pub fn get_cdr1(x: &TigData1, left: i64, right: i64) -> Option<String> {
    let left = left * 3;
    let right = right * 3;
    if x.cdr1_start.is_some()
        && x.fr2_start.is_some()
        && x.cdr1_start.unwrap() <= x.fr2_start.unwrap()
    {
        let mut dna = Vec::<u8>::new();
        if x.cdr1_start.unwrap() as i64 - left >= 0
            && x.cdr1_start.unwrap() as i64 - left < x.seq_del_amino.len() as i64
            && x.fr2_start.unwrap() as i64 + right > 0
            && x.fr2_start.unwrap() as i64 + right <= x.seq_del_amino.len() as i64
        {
            for p in x.cdr1_start.unwrap() as i64 - left..x.fr2_start.unwrap() as i64 + right {
                let p = p as usize;
                for j in 0..x.ins.len() {
                    if x.ins[j].0 == p {
                        let mut z = x.ins[j].1.clone();
                        dna.append(&mut z);
                    }
                }
                if x.seq_del_amino[p] != b'-' {
                    dna.push(x.seq_del_amino[p]);
                }
            }
            test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa).unwrap();
            return Some(stringme(&dna));
        }
    }
    None
}

pub fn get_cdr2(x: &TigData1, left: i64, right: i64) -> Option<String> {
    let left = left * 3;
    let right = right * 3;
    if x.cdr2_start.is_some()
        && x.fr3_start.is_some()
        && x.cdr2_start.unwrap() <= x.fr3_start.unwrap()
    {
        let mut dna = Vec::<u8>::new();
        if x.cdr2_start.unwrap() as i64 - left >= 0
            && x.cdr2_start.unwrap() as i64 - left < x.seq_del_amino.len() as i64
            && x.fr3_start.unwrap() as i64 + right > 0
            && x.fr3_start.unwrap() as i64 + right <= x.seq_del_amino.len() as i64
        {
            for p in x.cdr2_start.unwrap() as i64 - left..x.fr3_start.unwrap() as i64 + right {
                let p = p as usize;
                for j in 0..x.ins.len() {
                    if x.ins[j].0 == p {
                        let mut z = x.ins[j].1.clone();
                        dna.append(&mut z);
                    }
                }
                if x.seq_del_amino[p] != b'-' {
                    dna.push(x.seq_del_amino[p]);
                }
            }
            test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa).unwrap();
            return Some(stringme(&dna));
        }
    }
    None
}

pub fn get_cdr3(x: &TigData1, left: i64, right: i64) -> Option<String> {
    let left = left * 3;
    let right = right * 3;
    let mut dna = Vec::<u8>::new();
    if x.cdr3_start as i64 - left >= 0
        && x.cdr3_start as i64 - left < x.seq_del_amino.len() as i64
        && x.cdr3_start as i64 + 3 * x.cdr3_aa.len() as i64 + right > 0
        && x.cdr3_start as i64 + 3 * x.cdr3_aa.len() as i64 + right <= x.seq_del_amino.len() as i64
    {
        for p in
            x.cdr3_start as i64 - left..x.cdr3_start as i64 + 3 * x.cdr3_aa.len() as i64 + right
        {
            let p = p as usize;
            for j in 0..x.ins.len() {
                if x.ins[j].0 == p {
                    let mut z = x.ins[j].1.clone();
                    dna.append(&mut z);
                }
            }
            if x.seq_del_amino[p] != b'-' {
                dna.push(x.seq_del_amino[p]);
            }
        }
        test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa).unwrap();
        return Some(stringme(&dna));
    }
    None
}

pub fn get_fwr1(x: &TigData1) -> Option<String> {
    if x.cdr1_start.is_some() && x.fr1_start <= x.cdr1_start.unwrap() {
        let mut dna = Vec::<u8>::new();
        for p in x.fr1_start..x.cdr1_start.unwrap() {
            for j in 0..x.ins.len() {
                if x.ins[j].0 == p {
                    let mut z = x.ins[j].1.clone();
                    dna.append(&mut z);
                }
            }
            if x.seq_del_amino[p] != b'-' {
                dna.push(x.seq_del_amino[p]);
            }
        }
        test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa).unwrap();
        return Some(stringme(&dna));
    }
    None
}

pub fn get_fwr2(x: &TigData1) -> Option<String> {
    if x.fr2_start.unwrap() <= x.cdr2_start.unwrap() {
        let mut dna = Vec::<u8>::new();
        for p in x.fr2_start.unwrap()..x.cdr2_start.unwrap() {
            for j in 0..x.ins.len() {
                if x.ins[j].0 == p {
                    let mut z = x.ins[j].1.clone();
                    dna.append(&mut z);
                }
            }
            if x.seq_del_amino[p] != b'-' {
                dna.push(x.seq_del_amino[p]);
            }
        }
        test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa).unwrap();
        return Some(stringme(&dna));
    }
    None
}

pub fn get_fwr3(x: &TigData1) -> Option<String> {
    if x.fr3_start.is_some() && x.fr3_start.unwrap() <= x.cdr3_start - x.ins_len() {
        let mut dna = Vec::<u8>::new();
        for p in x.fr3_start.unwrap()..x.cdr3_start - x.ins_len() {
            for j in 0..x.ins.len() {
                if x.ins[j].0 == p {
                    let mut z = x.ins[j].1.clone();
                    dna.append(&mut z);
                }
            }
            if x.seq_del_amino[p] != b'-' {
                dna.push(x.seq_del_amino[p]);
            }
        }
        test_internal_error_seq(&x.seq, &dna, &x.cdr3_aa).unwrap();
        return Some(stringme(&dna));
    }
    None
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn compute_field_types(
    ctl: &EncloneControl,
    rsi: &ColInfo,
    show_aa: &Vec<Vec<usize>>,
) -> Vec<Vec<u8>> {
    let cols = rsi.mat.len();
    let mut field_types = vec![Vec::new(); cols];
    for cx in 0..cols {
        let mut ft = vec![0_u8; show_aa[cx].len()];
        let cs1 = rsi.cdr1_starts[cx];
        let cs2 = rsi.cdr2_starts[cx];
        let cs3 = rsi.cdr3_starts[cx];
        let n3 = rsi.cdr3_lens[cx];
        let fs1 = rsi.fr1_starts[cx];
        let fs2 = rsi.fr2_starts[cx];
        let fs3 = rsi.fr3_starts[cx];
        let show_cdr1 = cs1.is_some()
            && fs2.is_some()
            && cs1.unwrap() <= fs2.unwrap()
            && ctl.clono_print_opt.amino.contains(&"cdr1".to_string());
        let show_cdr2 = cs2.is_some()
            && fs3.is_some()
            && cs2.unwrap() <= fs3.unwrap()
            && ctl.clono_print_opt.amino.contains(&"cdr2".to_string());
        let show_cdr3 = ctl.clono_print_opt.amino.contains(&"cdr3".to_string());
        let show_fwr1 = cs1.is_some()
            && rsi.fr1_starts[cx] <= cs1.unwrap()
            && ctl.clono_print_opt.amino.contains(&"fwr1".to_string());
        let show_fwr2 = fs2.is_some()
            && cs2.is_some()
            && fs2.unwrap() <= cs2.unwrap()
            && ctl.clono_print_opt.amino.contains(&"fwr2".to_string());
        let show_fwr3 = fs3.is_some()
            && fs3.unwrap() <= rsi.cdr3_starts[cx]
            && ctl.clono_print_opt.amino.contains(&"fwr3".to_string());
        let show_fwr4 = ctl.clono_print_opt.amino.contains(&"fwr4".to_string());
        for (j, p) in show_aa[cx].iter().enumerate() {
            if show_cdr1 && *p >= cs1.unwrap() / 3 && *p < fs2.unwrap() / 3 {
                ft[j] = 1;
            } else if show_cdr2 && *p >= cs2.unwrap() / 3 && *p < fs3.unwrap() / 3 {
                ft[j] = 2;
            } else if show_cdr3 && *p >= cs3 / 3 && *p < cs3 / 3 + n3 {
                ft[j] = 3;
            } else if show_fwr1 && *p >= fs1 / 3 && *p < cs1.unwrap() / 3 {
                ft[j] = 4;
            } else if show_fwr2 && *p >= fs2.unwrap() / 3 && *p < cs2.unwrap() / 3 {
                ft[j] = 5;
            } else if show_fwr3 && *p >= fs3.unwrap() / 3 && *p < rsi.cdr3_starts[cx] / 3 {
                ft[j] = 6;
            } else if show_fwr4 && *p >= cs3 / 3 + n3 {
                ft[j] = 7;
            }
        }
        field_types[cx] = ft;
    }
    field_types
}

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

    let log0 = stringme(mlog);
    let mut log = String::new();
    if ctl.debug_table_printing {
        for i in 0..rows.len() {
            println!();
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
        println!();
    }
    print_tabular_vbox(&mut log, rows, 2, justify, ctl.debug_table_printing, false);
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
                *logz += &"[01m#".to_string();
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
            *logz += &"[01m".to_string();
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
    } else if i == 0 {
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

// Start to generate parseable output, warn about multi-donor clonotypes, and other things.

pub fn start_gen(
    ctl: &EncloneControl,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    out_data: &mut Vec<HashMap<String, String>>,
    mut mlog: &mut Vec<u8>,
    extra_args: &Vec<String>,
) {
    let pcols_sort = &ctl.parseable_opt.pcols_sort;
    macro_rules! speak {
        ($u:expr, $var:expr, $val:expr) => {
            if ctl.parseable_opt.pout.len() > 0 || extra_args.len() > 0 {
                if pcols_sort.is_empty()
                    || bin_member(&pcols_sort, &$var.to_string())
                    || bin_member(&extra_args, &$var.to_string())
                {
                    out_data[$u].insert($var.to_string(), $val);
                }
            }
        };
    }
    let nexacts = exacts.len();
    let mut n = 0;
    for u in 0..nexacts {
        n += exact_clonotypes[exacts[u]].ncells();
    }
    if !ctl.parseable_opt.pout.is_empty() || !extra_args.is_empty() {
        *out_data = vec![HashMap::<String, String>::new(); nexacts];
    }
    for u in 0..exacts.len() {
        let mut bc = Vec::<String>::new();
        for x in exact_clonotypes[exacts[u]].clones.iter() {
            bc.push(x[0].barcode.clone());
        }
        bc.sort();
        speak!(u, "barcodes", format!("{}", bc.iter().format(",")));
        for d in ctl.origin_info.dataset_list.iter() {
            if !d.is_empty() {
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
            speak!(u, "barcode", format!("{}", bc.iter().format(POUT_SEP)));
            for d in ctl.origin_info.dataset_list.iter() {
                if !d.is_empty() {
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
                        format!("{}", bc.iter().format(POUT_SEP))
                    );
                }
            }
        }
    }

    // Start to print the clonotype.

    let mut donors = Vec::<usize>::new();
    for u in 0..exacts.len() {
        let ex = &exact_clonotypes[exacts[u]];
        for m in 0..ex.clones.len() {
            if ex.clones[m][0].donor_index.is_some() {
                let d = ex.clones[m][0].donor_index.unwrap();
                if !ctl.origin_info.donor_list[d].is_empty() {
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
            emit_end_escape(mlog);
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
                    if ex.clones[l][0].donor_index.is_some()
                        && ex.clones[l][0].donor_index.unwrap() == donors[i]
                    {
                        datasets.push(
                            ctl.origin_info.dataset_id[ex.clones[l][0].dataset_index].clone(),
                        );
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
    ctl: &EncloneControl,
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
                    if rsi.cvars[cx][m] == *"amino" {
                        for p in show_aa[cx].iter() {
                            digits = max(digits, ndigits(*p));
                        }
                    } else if rsi.cvars[cx][m] == *"var" {
                        for p in vars[cx].iter() {
                            digits = max(digits, ndigits(*p));
                        }
                    }
                } else {
                    for i in 0..digits {
                        if rsi.cvars[cx][m] == *"amino" {
                            let mut ds = String::new();
                            for (j, p) in show_aa[cx].iter().enumerate() {
                                if j > 0
                                    && field_types[cx][j] != field_types[cx][j - 1]
                                    && !ctl.gen_opt.nospaces
                                {
                                    ds += " ";
                                }
                                print_digit(*p, i, digits, &mut ds);
                            }
                            drows[i].push(ds);
                        } else if rsi.cvars[cx][m] == *"var" {
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn color_codon(
    ctl: &EncloneControl,
    seq_amino: &Vec<u8>,
    ref_diff_pos: &Vec<Vec<Vec<usize>>>,
    x: &Vec<(usize, u8, u32)>,
    col: usize,
    mid: usize,
    p: usize,
    u: usize,
    last_color: &mut String,
    last: bool,
    cdr3_con: &Vec<Vec<u8>>,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
) -> Vec<u8> {
    let mut log = Vec::<u8>::new();
    let codon = &seq_amino[3 * p..3 * p + 3];
    let aa = codon_to_aa(codon);
    if ctl.gen_opt.color == *"codon" || ctl.gen_opt.color == *"codon-diffs" {
        let mut diff = false;
        if !ref_diff_pos.is_empty() && ctl.gen_opt.color == *"codon-diffs" {
            for j in 0..3 {
                if bin_member(&ref_diff_pos[col][u], &(3 * p + j)) {
                    diff = true;
                }
            }
            let cdr3_start = exact_clonotypes[exacts[u]].share[mid].cdr3_start;
            let cdr3 = &exact_clonotypes[exacts[u]].share[mid].cdr3_dna.as_bytes();
            if 3 * p >= cdr3_start && 3 * p < cdr3_start + cdr3.len() {
                let cdr3_con = &cdr3_con[col];
                for j in 0..3 {
                    let cp = 3 * p - cdr3_start + j;
                    if cdr3[cp] != cdr3_con[cp] {
                        diff = true;
                    }
                }
            }
        }
        if !ref_diff_pos.is_empty() && !diff && ctl.gen_opt.color == *"codon-diffs" {
            log.append(&mut b"[01m[38;5;254m".to_vec());
        } else {
            emit_codon_color_escape(codon, &mut log);
        }
        log.push(aa);
        emit_end_escape(&mut log);
    } else if ctl.gen_opt.color == *"property" {
        color_by_property(&[aa], &mut log);
    } else {
        let (low, high) = (lower_bound1_3(x, &p), upper_bound1_3(x, &p));
        let (mut total, mut this) = (0.0, 0.0);
        for u in low..high {
            total += x[u as usize].2 as f64;
            if x[u as usize].1 == aa {
                this = x[u as usize].2 as f64;
            }
        }
        let mut color = "black".to_string();
        if total > 0.0 && 100.0 * this / total <= ctl.gen_opt.color_by_rarity_pc {
            if this == 0.0 {
                color = "red".to_string();
            } else {
                color = "blue".to_string();
            }
        }
        if color != *last_color {
            if color == *"black" {
                emit_end_escape(&mut log);
            } else {
                if color == *"red" {
                    emit_red_escape(&mut log);
                } else {
                    emit_eight_bit_color_escape(&mut log, 6);
                }
                emit_bold_escape(&mut log);
            }
            *last_color = color;
        }
        fwrite!(log, "{}", aa as char);
    }
    if last && *last_color != "black".to_string() {
        emit_end_escape(&mut log);
    }
    log
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn aa_classes() -> Vec<(char, Vec<u8>)> {
    let mut classes = Vec::new();
    classes.push(('B', b"DN".to_vec()));
    classes.push(('Z', b"EQ".to_vec()));
    classes.push(('J', b"IL".to_vec()));
    classes.push(('-', b"DE".to_vec()));
    classes.push(('+', b"KHR".to_vec()));
    classes.push(('Ψ', b"ILMV".to_vec()));
    classes.push(('π', b"AGPS".to_vec()));
    classes.push(('Ω', b"FHWY".to_vec()));
    classes.push(('Φ', b"IFLMVWY".to_vec()));
    classes.push(('ζ', b"DEHKNQRST".to_vec()));
    classes.push(('X', b"ACDEFGHIKLMNPQRSTVWY".to_vec()));
    classes
}

pub fn cdr3_aa_con(
    style: &str,
    col: usize,
    exacts: &Vec<usize>,
    exact_clonotypes: &Vec<ExactClonotype>,
    rsi: &ColInfo,
) -> String {
    let mat = &rsi.mat;
    let mut cdr3s = Vec::<String>::new();
    for v in 0..exacts.len() {
        let m = mat[col][v];
        if m.is_some() {
            let ex = &exact_clonotypes[exacts[v]];
            cdr3s.push(ex.share[m.unwrap()].cdr3_aa.clone());
        }
    }
    let classes = aa_classes();
    let mut c = String::new();
    for i in 0..cdr3s[0].len() {
        let mut vals = Vec::<u8>::new();
        for j in 0..cdr3s.len() {
            vals.push(cdr3s[j].as_bytes()[i]);
        }
        unique_sort(&mut vals);
        if vals.solo() {
            c.push(vals[0] as char);
        } else if style == "x" {
            c.push('X');
        } else {
            for m in classes.iter() {
                if meet_size(&vals, &m.1) == vals.len() {
                    c.push(m.0);
                    break;
                }
            }
        }
    }
    c
}

pub fn get_gex_matrix_entry(
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    fid: usize,
    d_all: &Vec<Vec<u32>>,
    ind_all: &Vec<Vec<u32>>,
    li: usize,
    l: usize,
    p: usize,
    y: &str,
) -> f64 {
    let mut raw_count = 0 as f64;
    if gex_info.gex_matrices[li].initialized() {
        raw_count = gex_info.gex_matrices[li].value(p as usize, fid) as f64;
    } else {
        for j in 0..d_all[l].len() {
            if ind_all[l][j] == fid as u32 {
                raw_count = d_all[l][j] as f64;
                break;
            }
        }
    }
    let mult: f64;
    if y.ends_with("_g") {
        mult = gex_info.gex_mults[li];
    } else {
        mult = gex_info.fb_mults[li];
    }
    if !ctl.gen_opt.full_counts {
        raw_count *= mult;
    }
    raw_count
}

pub fn extra_args(ctl: &EncloneControl) -> Vec<String> {
    let mut extra_args = ctl.gen_opt.tree.clone();
    if !ctl.plot_opt.plot_xy_filename.is_empty() {
        extra_args.push(ctl.plot_opt.plot_xy_xvar.clone());
        extra_args.push(ctl.plot_opt.plot_xy_yvar.clone());
    }
    match ctl.plot_opt.cell_color {
        CellColor::ByVariableValue(ref x) => {
            extra_args.push(x.var.clone());
        }
        _ => {}
    };
    for i in 0..ctl.clono_filt_opt.bounds.len() {
        extra_args.append(&mut ctl.clono_filt_opt.bounds[i].var.clone());
    }
    if ctl.gen_opt.gene_scan_test.is_some() {
        extra_args.append(&mut ctl.gen_opt.gene_scan_test.as_ref().unwrap().var.clone());
        extra_args.append(&mut ctl.gen_opt.gene_scan_control.as_ref().unwrap().var.clone());
    }
    extra_args.append(&mut ctl.plot_opt.sim_mat_plot_vars.clone());
    for i in 0..ctl.gen_opt.var_def.len() {
        let x = &ctl.gen_opt.var_def[i].2;
        for v in vars_of_node(x).iter() {
            extra_args.push(decode_arith(v));
        }
    }
    unique_sort(&mut extra_args);
    extra_args
}
