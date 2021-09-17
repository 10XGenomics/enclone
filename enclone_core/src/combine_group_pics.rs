// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::defs::*;
use ansi_escape::*;
use io_utils::*;
use std::io::Write;
use string_utils::*;
use tables::*;

pub fn combine_group_pics(
    group_pics: &Vec<String>,
    last_widths: &Vec<u32>,
    parseable_stdouth: bool,
    noprint: bool,
    noprintx: bool,
    html: bool,
    ngroup: bool,
    pretty: bool,
) -> String {
    let mut glog = Vec::<u8>::new();
    let mut done = false;
    if noprint && parseable_stdouth && group_pics.len() > 0 {
        let mut rows = Vec::<Vec<String>>::new();
        for i in 0..group_pics.len() {
            let r: Vec<String> = group_pics[i].split('\n').map(str::to_owned).collect();
            for j in 0..r.len() - 1 {
                let s = r[j].split('\t').map(str::to_owned).collect();
                rows.push(s);
            }
        }
        let mut same = true;
        let n = rows[0].len();
        for i in 1..rows.len() {
            if rows[i].len() != n {
                same = false;
            }
        }
        if same {
            let mut justify = Vec::<u8>::new();
            for x in rows[0].iter() {
                justify.push(justification(&x));
            }
            print_tabular(&mut glog, &rows, 2, Some(justify));
            done = true;
        }
    }

    if !done {
        // Get the newlines right is tricky, so they're marked.

        for i in 0..group_pics.len() {
            if !noprint {
                if !html && !ngroup && (!noprintx || i > 0) {
                    fwriteln!(glog, ""); // NEWLINE 1
                }

                // If we just printed a clonotype box, output a bar.

                if i > 0 && last_widths[i - 1] > 0 {
                    if ngroup || html {
                        fwriteln!(glog, ""); // NEWLINE 2
                    }
                    if pretty {
                        let mut log = Vec::<u8>::new();
                        emit_eight_bit_color_escape(&mut log, 44);
                        fwrite!(glog, "{}", strme(&log));
                    }
                    fwrite!(glog, "╺{}╸", "━".repeat((last_widths[i - 1] - 2) as usize));
                    if !ngroup {
                        fwriteln!(glog, ""); // NEWLINE 3
                    }
                    fwriteln!(glog, ""); // NEWLINE 4
                    if pretty {
                        let mut log = Vec::<u8>::new();
                        emit_end_escape(&mut log);
                        fwrite!(glog, "{}", strme(&log));
                    }
                }
            }
            glog.append(&mut group_pics[i].as_bytes().to_vec());
        }
    }
    stringme(&glog)
}
