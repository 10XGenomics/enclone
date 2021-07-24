// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use ansi_escape::*;
use io_utils::*;
use std::io::Write;
use string_utils::*;

pub fn combine_group_pics(
    group_pics: &Vec<String>,
    last_widths: &Vec<u32>,
    noprint: bool,
    noprintx: bool,
    html: bool,
    ngroup: bool,
    pretty: bool,
) -> String {
    // Get the newlines right is tricky, so they're marked.

    let mut glog = Vec::<u8>::new();
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
    stringme(&glog)
}
