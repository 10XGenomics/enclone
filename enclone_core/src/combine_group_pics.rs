// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::defs::EncloneControl;
use ansi_escape::*;
use io_utils::*;
use std::io::Write;
use string_utils::*;

pub fn combine_group_pics(
    group_pics: &Vec<String>,
    last_widths: &Vec<usize>,
    ctl: &EncloneControl,
) -> String {
    // Get the newlines right is tricky, so they're marked.

    let mut glog = Vec::<u8>::new();
    for i in 0..group_pics.len() {
        if !ctl.gen_opt.noprint {
            if !ctl.gen_opt.html && !ctl.clono_group_opt.ngroup && !ctl.gen_opt.noprintx {
                fwriteln!(glog, ""); // NEWLINE 1
            }

            // If we just printed a clonotype box, output a bar.

            if i > 0 && last_widths[i - 1] > 0 {
                if ctl.clono_group_opt.ngroup || ctl.gen_opt.html {
                    fwriteln!(glog, ""); // NEWLINE 2
                }
                if ctl.pretty {
                    let mut log = Vec::<u8>::new();
                    emit_eight_bit_color_escape(&mut log, 44);
                    fwrite!(glog, "{}", strme(&log));
                }
                fwrite!(glog, "╺{}╸", "━".repeat(last_widths[i - 1] - 2));
                if !ctl.clono_group_opt.ngroup {
                    fwriteln!(glog, ""); // NEWLINE 3
                }
                fwriteln!(glog, ""); // NEWLINE 4
                if ctl.pretty {
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
