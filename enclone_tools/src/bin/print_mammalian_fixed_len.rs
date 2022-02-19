// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Print mammalian fixed len table as input to mammalian_fixed_len.rs.

use enclone_denovo::make_mammalian_fixed_len::*;
use io_utils::*;
use pretty_trace::*;
use std::io::Write;

fn main() {
    PrettyTrace::new().on();
    let x = make_mammalian_fixed_len();
    let mut log = open_for_write_new!["src/mammalian_fixed_len.table"];
    for i in 0..x.len() {
        let y = &x[i];
        fwrite!(log, "{},{},{},", y.0, y.1, y.2);
        for j in 0..y.3.len() {
            let z = &y.3[j];
            for k in 0..z.len() {
                let w = &z[k];
                fwrite!(log, "{}:{}", w.0, w.1 as char);
                if k < z.len() - 1 {
                    fwrite!(log, "/");
                }
            }
            if j < y.3.len() - 1 {
                fwrite!(log, "+");
            }
        }
        fwriteln!(log, "");
    }
}
