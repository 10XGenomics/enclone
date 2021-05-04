// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Process the SUBSET_JSON option.

use enclone_core::defs::*;
use io_utils::*;
use serde_json::Value;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use string_utils::*;
use vector_utils::*;

pub fn subset_json(
    ctl: &EncloneControl,
    exact_clonotypes: &Vec<ExactClonotype>,
    exacts: &Vec<Vec<usize>>,
    ann: &str,
) {
    if ctl.gen_opt.subset_json.len() > 0 {
        let mut barcodes = Vec::<String>::new();
        for l in 0..exacts.len() {
            for u in 0..exacts[l].len() {
                let ex = &exact_clonotypes[exacts[l][u]];
                for j in 0..ex.clones.len() {
                    barcodes.push(ex.clones[j][0].barcode.clone());
                }
            }
        }
        unique_sort(&mut barcodes);
        let mut g = open_for_write_new![&ctl.gen_opt.subset_json];
        fwriteln!(g, "[");
        let mut written = false;
        for li in 0..ctl.origin_info.dataset_path.len() {
            let json = format!("{}/{}", ctl.origin_info.dataset_path[li], ann);
            let mut jsonx = json.clone();
            if !path_exists(&json) {
                jsonx = format!("{}.lz4", json);
            }
            let mut xs = Vec::<Vec<u8>>::new();
            let mut f = BufReader::new(open_maybe_compressed(&jsonx));
            loop {
                match read_vector_entry_from_json(&mut f) {
                    None => break,
                    Some(x) => {
                        let v: Value = serde_json::from_str(strme(&x)).unwrap();
                        let barcode = &v["barcode"].to_string().between("\"", "\"").to_string();
                        if bin_member(&barcodes, &barcode) {
                            xs.push(x);
                        }
                    }
                }
            }
            for j in 0..xs.len() {
                if j == 0 && written {
                    fwriteln!(g, ",");
                }
                written = true;
                fwrite!(g, "{}", strme(&xs[j]));
                if j < xs.len() - 1 {
                    fwrite!(g, ",");
                    fwriteln!(g, "");
                }
            }
        }
        if written {
            fwriteln!(g, "");
        }
        fwriteln!(g, "]");
    }
}
