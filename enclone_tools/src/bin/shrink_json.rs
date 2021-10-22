// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Usage: shrink_json <input all_contig_annotations.json> <output all_contig_annotations.json>
//
// Read and remove certain fields that are not needed, at least if BUILT_IN and REPROD are used.

use io_utils::*;
use pretty_trace::PrettyTrace;
use serde_json::Value;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let json = std::fs::read_to_string(&args[1]).unwrap();
    let v: Value = serde_json::from_str(&json).unwrap();
    let entries = &mut v.as_array().unwrap().clone();
    for i in 0..entries.len() {
        let mut x = entries[i].as_object().unwrap().clone();
        x.remove("aa_sequence").unwrap();
        x.remove("annotations").unwrap();
        x.remove("clonotype").unwrap();
        x.remove("frame").unwrap();
        let _ = x.remove("full_length");
        x.remove("info").unwrap();
        x.remove("start_codon_pos").unwrap();
        x.remove("stop_codon_pos").unwrap();
        entries[i] = serde_json::Value::Object(x);
    }
    let mut f = open_for_write_new![&args[2]];
    fwriteln!(f, "[");
    for i in 0..entries.len() {
        fwrite!(f, "{}", entries[i]);
        if i < entries.len() - 1 {
            fwrite!(f, ",");
        }
        fwriteln!(f, "");
    }
    fwriteln!(f, "]");
}
