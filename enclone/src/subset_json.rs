// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Extract the entries in a given all_contig_annotations.json file that corrrespond to barcodes
// in a given sorted vector.

use io_utils::*;
use std::io::BufRead;
use string_utils::*;
use vector_utils::*;

pub fn subset_all_contig_annotations_json(filename: &str, barcodes: &Vec<String>) -> String {
    let mut x = "[\n".to_string();
    let mut lines = Vec::<String>::new();
    let f = open_userfile_for_read(filename);
    let mut keep = false;
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with('[') {
            continue;
        }
        if s == "]" {
            if keep {
                for i in 0..lines.len() {
                    x += &format!("{}\n", lines[i]);
                }
            }
            break;
        }
        lines.push(s.clone());
        if s.starts_with("        \"barcode\": \"") {
            let t = s.between("        \"barcode\": \"", "\"");
            keep = bin_member(barcodes, &t.to_string());
        } else if s.starts_with("    }") {
            if keep {
                for i in 0..lines.len() {
                    x += &format!("{}\n", lines[i]);
                }
            }
            lines.clear();
        }
    }
    x += "]\n";
    x
}
