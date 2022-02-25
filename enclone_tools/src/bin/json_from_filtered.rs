// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Make all_contig_annotations.json from filtered_contig.fasta and filtered_contig_annotations.csv.
//
// Limitations:
// 1. Fakes the quality scores as perfect because they're not there.
// 2. Doesn't annotate so you need to rerun enclone with BUILT_IN.

use io_utils::*;
use pretty_trace::PrettyTrace;
use std::collections::HashMap;
use std::io::BufRead;
use std::io::Write;
use string_utils::*;

pub struct Contig {
    pub barcode: String,
    pub contig_id: String,
    pub is_cell: bool,
    pub high_confidence: bool,
    pub productive: bool,
    pub reads: usize,
    pub umis: usize,
    pub seq: Vec<u8>,
}

fn main() {
    PrettyTrace::new().on();
    let fasta = open_for_read!["filtered_contig.fasta"];
    let csv = open_for_read!["filtered_contig_annotations.csv"];
    let mut to_field = HashMap::<String, usize>::new();
    let mut contigs = Vec::<Contig>::new();
    for (i, line) in csv.lines().enumerate() {
        let s = line.unwrap();
        let fields = parse_csv(&s);
        if i == 0 {
            for j in 0..fields.len() {
                to_field.insert(fields[j].to_string(), j);
            }
        } else {
            contigs.push(Contig {
                barcode: fields[to_field["barcode"]].to_string(),
                contig_id: fields[to_field["contig_id"]].to_string(),
                is_cell: fields[to_field["is_cell"]] == "True",
                high_confidence: fields[to_field["high_confidence"]] == "True",
                productive: fields[to_field["productive"]] == "True",
                reads: fields[to_field["reads"]].force_usize(),
                umis: fields[to_field["umis"]].force_usize(),
                seq: Vec::new(),
            });
        }
    }
    let mut to_tig = HashMap::<String, usize>::new();
    for i in 0..contigs.len() {
        to_tig.insert(contigs[i].contig_id.clone(), i);
    }
    let mut name = String::new();
    for (i, line) in fasta.lines().enumerate() {
        let s = line.unwrap();
        if i % 2 == 0 {
            name = s.after(">").to_string();
        } else {
            contigs[to_tig[&name]].seq = s.as_bytes().to_vec();
        }
    }
    let mut f = open_for_write_new!["all_contig_annotations.json"];
    fwriteln!(f, "[");
    for i in 0..contigs.len() {
        fwriteln!(f, "    {{");
        fwriteln!(f, "        \"barcode\": \"{}\",", contigs[i].barcode);
        fwriteln!(f, "        \"contig_name\": \"{}\",", contigs[i].contig_id);
        fwriteln!(f, "        \"is_cell\": {},", contigs[i].is_cell);
        fwriteln!(
            f,
            "        \"high_confidence\": {},",
            contigs[i].high_confidence
        );
        fwriteln!(f, "        \"productive\": {},", contigs[i].productive);
        fwriteln!(f, "        \"read_count\": {},", contigs[i].reads);
        fwriteln!(f, "        \"umi_count\": {},", contigs[i].umis);
        fwriteln!(f, "        \"sequence\": \"{}\",", strme(&contigs[i].seq));
        fwriteln!(
            f,
            "        \"quals\": \"{}\"",
            strme(&vec![b']'; contigs[i].seq.len()])
        );
        if i < contigs.len() - 1 {
            fwriteln!(f, "    }},");
        } else {
            fwriteln!(f, "    }}");
        }
    }
    fwriteln!(f, "]");
}
