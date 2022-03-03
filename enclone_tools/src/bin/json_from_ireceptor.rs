// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Make all_contig_annotations.json.lz4 from an IReceptor tsv file.
//
// Usage: json_from_ireceptor whatever.tsv
//
// Limitations:
// 1. Fakes the quality scores as perfect because they're not there.
// 2. Doesn't annotate so you need to rerun enclone with BUILT_IN.
// 3. Sets high_confidence and is_cell to the given value for productive.
// 4. Contig names are assigned arbitrarily.
// 5. Read counts are set to zero.

use io_utils::*;
use lz4::EncoderBuilder;
use pretty_trace::PrettyTrace;
use std::collections::HashMap;
use std::env;
use std::fs::{remove_file, File};
use std::io::BufRead;
use std::io::Write;
use string_utils::*;

fn compress(source: &str, destination: &str) {
    let mut input_file = File::open(source).unwrap();
    let output_file = File::create(destination).unwrap();
    let mut encoder = EncoderBuilder::new().level(4).build(output_file).unwrap();
    std::io::copy(&mut input_file, &mut encoder).unwrap();
    let _ = encoder.finish();
}

fn lz4_file(f: &str) {
    compress(f, &format!("{}.lz4", f));
    remove_file(&f).unwrap();
}

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
    let args: Vec<String> = env::args().collect();
    let tsv = open_for_read![&args[1]];
    let mut to_field = HashMap::<String, usize>::new();
    let mut contigs = Vec::<Contig>::new();
    for (i, line) in tsv.lines().enumerate() {
        let s = line.unwrap();
        let fields = s.split('\t').collect::<Vec<&str>>();
        if i == 0 {
            for j in 0..fields.len() {
                to_field.insert(fields[j].to_string(), j);
            }
        } else {
            let productive = fields[to_field["productive"]] == "T";
            contigs.push(Contig {
                barcode: fields[to_field["cell_id"]].to_string(),
                contig_id: format!("contig_{}", i),
                is_cell: productive,
                high_confidence: productive,
                productive: productive,
                reads: 0,
                umis: fields[to_field["consensus_count"]].force_usize(),
                seq: fields[to_field["sequence"]].to_string().as_bytes().to_vec(),
            });
        }
    }
    {
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
    lz4_file("all_contig_annotations.json");
}
