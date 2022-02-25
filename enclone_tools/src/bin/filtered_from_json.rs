// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Make filtered_contig.fasta and filtered_contig_annotations.csv from all_contig_annotations.json.

use io_utils::*;
use pretty_trace::PrettyTrace;
use serde_json::Value;
use std::io::{BufReader, Write};
use string_utils::*;

fn strip(x: &str) -> String {
    if x == "null" {
        String::new()
    } else {
        x.between("\"", "\"").to_string()
    }
}

fn main() {
    PrettyTrace::new().on();

    // Read in the json file and break it into entries.

    let mut f = BufReader::new(open_maybe_compressed("all_contig_annotations.json"));
    let mut xs = Vec::<Vec<u8>>::new();
    loop {
        match read_vector_entry_from_json(&mut f).unwrap() {
            None => break,
            Some(x) => {
                xs.push(x);
            }
        }
    }

    // Go through the json entries.

    let mut csv = open_for_write_new!["filtered_contig_annotations.csv"];
    let mut fasta = open_for_write_new!["filtered_contig.fasta"];
    fwriteln!(
        csv,
        "barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,\
        c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id"
    );
    for i in 0..xs.len() {
        let v: Result<Value, _> = serde_json::from_str(strme(&xs[i]));
        if v.is_err() {
            eprintln!(
                "\nInternal error, failed to parse a value from a string.  The string is:\n{}\n",
                strme(&xs[i])
            );
            std::process::exit(1);
        }
        let v = v.unwrap();
        let dataset = strip(&v["dataset"].to_string());
        let barcode = strip(&v["barcode"].to_string());
        let bc = barcode.before("-");
        let barcode = format!("{}-{}", bc, dataset);
        let is_cell = v["is_cell"].as_bool().unwrap_or(false);
        let mut contig_id = strip(&v["contig_name"].to_string());
        let contig = contig_id.rev_after("_");
        contig_id = format!("{}-{}_contig_{}", bc, dataset, contig);
        let high_confidence = v["high_confidence"].as_bool().unwrap_or(false);
        let seq = strip(&v["sequence"].to_string());
        let length = seq.len();
        let ann = v["annotations"].as_array().unwrap();
        let chain = strip(&ann[0]["feature"]["chain"].to_string());
        let mut v_gene = String::new();
        let mut d_gene = String::new();
        let mut j_gene = String::new();
        let mut c_gene = String::new();
        for j in 0..ann.len() {
            if ann[j]["feature"]["region_type"] == "L-REGION+V-REGION" {
                v_gene = strip(&ann[j]["feature"]["gene_name"].to_string());
            }
            if ann[j]["feature"]["region_type"] == "D-REGION" {
                d_gene = strip(&ann[j]["feature"]["gene_name"].to_string());
            }
            if ann[j]["feature"]["region_type"] == "J-REGION" {
                j_gene = strip(&ann[j]["feature"]["gene_name"].to_string());
            }
            if ann[j]["feature"]["region_type"] == "C-REGION" {
                c_gene = strip(&ann[j]["feature"]["gene_name"].to_string());
            }
        }
        let full_length = v["full_length"].as_bool().unwrap_or(false);
        let productive = v["productive"].as_bool().unwrap_or(false);
        let cdr3 = strip(&v["cdr3"].to_string());
        let cdr3_nt = strip(&v["cdr3_seq"].to_string());
        let reads = v["read_count"].as_i64().unwrap() as usize;
        let umis = v["umi_count"].as_i64().unwrap() as usize;
        let (raw_clonotype_id, raw_consensus_id) = (String::new(), String::new());
        fwriteln!(
            csv,
            "{barcode},{is_cell},{contig_id},{high_confidence},{length},{chain},\
            {v_gene},{d_gene},{j_gene},{c_gene},{full_length},{productive},{cdr3},{cdr3_nt},\
            {reads},{umis},{raw_clonotype_id},{raw_consensus_id}"
        );
        fwriteln!(fasta, ">{contig_id}\n{seq}");
    }
}
