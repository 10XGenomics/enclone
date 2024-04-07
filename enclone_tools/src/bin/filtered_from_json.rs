use enclone_main::subset::AnnotationWithDataset;
// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Make filtered_contig.fasta and filtered_contig_annotations.csv from all_contig_annotations.json.
use io_utils::{fwriteln, open_for_write_new, open_maybe_compressed};
use std::io::Write;
use string_utils::TextUtils;
use vdj_types::VdjRegion;

fn main() {
    // Read in the json file and break it into entries.

    let mut contents = String::new();
    open_maybe_compressed("all_contig_annotations.json")
        .read_to_string(&mut contents)
        .unwrap();

    // Go through the json entries.

    let mut csv = open_for_write_new!["filtered_contig_annotations.csv"];
    let mut fasta = open_for_write_new!["filtered_contig.fasta"];
    fwriteln!(
        csv,
        "barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,\
        c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id"
    );
    for ann in serde_json::Deserializer::from_str(&contents).into_iter::<AnnotationWithDataset>() {
        let ann = ann.unwrap();
        let dataset = ann.dataset.unwrap_or("null".to_string());
        let ann = ann.data;

        let v_gene = ann.get_gene_name(VdjRegion::V).cloned().unwrap_or_default();
        let d_gene = ann.get_gene_name(VdjRegion::D).cloned().unwrap_or_default();
        let j_gene = ann.get_gene_name(VdjRegion::J).cloned().unwrap_or_default();
        let c_gene = ann.get_gene_name(VdjRegion::C).cloned().unwrap_or_default();

        let bc = ann.barcode.before("-");
        let barcode = format!("{bc}-{dataset}");
        let is_cell = ann.is_cell;
        let mut contig_id = ann.contig_name;
        let contig = contig_id.rev_after("_");
        contig_id = format!("{bc}-{dataset}_contig_{contig}");

        let length = ann.sequence.len();

        let chain = ann.annotations[0].feature.chain;
        let high_confidence = ann.high_confidence;
        let full_length = ann.full_length.unwrap_or_default();
        let productive = ann.productive.unwrap_or_default();
        let cdr3 = ann.cdr3.unwrap_or_default();
        let cdr3_nt = ann.cdr3_seq.unwrap_or_default();
        let reads = ann.read_count;
        let umis = ann.umi_count;
        let seq = ann.sequence;
        fwriteln!(
            csv,
            "{barcode},{is_cell},{contig_id},{high_confidence},{length},{chain},\
            {v_gene},{d_gene},{j_gene},{c_gene},{full_length},{productive},{cdr3},{cdr3_nt},\
            {reads},{umis},,"
        );
        fwriteln!(fasta, ">{contig_id}\n{seq}");
    }
}
