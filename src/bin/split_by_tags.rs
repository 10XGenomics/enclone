// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Split a given all_contig_annotations.json file using multiplexing data and put the results
// in individual directories.
//
// Usage:
// split_by_tags vdj_lena_id full_path_of_tag_calls_file output_directory tag_scheme
// where tag_scheme is like "tag1,...,tagn;tag1',...tagn'".
//
// Test cases:
//
// 1. two PBMCs, BCR
//
//   split_by_tags 165809
//   /mnt/deck6/narek/multiplexing/umi_filtering/filtering_166824/outs/tag_calls_per_cell.csv
//   /mnt/assembly/vdj/splits/165809 "JWY119,JWY120,JWY121;JWY122,JWY123,JWY124"
//
//   enclone PRE=/mnt/assembly/vdj/splits/165809 BCR="165809_1;165809_2" FAIL_ONLY=true
//
//   enclone PRE=/mnt/assembly/vdj/splits/165809 BCR="165809_1;165809_2" CDR3=CAREGMYYDFWSLRSRVGMDVW
//
//   enclone PRE=/mnt/assembly/vdj BCR="splits/165809/165809_1,165807;splits/165809/165809_2,165808"
//   FAIL_ONLY=true
//
// 2. same two PBMCs, TCR
//
//   split_by_tags 165822
//   /mnt/deck6/narek/multiplexing/umi_filtering/filtering_166824/outs/tag_calls_per_cell.csv
//   /mnt/assembly/vdj/splits/165822 "JWY119,JWY120,JWY121;JWY122,JWY123,JWY124"
//
//   enclone PRE=/mnt/assembly/vdj/splits/165822 TCR="165822_1;165822_2" FAIL_ONLY=true
//
//   enclone PRE=/mnt/assembly/vdj TCR="splits/165822/165822_1,165820;splits/165822/165822_2,165821"
//   FAIL_ONLY=true

extern crate enclone;
extern crate io_utils;
extern crate marsoc;
extern crate pretty_trace;

use enclone::*;
use io_utils::*;
use marsoc::*;
use pretty_trace::*;
use std::env;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use subset_json::*;

pub fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let (lena, tags_file, out_dir, tag_scheme) = (&args[1], &args[2], &args[3], &args[4]);
    let groups = tag_scheme.split(';').collect::<Vec<&str>>();
    let mut tags = Vec::<Vec<&str>>::new();
    for i in 0..groups.len() {
        tags.push(groups[i].split(',').collect::<Vec<&str>>());
    }
    let json = format!("{}/all_contig_annotations.json", get_outs(&lena));
    for i in 0..groups.len() {
        let _ = fs::create_dir_all(&format!("{}/{}_{}/outs", out_dir, lena, i + 1));
        let f = open_for_read![&tags_file];
        let mut barcodes = Vec::<String>::new();
        for line in f.lines() {
            let s = line.unwrap();
            let fields = s.split(',').collect::<Vec<&str>>();
            if fields[0] == "cell_barcode" {
                continue;
            }
            let (barcode, assignment) = (fields[0], fields[2]);
            let mut keep = false;
            for j in 0..tags[i].len() {
                if tags[i][j] == assignment {
                    keep = true;
                    break;
                }
            }
            if keep {
                barcodes.push(barcode.to_string());
            }
        }
        barcodes.sort();
        println!("group {}, selected {} barcodes", i + 1, barcodes.len());
        let x = subset_all_contig_annotations_json(&json, &barcodes);
        let mut f = open_for_write_new![&format!(
            "{}/{}_{}/outs/all_contig_annotations.json",
            out_dir,
            lena,
            i + 1
        )];
        fwrite!(f, "{}", x);
    }
}
