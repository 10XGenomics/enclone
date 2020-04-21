// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Build html files by generating and inserting other html files.  To be expanded.

use enclone::html::insert_html;
use pretty_trace::*;
use std::fs::File;
use std::process::{Command, Stdio};

fn main() {
    PrettyTrace::new().on();

    let out_file = "pages/auto/clonotype_with_gex.html";
    let outputs = File::create(&out_file).unwrap();
    Command::new("target/release/enclone")
        .arg("BCR=123085")
        .arg("CDR3=CQQRSNWPPSITF")
        .arg("GEX=123749")
        .arg("LVARSP=gex,IGHV3-49_g,CD19_ab")
        .arg("HTML")
        .stdout(Stdio::from(outputs))
        .spawn()
        .unwrap()
        .wait_with_output()
        .unwrap();

    insert_html("pages/index.html.src", "index.html");
}
