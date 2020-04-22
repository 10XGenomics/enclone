// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Build html files by generating and inserting other html files.

use enclone::html::insert_html;
use enclone::testlist::SITE_EXAMPLES;
use pretty_trace::*;
use std::fs::File;
use std::process::{Command, Stdio};

fn main() {
    PrettyTrace::new().on();

    for i in 0..SITE_EXAMPLES.len() {
        let example_name = SITE_EXAMPLES[i].0;
        let test = SITE_EXAMPLES[i].1;
        let out_file = format!("pages/auto/{}.html", example_name);
        let args = test.split(' ').collect::<Vec<&str>>();
        let outputs = File::create(&out_file).unwrap();
        Command::new("target/release/enclone")
            .args(&args)
            .arg("HTML")
            .stdout(Stdio::from(outputs))
            .spawn()
            .unwrap()
            .wait_with_output()
            .unwrap();
    }

    insert_html("pages/index.html.src", "index.html");
    insert_html("pages/expanded.html.src", "pages/auto/expanded.html");
}
