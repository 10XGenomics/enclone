// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Build html files by inserting other html files.

use enclone::html::insert_html;
use pretty_trace::*;

fn main() {
    PrettyTrace::new().on();

    insert_html("pages/index.html.src", "index.html");
    insert_html("pages/expanded.html.src", "pages/auto/expanded.html");
}
