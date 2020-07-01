// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Build html files by inserting other html files.

use enclone::html::insert_html;
use pretty_trace::*;

fn main() {
    PrettyTrace::new().on();

    insert_html("pages/index.html.src", "index.html", false);
    insert_html("pages/expanded.html.src", "pages/auto/expanded.html", false);
    insert_html("pages/tree.html.src", "pages/auto/tree.html", false);
    insert_html("pages/plot.html.src", "pages/auto/plot.html", false);
    insert_html("pages/windows.html.src", "pages/auto/windows.html", false);
    insert_html("pages/history.html.src", "pages/auto/history.html", false);
    insert_html("pages/compile.html.src", "pages/auto/compile.html", false);
    insert_html(
        "pages/dang_i_cannot_install.html.src",
        "pages/auto/dang_i_cannot_install.html",
        false,
    );
    insert_html(
        "pages/installation_details.html.src",
        "pages/auto/installation_details.html",
        false,
    );
    insert_html(
        "pages/heuristics.html.src",
        "pages/auto/heuristics.html",
        false,
    );
}
