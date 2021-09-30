// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Sort and replace the vars file.

use enclone_vars::sort_vars;
use pretty_trace::PrettyTrace;
use std::io::Write;

fn main() {
    PrettyTrace::new().on();
    let old = std::fs::read_to_string("enclone_vars/src/vars").unwrap();
    let new = sort_vars(&old);
    if new != old {
        let mut f = std::fs::File::create("enclone_vars/src/vars").unwrap();
        f.write_all(new.as_bytes()).unwrap();
    }
}
