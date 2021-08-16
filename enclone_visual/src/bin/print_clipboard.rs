// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Print the contents of the clipboard, Mac only for now.

use enclone_visual::*;

fn main() {
    let copy = get_clipboard_content();
    println!("{}", copy);
}
