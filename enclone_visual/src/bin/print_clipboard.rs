// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Print the contents of the clipboard.

use clipboard::ClipboardContext;
use clipboard::ClipboardProvider;

fn main() {
    let mut ctx: ClipboardContext = ClipboardProvider::new().unwrap();
    // println!("{:?}", ctx.get_contents());
    println!("{}", ctx.get_contents().unwrap());
}
