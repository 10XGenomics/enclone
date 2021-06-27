// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Run some tests of enclone visual on a Mac.  As part of the tests, this opens a window.

use perf_stats::*;
use pretty_trace::*;
use std::fs::File;
use std::io::Read;
use std::process::Command;
use std::time::Instant;
use string_utils::*;

fn main() {
    let t = Instant::now();
    PrettyTrace::new().on();
    let o = Command::new("enclone")
        .arg(&"VIS")
        .arg(&"TEST")
        .output()
        .expect("failed to execute enclone visual test");
    if o.status.code() != Some(0) {
        eprintln!("\nnonzero exit code from enclone visual test\n");
        std::process::exit(1);
    }
    print!("{}", strme(&o.stdout));
    let mut fail = false;
    for i in 1..=3 {
        let (mut image_old, mut image_new) = (Vec::<u8>::new(), Vec::<u8>::new());
        let mut f = File::open(&format!("enclone_visual/regression_images/test{}.png", i)).unwrap();
        f.read_to_end(&mut image_old).unwrap();
        let mut f = File::open(&format!("enclone_visual/outputs/test{}.png", i)).unwrap();
        f.read_to_end(&mut image_new).unwrap();
        let (_, image_data_old) = png_decoder::decode(&image_old).unwrap();
        let (_, image_data_new) = png_decoder::decode(&image_new).unwrap();
        if image_data_old.len() != image_data_new.len() {
            eprintln!("\nimage size for test {} changed", i);
            std::process::exit(1);
        }
        const BIG: isize = 6;
        const MAX_DIFFS: usize = 400;
        let mut big_diffs = 0;
        for i in 0..image_data_old.len() {
            if ((image_data_old[i] as isize) - (image_data_new[i] as isize)).abs() > BIG {
                big_diffs += 1;
            }
        }
        if big_diffs > MAX_DIFFS {
            eprintln!("\nThere are {} big diffs for test {}.", big_diffs, i);
            fail = true;
        }
    }
    let state = if fail { "unsuccessfully" } else { "successfully" };
    println!("\nenclone visual tests completely {} in {:.1} seconds\n", state, elapsed(&t));
    if fail {
        std::process::exit(1);
    }
}
