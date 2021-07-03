// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Run some tests of enclone visual on a Mac.  As part of the tests, this opens a window.
//
// If you run with the single argument UPDATE, failing results will be replaced.

use enclone_visual::testsuite::TESTS;
use io_utils::*;
use perf_stats::*;
use pretty_trace::*;
use std::env;
use std::fs::{copy, File};
use std::io::Read;
use std::process::Command;
use std::time::Instant;
use string_utils::*;

fn main() {
    let t = Instant::now();
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let mut update = false;
    if args.len() >= 2 && args[1] == "UPDATE" {
        update = true;
    }
    let o = Command::new("enclone")
        .arg(&"VIS")
        .arg(&"TEST")
        .output()
        .expect("failed to execute enclone visual test");
    if o.status.code() != Some(0) {
        eprintln!("\nnonzero exit code from enclone visual test\n");
        eprintln!("stderr =\n{}", strme(&o.stderr));
        std::process::exit(1);
    }
    print!("{}", strme(&o.stdout));
    let mut fail = false;
    for i in 1..=TESTS.len() {
        let (mut image_old, mut image_new) = (Vec::<u8>::new(), Vec::<u8>::new());
        let old_file = format!("enclone_visual/regression_images/test{}.png", i);
        if !path_exists(&old_file) {
            eprintln!(
                "\nLooks like you've added a test.  Please look at outputs/test{}.png and\n\
                if it's right, copy to regression_tests and git add it.\n",
                i,
            );
            std::process::exit(1);
        }
        let mut f = File::open(&old_file).unwrap();
        f.read_to_end(&mut image_old).unwrap();
        let new_file = format!("enclone_visual/outputs/test{}.png", i);
        let mut f = File::open(&new_file).unwrap();
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
            if update {
                copy(&new_file, &old_file).unwrap();
            }
        }
    }
    let state = if fail {
        "unsuccessfully"
    } else {
        "successfully"
    };
    let used = elapsed(&t);
    const MAX_TIME: f64 = 11.0;
    if used > MAX_TIME {
        eprintln!(
            "\nUsed {:.1} seconds, exceeding max test time of {} seconds.",
            used, MAX_TIME,
        );
        eprintln!(
            "You might want to retry a second time.  Note that if you're disconnected from\n\
            the internet, then the Mac Gatekeeper could introduce very long delays.\n"
        );
        std::process::exit(1);
    }
    println!(
        "\nenclone visual tests completely {} in {:.1} seconds\n",
        state, used,
    );
    if fail {
        std::process::exit(1);
    }
}
