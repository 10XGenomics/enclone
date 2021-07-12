// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Run some tests of enclone visual on a Mac.  As part of the tests, this opens a window.
//
// If you run with the single argument UPDATE, failing results will be replaced.
//
// Argument: QUIET.
//
// See also show_diffs.

use enclone_visual::compare_images::*;
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
    PrettyTrace::new().on();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // PRETEST
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Run enclone once to get it in cache.  This doesn't totally make sense but seems to improve
    // the reproducibility of timing of the actual work.

    let o = Command::new("enclone")
        .arg(&"--version")
        .output()
        .expect("failed to execute enclone visual pretest");
    if o.status.code() != Some(0) {
        eprintln!("\nnonzero exit code from enclone visual pretest\n");
        eprintln!("stderr =\n{}", strme(&o.stderr));
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // MAIN TESTS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    let t = Instant::now();
    let args: Vec<String> = env::args().collect();
    let mut update = false;
    let mut quiet = false;
    if args.len() >= 2 && args[1] == "UPDATE" {
        update = true;
    }
    if args.len() >= 2 && args[1] == "QUIET" {
        quiet = true;
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
    if !quiet {
        print!("{}", strme(&o.stdout));
    }
    let mut fail = false;
    const BIG_DIFFS: usize = 5;
    for i in 1..=TESTS.len() {
        if TESTS[i - 1].2.len() == 0 {
            continue;
        }
        let (mut image_old, mut image_new) = (Vec::<u8>::new(), Vec::<u8>::new());
        let old_file = format!("enclone_visual/regression_images/{}.png", TESTS[i - 1].2);
        if !path_exists(&old_file) {
            eprintln!(
                "\nLooks like you've added a test.  Please look at \
                enclone_visual/outputs/{}.png and\n\
                if it's right, copy to regression_tests and git add it.\n",
                TESTS[i - 1].2,
            );
            std::process::exit(1);
        }
        let mut f = File::open(&old_file).unwrap();
        f.read_to_end(&mut image_old).unwrap();
        let new_file = format!("enclone_visual/outputs/{}.png", TESTS[i - 1].2);
        let mut f = File::open(&new_file).unwrap();
        f.read_to_end(&mut image_new).unwrap();
        let (header, image_data_old) = png_decoder::decode(&image_old).unwrap();
        let (width, height) = (header.width as usize, header.height as usize);
        let (_, image_data_new) = png_decoder::decode(&image_new).unwrap();
        if image_data_old.len() != image_data_new.len() {
            eprintln!("\nimage size for test {} changed", i);
            std::process::exit(1);
        }
        let big_diffs = compare_images(&image_data_old, &image_data_new, width, height, true);
        if big_diffs > BIG_DIFFS {
            eprintln!(
                "\nThere are {} big diffs for {}.",
                big_diffs,
                TESTS[i - 1].2
            );
            fail = true;
            if update {
                copy(&new_file, &old_file).unwrap();
            }
        }
    }
    let state = if fail {
        "UNSUCCESSFULLY"
    } else {
        "successfully"
    };

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // TEST COMPUTATIONAL PERFORMANCE
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    let used = elapsed(&t);
    const EXPECTED_TIME: f64 = 21.3; // this is supposed to be the lowest observed value
    const MAX_PERCENT_OVER: f64 = 4.0;
    let percent_over = 100.0 * (used - EXPECTED_TIME) / EXPECTED_TIME;
    if percent_over > MAX_PERCENT_OVER {
        eprintln!(
            "\nUsed {:.1} seconds, exceeding expected test time of {:.1} seconds by {:.1}%, \
                versus max allowed = {}%.",
            used, EXPECTED_TIME, percent_over, MAX_PERCENT_OVER,
        );
        eprintln!(
            "You might want to retry a second time.  Note that if you're disconnected from\n\
            the internet, then the Mac Gatekeeper could introduce very long delays.\n"
        );
        std::process::exit(1);
    }

    // Get peak memory.  This is sensitive to changes of a few percent or so.

    let maxrss_children;
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        let retval = libc::getrusage(libc::RUSAGE_CHILDREN, &mut rusage as *mut _);
        assert_eq!(retval, 0);
        maxrss_children = rusage.ru_maxrss;
    }
    let peak_mem_mb = maxrss_children as f64 / ((1024 * 1024) as f64);
    const MAX_PEAK_MEM: f64 = 131.6; // this is supposed to be the lowest observed value
    const MAX_PERCENT_OVER_MEM: f64 = 16.2;
    let percent_over = 100.0 * (peak_mem_mb - MAX_PEAK_MEM) / MAX_PEAK_MEM;

    eprintln!(
        "\nPeak mem {:.1} MB, exceeded expected peak mem of {:.1} seconds by {:.1}%, \
            versus max allowed = {}%.",
        peak_mem_mb, MAX_PEAK_MEM, percent_over, MAX_PERCENT_OVER_MEM,
    );
    if percent_over > MAX_PERCENT_OVER_MEM {
        eprintln!("That's too high.  This happens occasionally, so please retry.\n");
        std::process::exit(1);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // DONE, REPORT STATUS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    println!(
        "\nenclone visual tests completely {} in {:.1} seconds using {:.1} MB\n",
        state, used, peak_mem_mb,
    );
    if fail {
        std::process::exit(1);
    }
}
