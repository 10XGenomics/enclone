// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Run some tests of enclone visual on a Mac.  As part of the tests, this opens a window.
//
// If you run with the single argument UPDATE, failing results will be replaced.
//
// Argument: QUIET.
//
// This code works by comparing lowest resolution JPEG files.  We use that format to avoid
// having larger files in git.  A better solution would be to use lowest resolution
// JPEG2000 files, which would be even smaller.
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

use image::codecs::jpeg::JpegEncoder;
use image::ColorType::Rgba8;
use jpeg_decoder::Decoder;
use std::io::BufReader;
use std::io::BufWriter;

fn main() {
    PrettyTrace::new().on();
    let tall = Instant::now();
    if !path_exists("enclone_visual") {
        eprintln!("\nYou need to run this from the top level directory of the enclone repo.\n");
        std::process::exit(1);
    }

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
    const MAX_DIFFS: usize = 0;
    for i in 1..=TESTS.len() {
        if TESTS[i - 1].2.len() == 0 {
            continue;
        }
        let mut image_new = Vec::<u8>::new();
        let old_png_file = format!("enclone_visual/regression_images/{}.png", TESTS[i - 1].2);
        let new_png_file = format!("enclone_visual/outputs/{}.png", TESTS[i - 1].2);
        let mut f = File::open(&new_png_file).unwrap();
        f.read_to_end(&mut image_new).unwrap();
        let (header, image_data_new) = png_decoder::decode(&image_new).unwrap();
        let (width, height) = (header.width as usize, header.height as usize);

        // Convert the new png file to a jpg, and read that back in as a bit image.

        let new_jpg_file = format!("{}.jpg", new_png_file.rev_before(".png"));
        {
            let quality = 1 as u8; // lowest quality
                                   // This file removal is to circumvent a bug in Carbon Black.
            std::fs::remove_file(&new_jpg_file).unwrap();
            let mut f = open_for_write_new![&new_jpg_file];
            let mut buff = BufWriter::new(&mut f);
            let mut encoder = JpegEncoder::new_with_quality(&mut buff, quality);
            encoder
                .encode(&image_data_new, width as u32, height as u32, Rgba8)
                .unwrap();
        }
        let file = open_for_read![&new_jpg_file];
        let mut decoder = Decoder::new(BufReader::new(file));
        let image_data_new = decoder.decode().expect("failed to decode image");

        // Check for existence of old jpg file.

        let old_jpg_file = format!("{}.jpg", old_png_file.rev_before(".png"));
        if !path_exists(&old_jpg_file) {
            eprintln!(
                "\nLooks like you've added a test.  Please look at \
                enclone_visual/outputs/{}.png and\n\
                if it's right, copy it and the jpg file to regression_tests and git add the jpg.\n",
                TESTS[i - 1].2,
            );
            std::process::exit(1);
        }

        // Read in the old jpg file as a bit image.

        let file = open_for_read![&old_jpg_file];
        let mut decoder = Decoder::new(BufReader::new(file));
        let image_data_old = decoder.decode().expect("failed to decode image");

        // Test for differences.

        if image_data_old.len() != image_data_new.len() {
            eprintln!("\nimage size for test {} changed", i);
            std::process::exit(1);
        }
        let diffs = compare_images(&image_data_old, &image_data_new, width, height, false);
        if diffs > MAX_DIFFS {
            eprintln!("\nThere are {} diffs for {}.", diffs, TESTS[i - 1].2);
            fail = true;
            if update {
                copy(&new_png_file, &old_png_file).unwrap();
                copy(&new_jpg_file, &old_jpg_file).unwrap();
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
    const EXPECTED_TIME: f64 = 23.2; // this is supposed to be the lowest observed value
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
    // TEST FOR STRAY PROCESSES
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    let o = Command::new("ps")
        .arg(&"ux")
        .output()
        .expect("failed to execute ps");
    let out = strme(&o.stdout);
    for line in out.lines() {
        if line.contains("ssh -p") {
            eprintln!(
                "\nLooks like you may have a stray process running:\n{}\n\n\
                Treating this as an error.\n",
                line
            );
            std::process::exit(1);
        }
    }
    for (key, value) in env::vars() {
        if key == "ENCLONE_CONFIG" {
            let host = value.before(":");
            let o = Command::new("ssh")
                .arg(&host)
                .arg(&"ps x")
                .output()
                .expect("failed to execute ssh");
            let out = strme(&o.stdout);
            for line in out.lines() {
                if line.contains("enclone") {
                    eprintln!(
                        "\nLooks like you may have a stray process running on the \
                        remote server:\n{}\n\n\
                        Treating this as an error.\n",
                        line
                    );
                    std::process::exit(1);
                }
            }
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // DONE, REPORT STATUS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    println!(
        "\nenclone visual tests completely {} in {:.1} seconds using {:.1} MB\n",
        state, used, peak_mem_mb,
    );
    println!("actual total time = {:.1} seconds\n", elapsed(&tall));
    if fail {
        std::process::exit(1);
    }
}
