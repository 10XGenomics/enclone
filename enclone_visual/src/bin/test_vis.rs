// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Run some tests of enclone visual on a Mac.  As part of the tests, this opens a window.
//
// If you run with the single argument UPDATE, failing results will be replaced.
//
// Argument: QUIET.
// Argument: VERBOSE.
//
// This code works by comparing lowest resolution JPEG files.  We use that format to avoid
// having larger files in git.  A better solution would be to use lowest resolution
// JPEG2000 files, which would be even smaller.
//
// You need the following datasets to run this:
// dataset   notes
// 123085    public
// 123217    public, but we use our internal copy, which includes feature_barcode_matrix.bin
// 1145040   not public
// 1142282   not public
//
// See also show_diffs.

use enclone_visual::compare_images::*;
use enclone_visual::testsuite::TESTS;
use image::codecs::jpeg::JpegEncoder;
use image::ColorType::Rgba8;
use io_utils::*;
use jpeg_decoder::Decoder;
use perf_stats::*;
use pretty_trace::*;
use std::env;
use std::fs::{copy, File};
use std::io::{BufReader, BufWriter, Read};
use std::process::Command;
use std::time::Instant;
use string_utils::*;
use tilde_expand::*;

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
    let mut verbose = false;
    if args.len() >= 2 && args[1] == "UPDATE" {
        update = true;
    }
    if args.len() >= 2 && args[1] == "QUIET" {
        quiet = true;
    }
    if args.len() >= 2 && args[1] == "VERBOSE" {
        verbose = true;
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
    const MAX_DIFFS: usize = 150;
    for i in 1..=TESTS.len() {
        if TESTS[i - 1].2.len() == 0 {
            continue;
        }
        let mut image_new = Vec::<u8>::new();
        let old_png_file = format!("enclone_visual/regression_images/{}.png", TESTS[i - 1].2);
        let new_png_file = format!("enclone_visual/outputs/{}.png", TESTS[i - 1].2);
        let mut f = File::open(&new_png_file).unwrap();
        f.read_to_end(&mut image_new).unwrap();
        let (header, image_data_new0) = png_decoder::decode(&image_new).unwrap();
        let (width, height) = (header.width as usize, header.height as usize);

        // Convert the new png file to a jpg, and read that back in as a bit image.

        let new_jpg_file = format!("{}.jpg", new_png_file.rev_before(".png"));
        {
            let quality = 1 as u8; // lowest quality
            if path_exists(&new_jpg_file) {
                // This file removal is to circumvent a bug in Carbon Black.
                std::fs::remove_file(&new_jpg_file).unwrap();
            }
            let mut f = open_for_write_new![&new_jpg_file];
            let mut buff = BufWriter::new(&mut f);
            let mut encoder = JpegEncoder::new_with_quality(&mut buff, quality);
            encoder
                .encode(&image_data_new0, width as u32, height as u32, Rgba8)
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
        let diffs = compare_images(&image_data_old, &image_data_new, width, height, verbose);
        if diffs > MAX_DIFFS {
            eprintln!(
                "\nThere are {} diffs for {}.  Please open enclone_visual/outputs/\
                joint.{}.jpg.",
                diffs,
                TESTS[i - 1].2,
                TESTS[i - 1].2
            );

            // Create and save concatenated image.  Note that we're depending on the png
            // in regression_images as being current.  We could do the same thing with the jpg
            // versions, which would be guaranteed to be correct, but of lower quality.  We did
            // this in 8949a053552c478af5c952ee407416d0e52ab8a0 of dj/189, if you want to go back
            // to that.

            let mut f =
                File::open(&old_png_file).expect(&format!("\nCan't find {}.\n", old_png_file));
            let mut image_old = Vec::<u8>::new();
            f.read_to_end(&mut image_old).unwrap();
            let (_, image_data_old0) = png_decoder::decode(&image_old).unwrap();
            let mut joint = Vec::<u8>::new();
            for i in 0..height {
                let start = i * width * 4;
                let stop = (i + 1) * width * 4;
                joint.append(&mut image_data_old0[start..stop].to_vec());
                joint.append(&mut image_data_new0[start..stop].to_vec());
                for j in start..stop {
                    let diff = image_data_new0[j].wrapping_sub(image_data_old0[j]);
                    joint.push(255 - diff);
                }
            }
            let new_jpg_file = format!("enclone_visual/outputs/joint.{}.jpg", TESTS[i - 1].2);
            let quality = 80 as u8;
            let mut f = open_for_write_new![&new_jpg_file];
            let mut buff = BufWriter::new(&mut f);
            let mut encoder = JpegEncoder::new_with_quality(&mut buff, quality);
            encoder
                .encode(&joint, (width * 3) as u32, height as u32, Rgba8)
                .unwrap();

            // Keep going.

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
    const EXPECTED_TIME: f64 = 36.1; // this is supposed to be the lowest observed value
    const MAX_PERCENT_OVER: f64 = 4.2;
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
    const MAX_PEAK_MEM: f64 = 152.2; // this is supposed to be the lowest observed value
    const MAX_PERCENT_OVER_MEM: f64 = 17.6;
    let percent_over = 100.0 * (peak_mem_mb - MAX_PEAK_MEM) / MAX_PEAK_MEM;

    eprintln!(
        "\nPeak mem {:.1} MB, exceeded expected peak mem of {:.1} MB by {:.1}%, \
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
    // MAC TESTS THAT DON'T REALLY BELONG HERE BUT DON'T HAVE A BETTER HOME
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Test that a particular case of tilde expansion works on a Mac.

    let o = Command::new("enclone")
        .arg("BCR=123085")
        .arg("MIN_CELLS=10")
        .arg("PLOT_BY_ISOTYPE=gui")
        .arg("HONEY_OUT=~/enclone_temp_test_file")
        .arg("NOPRINT")
        .output()
        .expect("failed to execute enclone tilde test");
    if o.status.code() != Some(0) {
        eprintln!("\nnonzero exit code from enclone tilde test\n");
        eprintln!("stderr =\n{}", strme(&o.stderr));
        std::process::exit(1);
    }
    std::fs::remove_file(&strme(&tilde_expand("~/enclone_temp_test_file".as_bytes()))).unwrap();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // DONE, REPORT STATUS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    println!(
        "\nenclone visual tests completed {} in {:.1} seconds",
        state, used,
    );
    println!("actual total time = {:.1} seconds\n", elapsed(&tall));
    if fail {
        std::process::exit(1);
    }
}
