// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Show differences for one test, after having run test_vis.

use enclone_visual::testsuite::TESTS;
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::Read;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    for i in 1..=TESTS.len() {
        if TESTS[i - 1].2 == args[1] {
            let (mut image_old, mut image_new) = (Vec::<u8>::new(), Vec::<u8>::new());
            let old_file = format!("enclone_visual/regression_images/{}.png", TESTS[i - 1].2);
            let mut f = File::open(&old_file).unwrap();
            f.read_to_end(&mut image_old).unwrap();
            let new_file = format!("enclone_visual/outputs/{}.png", TESTS[i - 1].2);
            let mut f = File::open(&new_file).unwrap();
            f.read_to_end(&mut image_new).unwrap();
            let (header, image_data_old) = png_decoder::decode(&image_old).unwrap();
            let (width, height) = (header.width as usize, header.height as usize);
            println!("\nimage is {} x {}\n", width, height);
            let (_, image_data_new) = png_decoder::decode(&image_new).unwrap();
            if image_data_old.len() != image_data_new.len() {
                eprintln!("\nimage size for test {} changed", i);
                std::process::exit(1);
            }
            const BIG: isize = 0;
            const MAX_DIFFS: usize = 477;
            let mut big_diffs = 0;
            for i in 0..image_data_old.len() {
                let (vold, vnew) = (image_data_old[i] as isize, image_data_new[i] as isize);
                let x = i / (height * 4);
                let y = (i % (height * 4)) / 4;
                let c = i % 4;
                if c == 3 && vold <= 14 && vnew <= 14 {
                    continue;
                }
                if c == 3 && (vold - vnew).abs() <= 5 {
                    continue;
                }
                if (vold - vnew).abs() > BIG {
                    println!("x = {}, y = {}, c = {}, vold = {}, vnew = {}",
                        x, y, c, vold, vnew
                    );
                    big_diffs += 1;
                }
            }
            if big_diffs > MAX_DIFFS {
                eprintln!(
                    "\nThere are {} big diffs for {}.\n",
                    big_diffs,
                    TESTS[i - 1].2
                );
            }
        }
    }
}
