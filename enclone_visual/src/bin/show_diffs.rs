// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Show differences for one test, after having run test_vis.

use enclone_visual::testsuite::TESTS;
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::Read;
use vector_utils::*;

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
            const BIG: isize = 2;
            const MAX_GRAY_DIFF: isize = 80;
            let mut big_diffs = 0;
            for x in 0..width {
                for y in 0..height {

                    // Test for small grayscale differences, ignoring transparency.

                    let mut olds = Vec::<isize>::new();
                    let mut news = Vec::<isize>::new();
                    for c in 0..3 {
                        let i = x * (height * 4) + y * 4 + c;
                        let (vold, vnew) = (image_data_old[i] as isize, image_data_new[i] as isize);
                        olds.push(vold);
                        news.push(vnew);
                    }
                    unique_sort(&mut olds);
                    unique_sort(&mut news);
                    if olds.solo() && news.solo() && (olds[0] - news[0]).abs() <= MAX_GRAY_DIFF {
                        continue;
                    }

                    // Test for other differences.

                    for c in 0..4 {
                        let i = x * (height * 4) + y * 4 + c;
                        let (vold, vnew) = (image_data_old[i] as isize, image_data_new[i] as isize);
                        if (vold - vnew).abs() > BIG {
                            println!("x = {}, y = {}, c = {}, vold = {}, vnew = {}",
                                x, y, c, vold, vnew
                            );
                            big_diffs += 1;
                        }
                    }
                }
            }
            eprintln!(
                "\nThere are {} big diffs for {}.\n",
                big_diffs,
                TESTS[i - 1].2
            );
        }
    }
}
