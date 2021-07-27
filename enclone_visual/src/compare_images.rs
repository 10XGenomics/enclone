// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
//
// Compare two RGB8 images, returning the number of pixel-color positions at which they differ,
// excluding certain expected differences.

use io_utils::*;
use rayon::prelude::*;
use std::io::Write;
use string_utils::*;
use vector_utils::*;

pub fn compare_images(
    image_data_old: &Vec<u8>,
    image_data_new: &Vec<u8>,
    width: usize,
    height: usize,
    verbose: bool,
) -> usize {
    const MAX_DIFF: u8 = 61;
    const MAX_GRAY_DIFF: u8 = 134;
    let mut results = Vec::<(usize, usize, Vec<u8>)>::new();
    for x in 0..width {
        results.push((x, 0, Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let x = res.0;
        let (mut olds, mut news) = (Vec::<u8>::new(), Vec::<u8>::new());
        for y in 0..height {
            //
            // Test for relatively small differences.

            let mut close = true;
            for c in 0..3 {
                let i = x * (height * 3) + y * 3 + c;
                let (old, new) = (image_data_old[i], image_data_new[i]);
                let diff = if old <= new { new - old } else { old - new };
                if diff > MAX_DIFF {
                    close = false;
                }
            }
            if close {
                continue;
            }

            // Test for relatively small grayscale differences.

            olds.clear();
            news.clear();
            for c in 0..3 {
                let i = x * (height * 3) + y * 3 + c;
                olds.push(image_data_old[i]);
                news.push(image_data_new[i]);
            }
            unique_sort(&mut olds);
            unique_sort(&mut news);
            if olds.solo() && news.solo() {
                let diff = if olds[0] <= news[0] {
                    news[0] - olds[0]
                } else {
                    olds[0] - news[0]
                };
                if diff <= MAX_GRAY_DIFF {
                    continue;
                }
            }

            // Test for other differences.

            for c in 0..3 {
                let i = x * (height * 3) + y * 3 + c;
                let (vold, vnew) = (image_data_old[i], image_data_new[i]);
                if vnew != vold {
                    if verbose {
                        fwriteln!(
                            res.2,
                            "x = {}, y = {}, c = {}, vold = {}, vnew = {}",
                            x,
                            y,
                            c,
                            vold,
                            vnew
                        );
                    }
                    res.1 += 1;
                }
            }
        }
    });
    let mut diffs = 0;
    for x in 0..results.len() {
        if verbose {
            print!("{}", strme(&results[x].2));
        }
        diffs += results[x].1;
    }
    diffs
}
