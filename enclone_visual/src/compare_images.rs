// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
//
// Compare two images.

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
    const BIG: isize = 4;
    const MAX_GRAY_DIFF: isize = 80;
    let mut results = Vec::<(usize, usize, Vec<u8>)>::new();
    for x in 0..width {
        results.push((x, 0, Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let x = res.0;
        let mut olds = Vec::<isize>::new();
        let mut news = Vec::<isize>::new();
        for y in 0..height {
            // Test for small grayscale differences, ignoring transparency.

            olds.clear();
            news.clear();
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
    let mut big_diffs = 0;
    for x in 0..results.len() {
        print!("{}", strme(&results[x].2));
        big_diffs += results[x].1;
    }
    big_diffs
}
