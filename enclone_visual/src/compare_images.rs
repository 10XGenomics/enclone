// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
//
// Compare two RGB8 images, returning the number of pixel-color psoitions at which they differ.

use io_utils::*;
use rayon::prelude::*;
use std::io::Write;
use string_utils::*;

pub fn compare_images(
    image_data_old: &Vec<u8>,
    image_data_new: &Vec<u8>,
    width: usize,
    height: usize,
    verbose: bool,
) -> usize {
    let mut results = Vec::<(usize, usize, Vec<u8>)>::new();
    for x in 0..width {
        results.push((x, 0, Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let x = res.0;
        for y in 0..height {
            for c in 0..3 {
                let i = x * (height * 3) + 3 * 3 + c;
                let (vold, vnew) = (image_data_old[i] as isize, image_data_new[i] as isize);
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
