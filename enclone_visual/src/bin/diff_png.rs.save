// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use pretty_trace::*;

fn main() {
    PrettyTrace::new().on();
    let image_bytes = include_bytes!("../../../test1.png");
    let (header, image_data1) = png_decoder::decode(image_bytes).unwrap();

    println!("Header: {:#?}", header);
    println!("Image data size: {}x{}x4 {}", header.width, header.height, image_data1.len());

    let image_bytes = include_bytes!("../../../test1.png.save");
    let (header, image_data2) = png_decoder::decode(image_bytes).unwrap();

    println!("Header: {:#?}", header);
    println!("Image data size: {}x{}x4 {}", header.width, header.height, image_data2.len());

    assert!(image_data1.len() == image_data2.len());
    let mut diffs = 0;
    let mut z1 = 0;
    let mut z2 = 0;
    let mut delta = 0;
    let mut big_diffs = 0;
    let mut lt50 = 0;
    let mut gt100 = 0;
    for i in 0..image_data1.len() {
        if image_data1[i] != image_data2[i] {
            diffs += 1;
        }
        if image_data1[i] == 255 {
            z1 += 1;
        }
        if image_data2[i] == 255 {
            z2 += 1;
        }
        if image_data1[i] == 255 && image_data2[i] != 255 {
            delta += 1;
        }
        if ((image_data1[i] as isize) - (image_data2[i] as isize)).abs() > 5 {
            big_diffs += 1;
            // println!("{} ==> {} == {}", i, image_data1[i], image_data2[i]);
        }
        if image_data1[i] < 50 {
            lt50 += 1;
        }
        if image_data1[i] > 100 {
            gt100 += 1;
        }
    }
    println!("diffs = {} of {}", diffs, image_data1.len());
    println!("z1 = {}, z2 = {}", z1, z2);
    println!("delta = {}", delta);
    println!("big_diffs = {}", big_diffs);
    println!("lt50 = {}", lt50);
    println!("gt100 = {}", gt100);

}
