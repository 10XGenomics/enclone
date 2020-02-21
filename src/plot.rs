// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::defs::*;
use io_utils::*;
use std::fs::File;
use std::io::Write;
use std::io::*;

// For radius r and n = 0, 1, ..., consider a counterclockwise spiral of lattice-packed disks of 
// radius r, starting at the origin and going first to the right.  Return the coordinates of the 
// center of the nth disk.  See this picture:
// https://www.researchgate.net/profile/Guorui_Li4/publication/220270050/figure/fig1/
//         AS:393993713143808@1470946829076/The-hexagonal-coordinate-system.png
// There is no attempt at efficiency.

pub fn hex_coord(n: usize, r: f64) -> (f64,f64) {
    // Special case.
    if n == 0 { 
        return (0.0, 0.0);
    }
    // If the hexagons are numbered 0, 1, ... outward, which hexagon "hid" are we on and
    // which position "hpos" on that are we at?
    let mut hid = 1;
    let mut k = 6;
    let mut hpos = n - 1;
    loop {
        if hpos < k { // WAS <= *******************************************************************
            break;
        }
        hpos -= k;
        hid += 1;
        k += 6;
        // k *= 2;
    }
    // Find coordinates.
    let c = r * 3.0f64.sqrt() / 2.0; // center to center distance, divided by 2
    let mut x = hid as f64 * 2.0 * c;
    let mut y = 0.0;
    let mut p = hpos;
    println!(""); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    printme!( n, hid, hpos, p ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if p > 0 {
        // Traverse the six faces, as far as we have to go.
        for _ in 0..hid {
            x -= c;
            y += 1.5; // was c
            p -= 1;
            if p == 0 {
                println!("take 1"); // XXX
                break;
            }
        }
        if p > 0 {
            for _ in 0..hid {
                x -= 2.0*c;
                p -= 1;
                if p == 0 {
                    println!("take 2"); // XXX
                    break;
                }
            }
            if p > 0 {
                for _ in 0..hid {
                    x -= c;
                    y -= 1.5;
                    p -= 1;
                    if p == 0 {
                        println!("take 3"); // XXX
                        break;
                    }
                }
                if p > 0 {
                    for _ in 0..hid {
                        x += c;
                        y -= 1.5;
                        p -= 1;
                        if p == 0 {
                            println!("take 4"); // XXX
                            break;
                        }
                    }
                    if p > 0 {
                        for _ in 0..hid {
                            x += 2.0*c;
                            p -= 1;
                            if p == 0 {
                                break;
                            }
                        }
                        if p > 0 {
                            for _ in 0..hid - 1 {
                                x += c;
                                y += 1.5;
                                p -= 1;
                                if p == 0 {
                                    println!("take 5"); // XXX
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printme!(x, y); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    x *= 2.0 / 3.0f64.sqrt();
    y *= 2.0 / 3.0f64.sqrt();
    (x, y)
}

// Pack circles of given radii.  There is probably a literature on this, and this is probably
// a very crappy algorithm.  The answer is certainly not optimal.  The run time is O(n^2) where
// the constant includes a factor of 100.  Return centers for the circles.

pub fn pack_circles( r: &Vec<f64> ) -> Vec<(f64,f64)> {
    let mut c = Vec::<(f64,f64)>::new();
    if r.is_empty() {
        return c;
    }
    c.push( (0.0, 0.0) );
    let mut bigr = r[0];
    let mut rand = 0i64;
    // XXX -- supposed to be 100
    const SAMPLE : usize = 100000;
    const MUL : f64 = 1.5;
    for i in 1..r.len() {
        println!(""); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        printme!( i, bigr, r[i] ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        let mut q = Vec::<(f64,f64,f64)>::new();
        loop {
            for _ in 0..SAMPLE {
                // Get a random point in [-1,+1] x [-1,+1].
                let rand1 = 6_364_136_223_846_793_005i64
                    .wrapping_mul(rand)
                    .wrapping_add(1_442_695_040_888_963_407);
                let rand2 = 6_364_136_223_846_793_005i64
                    .wrapping_mul(rand1)
                    .wrapping_add(1_442_695_040_888_963_407);
                rand = rand2;
                let mut r1 = ( 2.0 * ( rand1 % 1_000_000i64 ) as f64 / 1_000_000.0 ) - 1.0;
                let mut r2 = ( 2.0 * ( rand2 % 1_000_000i64 ) as f64 / 1_000_000.0 ) - 1.0;
                // Make it bigger.
                r1 *= (bigr + r[i]) * MUL;
                r2 *= (bigr + r[i]) * MUL;
                // See if circle at (r1,r2) overlaps any of the existing circles.
                let mut ok = true;
                for k in 0..i {
                    let d = ( (c[k].0-r1)*(c[k].0-r1) + (c[k].1-r2)*(c[k].1-r2) ).sqrt();
                    if d < r[i] + r[k] {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    q.push( (r1*r1 + r2*r2, r1, r2) );
                }
            }
            q.sort_by(|a, b| a.partial_cmp(b).unwrap());
            if !q.is_empty() {
                break;
            }
        }
        printme!( q[0].1, q[0].2 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        c.push( (q[0].1, q[0].2) );
        bigr = bigr.max( r[i] + (c[i].0*c[i].0 + c[i].1*c[i].1).sqrt() );
    }
    c
}

// Given a collection of circles having specified colors, create an svg string that shows the
// circles on a canvas of fixed size.  The circles are moved and resized accordingly.

pub fn circles_to_svg( 
    center: &Vec<(f64,f64)>, 
    radius: &Vec<f64>, 
    color: &Vec<String>,
    width: usize,
    height: usize,
) -> String {
    let n = center.len();
    assert!( !center.is_empty() );
    assert!( radius.len() == n );
    assert!( color.len() == n );
    for i in 0..n {
        assert!( radius[i] > 0.0 );
    }
    let mut out = 
        format!("<svg version=\"1.1\"\n\
            baseProfile=\"full\"\n\
            width=\"{}\" height=\"{}\"\n\
            xmlns=\"http://www.w3.org/2000/svg\">\n", 
            width, height);
    let mut center = center.clone();
    let mut radius = radius.clone();
    let mut xmin = center[0].0;
    let mut xmax = center[0].0;
    let mut ymin = center[0].1;
    let mut ymax = center[0].1;
    for i in 0..n {
        xmin = xmin.min( center[i].0 - radius[i] );
        xmax = xmax.max( center[i].0 + radius[i] );
        ymin = ymin.min( center[i].1 - radius[i] );
        ymax = ymax.max( center[i].1 + radius[i] );
    }
    let scale = ((width as f64)/(xmax-xmin)).min( (height as f64)/(ymax-ymin) );
    for i in 0..n {
        center[i].0 -= xmin;
        center[i].1 -= ymin;
        center[i].0 *= scale;
        center[i].1 *= scale;
        radius[i] *= scale;
    }
    for i in 0..center.len() {
        out += &format!(
          "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"{}\" />\n", 
          center[i].0, center[i].1, radius[i], color[i] );
    }

    out += "</svg>\n";
    out
}

pub fn plot_clonotypes(
    ctl: &EncloneControl,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
) {
    if ctl.gen_opt.plot_file.is_empty() {
        return;
    }
    let mut clusters = Vec::<(Vec<String>,Vec<(f64,f64)>)>::new();
    let mut radii = Vec::<f64>::new();
    const SEP : f64 = 1.0; // separation between clusters
    for i in 0..exacts.len() {
        let mut colors = Vec::<String>::new();
        let mut coords = Vec::<(f64,f64)>::new();
        let mut n = 0;
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            for j in 0..ex.clones.len() {
                let mut color = "black".to_string();
                if ex.clones[j][0].sample_index.is_some() {
                    let s = &ctl.sample_info.sample_list[ex.clones[j][0].sample_index.unwrap()];
                    if ctl.gen_opt.sample_color_map.contains_key(&s.clone()) {
                        color = ctl.gen_opt.sample_color_map[s].clone();
                    }
                }
                colors.push(color);
                coords.push( hex_coord( n, 1.0 ) );
                n += 1;
            }
        }
        let mut radius = 0.0f64;
        for j in 0..coords.len() {
            radius 
                = radius.max( 1.0 + (coords[j].0*coords[j].0 + coords[j].1*coords[j].1).sqrt() );
        }
        radius += SEP;
        clusters.push( (colors, coords) );
        radii.push(radius);
    }
    let centers = pack_circles(&radii);

    // XXX:
    for i in 0..clusters.len() {
        println!( "\nCLUSTER {} ==> radius {} ==> {}, {}", i, radii[i], centers[i].0, centers[i].1 );
        for j in 0..clusters[i].1.len() {
            println!( "{} ==> {}, {}", j, clusters[i].1[j].0, clusters[i].1[j].1 );
        }
    }

    for i in 0..clusters.len() {
        for j in 0..clusters[i].1.len() {
            clusters[i].1[j].0 += centers[i].0;
            clusters[i].1[j].1 += centers[i].1;
        }
    }
    let mut center = Vec::<(f64,f64)>::new();
    let mut radius = Vec::<f64>::new();
    let mut color = Vec::<String>::new();
    for i in 0..clusters.len() {
        for j in 0..clusters[i].0.len() {
            color.push( clusters[i].0[j].clone() );
            center.push( ( clusters[i].1[j].0, clusters[i].1[j].1 ) );
            radius.push(1.0);
        }
    }
    const WIDTH : usize = 400;
    const HEIGHT : usize = 400;
    for i in 0..center.len() {
        center[i].1 = -center[i].1; // otherwise inverted, not sure why
    }
    let svg = circles_to_svg(&center, &radius, &color, WIDTH, HEIGHT );
    let mut f = open_for_write_new![ctl.gen_opt.plot_file];
    fwriteln!(f, "{}", svg);
}
