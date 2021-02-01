// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// The purpose of this file is the function plot_clonotypes.  It plots clonotypes as partial
// hexagonal closest packings.  This is visually kind of satisfying, but also a bit weird looking.
// In some cases, by eye, you can see rounder forms that could be created by relocating some of
// the cells.

use crate::polygon::*;
use crate::string_width::*;
use ansi_escape::*;
use enclone_core::defs::*;
use io_utils::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::io::*;
use std::time::Instant;
use string_utils::*;
use vdj_ann::refx::*;
use vector_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Change the width or height of an svg document.

fn set_svg_width(svg: &mut String, new_width: f64) {
    let (svg1, svg2) = (svg.before("width="), svg.after("width=\"").after("\""));
    *svg = format!("{}width=\"{}\"{}", svg1, new_width, svg2);
}

fn set_svg_height(svg: &mut String, new_height: f64) {
    let (svg1, svg2) = (svg.before("height="), svg.after("height=\"").after("\""));
    *svg = format!("{}height=\"{}\"{}", svg1, new_height, svg2);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// For radius r and n = 0, 1, ..., consider a counterclockwise spiral of lattice-packed disks of
// radius r, starting at the origin and going first to the right.  Return the coordinates of the
// center of the nth disk.  See this picture:
// https://www.researchgate.net/profile/Guorui_Li4/publication/220270050/figure/fig1/
//         AS:393993713143808@1470946829076/The-hexagonal-coordinate-system.png
// There is no attempt at efficiency.

fn hex_coord(n: usize, r: f64) -> (f64, f64) {
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
        if hpos < k {
            break;
        }
        hpos -= k;
        hid += 1;
        k += 6;
    }
    // Find coordinates.
    let c = r * 3.0f64.sqrt() / 2.0; // center to center distance, divided by 2
    let mut x = hid as f64 * 2.0 * c;
    let mut y = 0.0;
    let mut p = hpos;
    if p > 0 {
        // Traverse the six faces, as far as we have to go.
        for _ in 0..hid {
            x -= c;
            y += 1.5;
            p -= 1;
            if p == 0 {
                break;
            }
        }
        if p > 0 {
            for _ in 0..hid {
                x -= 2.0 * c;
                p -= 1;
                if p == 0 {
                    break;
                }
            }
            if p > 0 {
                for _ in 0..hid {
                    x -= c;
                    y -= 1.5;
                    p -= 1;
                    if p == 0 {
                        break;
                    }
                }
                if p > 0 {
                    for _ in 0..hid {
                        x += c;
                        y -= 1.5;
                        p -= 1;
                        if p == 0 {
                            break;
                        }
                    }
                    if p > 0 {
                        for _ in 0..hid {
                            x += 2.0 * c;
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
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    x *= 2.0 / 3.0f64.sqrt();
    y *= 2.0 / 3.0f64.sqrt();
    (x, y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Pack circles of given radii, which should be in descending order.  Return centers for the
// circles.  There is probably a literature on this, and this is probably a very crappy algorithm.
//
// Blacklisted polygons are avoided.

fn pack_circles(r: &Vec<f64>, blacklist: &Vec<Polygon>) -> Vec<(f64, f64)> {
    // Set up to track the centers of the circles.

    let mut c = Vec::<(f64, f64)>::new();
    if r.is_empty() {
        return c;
    }

    // Set up to track the x extents of the circles, which are intervals, and which we break into
    // subintervals equal to the extent of the smallest circle (and a possible fractional part).
    // We track the starting point of each subinterval, and the index of the circle.
    //
    // The overall code would presumably be faster if we did not store the data in a sorted
    // vector that needs to be sorted after each center is added.  However it's not actually
    // clear that much time is spent in the sort.

    let mut ints = Vec::<(f64, usize)>::new();
    let radius0 = *r.last().unwrap();

    // Define function to add a center.

    fn push_center(
        x: f64,
        y: f64,
        index: usize,
        radius: f64,
        radius0: f64,
        c: &mut Vec<(f64, f64)>,
        ints: &mut Vec<(f64, usize)>,
    ) {
        c.push((x, y));

        // need to cover [x-radius..x+radius] using segments of length radius0

        let mut pos = x - radius;
        while pos < x + radius {
            ints.push((pos, index));
            pos += radius0;
        }
        ints.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }

    // Define data for the first circle.

    if blacklist.is_empty() {
        push_center(0.0, 0.0, 0, r[0], radius0, &mut c, &mut ints);
    }

    // Compute the maximum distance "bigr" from the origin.

    let mut bigr = if blacklist.is_empty() { r[0] } else { 0.0 };
    for p in blacklist.iter() {
        for x in p.v.iter() {
            bigr = bigr.max(x.origin_dist());
        }
    }

    // Proceed.

    let mut rand = 0i64;
    // We use a ridiculously large sample.  Reducing it to 10,000 noticeably reduces symmetry.
    // Presumably as the number of clusters increases, the sample would need to be increased
    // (ideally) to increase symmetry.
    const SAMPLE: usize = 100000;
    const MUL: f64 = 1.5;
    let mut centers = vec![(0.0, 0.0); SAMPLE];
    let start = if blacklist.is_empty() { 1 } else { 0 };
    for i in start..r.len() {
        let mut found = false;
        let mut best_r1 = 0.0;
        let mut best_r2 = 0.0;
        let mut best_val = 0.0;

        // Loop until we find a placement for the circle.

        loop {
            for z in 0..SAMPLE {
                // Get a random point in [-1,+1] x [-1,+1].  Using a hand-rolled random number
                // generator (from the internet) for speed and reproducibility, although there
                // might be something better in the standard packages.

                let rand1 = 6_364_136_223_846_793_005i64
                    .wrapping_mul(rand)
                    .wrapping_add(1_442_695_040_888_963_407);
                let rand2 = 6_364_136_223_846_793_005i64
                    .wrapping_mul(rand1)
                    .wrapping_add(1_442_695_040_888_963_407);
                rand = rand2;
                let mut r1 = (2.0 * (rand1 % 1_000_000i64) as f64 / 1_000_000.0) - 1.0;
                let mut r2 = (2.0 * (rand2 % 1_000_000i64) as f64 / 1_000_000.0) - 1.0;

                // Make it bigger.

                r1 *= (bigr + r[i]) * MUL;
                r2 *= (bigr + r[i]) * MUL;
                centers[z] = (r1, r2);
            }
            let mut results = Vec::<(usize, f64, f64, bool)>::new();
            for z in 0..SAMPLE {
                results.push((z, 0.0, 0.0, false));
            }
            results.par_iter_mut().for_each(|res| {
                let z = res.0;
                let r1 = centers[z].0;
                let r2 = centers[z].1;

                // Test for overlap with a blacklisted polygon.

                let mut ok = true;
                let m = Point { x: r1, y: r2 };
                for p in blacklist.iter() {
                    if p.touches_disk(m, r[i]) {
                        ok = false;
                        break;
                    }
                }

                // See if circle at (r1,r2) with radius r[i] overlaps any of the existing circles.
                // This involves a binary search, which reduces the complexity of the algorithm.

                if ok {
                    let low = ints.binary_search_by(|v| {
                        v.partial_cmp(&(r1 - r[i] - radius0, 0))
                            .expect("Comparison failed.")
                    });
                    let low = match low {
                        Ok(i) => i,
                        Err(i) => i,
                    };
                    let high = ints.binary_search_by(|v| {
                        v.partial_cmp(&(r1 + r[i], 0)).expect("Comparison failed.")
                    });
                    let high = match high {
                        Ok(i) => i,
                        Err(i) => i,
                    };
                    for m in low..high {
                        if m > low && ints[m].1 == ints[m - 1].1 {
                            continue;
                        }
                        let k = ints[m].1;
                        let d = (c[k].0 - r1) * (c[k].0 - r1) + (c[k].1 - r2) * (c[k].1 - r2);
                        if d < (r[i] + r[k]) * (r[i] + r[k]) {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    res.1 = r1;
                    res.2 = r2;
                    res.3 = true;
                }
            });
            for z in 0..SAMPLE {
                if results[z].3 {
                    let (r1, r2) = (results[z].1, results[z].2);
                    let val;
                    if blacklist.is_empty() || i == 0 {
                        val = r1 * r1 + r2 * r2;
                    } else {
                        val = (r1 - c[0].0) * (r1 - c[0].0) + (r2 - c[0].1) * (r2 - c[0].1);
                    }
                    if !found || val < best_val {
                        best_r1 = r1;
                        best_r2 = r2;
                        best_val = val;
                    }
                    found = true;
                }
            }
            if found {
                break;
            }
        }
        push_center(best_r1, best_r2, i, r[i], radius0, &mut c, &mut ints);
        bigr = bigr.max(r[i] + (c[i].0 * c[i].0 + c[i].1 * c[i].1).sqrt());
    }
    c
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given a collection of circles having specified colors, create an svg string that shows the
// circles on a canvas of fixed size.  The circles are moved and resized accordingly.
// Also shades smoothed polygons.

fn circles_to_svg(
    center: &Vec<(f64, f64)>,
    radius: &Vec<f64>,
    color: &Vec<String>,
    shades: &Vec<Polygon>,
    shade_colors: &Vec<String>,
    shade_enclosures: &Vec<Polygon>,
    width: usize,
    height: usize,
    boundary: usize,
) -> String {
    let n = center.len();
    assert!(!center.is_empty());
    assert!(radius.len() == n);
    assert!(color.len() == n);
    assert!(boundary < width);
    assert!(boundary < height);
    for i in 0..n {
        assert!(radius[i] > 0.0);
    }
    let mut out = format!(
        "<svg version=\"1.1\"\n\
         baseProfile=\"full\"\n\
         width=\"{}\" height=\"{}\"\n\
         xmlns=\"http://www.w3.org/2000/svg\">\n",
        width, height
    );
    let mut center = center.clone();
    let mut radius = radius.clone();
    let mut xmin = center[0].0;
    let mut xmax = center[0].0;
    let mut ymin = center[0].1;
    let mut ymax = center[0].1;
    for i in 0..n {
        xmin = xmin.min(center[i].0 - radius[i]);
        xmax = xmax.max(center[i].0 + radius[i]);
        ymin = ymin.min(center[i].1 - radius[i]);
        ymax = ymax.max(center[i].1 + radius[i]);
    }
    for i in 0..shades.len() {
        for j in 0..shades[i].v.len() {
            xmin = xmin.min(shade_enclosures[i].v[j].x);
            xmax = xmax.max(shade_enclosures[i].v[j].x);
            ymin = ymin.min(shade_enclosures[i].v[j].y);
            ymax = ymax.max(shade_enclosures[i].v[j].y);
        }
    }
    let width = width - boundary;
    let height = height - boundary;
    let scale = ((width as f64) / (xmax - xmin)).min((height as f64) / (ymax - ymin));
    for i in 0..n {
        center[i].0 -= xmin;
        center[i].1 -= ymin;
        center[i].0 *= scale;
        center[i].1 *= scale;
        radius[i] *= scale;
        center[i].0 += boundary as f64;
        center[i].1 += boundary as f64;
    }
    let mut shades = shades.clone();
    for i in 0..shades.len() {
        for j in 0..shades[i].v.len() {
            shades[i].v[j].x -= xmin;
            shades[i].v[j].y -= ymin;
            shades[i].v[j].x *= scale;
            shades[i].v[j].y *= scale;
            shades[i].v[j].x += boundary as f64;
            shades[i].v[j].y += boundary as f64;
        }
    }
    for (g, p) in shades.iter().enumerate() {
        out += "<path d=\"";
        out += &format!("{}", p.catmull_bezier_svg());
        out += "\" ";
        out += &format!("fill=\"{}\"\n", shade_colors[g]);
        out += " stroke=\"rgb(150,150,150)\"";
        out += " stroke-width=\"0.2\"";
        out += "/>\n";
    }
    for i in 0..center.len() {
        out += &format!(
            "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"{}\" />\n",
            center[i].0, center[i].1, radius[i], color[i]
        );
    }
    out += "</svg>\n";
    out
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Here, and in "enclone help color", we swap the order of colors, placing the last three before
// the first three.  This is because the last three seem to make a better three-color palette.

fn substitute_enclone_color(color: &mut String) {
    if *color == "@1".to_string() {
        *color = "rgb(0,95,175)".to_string();
    } else if *color == "@2".to_string() {
        *color = "rgb(215,135,175)".to_string();
    } else if *color == "@3".to_string() {
        *color = "rgb(0,175,135)".to_string();
    } else if *color == "@4".to_string() {
        *color = "rgb(215,95,0)".to_string();
    } else if *color == "@5".to_string() {
        *color = "rgb(95,175,255)".to_string();
    } else if *color == "@6".to_string() {
        *color = "rgb(215,175,0)".to_string();
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn plot_clonotypes(
    ctl: &EncloneControl,
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    svg: &mut String,
) {
    let t = Instant::now();
    if ctl.gen_opt.plot_file.is_empty() {
        return;
    }
    if exacts.is_empty() {
        eprintln!("\nThere are no clonotypes to plot, giving up.\n");
        std::process::exit(1);
    }

    let mut const_names = Vec::<String>::new();
    for id in refdata.cs.iter() {
        if refdata.rtype[*id] == 0 {
            const_names.push(refdata.name[*id].clone());
        }
    }
    unique_sort(&mut const_names);
    if ctl.gen_opt.plot_by_isotype_color.len() > 0 {
        if const_names.len() + 1 > ctl.gen_opt.plot_by_isotype_color.len() {
            eprintln!(
                "\nUsing the PLOT_BY_ISOTYPE_COLOR argument, you specified {} colors, \
                but there are {} constant region\nnames, and one more color is needed for the \
                \"undetermined\" case.  Please add more colors.\n",
                ctl.gen_opt.plot_by_isotype_color.len(),
                const_names.len()
            );
            std::process::exit(1);
        }
    } else if ctl.gen_opt.plot_by_isotype && const_names.len() > 12 {
        eprintln!(
            "\nCurrently PLOT_BY_ISOTYPE only works if there are at most 12 constant \
            region names.  If this is a problem, please let us know and we will generalize it.\n"
        );
        std::process::exit(1);
    }
    let mut clusters = Vec::<(Vec<String>, Vec<(f64, f64)>)>::new();
    let mut radii = Vec::<f64>::new();
    const SEP: f64 = 1.0; // separation between clusters
    let mut origins = Vec::<String>::new();

    // Go through the clonotypes.

    for i in 0..exacts.len() {
        let mut colors = Vec::<String>::new();
        let mut coords = Vec::<(f64, f64)>::new();
        let mut n = 0;

        // For PLOT_BY_MARK, find the dataset having the largest number of cells.

        let mut dsx = 0;
        if ctl.gen_opt.plot_by_mark {
            let mut ds = Vec::<usize>::new();
            for j in 0..exacts[i].len() {
                let ex = &exact_clonotypes[exacts[i][j]];
                for j in 0..ex.clones.len() {
                    ds.push(ex.clones[j][0].dataset_index);
                }
            }
            ds.sort();
            let mut freq = Vec::<(u32, usize)>::new();
            make_freq(&ds, &mut freq);
            dsx = freq[0].1;
        }

        // Go through the exact subclonotypes in a clonotype.

        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            for j in 0..ex.clones.len() {
                let mut color = "black".to_string();

                // Determine color for PLOT_BY_ISOTYPE.

                if ctl.gen_opt.plot_by_isotype {
                    let mut crefs = Vec::<Option<usize>>::new();
                    for l in 0..ex.share.len() {
                        if ex.share[l].left {
                            crefs.push(ex.share[l].c_ref_id);
                        }
                    }
                    unique_sort(&mut crefs);
                    let mut color_id = 0;
                    if crefs.solo() && crefs[0].is_some() {
                        let c = &refdata.name[crefs[0].unwrap()];
                        let p = bin_position(&const_names, &c) as usize;
                        color_id = 1 + p;
                    }
                    if ctl.gen_opt.plot_by_isotype_color.is_empty() {
                        let x = print_color13(color_id);
                        color = format!("rgb({},{},{})", x.0, x.1, x.2);
                    } else {
                        color = ctl.gen_opt.plot_by_isotype_color[color_id].clone();
                    }

                // Determine color for PLOT_BY_MARK.
                } else if ctl.gen_opt.plot_by_mark {
                    let dom = ex.clones[j][0].dataset_index == dsx;
                    let marked = ex.clones[j][0].marked;
                    if dom {
                        if !marked {
                            color = "red".to_string();
                        } else {
                            color = "rgb(255,200,200)".to_string();
                        }
                    } else {
                        if !marked {
                            color = "blue".to_string();
                        } else {
                            color = "rgb(200,200,255)".to_string();
                        }
                    }

                // Determine color in other cases.
                } else {
                    if ex.clones[j][0].origin_index.is_some() {
                        let s = &ctl.origin_info.origin_list[ex.clones[j][0].origin_index.unwrap()];
                        origins.push(s.clone());
                        if ctl.gen_opt.origin_color_map.contains_key(&s.clone()) {
                            color = ctl.gen_opt.origin_color_map[s].clone();
                        }
                    }
                    if ctl.gen_opt.origin_color_map.is_empty() {
                        let mut dataset_colors = false;
                        for c in ctl.origin_info.color.iter() {
                            if !c.is_empty() {
                                dataset_colors = true;
                            }
                        }
                        let di = ex.clones[j][0].dataset_index;
                        if dataset_colors {
                            color = ctl.origin_info.color[di].clone();
                        } else {
                            let bc = &ex.clones[j][0].barcode;
                            if ctl.origin_info.barcode_color[di].contains_key(bc) {
                                color = ctl.origin_info.barcode_color[di][bc].clone();
                            }
                        }
                    }
                }
                colors.push(color);
                coords.push(hex_coord(n, 1.0));
                n += 1;
            }
        }
        unique_sort(&mut origins);

        // Move the colors around to get vertical separation, e.g. blues on the left, reds
        // on the right.

        colors.sort();
        coords.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Substitute enclone colors.

        for j in 0..colors.len() {
            substitute_enclone_color(&mut colors[j]);
        }

        // Save.

        let mut radius = 0.0f64;
        for j in 0..coords.len() {
            radius =
                radius.max(1.0 + (coords[j].0 * coords[j].0 + coords[j].1 * coords[j].1).sqrt());
        }
        radius += SEP;
        clusters.push((colors, coords));
        radii.push(radius);
    }

    // Set group specification.

    let mut group_id = vec![0; radii.len()];
    let mut group_color = vec!["".to_string()];
    let mut group_name = vec!["".to_string()];
    if ctl.gen_opt.clonotype_group_names.is_some() {
        let f = open_for_read![&ctl.gen_opt.clonotype_group_names.as_ref().unwrap()];
        let mut first = true;
        let mut group_id_field = 0;
        let mut new_group_name_field = 0;
        let mut new_group_names = vec![None; radii.len()];
        for line in f.lines() {
            let s = line.unwrap();
            for c in s.chars() {
                if c.is_control() || c == '\u{FEFF}' {
                    eprintln!(
                        "\nThe first line in your CLONOTYPE_GROUP_NAMES file contains a \
                        nonprinting character.\n"
                    );
                    std::process::exit(1);
                }
            }
            let fields = s.split(',').collect::<Vec<&str>>();
            if first {
                let p = position(&fields, &"group_id");
                if p < 0 {
                    eprintln!("\nThe CLONOTYPE_GROUP_NAMES file does not have a group_id field.\n");
                    std::process::exit(1);
                }
                group_id_field = p as usize;
                let p = position(&fields, &"new_group_name");
                if p < 0 {
                    eprintln!(
                        "\nThe CLONOTYPE_GROUP_NAMES file does not have a \
                         new_group_name field.\n"
                    );
                    std::process::exit(1);
                }
                new_group_name_field = p as usize;
                first = false;
            } else {
                let group_id = &fields[group_id_field];
                if !group_id.parse::<usize>().is_ok() || group_id.force_usize() == 0 {
                    eprintln!(
                        "\nThe group_id {} in your CLONOTYPE_GROUP_NAMES file is not a \
                        positive integer.\n",
                        group_id
                    );
                    std::process::exit(1);
                }
                let group_id = group_id.force_usize() - 1;
                if group_id > radii.len() {
                    eprintln!(
                        "\nThe group_id {} in your CLONOTYPE_GROUP_NAMES file is larger \
                        than the number of clonotypes, which is {}.\n",
                        group_id,
                        radii.len()
                    );
                    std::process::exit(1);
                }
                let new_group_name = fields[new_group_name_field].to_string();
                if new_group_names[group_id].is_none() {
                    new_group_names[group_id] = Some(new_group_name.clone());
                } else if *new_group_names[group_id].as_ref().unwrap() != new_group_name {
                    eprintln!(
                        "\nThe group_id {} in your CLONOTYPE_GROUP_NAMES file is assigned \
                        the different new_group_names {} and {}.\n",
                        group_id,
                        new_group_names[group_id].as_ref().unwrap(),
                        new_group_name,
                    );
                    std::process::exit(1);
                }
            }
        }

        // Reverse sort by total number of cells associated to a name.  This defines group names.

        {
            let mut nx = Vec::<(String, usize)>::new();
            for i in 0..new_group_names.len() {
                if new_group_names[i].is_some() {
                    nx.push((
                        new_group_names[i].as_ref().unwrap().to_string(),
                        clusters[i].1.len(),
                    ));
                }
            }
            nx.sort();
            let mut ny = Vec::<(usize, String)>::new();
            let mut i = 0;
            while i < nx.len() {
                let j = next_diff1_2(&nx, i as i32) as usize;
                let mut n = 0;
                for k in i..j {
                    n += nx[k].1;
                }
                ny.push((n, nx[i].0.clone()));
                i = j;
            }
            reverse_sort(&mut ny);
            group_name.clear();
            for i in 0..ny.len() {
                group_name.push(ny[i].1.clone());
            }
        }

        // Define group ids.

        group_id.clear();
        for i in 0..new_group_names.len() {
            if new_group_names[i].is_some() {
                let p = position(&group_name, &new_group_names[i].as_ref().unwrap());
                group_id.push(p as usize);
            }
        }

        // Build colors.

        let sum = 40; // a+b+c, where color is rgb(255-a, 255-b, 255-c); smaller is closer to white
        let mut rand = 0i64;
        let mut points = Vec::<(f64, f64, f64)>::new();
        while points.len() < 10_000 {
            let rand1 = 6_364_136_223_846_793_005i64
                .wrapping_mul(rand)
                .wrapping_add(1_442_695_040_888_963_407);
            let rand2 = 6_364_136_223_846_793_005i64
                .wrapping_mul(rand1)
                .wrapping_add(1_442_695_040_888_963_407);
            let rand3 = 6_364_136_223_846_793_005i64
                .wrapping_mul(rand2)
                .wrapping_add(1_442_695_040_888_963_407);
            rand = rand3;
            let mut r1 = (rand1 % 1_000_000i64) as f64 / 1_000_000.0;
            let mut r2 = (rand2 % 1_000_000i64) as f64 / 1_000_000.0;
            let mut r3 = (rand3 % 1_000_000i64) as f64 / 1_000_000.0;
            r1 = (r1 + 1.0) / 2.0;
            r2 = (r2 + 1.0) / 2.0;
            r3 = (r3 + 1.0) / 2.0;
            let r = r1 + r2 + r3;
            if r > 0.0 {
                r1 /= r;
                r2 /= r;
                r3 /= r;
                points.push((r1, r2, r3));
            }
        }

        // Prepopulate colors with what appear to be an optimal sequence.  This is as measured by
        // distance from each other.  These colors are not optimal relative to color blindness.

        let mut fracs = vec![
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
            (0.0, 1.0 / 3.0, 2.0 / 3.0),
            (1.0 / 3.0, 2.0 / 3.0, 0.0),
            (2.0 / 3.0, 1.0 / 3.0, 0.0),
            (2.0 / 3.0, 0.0, 1.0 / 3.0),
            (0.0, 2.0 / 3.0, 1.0 / 3.0),
            (1.0 / 3.0, 0.0, 2.0 / 3.0),
            (5.0 / 9.0, 2.0 / 9.0, 2.0 / 9.0),
            (4.0 / 9.0, 1.0 / 9.0, 4.0 / 9.0),
            (7.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0),
            (1.0 / 9.0, 4.0 / 9.0, 4.0 / 9.0),
            (2.0 / 9.0, 5.0 / 9.0, 2.0 / 9.0),
            (4.0 / 9.0, 4.0 / 9.0, 1.0 / 9.0),
            (2.0 / 9.0, 2.0 / 9.0, 5.0 / 9.0),
            (1.0 / 9.0, 1.0 / 9.0, 7.0 / 9.0),
            (1.0 / 9.0, 7.0 / 9.0, 1.0 / 9.0),
            (1.0 / 9.0, 2.0 / 9.0, 6.0 / 9.0),
            (3.0 / 9.0, 1.0 / 9.0, 5.0 / 9.0),
            (3.0 / 9.0, 2.0 / 9.0, 4.0 / 9.0),
            (1.0 / 9.0, 8.0 / 9.0, 0.0 / 9.0),
            (3.0 / 9.0, 4.0 / 9.0, 2.0 / 9.0),
            (5.0 / 9.0, 3.0 / 9.0, 1.0 / 9.0),
        ];
        if fracs.len() >= group_name.len() {
            fracs.truncate(group_name.len());

        // Add more colors if needed.
        } else {
            for _ in fracs.len()..group_name.len() {
                let mut max_dist = 0.0;
                let mut best = (0.0, 0.0, 0.0);
                for p in points.iter() {
                    let mut min_dist = 100.0_f64;
                    for q in fracs.iter() {
                        let d = (p.0 - q.0) * (p.0 - q.0)
                            + (p.1 - q.1) * (p.1 - q.1)
                            + (p.2 - q.2) * (p.2 - q.2);
                        min_dist = min_dist.min(d);
                    }
                    if min_dist > max_dist {
                        max_dist = min_dist;
                        best = *p;
                    }
                }
                fracs.push(best);
            }
        }
        group_color.clear();
        for i in 0..fracs.len() {
            let r = 255 - (fracs[i].0 * sum as f64).round() as usize;
            let g = 255 - (fracs[i].1 * sum as f64).round() as usize;
            let b = 255 - (fracs[i].2 * sum as f64).round() as usize;
            group_color.push(format!("rgb({},{},{})", r, g, b));
        }
    }
    let ngroups = group_color.len();
    ctl.perf_stats(&t, "in preamble to plotting clonotypes");

    // Traverse the groups.

    let t = Instant::now();
    let using_shading = ngroups > 1 || group_color[0].len() > 0;
    let mut blacklist = Vec::<Polygon>::new();
    let mut shades = Vec::<Polygon>::new();
    let mut shade_colors = Vec::<String>::new();
    let mut shade_enclosures = Vec::<Polygon>::new();
    let mut centers = vec![(0.0, 0.0); radii.len()];
    for g in 0..ngroups {
        // Gather the group.

        let mut ids = Vec::<usize>::new();
        let mut radiix = Vec::<f64>::new();
        for i in 0..group_id.len() {
            if group_id[i] == g {
                ids.push(i);
                radiix.push(radii[i]);
            }
        }

        // Find circle centers.

        let centersx = pack_circles(&radiix, &blacklist);
        for i in 0..ids.len() {
            centers[ids[i]] = centersx[i];
        }

        // Find polygon around the group.

        if using_shading {
            let mut z = Vec::<(f64, f64, f64)>::new();
            for i in 0..centersx.len() {
                z.push((radiix[i], centersx[i].0, centersx[i].1));
            }
            let d = 5.0; // distance of polygon from the circles
            let n = 35; // number of vertices on polygon
            let mut p = enclosing_polygon(&z, d, n);
            shades.push(p.clone());
            shade_colors.push(group_color[g].clone());

            // Build an enlarged polygon that would hopefully include the smoothed polygonal
            // curve.  We should actually calculate to enforce this.

            p.enlarge(25.0);
            shade_enclosures.push(p.clone());
            blacklist.push(p);
        }

        // Reorganize constant-color clusters so that like-colored clusters are proximate,
        // We got this idea from Ganesh Phad, who showed us a picture!  The primary effect is on
        // single-cell clonotypes.

        let mut ccc = Vec::<(usize, String, usize)>::new(); // (cluster size, color, index)
        for i in 0..ids.len() {
            let id = ids[i];
            let mut c = clusters[id].0.clone();
            unique_sort(&mut c);
            if c.solo() {
                ccc.push((clusters[i].0.len(), c[0].clone(), i));
            }
        }
        ccc.sort();
        let mut i = 0;
        while i < ccc.len() {
            let j = next_diff1_3(&ccc, i as i32) as usize;
            let mut angle = vec![(0.0, 0); j - i];
            for k in i..j {
                let id = ccc[k].2;
                angle[k - i] = (centersx[id].1.atan2(centersx[id].0), id);
            }
            angle.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for k in i..j {
                let new_id = angle[k - i].1;
                for u in 0..clusters[ids[new_id]].0.len() {
                    clusters[ids[new_id]].0[u] = ccc[k].1.clone();
                }
            }
            i = j;
        }
    }
    ctl.perf_stats(&t, "plotting clonotypes");

    // Build the svg file.

    let t = Instant::now();
    for i in 0..clusters.len() {
        for j in 0..clusters[i].1.len() {
            clusters[i].1[j].0 += centers[i].0;
            clusters[i].1[j].1 += centers[i].1;
        }
    }
    let mut center = Vec::<(f64, f64)>::new();
    let mut radius = Vec::<f64>::new();
    let mut color = Vec::<String>::new();
    for i in 0..clusters.len() {
        for j in 0..clusters[i].0.len() {
            color.push(clusters[i].0[j].clone());
            center.push((clusters[i].1[j].0, clusters[i].1[j].1));
            radius.push(1.0);
        }
    }
    const WIDTH: usize = 400;
    const HEIGHT: usize = 400;
    const BOUNDARY: usize = 10;

    // Negate y coordinates, as otherwise images are inverted, not sure why.

    for i in 0..center.len() {
        center[i].1 = -center[i].1;
    }
    for i in 0..shades.len() {
        for j in 0..shades[i].v.len() {
            shades[i].v[j].y = -shades[i].v[j].y;
        }
    }
    for i in 0..shade_enclosures.len() {
        for j in 0..shade_enclosures[i].v.len() {
            shade_enclosures[i].v[j].y = -shade_enclosures[i].v[j].y;
        }
    }

    // Generate svg.

    *svg = circles_to_svg(
        &center,
        &radius,
        &color,
        &shades,
        &shade_colors,
        &shade_enclosures,
        WIDTH,
        HEIGHT,
        BOUNDARY,
    );

    // Calculate the actual height and width of the svg.

    let mut actual_height = 0.0f64;
    let fields = svg.split(' ').collect::<Vec<&str>>();
    let mut y = 0.0;
    for i in 0..fields.len() {
        if fields[i].starts_with("cy=") {
            y = fields[i].between("\"", "\"").force_f64();
        }
        if fields[i].starts_with("r=") {
            let r = fields[i].between("\"", "\"").force_f64();
            actual_height = actual_height.max(y + r);
        }
    }
    let mut actual_width = 0.0f64;
    let fields = svg.split(' ').collect::<Vec<&str>>();
    let mut x = 0.0;
    for i in 0..fields.len() {
        if fields[i].starts_with("cx=") {
            x = fields[i].between("\"", "\"").force_f64();
        }
        if fields[i].starts_with("r=") {
            let r = fields[i].between("\"", "\"").force_f64();
            actual_width = actual_width.max(x + r);
        }
    }

    // Add legend for shading.

    let mut font_size = 20;
    const LEGEND_BOX_STROKE_WIDTH: usize = 2;
    let mut legend_xstop_shading = 0.0;
    if using_shading {
        font_size = 16;
        let n = ngroups;
        let mut max_string_width = 0.0f64;
        for s in group_name.iter() {
            max_string_width = max_string_width.max(arial_width(s, font_size));
        }
        let color_bar_width = 50.0;
        let vsep = 3.0;
        let legend_height = ((font_size + BOUNDARY / 2) * n + BOUNDARY) as f64 + n as f64 * vsep;
        let legend_width = BOUNDARY as f64 * 2.5 + color_bar_width + max_string_width + 5.0;
        let legend_xstart = actual_width + BOUNDARY as f64 + 10.0;
        let legend_ystart = 50.0;
        *svg = svg.rev_before("<").to_string();
        *svg += &format!(
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" \
             style=\"fill:white;stroke:black;stroke-width:{}\" />\n",
            legend_xstart, legend_ystart, legend_width, legend_height, LEGEND_BOX_STROKE_WIDTH
        );
        legend_xstop_shading = legend_xstart + legend_width;
        for i in 0..n {
            // Determine y start.
            let y = legend_ystart as f64
                + BOUNDARY as f64 * 2.5
                + ((font_size + BOUNDARY / 2) * i) as f64
                + i as f64 * vsep;
            // Add group name.
            *svg += &format!(
                "<text x=\"{}\" y=\"{}\" font-family=\"Arial\" \
                 font-size=\"{}\">{}</text>\n",
                legend_xstart + color_bar_width + BOUNDARY as f64 * 2.0,
                y - BOUNDARY as f64 * 0.5,
                font_size,
                group_name[i]
            );
            // Add color bar.
            *svg += &format!(
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\" />\n",
                legend_xstart + (BOUNDARY * 1) as f64,
                y - BOUNDARY as f64 * 2.0,
                color_bar_width,
                font_size + BOUNDARY / 2,
                group_color[i]
            );
        }
        let new_width = legend_xstart + legend_width + 5.0;
        set_svg_width(svg, new_width);
        *svg += "</svg>";
    }

    // Add main legend.

    if ctl.gen_opt.use_legend || ctl.gen_opt.plot_by_isotype || ctl.gen_opt.plot_by_mark {
        let (mut colors, mut labels) = (Vec::<String>::new(), Vec::<String>::new());
        let mut max_string_width = 0.0f64;
        if ctl.gen_opt.plot_by_isotype {
            for i in 0..const_names.len() {
                labels.push(const_names[i].clone());
                let color_id = i + 1;
                if ctl.gen_opt.plot_by_isotype_color.is_empty() {
                    let x = print_color13(color_id);
                    let color = format!("rgb({},{},{})", x.0, x.1, x.2);
                    colors.push(color);
                } else {
                    let color = ctl.gen_opt.plot_by_isotype_color[color_id].clone();
                    colors.push(color);
                }
            }
            labels.push("undetermined".to_string());
            let color_id = 0;
            if ctl.gen_opt.plot_by_isotype_color.is_empty() {
                let x = print_color13(color_id);
                let color = format!("rgb({},{},{})", x.0, x.1, x.2);
                colors.push(color);
            } else {
                let color = ctl.gen_opt.plot_by_isotype_color[color_id].clone();
                colors.push(color);
            }
        } else if ctl.gen_opt.plot_by_mark {
            colors.push("red".to_string());
            labels.push("in most common dataset, !marked".to_string());
            colors.push("rgb(255,200,200)".to_string());
            labels.push("in most common dataset, marked".to_string());
            colors.push("blue".to_string());
            labels.push("not in most common dataset, !marked".to_string());
            colors.push("rgb(200,200,255)".to_string());
            labels.push("not in most common dataset, marked".to_string());
        } else {
            if ctl.gen_opt.legend.len() == 0 {
                for s in origins.iter() {
                    let mut color = "black".to_string();
                    if ctl.gen_opt.origin_color_map.contains_key(&s.clone()) {
                        color = ctl.gen_opt.origin_color_map[s].clone();
                    }
                    colors.push(color);
                }
            } else {
                origins.clear();
                for i in 0..ctl.gen_opt.legend.len() {
                    colors.push(ctl.gen_opt.legend[i].0.clone());
                    origins.push(ctl.gen_opt.legend[i].1.clone());
                }
            }
            for i in 0..colors.len() {
                substitute_enclone_color(&mut colors[i]);
            }
            labels = origins.clone();
        }
        for s in labels.iter() {
            max_string_width = max_string_width.max(arial_width(s, font_size));
        }

        // Build the legend.

        let n = labels.len();
        const LEGEND_CIRCLE_RADIUS: usize = 4;
        let legend_height = (font_size + BOUNDARY / 2) * n + BOUNDARY;
        let legend_width = BOUNDARY as f64 * 2.5 + max_string_width;
        let mut legend_xstart = BOUNDARY as f64;
        let mut legend_ystart = actual_height + (BOUNDARY as f64) * 1.5;
        if using_shading {
            legend_xstart = legend_xstop_shading + 10.0;
            legend_ystart = 50.0;
        }
        *svg = svg.rev_before("<").to_string();
        *svg += &format!(
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" \
             style=\"fill:white;stroke:black;stroke-width:{}\" />\n",
            legend_xstart, legend_ystart, legend_width, legend_height, LEGEND_BOX_STROKE_WIDTH
        );
        for i in 0..labels.len() {
            let y = legend_ystart as f64
                + BOUNDARY as f64 * 2.5
                + ((font_size + BOUNDARY / 2) * i) as f64;
            *svg += &format!(
                "<text x=\"{}\" y=\"{}\" font-family=\"Arial\" \
                 font-size=\"{}\">{}</text>\n",
                legend_xstart + BOUNDARY as f64 * 2.0,
                y,
                font_size,
                labels[i]
            );
            *svg += &format!(
                "<circle cx=\"{}\" cy=\"{}\" r=\"{}\" fill=\"{}\" />\n",
                legend_xstart + BOUNDARY as f64,
                y - BOUNDARY as f64 / 2.0,
                LEGEND_CIRCLE_RADIUS,
                colors[i]
            );
        }
        if !using_shading {
            let new_height = legend_ystart + (legend_height + LEGEND_BOX_STROKE_WIDTH) as f64;
            set_svg_height(svg, new_height);
        } else {
            let new_width = legend_xstart + legend_width + LEGEND_BOX_STROKE_WIDTH as f64;
            set_svg_width(svg, new_width);
        }
        *svg += "</svg>";
    }

    // Output the svg file.

    if ctl.gen_opt.plot_file != "stdout".to_string() {
        let f = File::create(&ctl.gen_opt.plot_file);
        if f.is_err() {
            eprintln!(
                "\nThe file {} in your PLOT argument could not be created.\n",
                ctl.gen_opt.plot_file
            );
            std::process::exit(1);
        }
        let mut f = BufWriter::new(f.unwrap());
        fwriteln!(f, "{}", svg);
    }
    ctl.perf_stats(&t, "building svg file");
}
