// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// The purpose of this file is the function plot_clonotypes.  It plots clonotypes as partial
// hexagonal closest packings.  This is visually kind of satisfying, but also a bit weird looking.
// In some cases, by eye, you can see rounder forms that could be created by relocating some of
// the cells.

use crate::hex::*;
use crate::pack_circles::*;
use crate::polygon::*;
use crate::string_width::*;
use ansi_escape::*;
use enclone_core::defs::*;
use io_utils::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
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

fn get_svg_height(svg: &String) -> f64 {
    svg.between("height=\"", "\"").force_f64()
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given a collection of circles having specified colors, create an svg string that shows the
// circles on a canvas of fixed size.  The circles are moved and resized accordingly.
// Also shades smoothed polygons.  Also add tooltip notes if requested.

fn circles_to_svg(
    center: &Vec<(f64, f64)>,
    radius: &Vec<f64>,
    color: &Vec<String>,
    shades: &Vec<Polygon>,
    shade_colors: &Vec<String>,
    shade_enclosures: &Vec<Polygon>,
    group_index2: &Vec<usize>,
    clonotype_index2: &Vec<usize>,
    width: usize,
    height: usize,
    boundary: usize,
    tooltip: bool,
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
        const BOUNDING_CURVE_BOUND: f64 = 15.0; // must be smaller than POLYGON_ENLARGEMENT
        out += &format!("{}", p.catmull_bezier_bounded_svg(BOUNDING_CURVE_BOUND));
        out += "\" ";
        out += &format!("fill=\"{}\"\n", shade_colors[g]);
        out += " stroke=\"rgb(150,150,150)\"";
        out += " stroke-width=\"0.2\"";
        out += "/>\n";
    }
    for i in 0..center.len() {
        let mut tooltipx = String::new();
        if tooltip {
            tooltipx = format!(
                " tooltip=\"group_id={},clonotype_id={}\"",
                group_index2[i] + 1,
                clonotype_index2[i] + 1,
            );
        }
        out += &format!(
            "<circle{} cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"{}\" />\n",
            tooltipx, center[i].0, center[i].1, radius[i], color[i]
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
    plot_opt: &PlotOpt, // overrides ctl
    refdata: &RefData,
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    groups: &Vec<Vec<(i32, String)>>,
    svg: &mut String,
) -> Result<(), String> {
    let t = Instant::now();
    if plot_opt.plot_file.is_empty() {
        return Ok(());
    }
    if exacts.is_empty() {
        return Err(format!("\nThere are no clonotypes to plot, giving up.\n"));
    }

    let mut const_names = Vec::<String>::new();
    for id in refdata.cs.iter() {
        if refdata.rtype[*id] == 0 {
            const_names.push(refdata.name[*id].clone());
        }
    }
    unique_sort(&mut const_names);
    if plot_opt.plot_by_isotype_color.len() > 0 {
        if const_names.len() + 1 > plot_opt.plot_by_isotype_color.len() {
            return Err(format!(
                "\nUsing the PLOT_BY_ISOTYPE_COLOR argument, you specified {} colors, \
                but there are {} constant region\nnames, and one more color is needed for the \
                \"undetermined\" case.  Please add more colors.\n",
                plot_opt.plot_by_isotype_color.len(),
                const_names.len()
            ));
        }
    } else if plot_opt.plot_by_isotype && const_names.len() > 12 {
        return Err(format!(
            "\nCurrently PLOT_BY_ISOTYPE only works if there are at most 12 constant \
            region names.  If this is a problem, please let us know and we will generalize it.\n"
        ));
    }
    let mut clusters = Vec::<(Vec<String>, Vec<(f64, f64)>, usize)>::new();
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
        if plot_opt.plot_by_mark {
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

                if plot_opt.plot_by_isotype {
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
                        // Note that is possible for p to be -1 in the following.  This is known
                        // to happen if a heavy chain V gene is on the same contig as a light
                        // chain C gene (which may be an artifact).  There is an example in
                        // enclone_main/testx/inputs/flaky.
                        let p = bin_position(&const_names, &c);
                        color_id = (1 + p) as usize;
                    }
                    if plot_opt.plot_by_isotype_color.is_empty() {
                        let x = print_color13(color_id);
                        color = format!("rgb({},{},{})", x.0, x.1, x.2);
                    } else {
                        color = plot_opt.plot_by_isotype_color[color_id].clone();
                    }

                // Determine color for PLOT_BY_MARK.
                } else if plot_opt.plot_by_mark {
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
        clusters.push((colors, coords, i));
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
                    return Err(format!(
                        "\nThe first line in your CLONOTYPE_GROUP_NAMES file contains a \
                        nonprinting character.\n"
                    ));
                }
            }
            let fields = s.split(',').collect::<Vec<&str>>();
            if first {
                let p = position(&fields, &"group_id");
                if p < 0 {
                    return Err(format!(
                        "\nThe CLONOTYPE_GROUP_NAMES file does not have a group_id field.\n"
                    ));
                }
                group_id_field = p as usize;
                let p = position(&fields, &"new_group_name");
                if p < 0 {
                    return Err(format!(
                        "\nThe CLONOTYPE_GROUP_NAMES file does not have a \
                         new_group_name field.\n"
                    ));
                }
                new_group_name_field = p as usize;
                first = false;
            } else {
                let group_id = &fields[group_id_field];
                if !group_id.parse::<usize>().is_ok() || group_id.force_usize() == 0 {
                    return Err(format!(
                        "\nThe group_id {} in your CLONOTYPE_GROUP_NAMES file is not a \
                        positive integer.\n",
                        group_id
                    ));
                }
                let group_id = group_id.force_usize() - 1;
                if group_id > radii.len() {
                    return Err(format!(
                        "\nThe group_id {} in your CLONOTYPE_GROUP_NAMES file is larger \
                        than the number of clonotypes, which is {}.\n",
                        group_id,
                        radii.len()
                    ));
                }
                let new_group_name = fields[new_group_name_field].to_string();
                if new_group_names[group_id].is_none() {
                    new_group_names[group_id] = Some(new_group_name.clone());
                } else if *new_group_names[group_id].as_ref().unwrap() != new_group_name {
                    return Err(format!(
                        "\nThe group_id {} in your CLONOTYPE_GROUP_NAMES file is assigned \
                        the different new_group_names {} and {}.\n",
                        group_id,
                        new_group_names[group_id].as_ref().unwrap(),
                        new_group_name,
                    ));
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

        let centersx = pack_circles(&radiix, &blacklist, plot_opt.plot_quad);
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

            // Build an enlarged polygon that includes the smoothed polygonal curve.

            const POLYGON_ENLARGEMENT: f64 = 22.5; // must be larger than BOUNDING_CURVE_BOUND
            p.enlarge(POLYGON_ENLARGEMENT);
            p.precompute();
            shade_enclosures.push(p.clone());
            blacklist.push(p);
        }

        // Reorganize constant-color clusters so that like-colored clusters are proximate,
        // We got this idea from Ganesh Phad, who showed us a picture!  The primary effect is on
        // single-cell clonotypes.

        let mut honey_map_in = Vec::<(usize, usize)>::new();
        if ctl.plot_opt.honey_in.is_some() {
            let f = open_for_read![&ctl.plot_opt.honey_in.as_ref().unwrap()];
            for line in f.lines() {
                let s = line.unwrap();
                if !s.contains(",")
                    || !s.before(",").parse::<usize>().is_ok()
                    || !s.after(",").parse::<usize>().is_ok()
                {
                    return Err(format!("\nHONEY_IN file incorrectly formatted.\n"));
                }
                honey_map_in.push((s.before(",").force_usize(), s.after(",").force_usize()));
            }
        }
        let mut honey_map_out = Vec::<(usize, usize)>::new();
        if ctl.plot_opt.honey_in.is_none() {
            let mut ccc = Vec::<(usize, String, usize)>::new(); // (cluster size, color, index)
            for i in 0..ids.len() {
                let id = ids[i];
                let mut c = clusters[id].0.clone();
                unique_sort(&mut c);
                if c.solo() {
                    ccc.push((clusters[i].0.len(), c[0].clone(), i));
                } else {
                    honey_map_out.push((i, i));
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
                    honey_map_out.push((ids[new_id], ccc[k].2));
                    for u in 0..clusters[ids[new_id]].0.len() {
                        clusters[ids[new_id]].0[u] = ccc[k].1.clone();
                    }
                }
                i = j;
            }
            if ctl.plot_opt.honey_out.len() > 0 {
                let mut f = open_for_write_new![&ctl.plot_opt.honey_out];
                for i in 0..honey_map_out.len() {
                    fwriteln!(f, "{},{}", honey_map_out[i].0, honey_map_out[i].1);
                }
            }
        } else {
            if honey_map_in.len() != clusters.len() {
                return Err(format!(
                    "\nHONEY_IN file appears to come from data having {} clusters, \
                    whereas the current data have {} clusters.\n",
                    honey_map_in.len(),
                    clusters.len(),
                ));
            }
            let mut clusters2 = clusters.clone();
            for i in 0..clusters.len() {
                clusters2[honey_map_in[i].0].0 = clusters[honey_map_in[i].1].0.clone();
            }
            clusters = clusters2;
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
    let mut group_index = Vec::<usize>::new();
    let mut clonotype_index = Vec::<usize>::new();
    for i in 0..groups.len() {
        for j in 0..groups[i].len() {
            group_index.push(i);
            clonotype_index.push(j);
        }
    }
    let mut group_index2 = Vec::<usize>::new();
    let mut clonotype_index2 = Vec::<usize>::new();
    for i in 0..clusters.len() {
        for j in 0..clusters[i].0.len() {
            color.push(clusters[i].0[j].clone());
            center.push((clusters[i].1[j].0, clusters[i].1[j].1));
            radius.push(1.0);
            let ind = clusters[i].2;
            group_index2.push(group_index[ind]);
            clonotype_index2.push(clonotype_index[ind]);
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
        &group_index2,
        &clonotype_index2,
        WIDTH,
        HEIGHT,
        BOUNDARY,
        ctl.plot_opt.plot_file == "gui" || ctl.plot_opt.plot_file == "gui_stdout",
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
    set_svg_width(svg, actual_width + BOUNDARY as f64);

    // Add legend for shading.

    let mut font_size = 20;
    const LEGEND_BOX_STROKE_WIDTH: usize = 2;
    let mut legend_xstop_shading = 0.0;
    if using_shading {
        font_size = 16;
        let n = ngroups;
        let mut max_string_width = 0.0f64;
        for s in group_name.iter() {
            max_string_width = max_string_width.max(arial_width(s, font_size as f64));
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
        let legend_height_plus = legend_height + vsep + 15.0;
        if legend_height_plus > get_svg_height(&svg) {
            set_svg_height(svg, legend_height_plus);
        }
        *svg += "</svg>";
    }

    // Add main legend.

    if plot_opt.use_legend
        || (plot_opt.plot_by_isotype && !plot_opt.plot_by_isotype_nolegend)
        || plot_opt.plot_by_mark
    {
        let (mut colors, mut labels) = (Vec::<String>::new(), Vec::<String>::new());
        let mut max_string_width = 0.0f64;
        if plot_opt.plot_by_isotype {
            for i in 0..const_names.len() {
                labels.push(const_names[i].clone());
                let color_id = i + 1;
                if plot_opt.plot_by_isotype_color.is_empty() {
                    let x = print_color13(color_id);
                    let color = format!("rgb({},{},{})", x.0, x.1, x.2);
                    colors.push(color);
                } else {
                    let color = plot_opt.plot_by_isotype_color[color_id].clone();
                    colors.push(color);
                }
            }
            labels.push("undetermined".to_string());
            let color_id = 0;
            if plot_opt.plot_by_isotype_color.is_empty() {
                let x = print_color13(color_id);
                let color = format!("rgb({},{},{})", x.0, x.1, x.2);
                colors.push(color);
            } else {
                let color = plot_opt.plot_by_isotype_color[color_id].clone();
                colors.push(color);
            }
        } else if plot_opt.plot_by_mark {
            colors.push("red".to_string());
            labels.push("in most common dataset, !marked".to_string());
            colors.push("rgb(255,200,200)".to_string());
            labels.push("in most common dataset, marked".to_string());
            colors.push("blue".to_string());
            labels.push("not in most common dataset, !marked".to_string());
            colors.push("rgb(200,200,255)".to_string());
            labels.push("not in most common dataset, marked".to_string());
        } else {
            if plot_opt.legend.len() == 0 {
                for s in origins.iter() {
                    let mut color = "black".to_string();
                    if ctl.gen_opt.origin_color_map.contains_key(&s.clone()) {
                        color = ctl.gen_opt.origin_color_map[s].clone();
                    }
                    colors.push(color);
                }
            } else {
                origins.clear();
                for i in 0..plot_opt.legend.len() {
                    colors.push(plot_opt.legend[i].0.clone());
                    origins.push(plot_opt.legend[i].1.clone());
                }
            }
            for i in 0..colors.len() {
                substitute_enclone_color(&mut colors[i]);
            }
            labels = origins.clone();
        }
        for s in labels.iter() {
            max_string_width = max_string_width.max(arial_width(s, font_size as f64));
        }

        // Build the legend.

        let n = labels.len();
        const LEGEND_CIRCLE_RADIUS: usize = 4;
        let legend_height = (font_size + BOUNDARY / 2) * n + BOUNDARY;
        let legend_width = BOUNDARY as f64 * 2.5 + max_string_width;
        let mut legend_xstart = actual_width + 20.0;
        let mut legend_ystart = BOUNDARY as f64;
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
        let new_height = actual_height.max(legend_height as f64) + BOUNDARY as f64;
        let new_width = actual_width + legend_width as f64 + 20.0 + BOUNDARY as f64;
        if !using_shading {
            set_svg_height(svg, new_height);
            set_svg_width(svg, new_width);
        } else {
            if new_height > get_svg_height(&svg) {
                set_svg_height(svg, new_height);
            }
            set_svg_width(svg, new_width);
        }
        *svg += "</svg>";
    }

    // Output the svg file.

    if plot_opt.plot_file != "stdout"
        && plot_opt.plot_file != "gui"
        && plot_opt.plot_file != "gui_stdout"
    {
        let f = File::create(&plot_opt.plot_file);
        if f.is_err() {
            return Err(format!(
                "\nThe file {} in your PLOT argument could not be created.\n",
                plot_opt.plot_file
            ));
        }
        let mut f = BufWriter::new(f.unwrap());
        fwriteln!(f, "{}", svg);
    }
    ctl.perf_stats(&t, "building svg file");
    Ok(())
}
