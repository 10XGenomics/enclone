// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// The purpose of this file is the function plot_clonotypes.  It plots clonotypes as partial
// hexagonal closest packings.  This is visually kind of satisfying, but also a bit weird looking.
// In some cases, by eye, you can see rounder forms that could be created by relocating some of
// the cells.

use crate::assign_cell_color::{VAR_HIGH, VAR_LOW};
use crate::circles_to_svg::circles_to_svg;
use crate::colors::{turbo_color_names, TURBO_SRGB_BYTES};
use crate::group_colors::make_group_colors;
use crate::legend::add_legend_for_color_by_variable;
use crate::pack_circles::pack_circles;
use crate::plot_utils::build_clusters;
use crate::polygon::{enclosing_polygon, Polygon};
use crate::string_width::arial_width;
use crate::{get_svg_height, set_svg_height, set_svg_width, substitute_enclone_color, BOUNDARY};
use ansi_escape::print_color13;
use enclone_core::cell_color::CellColor;
use enclone_core::convert_svg_to_png::convert_svg_to_png;
use enclone_core::defs::{EncloneControl, ExactClonotype, PlotOpt, POUT_SEP};
use io_utils::{fwriteln, open_for_read, open_for_write_new};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;
use string_utils::TextUtils;
use vdj_ann::refx::RefData;
use vector_utils::{
    bin_position, next_diff1_2, next_diff1_3, position, reverse_sort, unique_sort, VecUtils,
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn plot_clonotypes(
    ctl: &EncloneControl,
    plot_opt: &PlotOpt, // overrides ctl
    refdata: &RefData,
    // exacts: One entry for each clonotype.
    exacts: &Vec<Vec<usize>>,
    exact_clonotypes: &Vec<ExactClonotype>,
    out_datas: &Vec<Vec<HashMap<String, String>>>,
    // groups: There is one entry for each group of clonotypes.  The first entries of the inner
    // vectors indexes into exacts, and the second entry (String) is not used here.
    groups: &Vec<Vec<(i32, String)>>,
    svg: &mut String,
) -> Result<(), String> {
    let t = Instant::now();
    if plot_opt.plot_file.is_empty() {
        return Ok(());
    }
    if exacts.is_empty() {
        return Err("\nThere are no clonotypes to plot, giving up.\n".to_string());
    }

    let mut const_names = Vec::<String>::new();
    for id in refdata.cs.iter() {
        if refdata.rtype[*id] == 0 {
            const_names.push(refdata.name[*id].clone());
        }
    }
    unique_sort(&mut const_names);
    if !plot_opt.plot_by_isotype_color.is_empty() {
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
        return Err(
            "\nCurrently PLOT_BY_ISOTYPE only works if there are at most 12 constant \
            region names.  If this is a problem, please let us know and we will generalize it.\n"
                .to_string(),
        );
    }

    // Get origins.

    let mut origins = Vec::<String>::new();
    for i in 0..exacts.len() {
        for j in 0..exacts[i].len() {
            let ex = &exact_clonotypes[exacts[i][j]];
            for k in 0..ex.clones.len() {
                if ex.clones[k][0].origin_index.is_some() {
                    let s = &ctl.origin_info.origin_list[ex.clones[k][0].origin_index.unwrap()];
                    origins.push(s.clone());
                }
            }
        }
    }
    unique_sort(&mut origins);

    // Build one cluster for each clonotype.

    let mut clusters = build_clusters(
        ctl,
        plot_opt,
        refdata,
        exacts,
        exact_clonotypes,
        out_datas,
        &const_names,
    );
    let mut radii = Vec::<f64>::new();
    for i in 0..clusters.len() {
        radii.push(clusters[i].radius);
    }

    // Set group specification, if CLONOTYPE_GROUP_NAMES was specified.
    // Note that CLONOTYPE_GROUP_NAMES is probably broken now.

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
                    return Err(
                        "\nThe first line in your CLONOTYPE_GROUP_NAMES file contains a \
                        nonprinting character.\n"
                            .to_string(),
                    );
                }
            }
            let fields = s.split(',').collect::<Vec<&str>>();
            if first {
                let p = position(&fields, &"group_id");
                if p < 0 {
                    return Err(
                        "\nThe CLONOTYPE_GROUP_NAMES file does not have a group_id field.\n"
                            .to_string(),
                    );
                }
                group_id_field = p as usize;
                let p = position(&fields, &"new_group_name");
                if p < 0 {
                    return Err("\nThe CLONOTYPE_GROUP_NAMES file does not have a \
                         new_group_name field.\n"
                        .to_string());
                }
                new_group_name_field = p as usize;
                first = false;
            } else {
                let group_id = &fields[group_id_field];
                if group_id.parse::<usize>().is_err() || group_id.force_usize() == 0 {
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
                        clusters[i].coords.len(),
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

        // Define group ids and colors.

        group_id.clear();
        for i in 0..new_group_names.len() {
            if new_group_names[i].is_some() {
                let p = position(&group_name, new_group_names[i].as_ref().unwrap());
                group_id.push(p as usize);
            }
        }
        group_color = make_group_colors(group_name.len());
    }
    let ngroups = group_color.len(); // THESE ARE SHADING GROUPS!
    ctl.perf_stats(&t, "in preamble to plotting clonotypes");

    // Traverse the shading groups.  In the default case, there is just one!!!!!!!!!!!!!!!!!!!!!!!!

    let t = Instant::now();
    let using_shading = ngroups > 1 || !group_color[0].is_empty();
    let mut blacklist = Vec::<Polygon>::new();
    let mut shades = Vec::<Polygon>::new();
    let mut shade_colors = Vec::<String>::new();
    let mut shade_enclosures = Vec::<Polygon>::new();
    let mut centers = vec![(0.0, 0.0); radii.len()];
    for g in 0..ngroups {
        // Gather the group.  In the default case, ids = 0..(number of clusters).

        let mut ids = Vec::<usize>::new();
        let mut radiix = Vec::<f64>::new();
        for i in 0..group_id.len() {
            if group_id[i] == g {
                ids.push(i);
                radiix.push(radii[i]);
            }
        }

        // Find circle centers.

        let mut centersx = vec![(0.0, 0.0); radiix.len()];
        let mut xshift = vec![0.0; radiix.len()];
        if plot_opt.split_plot_by_origin {
            let passes = ctl.origin_info.origin_list.len();
            let mut xstart = 0.0;
            for pass in 0..passes {
                let mut radiiy = Vec::<f64>::new();
                let mut indices = Vec::<usize>::new();
                for i in 0..radiix.len() {
                    let li = clusters[i].barcodes[0].0;
                    let p =
                        bin_position(&ctl.origin_info.origin_list, &ctl.origin_info.origin_id[li]);
                    if pass != p as usize {
                        continue;
                    }
                    radiiy.push(radiix[i]);
                    indices.push(i);
                }
                let centersy = pack_circles(&radiiy, &blacklist, plot_opt.plot_quad);
                let mut left = 0.0_f64;
                for j in 0..centersy.len() {
                    left = left.max(-centersy[j].0 + radiiy[j]);
                }
                xstart += left;
                for j in 0..centersy.len() {
                    centersx[indices[j]] = (centersy[j].0 + xstart, centersy[j].1);
                    xshift[indices[j]] = xstart;
                }
                let mut right = 0.0_f64;
                for j in 0..centersy.len() {
                    right = right.max(centersy[j].0 + radiiy[j]);
                }
                const HSEP: f64 = 10.0;
                xstart += right + HSEP;
            }
        } else {
            centersx = pack_circles(&radiix, &blacklist, plot_opt.plot_quad);
        }
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
        //
        // We do the split_plot_by_origin case second.  It is a more complicated version of the
        // same algorithm.  The second part definitely does not work with grouping.

        if !plot_opt.split_plot_by_origin {
            let mut ccc = Vec::<(usize, String, usize)>::new(); // (cluster size, color, index)
            let mut clusters2 = clusters.clone();
            for i in 0..ids.len() {
                let id = ids[i];
                let mut c = clusters[id].colors.clone();
                unique_sort(&mut c);
                if c.solo() {
                    // Note confusion here between the last argument, i, and clusters[i].2:
                    // We had clusters[i].colors.len() below and that appears to have been a bug.
                    ccc.push((clusters[id].colors.len(), c[0].clone(), i));
                }
            }
            ccc.sort();
            let mut i = 0;
            while i < ccc.len() {
                // On a given iteration of the while loop, we process all the constant-color
                // clusters that have the same size.  First we do the clusters that contains
                // just one cell, and so forth.
                let j = next_diff1_3(&ccc, i as i32) as usize;
                let mut angle = vec![(0.0, 0); j - i];
                for k in i..j {
                    let id = ccc[k].2;
                    angle[k - i] = (centersx[id].1.atan2(centersx[id].0), id);
                }
                angle.sort_by(|a, b| a.partial_cmp(b).unwrap());
                for k in i..j {
                    let new_id = angle[k - i].1;
                    let id = ccc[k].2;
                    clusters2[ids[new_id]].colors = clusters[id].colors.clone();
                    clusters2[ids[new_id]].clonotype_index = clusters[id].clonotype_index;
                    clusters2[ids[new_id]].barcodes = clusters[id].barcodes.clone();
                }
                i = j;
            }
            clusters = clusters2;
        } else {
            // WORK IN PROGRESS

            let mut clusters2 = clusters.clone();
            let mut centersp = centersx.clone();
            for i in 0..centersp.len() {
                centersp[i].0 -= xshift[i];
            }
            let passes = ctl.origin_info.origin_list.len();
            for pass in 0..passes {
                let mut ccc = Vec::<(usize, String, usize)>::new();
                let mut this = Vec::<usize>::new();
                for i in 0..radiix.len() {
                    let li = clusters[i].barcodes[0].0;
                    let p =
                        bin_position(&ctl.origin_info.origin_list, &ctl.origin_info.origin_id[li]);
                    if pass != p as usize {
                        continue;
                    }
                    let mut c = clusters[i].colors.clone();
                    unique_sort(&mut c);
                    if c.solo() {
                        ccc.push((clusters[i].colors.len(), c[0].clone(), i));
                        this.push(i);
                    }
                }
                ccc.sort();
                let mut i = 0;
                while i < ccc.len() {
                    let j = next_diff1_3(&ccc, i as i32) as usize;
                    let mut angle = vec![(0.0, 0); j - i];
                    for k in i..j {
                        let id = ccc[k].2;
                        angle[k - i] = (centersp[id].1.atan2(centersp[id].0), id);
                    }
                    angle.sort_by(|a, b| a.partial_cmp(b).unwrap());

                    for k in i..j {
                        let new_id = angle[k - i].1;
                        let id = ccc[k].2;
                        clusters2[new_id].colors = clusters[id].colors.clone();
                        clusters2[new_id].clonotype_index = clusters[id].clonotype_index;
                        clusters2[new_id].barcodes = clusters[id].barcodes.clone();
                    }
                    i = j;
                }
            }
            clusters = clusters2;
        }

        // Finish turbo color translation.

        let tcn = turbo_color_names();
        for i in 0..clusters.len() {
            for j in 0..clusters[i].colors.len() {
                if clusters[i].colors[j].starts_with("turbo-") {
                    let n = bin_position(&tcn, &clusters[i].colors[j]);
                    let c = &TURBO_SRGB_BYTES[n as usize];
                    let color = format!("rgb({},{},{})", c[0], c[1], c[2]);
                    clusters[i].colors[j] = color;
                }
            }
        }
    }
    ctl.perf_stats(&t, "plotting clonotypes");

    // Build the svg file.

    let t = Instant::now();
    for i in 0..clusters.len() {
        for j in 0..clusters[i].coords.len() {
            clusters[i].coords[j].0 += centers[i].0;
            clusters[i].coords[j].1 += centers[i].1;
        }
    }
    let mut center = Vec::<(f64, f64)>::new();
    let mut radius = Vec::<f64>::new();
    let mut color = Vec::<String>::new();
    let mut barcodes = Vec::<(usize, String)>::new();
    let mut group_index = HashMap::<usize, usize>::new();
    let mut clonotype_index = Vec::<usize>::new();
    for i in 0..groups.len() {
        for j in 0..groups[i].len() {
            group_index.insert(groups[i][j].0 as usize, i);
            clonotype_index.push(j);
        }
    }
    let mut group_index2 = Vec::<usize>::new();
    let mut clonotype_index2 = Vec::<usize>::new();
    for i in 0..clusters.len() {
        for j in 0..clusters[i].colors.len() {
            color.push(clusters[i].colors[j].clone());
            barcodes.push(clusters[i].barcodes[j].clone());
            center.push((clusters[i].coords[j].0, clusters[i].coords[j].1));
            radius.push(1.0);
            let ind = clusters[i].clonotype_index;
            group_index2.push(group_index[&ind]);
            clonotype_index2.push(clonotype_index[ind]);
        }
    }
    const WIDTH: usize = 400;
    const HEIGHT: usize = 400;

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

    // Implement HONEY_IN.

    if plot_opt.honey_in.is_some() {
        let mut honey_map = HashMap::<(usize, String), (f64, f64)>::new();
        let f = open_for_read![&plot_opt.honey_in.as_ref().unwrap()];
        for line in f.lines() {
            let s = line.unwrap();
            let fields = s.split(',').collect::<Vec<&str>>();
            if fields.len() != 4 {
                return Err("\nHONEY_IN file incorrectly formatted.\n".to_string());
            }
            if fields[0].parse::<usize>().is_err() {
                return Err("\nHONEY_IN file incorrectly formatted.\n".to_string());
            }
            if fields[2].parse::<f64>().is_err() || fields[3].parse::<f64>().is_err() {
                return Err("\nHONEY_IN file incorrectly formatted.\n".to_string());
            }
            honey_map.insert(
                (fields[0].force_usize(), fields[1].to_string()),
                (fields[2].force_f64(), fields[3].force_f64()),
            );
        }
        if honey_map.len() != barcodes.len() {
            return Err(format!(
                "\nHONEY_IN file appears to have come from data having {} \
                cells, whereas the current data have {} cells.\n",
                honey_map.len(),
                barcodes.len(),
            ));
        }
        for i in 0..barcodes.len() {
            if !honey_map.contains_key(&barcodes[i]) {
                return Err(
                    "\nHONEY_IN file appears to have come from different data.\n".to_string(),
                );
            }
            center[i] = honey_map[&barcodes[i]];
        }
    }

    // Implement HONEY_OUT.

    if !plot_opt.honey_out.is_empty() {
        let mut honey_map = Vec::<((usize, String), f64, f64)>::new();
        for i in 0..barcodes.len() {
            honey_map.push((barcodes[i].clone(), center[i].0, center[i].1));
        }
        honey_map.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut f = open_for_write_new![&plot_opt.honey_out];
        for x in honey_map.iter() {
            fwriteln!(f, "{},{},{},{}", x.0 .0, x.0 .1, x.1, x.2);
        }
    }

    // Determine if we are coloring cells by variable value.

    let mut by_var = false;
    let mut var = String::new();
    match plot_opt.cell_color {
        CellColor::ByVariableValue(ref x) => {
            by_var = true;
            var = x.var.clone();
        }
        _ => {}
    };

    // If so, look up the variable value.

    let mut barcode_to_var_value = HashMap::<(usize, String), String>::new();
    if by_var {
        for i in 0..out_datas.len() {
            for j in 0..out_datas[i].len() {
                if out_datas[i][j].contains_key(&var) {
                    let val_list = &out_datas[i][j][&var];
                    let vals = val_list.split(POUT_SEP).collect::<Vec<&str>>();
                    let ex = &exact_clonotypes[exacts[i][j]];
                    for k in 0..ex.ncells() {
                        let val = if vals.len() > 1 { &vals[k] } else { &vals[0] };
                        let li = ex.clones[k][0].dataset_index;
                        let bc = &ex.clones[k][0].barcode;
                        barcode_to_var_value.insert((li, bc.clone()), val.to_string());
                    }
                }
            }
        }
    }

    // Generate svg.

    *svg = circles_to_svg(
        &center,
        &radius,
        &color,
        &barcodes,
        &shades,
        &shade_colors,
        &shade_enclosures,
        &group_index2,
        &clonotype_index2,
        WIDTH,
        HEIGHT,
        BOUNDARY,
        plot_opt.plot_file == "gui" || plot_opt.plot_file == "gui_stdout",
        by_var,
        var,
        &barcode_to_var_value,
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
    set_svg_height(svg, actual_height + BOUNDARY as f64);

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
        let height = get_svg_height(svg).max(legend_height_plus);
        set_svg_height(svg, height + BOUNDARY as f64);
        *svg += "</svg>";
    }

    // Add legend for color by variable.

    if by_var && plot_opt.use_legend {
        add_legend_for_color_by_variable(plot_opt, svg, &color, actual_width, actual_height);
    } else if plot_opt.use_legend
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
            if plot_opt.legend.is_empty() {
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
            if new_height > get_svg_height(svg) {
                set_svg_height(svg, new_height);
            }
            set_svg_width(svg, new_width);
        }
        *svg += "</svg>";
    }
    VAR_LOW.lock().unwrap().clear();
    VAR_HIGH.lock().unwrap().clear();

    // Output the svg or png file.

    if plot_opt.plot_file == "stdout.png" {
        let png = convert_svg_to_png(svg.as_bytes(), 2000);
        std::io::stdout().write_all(&png).unwrap();
    } else if plot_opt.plot_file != "stdout"
        && plot_opt.plot_file != "gui"
        && plot_opt.plot_file != "gui_stdout"
    {
        let f = File::create(&plot_opt.plot_file);
        if f.is_err() {
            return Err(format!(
                "\nThe file {} in your PLOT or HONEY argument could not be created.\n",
                plot_opt.plot_file
            ));
        }
        let mut f = f.unwrap();
        if !plot_opt.plot_file.ends_with(".png") {
            let mut f = BufWriter::new(f);
            fwriteln!(f, "{}", svg);
        } else {
            let mut width = 2000;
            if plot_opt.png_width.is_some() {
                width = plot_opt.png_width.unwrap();
            }
            let png = convert_svg_to_png(svg.as_bytes(), width as u32);
            f.write_all(&png).unwrap();
        }
    }
    ctl.perf_stats(&t, "building svg file");
    Ok(())
}
