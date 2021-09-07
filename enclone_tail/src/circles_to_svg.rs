// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Given a collection of circles having specified colors, create an svg string that shows the
// circles on a canvas of fixed size.  The circles are moved and resized accordingly.
// Also shades smoothed polygons.  Also add tooltip notes if requested.

use crate::polygon::*;

pub fn circles_to_svg(
    center: &Vec<(f64, f64)>,
    radius: &Vec<f64>,
    color: &Vec<String>,
    barcodes: &Vec<(usize, String)>,
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
                " tooltip=\"group_id={},clonotype_id={},barcode={}\"",
                group_index2[i] + 1,
                clonotype_index2[i] + 1,
                barcodes[i].1,
            );
        }
        if color[i] != "undefined" {
            out += &format!(
                "<circle{} cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"{}\" />\n",
                tooltipx, center[i].0, center[i].1, radius[i], color[i]
            );
        } else {
            out += &format!(
                "<circle{} cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" stroke=\"red\" \
                 stroke-width=\"0.5\" fill=\"white\" />\n",
                tooltipx, center[i].0, center[i].1, radius[i]
            );
        }
    }
    out += "</svg>\n";
    out
}
