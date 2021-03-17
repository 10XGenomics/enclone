// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use plotters::prelude::*;

pub fn plot_points(points: &Vec<(f32, f32)>, xvar: &str, yvar: &str, svg_filename: &str) {
    // Requirements.

    assert!(!points.is_empty());

    // Define parameters of the plot.

    let title = format!("{} versus {}", xvar, yvar);

    // Possibly universal constants.

    let title_font_size = 30;
    let font = "sans-serif";
    let tic_font_size = 20;
    let axis_tics = 5;
    let point_size = 4;
    let margin = 25;
    let xsize = 800;
    let ysize = 600;
    let point_color = RED;
    let range_ext = 0.05;

    // Determine the plot ranges using the extreme values of the points, extended a little bit.

    let mut xlow = points[0].0;
    let mut xhigh = points[0].0;
    let mut ylow = points[0].1;
    let mut yhigh = points[0].1;
    for i in 1..points.len() {
        xlow = xlow.min(points[i].0);
        xhigh = xhigh.max(points[i].0);
        ylow = ylow.min(points[i].1);
        yhigh = yhigh.max(points[i].1);
    }
    if xlow > 0.0 && xlow / (xhigh - xlow) < range_ext {
        xlow = 0.0;
    } else {
        xlow *= 1.0 - range_ext;
    }
    xhigh *= 1.0 + range_ext;
    if ylow > 0.0 && ylow / (yhigh - ylow) < range_ext {
        ylow = 0.0;
    } else {
        ylow *= 1.0 - range_ext;
    }
    yhigh *= 1.0 + range_ext;

    // Determine the precision for the x axis tics.

    let mut x_precision = 0;
    let m = xhigh.abs();
    if m > 0.0 {
        let extra = 2;
        x_precision = -m.log10() as usize + extra;
    }

    // Determine the precision for the y axis tics.

    let mut y_precision = 0;
    let m = yhigh.abs();
    if m > 0.0 {
        let extra = 2;
        y_precision = -m.log10() as usize + extra;
    }

    // Determine the area size for the x label.

    let x_label_area_size = (2.5 * tic_font_size as f32).round() as u32;

    // Determine the area size for the y label.

    let mut y_digits_to_left = 1;
    if yhigh.abs() > 1.0 {
        y_digits_to_left = yhigh.abs().log10().floor() as usize + 1
    }
    let y_len = y_digits_to_left + y_precision;
    let y_label_area_size = y_len as f64 * (tic_font_size as f64 * 0.6)
        + tic_font_size as f64 * 0.2
        + tic_font_size as f64;
    let y_label_area_size = y_label_area_size as u32;

    // Make the plot.

    let root = SVGBackend::new(&svg_filename, (xsize, ysize)).into_drawing_area();
    let root = root.margin(margin, margin, margin, margin);
    let mut chart = ChartBuilder::on(&root)
        .caption(&title, (font, title_font_size).into_font())
        .x_label_area_size(x_label_area_size)
        .y_label_area_size(y_label_area_size)
        .build_cartesian_2d(xlow..xhigh, ylow..yhigh)
        .unwrap();
    chart
        .configure_mesh()
        .label_style((font, tic_font_size).into_font())
        .x_labels(axis_tics)
        .y_labels(axis_tics)
        .x_label_formatter(&|x| format!("{:.1$}", x, x_precision))
        .y_label_formatter(&|x| format!("{:.1$}", x, y_precision))
        .x_desc(xvar)
        .y_desc(yvar)
        .draw()
        .unwrap();
    chart
        .draw_series(PointSeries::of_element(
            points.clone(),
            point_size,
            &point_color,
            &|c, s, st| {
                return EmptyElement::at(c) + Circle::new((0, 0), s, st.filled());
            },
        ))
        .unwrap();
}
