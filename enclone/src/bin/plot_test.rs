// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Plot experiment.

use plotters::prelude::*;

fn main() {

    // Define parameters of the plot.

    let points = vec![(0.0, 0.0), (5.0, 5.0), (8.0, 20.0)];
    let title = "This is our first plot";
    let font = "sans-serif";
    let title_font_size = 40;
    let tic_font_size = 20;
    let plotfile = "/mnt/home/david.jaffe/public_html/plotz.svg";
    let axis_tics = 5;
    let point_size = 5;
    let xlow = 0 as f32;
    let xhigh = 10 as f32;
    let ylow = 0 as f32;
    let yhigh = 20 as f32;
    let margin = 25;
    let xsize = 800;
    let ysize = 600;
    let x_label_area_size = 40;
    let y_label_area_size = 60;
    let point_color = RED;
    let x_precision = 2;
    let y_precision = 1;
    let xlabel = "Count";
    let ylabel = "Bucket";

    // Make the plot.

    let root = SVGBackend::new(&plotfile, (xsize, ysize)).into_drawing_area();
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
        .y_desc(xlabel)
        .x_desc(ylabel)
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
