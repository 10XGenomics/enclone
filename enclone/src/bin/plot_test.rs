// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Plot experiment.

use plotters::prelude::*;

fn main() {

    let points = vec![(0.0, 0.0), (5.0, 5.0), (8.0, 20.0)];
    let title = "This is our first plot";
    let title_font_size = 40;
    let plotfile = "/mnt/home/david.jaffe/public_html/plotz.svg";
    let axis_tics = 5;
    let point_size = 5;
    let xlow = 0 as f32;
    let xhigh = 10 as f32;
    let ylow = 0 as f32;
    let yhigh = 20 as f32;

    let root = SVGBackend::new(&plotfile, (800, 600)).into_drawing_area();
    let root = root.margin(10, 10, 10, 10);
    let mut chart = ChartBuilder::on(&root)
        .caption(&title, ("sans-serif", title_font_size).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        .build_cartesian_2d(xlow..xhigh, ylow..yhigh).unwrap();

    chart
        .configure_mesh()
        .x_labels(axis_tics)
        .y_labels(axis_tics)
        // We can also change the format of the label text
        .y_label_formatter(&|x| format!("{:.3}", x))
        .draw().unwrap();

    chart.draw_series(PointSeries::of_element(
        points.clone(),
        point_size,
        &RED,
        &|c, s, st| {
            return EmptyElement::at(c)
            + Circle::new((0,0),s,st.filled());
        },
    )).unwrap();
}
