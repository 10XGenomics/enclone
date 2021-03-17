// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Plot experiment.

use plotters::prelude::*;
use string_utils::*;

fn main() {
    // Define parameters of the plot.

    let title = "This is our first plot";
    let font = "sans-serif";
    let title_font_size = 40;
    let tic_font_size = 20;
    let plotfile = "/mnt/home/david.jaffe/public_html/plotz.svg";
    let axis_tics = 5;
    let point_size = 5;
    let margin = 25;
    let xsize = 800;
    let ysize = 600;
    let y_label_area_size = 60;
    let point_color = RED;
    let x_precision = 2;
    let y_precision = 1;
    let xlabel = "Count";
    let ylabel = "Bucket";

    // Load data, which are at the end of this file.

    // let points: Vec<(f32, f32)> = vec![(0.0, -2.111), (7.0, 2.0), (4.0, 5.0), (8.0, 21.3)];
    let mut points = Vec::<(f32, f32)>::new();
    let f = include_str!["plot_test.rs"];
    let mut in_data = false;
    for line in f.lines() {
        if line == "// DATA" {
            in_data = true;
        } else if in_data {
            points.push((
                line.between(" ", ",").force_f64() as f32,
                line.after(",").force_f64() as f32,
            ));
        }
    }

    // Requirements.

    assert!(!points.is_empty());

    // Determine the area size for the x label.

    let x_label_area_size = (2.5 * tic_font_size as f32).round() as u32;

    // Determine the plot ranges using the extreme values of the points.

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

// DATA
// 0.00001,12
// 0.00001,11
// 0.000007,9
// 0.00001,22
// 0.00001,7
// 1.07042E-05,11
// 0.00001,10
// 0.00001,11
// 0.000472314,7
// 0.000261191,31
// 0.00001,6
// 0.00001,11
// 0.00106997,17
// 0.00001,3
// 0.00001,23
// 0.00001,21
// 4.54121E-05,7
// 0.00001,14
// 0.00001,21
// 0.000210643,15
// 6.08361E-05,26
// 0.00001,41
// 0.000254493,3
// 0.000325303,4
// 1.51725E-05,13
// 1.51725E-05,5
// 0.00001,14
// 0.00001,11
// 3.82693E-05,27
// 0.000228595,17
// 3.58783E-05,18
// 1.36806E-05,27
// 3.58783E-05,3
// 0.00001,16
// 3.17014E-05,18
// 0.000556859,2
// 0.000713161,10
// 0.000341469,18
// 0.001579333,10
// 0.000929721,3
// 2.57335E-05,10
// 0.000755029,7
// 0.00001,19
// 1.61381E-05,11
// 0.000618834,13
// 0.001542719,0
// 0.000417928,9
// 0.00019921,10
// 0.000386902,15
// 0.00001,16
// 0.000282927,13
// 0.00001,34
// 0.000238067,18
// 0.000431399,17
// 0.00001,22
// 0.00001,2
// 0.00078585,16
// 2.4144E-05,8
// 1.45446E-05,11
// 0.00001,18
// 3.50504E-05,13
// 0.000184403,11
// 0.00001,24
