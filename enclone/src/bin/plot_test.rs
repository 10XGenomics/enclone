// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Plot experiment.

use plotters::prelude::*;
use string_utils::*;

fn main() {
    // Define parameters of the plot.

    let title = "wt_koff versus CD27_ab";
    let title_font_size = 30;
    let font = "sans-serif";
    let tic_font_size = 20;
    let plotfile = "/mnt/home/david.jaffe/public_html/plotz.svg";
    let axis_tics = 5;
    let point_size = 5;
    let margin = 25;
    let xsize = 800;
    let ysize = 600;
    let y_label_area_size = 60;
    let point_color = RED;
    let y_precision = 1;
    let xlabel = "wt_koff";
    let ylabel = "CD27_ab";

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

    // Determine the precision for the x axis tics.

    let mut x_precision = 0;
    let m = xhigh.abs();
    if m > 0.0 {
        let extra = 2;
        x_precision = -m.log10() as usize + extra;
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

// enclone
// SOURCE=~/repos/sendai_private/dataset_commands/1a_base NOPRINT POUT=stdout PCOLS=wt_koff,CD27_ab

// DATA
// 0.00001,221
// 0.00001,186
// 0.000007,247
// 0.00001,390
// 0.00001,49
// 1.07042E-05,272
// 0.00001,116
// 0.00001,423
// 0.000472314,143
// 0.000261191,163
// 0.00001,130
// 0.00001,128
// 0.00106997,216
// 0.00001,171
// 0.00001,99
// 0.00001,249
// 4.54121E-05,78
// 0.00001,111
// 0.00001,63
// 0.000210643,17
// 6.08361E-05,202
// 0.00001,133
// 0.000254493,285
// 0.000325303,352
// 1.51725E-05,45
// 1.51725E-05,24
// 0.00001,187
// 0.00001,59
// 3.82693E-05,373
// 0.000228595,96
// 3.58783E-05,248
// 1.36806E-05,70
// 3.58783E-05,157
// 0.00001,6
// 3.17014E-05,453
// 0.000556859,279
// 0.000713161,86
// 0.000341469,111
// 0.001579333,7
// 0.000929721,33
// 2.57335E-05,132
// 0.000755029,12
// 0.00001,225
// 1.61381E-05,185
// 0.000618834,117
// 0.001542719,124
// 0.000417928,146
// 0.00019921,298
// 0.000386902,371
// 0.00001,31
// 0.000282927,32
// 0.00001,365
// 0.000238067,4
// 0.000431399,118
// 0.00001,17
// 0.00001,73
// 0.00078585,5
// 2.4144E-05,348
// 1.45446E-05,157
// 0.00001,184
// 3.50504E-05,225
// 0.000184403,29
// 0.00001,29
