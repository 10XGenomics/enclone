// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Plot a given set of points.  This provides a function plot_points, built on top of the crate
// plotters 0.3.0.  That crate provides some useful capabilities, but does not include the
// full fledged ability to plot points.  The devil is in the details, and specifically in how
// exactly one sets the precision of tick mark labels, and how one positions axis labels so as
// to not overlap the tick mark labels and be the right distance from them.
//
// At some point we may wish to switch to using a different plotting crate or build our own.
// The string_width and ticks crates we have now might be part of this.
//

// If symmetric = true, produce a square plot having the same range and tic marks on both axes.

use enclone_tail::string_width::arial_width;
use enclone_tail::ticks::ticks;
use plotters::prelude::*;
use std::cmp::max;

pub fn plot_points(
    points: &Vec<(f32, f32)>,
    xvar: &str,
    yvar: &str,
    svg: &mut String,
    symmetric: bool,
) -> Result<(), String> {

    // Requirements.

    if points.is_empty() {
        return Err(format!(
            "\nPlot of {} versus {} can't be carried out because there are no data \
            points\n(consisting of pairs of numbers).\n",
            xvar, yvar
        ));
    }

    // Define parameters of the plot.

    // let title = format!("{} versus {}", xvar, yvar);
    let title = "clonotyping accuracy";

    // Possibly universal constants.

    let title_font_size = 30;
    let font = "arial";
    let tic_font_size = 20;
    let axis_ticks = 5;
    let point_size = 5;
    let margin = 35;
    let mut xsize = 800;
    let ysize = 600;
    if symmetric {
        xsize = ysize;
    }
    let point_color = RGBColor(255,153,51);

    // Determine the plot ranges using the extreme values of the points, extended a little bit.

    let (mut xlow, mut xhigh) = (points[0].0, points[0].0);
    let (mut ylow, mut yhigh) = (points[0].1, points[0].1);
    for i in 1..points.len() {
        xlow = xlow.min(points[i].0);
        xhigh = xhigh.max(points[i].0);
        ylow = ylow.min(points[i].1);
        yhigh = yhigh.max(points[i].1);
    }
    xlow = 0.0;
    xhigh = 1500000.0;
    ylow = 0.0;
    yhigh = 2.0;
    if symmetric {
        xlow = xlow.min(ylow);
        ylow = xlow;
        xhigh = xhigh.max(yhigh);
        yhigh = xhigh;
    }

    // Get the tick mark labels.  Note that these are not the actual labels used by plotters.
    // (We could not figure out exactly how plotters determines the labels.)
    // Rather, these are our best approximation to them, and these enable us to set the precisions
    // of the tick lables and to appropriately position the  axis labels relative to them.

    let x_ticks = ticks(xlow, xhigh, axis_ticks, false);
    let y_ticks = ticks(ylow, yhigh, axis_ticks, false);
    assert!(!x_ticks.is_empty());
    assert!(!y_ticks.is_empty());

    // Determine precision for axes ticks.

    let (mut x_precision, mut y_precision) = (0, 0);
    let x_dot = x_ticks[0].find('.');
    if x_dot.is_some() {
        x_precision = x_ticks[0].len() - x_dot.unwrap() - 1;
    }
    let y_dot = y_ticks[0].find('.');
    if y_dot.is_some() {
        y_precision = y_ticks[0].len() - y_dot.unwrap() - 1;
    }

    // Determine the area size for the x label.

    let x_label_area_size = (2.5 * tic_font_size as f32).round() as u32;

    // Determine the area size for the y label.

    let mut max_ytick_width = 0;
    for t in y_ticks.iter() {
        max_ytick_width = max(
            max_ytick_width,
            arial_width(&*t, tic_font_size as f64).ceil() as usize,
        );
    }
    let extra = (tic_font_size as f32 * 1.5) as usize;
    let y_label_area_size = (max_ytick_width + extra) as u32;

    // Make the plot.

    let root = SVGBackend::with_string(svg, (xsize, ysize)).into_drawing_area();
    let root = root.margin(margin, margin, margin, margin);
    let mut chart = ChartBuilder::on(&root)
        .caption(&title, (font, title_font_size).into_font())
        .x_label_area_size(x_label_area_size)
        .y_label_area_size(y_label_area_size)
        .build_cartesian_2d(xlow..xhigh, ylow..yhigh)
        .unwrap();
    chart
        .configure_mesh()
        .label_style((font, tic_font_size as u32).into_font())
        .x_labels(axis_ticks)
        .y_labels(axis_ticks)
        .x_label_formatter(&|x| format!("{:.1$}", x, x_precision))
        .y_label_formatter(&|x| format!("{:.1$}", x, y_precision))
        .x_desc(xvar)
        .y_desc(yvar)
        .draw()
        .unwrap();
    chart
        .draw_series(PointSeries::of_element(
            points.clone(),
            // point_size,
            2,
            &point_color,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
        ))
        .unwrap();

    let point_color2 = BLUE;
    let points2 = vec![
        (652537.0, 0.66), 
        (672653.0, 1.23),
        (405168.0, 0.75),
        (245629.0, 0.43),
        (301997.0, 0.68),
    ];
    chart
        .draw_series(PointSeries::of_element(
            points2.clone(),
            point_size,
            &point_color2,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
        ))
        .unwrap();

    let point_color3 = RED;
    let points3 = vec![(1325190.0, 0.90)];
    chart
        .draw_series(PointSeries::of_element(
            points3.clone(),
            point_size,
            &point_color3,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
        ))
        .unwrap();

    let point_color4 = RGBColor(255,153,51);
    let sr = [
        (0.1, 0.25),
        (0.2, 0.42),
        (0.3, 0.45),
        (0.4, 0.60),
        (0.5, 0.57),
        (0.6, 0.72),
        (0.7, 0.81),
        (0.8, 0.76), // preliminary
    ];
    let cells = 1325190.0 as f32;
    let mut points4 = Vec::<(f32, f32)>::new();
    for i in 0..sr.len() {
        points4.push((sr[i].0 as f32 * cells, sr[i].1));
    }
    chart
        .draw_series(PointSeries::of_element(
            points4.clone(),
            point_size,
            &point_color4,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
        ))
        .unwrap();

    Ok(())
}
