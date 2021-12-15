// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Plot a given set of points, which come with sizes and colors.  This provides a function
// plot_points, built on top of the crate
// plotters 0.3.0.  That crate provides some useful capabilities, but does not include the
// full fledged ability to plot points.  The devil is in the details, and specifically in how
// exactly one sets the precision of tick mark labels, and how one positions axis labels so as
// to not overlap the tick mark labels and be the right distance from them.
//
// At some point we may wish to switch to using a different plotting crate or build our own.
// The string_width and ticks crates we have now might be part of this.
//
// If symmetric = true, produce a square plot having the same range and tic marks on both axes.
//
// points = {(point size, point color, x, y)}

use crate::string_width::arial_width;
use crate::ticks::ticks;
use plotters::prelude::*;
use std::cmp::max;

pub fn plot_points(
    points: &Vec<(u32, (u8, u8, u8), f32, f32)>,
    xvar: &str,
    yvar: &str,
    svg: &mut String,
    symmetric: bool,
    // title may be specified:
    title: Option<String>,
    // plot boundaries may be specified:
    xlow: Option<f32>,
    xhigh: Option<f32>,
    ylow: Option<f32>,
    yhigh: Option<f32>,
    // optional margin:
    // It is a bug that this has to be passed sometimes.
    margin: Option<u32>,
) -> Result<(), String> {
    // Requirements.

    if points.is_empty() {
        return Err(format!(
            "\nPlot of {} versus {} can't be carried out because there are no data \
            points.\n",
            xvar, yvar
        ));
    }

    // Define parameters of the plot.

    let mut titlex = format!("{} versus {}", xvar, yvar);
    if title.is_some() {
        titlex = title.unwrap().clone();
    }
    let mut marginx = 25;
    if margin.is_some() {
        marginx = margin.unwrap();
    }

    // Possibly universal constants.

    let titlex_font_size = 30;
    let font = "arial";
    let tic_font_size = 20;
    let axis_ticks = 5;
    let mut xsize = 800;
    let ysize = 600;
    if symmetric {
        xsize = ysize;
    }
    let range_ext = 0.02;

    // Determine the plot ranges using the extreme values of the points, extended a little bit.

    let (mut xlowx, mut xhighx) = (points[0].2, points[0].2);
    let (mut ylowx, mut yhighx) = (points[0].3, points[0].3);
    for i in 1..points.len() {
        xlowx = xlowx.min(points[i].2);
        xhighx = xhighx.max(points[i].2);
        ylowx = ylowx.min(points[i].3);
        yhighx = yhighx.max(points[i].3);
    }
    if xlowx > 0.0 && xlowx / (xhighx - xlowx) < range_ext {
        xlowx = 0.0;
    } else if xlowx > 0.0 {
        xlowx *= 1.0 - range_ext;
    } else {
        xlowx *= 1.0 + range_ext;
    }
    if xhighx > 0.0 {
        xhighx *= 1.0 + range_ext;
    } else {
        xhighx *= 1.0 - range_ext;
    }
    if ylowx > 0.0 && ylowx / (yhighx - ylowx) < range_ext {
        ylowx = 0.0;
    } else if ylowx > 0.0 {
        ylowx *= 1.0 - range_ext;
    } else if ylowx < 0.0 {
        ylowx *= 1.0 + range_ext;
    }
    if yhighx > 0.0 {
        yhighx *= 1.0 + range_ext;
    } else {
        yhighx *= 1.0 - range_ext;
    }
    if symmetric {
        xlowx = xlowx.min(ylowx);
        ylowx = xlowx;
        xhighx = xhighx.max(yhighx);
        yhighx = xhighx;
    }
    if xlow.is_some() {
        xlowx = xlow.unwrap();
    }
    if xhigh.is_some() {
        xhighx = xhigh.unwrap();
    }
    if ylow.is_some() {
        ylowx = ylow.unwrap();
    }
    if xhigh.is_some() {
        yhighx = yhigh.unwrap();
    }

    // Get the tick mark labels.  Note that these are not the actual labels used by plotters.
    // (We could not figure out exactly how plotters determines the labels.)
    // Rather, these are our best approximation to them, and these enable us to set the precisions
    // of the tick lables and to appropriately position the  axis labels relative to them.

    let x_ticks = ticks(xlowx, xhighx, axis_ticks, false);
    let y_ticks = ticks(ylowx, yhighx, axis_ticks, false);
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
    let root = root.margin(marginx, marginx, marginx, marginx);
    let mut chart = ChartBuilder::on(&root)
        .caption(&titlex, (font, titlex_font_size).into_font())
        .x_label_area_size(x_label_area_size)
        .y_label_area_size(y_label_area_size)
        .build_cartesian_2d(xlowx..xhighx, ylowx..yhighx)
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

    let mut i = 0;
    while i < points.len() {
        let mut j = i + 1;
        while j < points.len() {
            if points[j].0 != points[i].0 || points[j].1 != points[i].1 {
                break;
            }
            j += 1;
        }
        let mut pts = Vec::<(f32, f32)>::new();
        for k in i..j {
            pts.push((points[k].2, points[k].3));
        }
        chart
            .draw_series(PointSeries::of_element(
                pts,
                points[i].0,
                &RGBColor(points[i].1 .0, points[i].1 .1, points[i].1 .2),
                &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
            ))
            .unwrap();
        i = j;
    }
    Ok(())
}
