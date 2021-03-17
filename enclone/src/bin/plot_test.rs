// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Plot experiment.

use enclone::plot_points::*;
use pretty_trace::*;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();

    // Names of stuff.

    let svg_filename = "/mnt/home/david.jaffe/public_html/plotz.svg";
    let xvar = "wt_koff";
    let yvar = "CD27_ab";

    // Load plot data, which are at the end of this file.

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

    // Make the plot.

    plot_points(&points, &xvar, &yvar, &svg_filename);
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
