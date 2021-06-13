// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Execute SIM_MAT_PLOT.

use enclone_core::defs::*;
use io_utils::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use string_utils::*;

pub fn sim_mat_plot(
    ctl: &EncloneControl,
    groups: &Vec<Vec<(i32, String)>>,
    out_datas: &Vec<Vec<HashMap<String, String>>>,
    svgs: &mut Vec<String>,
) {

    if ctl.plot_opt.sim_mat_plot_file.len() > 0 {
        let filename = &ctl.plot_opt.sim_mat_plot_file;
        let vars = &ctl.plot_opt.sim_mat_plot_vars;
        let n = vars.len();
        let mut mat = Vec::<Vec<f64>>::new();
        for i in 0..groups.len() {
            let mut o = Vec::<i32>::new();
            for j in 0..groups[i].len() {
                o.push(groups[i][j].0);
            }
            for j in 0..o.len() {
                let oo = o[j] as usize;
                let mut v = Vec::<f64>::new();
                for k in 0..vars.len() {
                    let mut x = 0.0;
                    for i in 0..out_datas[oo].len() {
                        let p = &out_datas[oo][i];
                        if p.contains_key(&vars[k].clone()) {
                            let z = &p[&vars[k].clone()];
                            if z.parse::<f64>().is_ok() {
                                x = z.force_f64();
                            }
                        }
                    }
                    v.push(x);
                }
                mat.push(v);
            }
        }
        let mut lens = vec![0.0; n];
        for i in 0..n {
            for j in 0..mat.len() {
                lens[i] += mat[j][i] * mat[j][i];
            }
            lens[i] = lens[i].sqrt();
        }
        let mut cos = vec![vec![0.0; n]; n];
        for i1 in 0..n {
            for i2 in 0..n {
                let mut dot = 0.0;
                for j in 0..mat.len() {
                    dot += mat[j][i1] * mat[j][i2];
                }
                if lens[i1] == 0.0 || lens[i2] == 0.0 {
                    cos[i1][i2] = 0.0;
                } else {
                    cos[i1][i2] = dot / (lens[i1] * lens[i2]);
                }
            }
        }
        let dim = 400;
        let (width, height) = (dim, dim);
        let mut svg = format!(
            "<svg version=\"1.1\"\n\
             baseProfile=\"full\"\n\
             width=\"{}\" height=\"{}\"\n\
             xmlns=\"http://www.w3.org/2000/svg\">\n",
            width, height
        );
        let font_size = 80.0 / n as f64;
        let dimn = dim as f64 / n as f64;
        for i1 in 0..n {
            for i2 in 0..n {
                let x = (i1 as f64) * (dim as f64 / n as f64);
                let y = (i2 as f64) * (dim as f64 / n as f64);
                let gray = 255 as f64 * (1.0 - cos[i1][i2]);
                svg += &mut format!("<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" \
                    style=\"fill:rgb({},{},{});stroke:black;stroke-width:1\" />\n",
                    x,
                    y,
                    dimn,
                    dimn,
                    gray,
                    gray,
                    gray,
                );
                svg += &mut format!("<text x=\"{}\" y=\"{}\" font-family=\"Arial\" \
                    font-size=\"{}\" text-anchor=\"middle\" fill=\"red\" \
                    dominant-baseline=\"middle\">{}</text>\n",
                    x + dimn/2.0,
                    y + dimn/2.0,
                    font_size,
                    format!("{:.2}", cos[i1][i2]),
                );
            }
        }
        svg += "</svg>";
        if filename == "stdout" || filename == "gui_stdout" {
            for line in svg.lines() {
                println!("{}", line);
            }
        } else if filename == "gui" {
            svgs.push(svg);
        } else {
            let mut f = open_for_write_new![&filename];
            fwrite!(f, "{}", svg);
        }
    }
}
