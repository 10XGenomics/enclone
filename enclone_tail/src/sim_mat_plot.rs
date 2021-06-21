// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Execute SIM_MAT_PLOT.

use enclone_core::defs::*;
use io_utils::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use string_utils::*;
use tables::*;

fn hex(c: u8) -> String {
    let x1 = c / 16;
    let x2 = c % 16;
    let c1 = if x1 < 10 { b'0' + x1 } else { b'A' + x1 - 10 };
    let c2 = if x2 < 10 { b'0' + x2 } else { b'A' + x2 - 10 };
    stringme(&[c1, c2])
}

fn hex_color(r: u8, g: u8, b: u8) -> String {
    format!("#{}{}{}", hex(r), hex(g), hex(b))
}

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
                for k in 0..n {
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
        let ncells = mat.len();
        let mut lens = vec![0.0; n];
        for i in 0..n {
            for j in 0..ncells {
                lens[i] += mat[j][i] * mat[j][i];
            }
            lens[i] = lens[i].sqrt();
        }
        let mut means = Vec::<f64>::new();
        for i in 0..n {
            let mut sum = 0.0;
            for j in 0..ncells {
                sum += mat[j][i];
            }
            means.push(sum / ncells as f64);
        }
        let mut cos = vec![vec![0.0; n]; n];
        for i1 in 0..n {
            for i2 in 0..n {
                let mut dot = 0.0;
                for j in 0..ncells {
                    dot += mat[j][i1] * mat[j][i2];
                }
                if lens[i1] == 0.0 || lens[i2] == 0.0 {
                    cos[i1][i2] = 0.0;
                } else {
                    cos[i1][i2] = dot / (lens[i1] * lens[i2]);
                }
            }
        }
        let dim = 500;
        let (width, height) = (dim, dim);
        let font_size = 130.0 / n as f64;
        let dimn = dim as f64 / n as f64;

        // Define the row text matrix.

        let mut rtm = Vec::<Vec<String>>::new();
        rtm.push(vec![
            "variable".to_string(),
            "mean".to_string(),
            "#".to_string(),
        ]);
        for i in 0..n {
            rtm.push(vec![
                vars[i].clone(),
                format!("{:.1}", means[i]),
                format!("{}", i + 1),
            ]);
        }
        let mut log = Vec::<u8>::new();
        print_tabular(&mut log, &rtm, 2, Some(b"lrl".to_vec()));
        let mut slong = stringme(&log);
        slong = slong.replace(" ", "\u{00A0}"); // convert spaces to non-breaking spaces
        let mut lines = Vec::<String>::new();
        for line in slong.lines() {
            lines.push(line.to_string());
        }
        const DEJA_SANS_MONO_WIDTH_HEIGHT_RATIO: f64 = 0.42; // guess

        // Define row titles.

        let mut max_title_width = 0.0 as f64;
        for i in 0..n {
            max_title_width = max_title_width
                .max(lines[i].len() as f64 * font_size * DEJA_SANS_MONO_WIDTH_HEIGHT_RATIO);
        }
        let sep = 10.0;
        let x0 = max_title_width + sep * 2.0;

        // Start making SVG.

        let mut svg = format!(
            "<svg version=\"1.1\"\n\
             baseProfile=\"full\"\n\
             width=\"{}\" height=\"{}\"\n\
             xmlns=\"http://www.w3.org/2000/svg\">\n",
            x0 + width as f64 + sep,
            sep + height as f64 + sep * 2.0 + font_size
        );

        // Load font if not GUI.

        if ctl.plot_opt.sim_mat_plot_file != "gui" {
            svg += "\
              <defs>\n\
                <style type=\"text/css\">\n\
                  @font-face {\n\
                    font-family: \"DejaVu LGC Sans Mono\";\n\
                    src: url('https://cdn.jsdelivr.net/npm/@deathbeds/\
                        jupyterlab-font-dejavu-sans-mono@1.0.0/style/fonts/DejaVuSansMono.woff2')\n\
                    format('woff2'),\n\
                    url('https://cdn.jsdelivr.net/npm/dejavu-fonts-ttf@2.37.3/ttf/\
                        DejaVuSansMono.ttf')\n\
                    format('truetype');\n\
                  }\n\
                </style>\n\
              </defs>\n";
        }

        // Add row titles.

        svg += &mut format!(
            "<text x=\"{}\" y=\"{}\" font-family=\"DejaVu LGC Sans Mono\" \
            font-size=\"{}\" text-anchor=\"start\" fill=\"black\">{}</text>\n",
            sep, font_size, font_size, lines[0],
        );
        for i in 0..n {
            let y = sep + (i as f64) * dimn;
            svg += &mut format!(
                "<text x=\"{}\" y=\"{}\" font-family=\"DejaVu LGC Sans Mono\" \
                font-size=\"{}\" text-anchor=\"start\" fill=\"black\">{}</text>\n",
                sep,
                y + dimn / 2.0 + font_size / 2.0,
                font_size,
                lines[i + 1],
            );
        }

        // Print the variable numbers at the bottom.

        for i in 0..n {
            let x = x0 + (i as f64) * dimn;
            svg += &mut format!(
                "<text x=\"{}\" y=\"{}\" font-family=\"DejaVu LGC Sans Mono\" \
                font-size=\"{}\" text-anchor=\"middle\" fill=\"black\">{}</text>\n",
                x + dimn / 2.0,
                dim as f64 + sep * 2.0 + font_size,
                font_size,
                format!("{}", i + 1),
            );
        }

        // Print the matrix.

        for i1 in 0..n {
            for i2 in 0..n {
                let x = x0 + (i1 as f64) * dimn;
                let y = sep + (i2 as f64) * dimn;
                let gray = (255 as f64 * (1.0 - cos[i1][i2])).round() as u8;
                svg += &mut format!(
                    "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" \
                    style=\"fill:{};stroke:black;stroke-width:1\" />\n",
                    x,
                    y,
                    dimn,
                    dimn,
                    hex_color(gray, gray, gray),
                );
                svg += &mut format!(
                    "<text x=\"{}\" y=\"{}\" font-family=\"DejaVu LGC Sans Mono\" \
                    font-size=\"{}\" text-anchor=\"middle\" fill=\"red\">{}</text>\n",
                    x + dimn / 2.0,
                    y + dimn / 2.0 + font_size / 2.0,
                    font_size,
                    format!("{:.2}", cos[i1][i2]),
                );
            }
        }

        // Finish.

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
