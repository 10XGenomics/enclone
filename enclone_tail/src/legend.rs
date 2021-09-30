// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::assign_cell_color::*;
use crate::colors::*;
use crate::string_width::*;
use crate::ticks::*;
use crate::*;
use enclone_core::cell_color::*;
use enclone_core::defs::*;
use string_utils::*;

pub fn add_legend_for_color_by_variable(
    plot_opt: &PlotOpt,
    svg: &mut String,
    color: &Vec<String>,
    actual_width: f64,
    actual_height: f64,
) {
    let mut var = String::new();
    let mut display_var = String::new();
    let mut xmin = None;
    let mut xmax = None;
    match plot_opt.cell_color {
        CellColor::ByVariableValue(ref x) => {
            var = x.var.clone();
            display_var = x.display_var.clone();
            xmin = x.min;
            xmax = x.max;
        }
        _ => {}
    };
    let mut defined = false;
    let mut have_undefined = false;
    for i in 0..color.len() {
        if color[i] == "undefined" {
            have_undefined = true;
        }
    }

    // Get the actual low and high values for the variable.

    let (mut low, mut high) = (0.0, 0.0);
    let n = VAR_LOW.lock().unwrap().len();
    for i in 0..n {
        if VAR_LOW.lock().unwrap()[i].0 == var {
            low = VAR_LOW.lock().unwrap()[i].1;
            high = VAR_HIGH.lock().unwrap()[i].1;
            defined = true;
        }
    }

    // Print the variable name.

    *svg = svg.rev_before("<").to_string();
    let font_size = 20;
    let name_bar_height = font_size as f64 + font_size as f64 / 2.0;
    let legend_xstart = actual_width + 20.0;
    let legend_ystart = BOUNDARY as f64 + name_bar_height;
    let band_width = 100.0;
    *svg += &format!(
        "<text text-anchor=\"start\" x=\"{}\" y=\"{}\" font-family=\"Arial\" \
         font-size=\"{}\">{}</text>\n",
        legend_xstart,
        BOUNDARY as f64 + font_size as f64 / 2.0,
        font_size,
        display_var,
    );

    // Handle the special case where all points are undefined.

    if !defined {
        let fail_text = "The variable is undefined for all points.";
        *svg += &format!(
            "<text text-anchor=\"start\" x=\"{:.2}\" y=\"{:.2}\" font-family=\"Arial\" \
             font-size=\"{}\">{}</text>\n",
            legend_xstart,
            BOUNDARY as f64 + 2.0 * font_size as f64,
            font_size,
            fail_text,
        );
        let mut max_text_width = arial_width(&display_var, font_size as f64);
        max_text_width = max_text_width.max(arial_width(fail_text, font_size as f64));
        let width = legend_xstart + max_text_width + BOUNDARY as f64;
        set_svg_width(svg, width);
        *svg += "</svg>";

    // Handle the case where there is at least one defined point.
    } else {
        let mut height_for_undefined = 0.0;
        if have_undefined {
            height_for_undefined = 20.0;
        }
        let available = actual_height - name_bar_height - height_for_undefined - BOUNDARY as f64;
        *svg += &format!(
            "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" \
             style=\"fill:white;stroke:black;stroke-width:1\" />\n",
            legend_xstart, legend_ystart, band_width, available,
        );

        // Make the color bar.  It would make sense to have 256 bars that abut exactly,
        // however rendering appears to be better if the bars overlap slightly.  Since when
        // there are overlapping rectangles, the later rectangle should dominate, this should
        // not affect the appearance at all.  But it does.

        let band_height = available / 256.0;
        for i in 0..256 {
            let ystart = legend_ystart + i as f64 * band_height;
            let c = &TURBO_SRGB_BYTES[i];
            let color = format!("rgb({},{},{})", c[0], c[1], c[2]);
            let mut add = 0.0;
            if i < 255 {
                add = band_height / 10.0;
            }
            *svg += &format!(
                "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" \
                 style=\"fill:{}\" />\n",
                legend_xstart,
                ystart,
                band_width,
                band_height + add,
                color,
            );
        }

        // Define the tick marks.

        const MAX_TICKS: usize = 5;
        let mut zlow = low;
        if xmin.is_some() {
            zlow = xmin.unwrap();
        }
        let mut zhigh = high;
        if xmax.is_some() {
            zhigh = xmax.unwrap();
        }
        let mut ticks = ticks(zlow as f32, zhigh as f32, MAX_TICKS, false);
        ticks.insert(0, format!("{}", zlow));
        ticks.push(format!("{}", zhigh));

        // Add the ticks.

        let mut max_text_width: f64 = 0.0;
        let sep_to_text = 10.0;
        let text_xstart = legend_xstart + band_width + sep_to_text;
        let mut text_ystarts = Vec::<f64>::new();
        for (i, text) in ticks.iter().enumerate() {
            // Define vertical shift for value text.  We vertically center the text at the
            // correct point, adding font_size/4 to get this to happen.  We don't understand
            // why four makes sense.  Also, we treat the first and last labels differently,
            // because it is aesthetically displeasing to have the text outside the boundaries
            // of the color box.

            let vshift;
            if i == 0 {
                vshift = font_size as f64 / 2.0 + 1.0;
            } else if i == ticks.len() - 1 {
                vshift = 0.0;
            } else {
                vshift = font_size as f64 / 4.0;
            }
            let ystart = legend_ystart + available * (text.force_f64() - zlow) / (zhigh - zlow);
            let text_ystart = ystart + vshift;
            text_ystarts.push(text_ystart);
        }
        for (i, text) in ticks.iter().enumerate() {
            // Generate the text.

            let ystart = legend_ystart + available * (text.force_f64() - zlow) / (zhigh - zlow);
            let text_ystart = text_ystarts[i];
            if i == 1 && text_ystart - text_ystarts[0] < font_size as f64 {
                continue;
            }
            if ticks.len() >= 2
                && i == ticks.len() - 2
                && text_ystarts[ticks.len() - 1] - text_ystart < font_size as f64
            {
                continue;
            }
            let mut textp = text.clone();
            if i == 0 && zlow > low {
                textp = format!("≤ {}", text);
            }
            if i == ticks.len() - 1 && zhigh < high {
                textp = format!("≥ {}", text);
            }
            *svg += &format!(
                "<text text-anchor=\"start\" x=\"{:.2}\" y=\"{:.2}\" font-family=\"Arial\" \
                 font-size=\"{}\">{}</text>\n",
                text_xstart, text_ystart, font_size, textp,
            );
            max_text_width = max_text_width.max(arial_width(&textp, font_size as f64));

            // Add tick lines.

            if i > 0 && i < ticks.len() - 1 {
                *svg += &format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
                     stroke=\"#000000\" stroke-width=\"0.5\"/>\n",
                    legend_xstart,
                    ystart,
                    legend_xstart + band_width,
                    ystart,
                );
            }
        }

        // Add legend for undefined points.

        if have_undefined {
            let r = 5.0;
            let vsep = 10.0;
            let y = legend_ystart + available + vsep + r;
            *svg += &format!(
                "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" stroke=\"red\" \
                 stroke-width=\"0.5\" fill=\"white\" />\n",
                legend_xstart + band_width - r,
                legend_ystart + available + vsep + r,
                r,
            );
            *svg += &format!(
                "<text text-anchor=\"start\" x=\"{:.2}\" y=\"{:.2}\" font-family=\"Arial\" \
                 font-size=\"{}\">{}</text>\n",
                text_xstart,
                y + font_size as f64 / 4.0 - 1.0,
                font_size,
                "undefined",
            );
            max_text_width = max_text_width.max(arial_width("undefined", font_size as f64));
        }

        // Finish.

        let mut legend_width = band_width + sep_to_text + max_text_width;
        legend_width = legend_width.max(arial_width(&display_var, font_size as f64));
        let mut width = legend_xstart + legend_width;
        width += BOUNDARY as f64;
        set_svg_width(svg, width + BOUNDARY as f64);
        set_svg_height(svg, actual_height + BOUNDARY as f64);
        *svg += "</svg>";
    }
}
