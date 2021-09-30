// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod align_n;
pub mod alluvial_fb;
pub mod assign_cell_color;
pub mod circles_to_svg;
pub mod clustal;
pub mod colors;
pub mod display_tree;
pub mod fasta;
pub mod group;
pub mod group_colors;
pub mod grouper;
pub mod hex;
pub mod legend;
pub mod neighbor;
pub mod newick;
pub mod pack_circles;
pub mod phylip;
pub mod plot;
pub mod plot_points;
pub mod plot_utils;
pub mod polygon;
pub mod print_stats;
pub mod requirements;
pub mod sim_mat_plot;
pub mod string_width;
pub mod tail;
pub mod ticks;
pub mod tree;

use string_utils::*;

const BOUNDARY: usize = 10;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Change the width or height of an svg document.

pub fn set_svg_width(svg: &mut String, new_width: f64) {
    let (svg1, svg2) = (svg.before("width="), svg.after("width=\"").after("\""));
    *svg = format!("{}width=\"{}\"{}", svg1, new_width, svg2);
}

fn set_svg_height(svg: &mut String, new_height: f64) {
    let (svg1, svg2) = (svg.before("height="), svg.after("height=\"").after("\""));
    *svg = format!("{}height=\"{}\"{}", svg1, new_height, svg2);
}

fn get_svg_height(svg: &String) -> f64 {
    svg.between("height=\"", "\"").force_f64()
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Here, and in "enclone help color", we swap the order of colors, placing the last three before
// the first three.  This is because the last three seem to make a better three-color palette.

pub fn substitute_enclone_color(color: &mut String) {
    if *color == "@1".to_string() {
        *color = "rgb(0,95,175)".to_string();
    } else if *color == "@2".to_string() {
        *color = "rgb(215,135,175)".to_string();
    } else if *color == "@3".to_string() {
        *color = "rgb(0,175,135)".to_string();
    } else if *color == "@4".to_string() {
        *color = "rgb(215,95,0)".to_string();
    } else if *color == "@5".to_string() {
        *color = "rgb(95,175,255)".to_string();
    } else if *color == "@6".to_string() {
        *color = "rgb(215,175,0)".to_string();
    }
}
