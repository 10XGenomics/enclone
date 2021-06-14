// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use pretty_trace::*;
use usvg::SystemFontDB;

fn main() {
    PrettyTrace::new().on();
    let svg = std::fs::read("/Users/david.jaffe/plot.svg").unwrap();
    let png = convert_svg_to_png(&svg);
    std::fs::write("/Users/david.jaffe/plot.png", &png).unwrap();
}

pub fn convert_svg_to_png(svg: &[u8]) -> Vec<u8> {
    let fontdb = load_fonts();
    let current_dir = std::env::current_dir().unwrap();
    let usvg = usvg::Options {
        resources_dir: Some(current_dir),
        dpi: 96.0,
        font_family: "Times New Roman".to_string(),
        font_size: 12.0, // don't think this is used in our applications
        languages: vec!["en".to_string()],
        shape_rendering: usvg::ShapeRendering::default(),
        text_rendering: usvg::TextRendering::default(),
        image_rendering: usvg::ImageRendering::default(),
        keep_named_groups: false,
        fontdb,
    };
    let tree = usvg::Tree::from_data(&svg, &usvg).unwrap();
    let fit_to = usvg::FitTo::Original;
    let size = fit_to
        .fit_to(tree.svg_node().size.to_screen_size())
        .unwrap();
    let mut pixmap = tiny_skia::Pixmap::new(size.width(), size.height()).unwrap();
    pixmap.fill(tiny_skia::Color::from_rgba8(255, 255, 255, 255));
    resvg::render(&tree, fit_to, pixmap.as_mut());
    pixmap.encode_png().unwrap()
}

fn load_fonts() -> usvg::fontdb::Database {
    let mut fontdb = usvg::fontdb::Database::new();
    fontdb.load_system_fonts();
    fontdb.set_generic_families();
    fontdb.set_serif_family("Times New Roman");
    fontdb.set_sans_serif_family("Arial");
    fontdb.set_cursive_family("Comic Sans MS");
    fontdb.set_fantasy_family("Impact");
    fontdb.set_monospace_family("Courier New");
    fontdb
}
