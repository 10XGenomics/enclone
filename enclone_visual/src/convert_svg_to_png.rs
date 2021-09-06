// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// This code is pretty much copied from the resvg crate,
// rev = 6b29007311edc5022635362fe56f6e5c0318fdeb, done June 14, 2021.

use string_utils::*;

pub fn convert_svg_to_png(svg: &[u8]) -> Vec<u8> {
    let fontdb = load_fonts();
    let usvg = usvg::OptionsRef {
        resources_dir: None,
        dpi: 96.0,
        default_size: usvg::Size::new(400.0, 400.0).unwrap(),
        font_family: "Arial",
        font_size: 12.0, // don't think this is used in our applications
        languages: &["en".to_string()],
        shape_rendering: usvg::ShapeRendering::default(),
        text_rendering: usvg::TextRendering::default(),
        image_rendering: usvg::ImageRendering::default(),
        keep_named_groups: false,
        fontdb: &fontdb,
    };
    let mut svg = stringme(&svg);
    svg = svg.replace("Arial", "Liberation Sans");
    let tree = usvg::Tree::from_data(&svg.as_bytes(), &usvg);
    if tree.is_err() {
        panic!(
            "svg conversion failed with message {} on\n{}\n",
            tree.err().unwrap(),
            svg
        );
    }
    let tree = tree.unwrap();
    let fit_to = usvg::FitTo::Original;
    let size = fit_to
        .fit_to(tree.svg_node().size.to_screen_size())
        .unwrap();
    let mut pixmap = tiny_skia::Pixmap::new(size.width(), size.height()).unwrap();
    pixmap.fill(tiny_skia::Color::from_rgba8(255, 255, 255, 255));
    resvg::render(&tree, fit_to, pixmap.as_mut());
    pixmap.encode_png().unwrap()
}

fn load_fonts() -> fontdb::Database {
    let mut fontdb = fontdb::Database::new();
    let deja = include_bytes!("../../fonts/DejaVuLGCSansMono.ttf").to_vec();
    fontdb.load_font_data(deja);
    let liberation_sans = include_bytes!("../../fonts/LiberationSans-Regular.ttf").to_vec();
    fontdb.load_font_data(liberation_sans);
    fontdb
}
