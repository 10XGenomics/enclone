// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crc::*;
use string_utils::*;

// Modify a given PNG file by changing the pixels per meter to the given value.  This adds or
// replaces a preexisting pHYs chunk in the file.  Note that presence of an iDOT chunk might
// prevent the pHYs change from having an effect.  The iDOT chunk is an undocumented apple-ism.
// Googling reveals partial descriptions of it obtained by reverse engineering.
//
// The only tested value for pixels per meter is 5669, which corresponds to 144 DPI.  Using this
// value on a Mac causes the image to appear at half the size and hence at higher resolution,
// closer to the true Mac resolution.  This was tested on one Mac and might not be fully
// generalizable.  In fact, see last reference, 5669 is a special value.
//
// References:
// 1. https://en.wikipedia.org/wiki/Portable_Network_Graphics
// 2. https://www.w3.org/TR/2003/REC-PNG-20031110/#11pHYs
// 3. https://stackoverflow.com/questions/33894790/what-is-the-idot-chunk
// 4. https://www.hackerfactor.com/blog/index.php?/archives/895-Connecting-the-iDOTs.html
// 5. https://forums.ldraw.org/thread-23525-post-33020.html

pub fn set_pixels_per_meter(png: &mut Vec<u8>, pixels_per_meter: u32) {
    // Form the pHYs chunk.

    let mut bytes = Vec::<u8>::new();
    {
        let mut data = Vec::<u8>::new();
        data.append(&mut pixels_per_meter.to_be_bytes().to_vec());
        data.append(&mut pixels_per_meter.to_be_bytes().to_vec());
        data.append(&mut vec![1 as u8]);
        let len = 9 as u32;
        bytes.append(&mut len.to_be_bytes().to_vec());
        bytes.append(&mut b"pHYs".to_vec());
        bytes.append(&mut data);
        let crc = Crc::<u32>::new(&CRC_32_ISO_HDLC);
        let cs = crc.checksum(&bytes[4..bytes.len()]);
        bytes.append(&mut cs.to_be_bytes().to_vec());
    }

    // Copy the png, inserting a new pHYs chunk, and deleting the existing one.  We put the
    // new one right after the IHDR chunk.  According to the spec, it needs to go after IHDR and
    // before the first IDAT.

    let mut png2 = png[0..8].to_vec();
    let mut first = true;
    let mut pos = 8;
    while pos < png.len() {
        let mut x = png[pos..pos + 4].to_vec();
        let len = u32::from_be_bytes([x[0], x[1], x[2], x[3]]);
        let total = len as usize + 12;
        x = png[pos + 4..pos + 8].to_vec();
        if x != b"pHYs" {
            x = png[pos..pos + total].to_vec();
            png2.append(&mut x);
        }
        if first {
            png2.append(&mut bytes);
            first = false;
        }
        pos += total;
    }
    *png = png2;
}

// The following code is pretty much copied from the resvg crate,
// rev = 6b29007311edc5022635362fe56f6e5c0318fdeb, done June 14, 2021.

pub fn convert_svg_to_png(svg: &[u8], width: u32) -> Vec<u8> {
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
    svg = svg.replace("arial", "Liberation Sans");
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

    // Proceed.

    let fit_to = usvg::FitTo::Width(width);
    let size = fit_to
        .fit_to(tree.svg_node().size.to_screen_size())
        .unwrap();
    let mut pixmap = tiny_skia::Pixmap::new(size.width(), size.height()).unwrap();
    pixmap.fill(tiny_skia::Color::from_rgba8(255, 255, 255, 255));
    resvg::render(&tree, fit_to, pixmap.as_mut());
    let mut png = pixmap.encode_png().unwrap();
    const PIXELS_PER_METER: u32 = 5669;
    set_pixels_per_meter(&mut png, PIXELS_PER_METER);
    png
}

fn load_fonts() -> fontdb::Database {
    let mut fontdb = fontdb::Database::new();
    let deja = include_bytes!("../../fonts/DejaVuLGCSansMono.ttf").to_vec();
    fontdb.load_font_data(deja);
    let liberation_sans = include_bytes!("../../fonts/LiberationSans-Regular.ttf").to_vec();
    fontdb.load_font_data(liberation_sans);
    fontdb
}
