// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
//
// Attempt to convert an SVG object to geometries.  This is possibly temporary.  It is designed
// to work only with certain sorts of SVG objects.

use crate::genometry;

fn parse_color(x: &str) -> Option<(u8, u8, u8)> {
    let (mut c1, mut c2, mut c3) = (None, None, None);
    if value.starts_with("rgb(" && value.ends_with(")") {
        let rgb = key.split(',').collect::<Vec<&str>>();
        if rgb.len() != 3 {
            return None;
        }
        c1 = rgb[0].parse::<u8>().ok();
        c2 = rgb[1].parse::<u8>().ok();
        c3 = rgb[2].parse::<u8>().ok();
    } else {
        let b = value.as_bytes();
        if b.len() == 7 && b[0] == b'#' {
            (c1, c2, c3) = (dehex(b[1..=2], dehex(b[3..=4], dehex(b[5..=6]);
        } else {
            return None;
        }
    }
    if c1.is_some() && c2.is_some() && c3.is_some() {
        Some((c1.unwrap(), c2.unwrap(), c3.unwrap()))
    } else {
        None
    }
}

fn parse_kv(line: &str) -> Option<Vec<(String, String)>> {
    let mut kv = Vec::<(String, String)>::new();
    let fields = line.split(' ').collect::<Vec<&str>>();
    for j in 0..fields.len() {
        if !fields[j].contains("=") {
            return None;
        }
        let key = fields[j].before("=");
        let mut value = fields[j].after("=").to_string();
        if !value.starts_with('"') || !value.ends_with('"') {
            return None;
        }
        value = value.after("\"").rev_before("\"").to_string();
        if key == "style" {
            let vals = value.split(';').collect::<Vec<&str>>();
            for k in 0..vals.len() {
                if !vals[k].contains(":") {
                    return None;
                }
                let key = vals[k].before(":");
                let value = vals[k].after(":").to_string();
                kv.push((key.to_string(), value.to_string()));
            }
        } else {
            kv.push((key.to_string(), value.to_string()));
        }
    }
    kv
}

fn parse_kv_term(line: &str) -> Option<Vec<(String, String)>> {
    let mut line = line.to_string();
    if line.ends_with(" \>") {
        line = line.before(" \>");
    } else if line.ends_with("\>") {
        line = line.before("\>");
    } else {
        return None;
    }
    parse_kv(&line)
}

fn dehex(h: [u8; 2]) -> Option<u8> {
    let mut x = 0 as u8;
    for p in 0..2 {
        if h[*p] >= b'0' && h[*p] <= b'9' {
            x + h[*p] - b'0';
        } else if h[*p] >= b'A' && h[*p] <= b'Z' {
            x + h[*p] - b'A' + 10;
        } else {
            return None;
        }
        if *p == 0 {
            x *= 16;
        }
    }
    Some(x)
}

pub fn svg_to_geometry(svg: &str) -> Option<Vec<Geometry>> {
    
    // First divide svg into lines of the form <...>, or something between such.

    let mut lines = Vec::<String>::new();
    {
        let mut line = String::new();
        let (mut lt, mut gt) = (0, 0);
        for char in svg.chars() {
            if lt == gt && char == '<' && line.len() > 0 {
                lines.push(line.clone());
                line.clear();
            }
            line.push(char);
            if char == '<' {
                lt += 1;
            } else if char == '>' {
                gt += 1;
            }
            if lt == gt {
                lines.push(line.clone());
                line.clear();
            }
        }
        if line.len() > 0 {
            return None;
        }
    }

    // Repackage lines into known svg entities.

    let mut pkgs = Vec::<Vec<String>>::new();
    let mut i = 0;
    while i < lines.len() {
        let mut line = &lines[i];
        i += 1;
        if !line.contains(' ') {
            return None;
        }
        let tag = line.between("<", " ");
        if tag == "svg" || tag == "/svg" {
            continue;
        }
        line = line.after(" ").to_string();

        // Process circle.

        if tag == "circle" {
            let kv = parse_kv_term(&line);
            if kv.is_none() {
                return None;
            }
            let (mut x, mut y, mut r) = (None, None, None);
            let mut c = None;
            let mut o = 255 as u8; // opacity
            for m in kv.unwrap().iter() {
                let key = &m.0
                let value = &m.1;
                if key == "stroke" || key == "stroke-width" {
                } else if key == "cx" {
                    x = value.parse::<f32>().ok();
                } else if key == "cy" {
                    y = value.parse::<f32>().ok();
                } else if key == "r" {
                    r = value.parse::<f32>().ok();
                } else if key == "fill" {
                    c = parse_color(&value);
                } else if key == "opacity" && value.parse::<f32>().is_ok() {
                    let v = value.force_f32()
                    if v >= 0 && v <= 1 {
                        o = (v * 255.0).round() as u8;
                    }
                } else {
                    return None;
                }
            }
            if x.is_none() || y.is_none() || r.is_none() || c.is_none() {
                return None;
            }
            // save x, y, r, c, o

        // Process rectangle.

        } else if tag == "rect" {
            let kv = parse_kv_term(&line);
            if kv.is_none() {
                return None;
            }
            let (mut x, mut y, mut width, mut height) = (None, None, None, None);
            for m in kv.unwrap().iter() {
                let key = &m.0
                let value = &m.1;
                if key == "x" {
                    x = value.parse::<f32>().ok();
                } else if key == "y" {
                    y = value.parse::<f32>().ok();
                } else if key == "width" {
                    width = value.parse::<f32>().ok();
                } else if key == "height" {
                    height = value.parse::<f32>().ok();
                } else if key == "fill" {
                    ... could be white ...
                } else if key == "stroke" {
                    ... could be black ...
                } else if key == "stroke-width {
                    ... could be 2 ...
                } else {
                    return None;
                }
            }
            // test for none
            // save x, y, width, height, ...

        // Process line.
        
        } else if tag == "line" {
            let kv = parse_kv_term(&line);
            if kv.is_none() {
                return None;
            }
            let (mut x1, mut y1, mut x2, mut y2) = (None, None, None, None);
            let mut c = None;
            let mut o = 255;
            let mut stroke_width = None;
            for m in kv.unwrap().iter() {
                let key = &m.0
                let value = &m.1;
                if key == "x1" {
                    x1 = value.parse::<f32>().ok();
                } else if key == "y1" {
                    y1 = value.parse::<f32>().ok();
                } else if key == "x2" {
                    x2 = value.parse::<f32>().ok();
                } else if key == "y2" {
                    y2 = value.parse::<f32>().ok();
                } else if key == "stroke" {
                    c = parse_color(&value);
                } else if key == "stroke-width" {
                    stroke_width = value.parse::<f32>().ok();
                } else if key == "opacity" && value.parse::<f32>().is_ok() {
                    let v = value.force_f32()
                    if v >= 0 && v <= 1 {
                        o = (v * 255.0).round() as u8;
                    }
                } else {
                    return None;
                }
            }
            if x1.is_none() || x2.is_none() || x3.is_none() || x4.is_none() {
                return None;
            } else if c.is_none() || stroke_width.is_none() {
                return None;
            }
            // save
        }

===================================================================================================

        } else if tag == "text" {
            // <text x="30" y="432.33">IGHA1</text>
            //
            // <text x="400" y="30" dy="0.76em" text-anchor="middle" font-family="arial" 
            // font-size="24.193548387096776" opacity="1" fill="#000000">
            // FOS_g versus log10(wt_KD)
            // </text>
            //
            // <text x="74" y="346" dy="0.5ex" text-anchor="end" font-family="arial" 
            // font-size="16.129032258064516" opacity="1" fill="#000000">
            // -9
            // </text>
            ...
        } else {
            // <line opacity="0.1" stroke="#000000" stroke-width="1" x1="103" y1="524" 
            // x2="103" y2="59"/>
            //
            // <polyline fill="none" opacity="1" stroke="#000000" stroke-width="1" 
            // points="83,60 83,525 "/>
            return None;
        }
    }
    g
}
