// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
//
// Attempt to convert an SVG object to geometries.  This is possibly temporary.  It is designed
// to work only with certain sorts of SVG objects.

use crate::geometry::*;
use crate::svg_to_geometry::HorizontalAlignment::*;
use string_utils::*;

fn get_opacity(key: &str, value: &str, o: &mut u8) -> bool {
    if key == "opacity" && value.parse::<f32>().is_ok() {
        let v = value.parse::<f32>().unwrap();
        if v >= 0.0 && v <= 1.0 {
            *o = (v * 255.0).round() as u8;
            return true;
        } else {
            return false;
        }
    }
    return false;
}

fn numeric(x: &str) -> Option<f32> {
    x.parse::<f32>().ok()
}

fn get_numeric(key: &str, value: &str, var: &str, x: &mut Option<f32>) -> bool {
    if key == var {
        *x = numeric(&value);
        x.is_some()
    } else {
        false
    }
}

fn dehex(h: &[u8]) -> Option<u8> {
    let mut x = 0 as u8;
    for p in 0..2 {
        if h[p] >= b'0' && h[p] <= b'9' {
            x += h[p] - b'0';
        } else if h[p] >= b'A' && h[p] <= b'Z' {
            x += h[p] - b'A' + 10;
        } else {
            return None;
        }
        if p == 0 {
            x *= 16;
        }
    }
    Some(x)
}

fn parse_color(x: &str) -> Option<(u8, u8, u8)> {
    let (c1, c2, c3);
    if x.starts_with("rgb(") && x.ends_with(")") {
        let rgb = x
            .after("rgb(")
            .rev_before(")")
            .split(',')
            .collect::<Vec<&str>>();
        if rgb.len() != 3 {
            return None;
        }
        c1 = rgb[0].parse::<u8>().ok();
        c2 = rgb[1].parse::<u8>().ok();
        c3 = rgb[2].parse::<u8>().ok();
    } else if x == "white" {
        c1 = Some(255);
        c2 = Some(255);
        c3 = Some(255);
    } else if x == "black" {
        c1 = Some(0);
        c2 = Some(0);
        c3 = Some(0);
    } else {
        let b = x.as_bytes();
        if b.len() == 7 && b[0] == b'#' {
            c1 = dehex(&b[1..=2]);
            c2 = dehex(&b[3..=4]);
            c3 = dehex(&b[5..=6]);
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
    Some(kv)
}

fn parse_kv_term(line: &str) -> Option<Vec<(String, String)>> {
    let mut line = line.to_string();
    if line.ends_with(" \\>") {
        line = line.before(" \\>").to_string();
    } else if line.ends_with("\\>") {
        line = line.before("\\>").to_string();
    } else {
        return None;
    }
    parse_kv(&line)
}

pub fn svg_to_geometry(svg: &str) -> Option<Vec<Thing>> {
    // First divide svg into lines of the form <...>, or something between such.

    let mut lines = Vec::<String>::new();
    {
        let mut line = String::new();
        let (mut lt, mut gt) = (0, 0);
        for char in svg.chars() {
            if lt == gt && char == '<' {
                if line == "\n" {
                    line.clear();
                } else if line.len() > 0 {
                    lines.push(line.clone());
                    line.clear();
                }
            }
            line.push(char);
            if char == '<' {
                lt += 1;
            } else if char == '>' {
                gt += 1;
            }
            if lt == gt {
                if line != "\n" {
                    lines.push(line.clone());
                }
                line.clear();
            }
        }
        if line.len() > 0 {
            return None;
        }
    }

    // Repackage lines into known svg entities.

    let mut geom = Vec::<Thing>::new();
    let mut i = 0;
    while i < lines.len() {
        let mut line = lines[i].clone();
        println!("line = {} = ${}$", lines[i].len(), lines[i]); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        i += 1;
        if line == "/svg" {
            break;
        }
        if !line.contains(' ') {
            return None;
        }
        let tag = line.between("<", " ").to_string();
        if tag == "svg" {
            continue;
        }
        line = line.after(" ").to_string();

        // Process circle.

        if tag == "circle" {
            println!("processing circle"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            let kv = parse_kv_term(&line);
            if kv.is_none() {
                return None;
            }
            let (mut x, mut y, mut r) = (None, None, None);
            let mut c = None;
            let mut o = 255 as u8; // opacity
            for m in kv.unwrap().iter() {
                let key = &m.0;
                let value = &m.1;
                if key == "stroke" || key == "stroke-width" {
                } else if get_numeric(&key, &value, "cx", &mut x) {
                } else if get_numeric(&key, &value, "cy", &mut y) {
                } else if get_numeric(&key, &value, "r", &mut r) {
                } else if key == "fill" {
                    c = parse_color(&value);
                } else if get_opacity(&key, &value, &mut o) {
                } else {
                    return None;
                }
            }
            if x.is_none() || y.is_none() || r.is_none() || c.is_none() {
                return None;
            }
            geom.push(Thing::Circle(Circle {
                p: Point::new(x.unwrap(), y.unwrap()),
                r: r.unwrap(),
                c: Color::new(c.unwrap().0, c.unwrap().1, c.unwrap().2, o),
            }));

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
                let key = &m.0;
                let value = &m.1;
                if get_numeric(&key, &value, "x1", &mut x1) {
                } else if get_numeric(&key, &value, "y1", &mut y1) {
                } else if get_numeric(&key, &value, "x2", &mut x2) {
                } else if get_numeric(&key, &value, "y2", &mut y2) {
                } else if key == "stroke" {
                    c = parse_color(&value);
                } else if key == "stroke-width" {
                    stroke_width = value.parse::<f32>().ok();
                } else if get_opacity(&key, &value, &mut o) {
                } else {
                    return None;
                }
            }
            if x1.is_none() || x2.is_none() || x2.is_none() || y2.is_none() {
                return None;
            } else if c.is_none() || stroke_width.is_none() {
                return None;
            }
            geom.push(Thing::Segment(Segment {
                p1: Point::new(x1.unwrap(), y1.unwrap()),
                p2: Point::new(x2.unwrap(), y2.unwrap()),
                w: stroke_width.unwrap(),
                c: Color::new(c.unwrap().0, c.unwrap().1, c.unwrap().2, o),
            }));

        // Process polyline.
        } else if tag == "polyline" {
            let kv = parse_kv_term(&line);
            if kv.is_none() {
                return None;
            }
            let mut p = Option::<Vec<Point>>::default();
            let mut c = None;
            let mut o = 255;
            let mut stroke_width = None;
            for m in kv.unwrap().iter() {
                let key = &m.0;
                let value = &m.1;
                if key == "points" {
                    let mut points = Vec::<Point>::new();
                    let ps = value.split(' ').collect::<Vec<&str>>();
                    for x in ps.iter() {
                        let z = x.split(',').collect::<Vec<&str>>();
                        if z.len() != 2 {
                            return None;
                        }
                        let mut p = vec![0.0; 2];
                        for j in 0..z.len() {
                            if z[j].parse::<f32>().is_err() {
                                return None;
                            } else {
                                p[j] = z[j].parse::<f32>().unwrap();
                            }
                        }
                        points.push(Point { x: p[0], y: p[1] });
                    }
                    p = Some(points);
                } else if key == "stroke" {
                    c = parse_color(&value);
                } else if key == "stroke-width" {
                    stroke_width = value.parse::<f32>().ok();
                } else if get_opacity(&key, &value, &mut o) {
                } else {
                    return None;
                }
            }
            if p.is_none() || c.is_none() || stroke_width.is_none() {
                return None;
            }
            geom.push(Thing::PolySegment(PolySegment {
                p: p.unwrap(),
                w: stroke_width.unwrap(),
                c: Color::new(c.unwrap().0, c.unwrap().1, c.unwrap().2, o),
            }));

        // Process rectangle.
        } else if tag == "rect" {
            let kv = parse_kv_term(&line);
            if kv.is_none() {
                return None;
            }
            let (mut x, mut y, mut width, mut height) = (None, None, None, None);
            let (mut fill_color, mut stroke_color) = (None, None);
            let mut stroke_width = None;
            for m in kv.unwrap().iter() {
                let key = &m.0;
                let value = &m.1;
                if get_numeric(&key, &value, "x", &mut x) {
                } else if get_numeric(&key, &value, "y", &mut y) {
                } else if get_numeric(&key, &value, "width", &mut width) {
                } else if get_numeric(&key, &value, "y", &mut y) {
                } else if get_numeric(&key, &value, "width", &mut width) {
                } else if get_numeric(&key, &value, "height", &mut height) {
                } else if key == "fill" {
                    fill_color = parse_color(&value);
                } else if key == "stroke" {
                    stroke_color = parse_color(&value);
                } else if key == "stroke-width" {
                    stroke_width = value.parse::<f32>().ok();
                } else {
                    return None;
                }
            }
            if x.is_none() || y.is_none() || width.is_none() || height.is_none() {
                return None;
            }
            if stroke_width.is_none() ^ stroke_color.is_none() {
                return None;
            }
            if fill_color.is_none() {
                return None;
            }
            if stroke_width.is_none() {
                geom.push(Thing::Rectangle(Rectangle {
                    p: Point::new(x.unwrap(), y.unwrap()),
                    width: width.unwrap(),
                    height: height.unwrap(),
                    fill_color: Color::from_tuple_plus(fill_color.unwrap(), 255),
                    stroke_width: 0.0,
                    stroke_color: Color::new(0, 0, 0, 0),
                }));
            } else {
                geom.push(Thing::Rectangle(Rectangle {
                    p: Point::new(x.unwrap(), y.unwrap()),
                    width: width.unwrap(),
                    height: height.unwrap(),
                    fill_color: Color::from_tuple_plus(fill_color.unwrap(), 255),
                    stroke_width: stroke_width.unwrap(),
                    stroke_color: Color::from_tuple_plus(stroke_color.unwrap(), 255),
                }));
            }

        // Process text.
        } else if tag == "text" && i + 1 < lines.len() && lines[i + 1] == "</text>" {
            let text = lines[i].to_string();
            let mut font_size = None;
            let (mut x, mut y) = (None, None);
            let mut c = Some((0, 0, 0));
            let mut o = 255;
            let mut text_anchor = "left".to_string();
            i += 1;
            let kv = parse_kv(&line);
            if kv.is_none() {
                return None;
            }
            for m in kv.unwrap().iter() {
                let key = &m.0;
                let value = &m.1;
                if get_numeric(&key, &value, "x", &mut x) {
                } else if get_numeric(&key, &value, "y", &mut y) {
                } else if get_numeric(&key, &value, "font-size", &mut font_size) {
                } else if key == "dx" || key == "dy" {
                } else if key == "font-family" && (value == "arial" || value == "Arial") {
                } else if key == "fill" {
                    c = parse_color(&value);
                } else if get_opacity(&key, &value, &mut o) {
                } else if key == "text-anchor" && value == "start" {
                    text_anchor = "start".to_string();
                } else if key == "text-anchor" && value == "middle" {
                    text_anchor = "middle".to_string();
                } else if key == "text-anchor" && value == "end" {
                    text_anchor = "end".to_string();
                } else {
                    return None;
                }
            }
            if font_size.is_none() || x.is_none() || y.is_none() {
                return None;
            }
            let halign;
            if text_anchor == "start" {
                halign = Left;
            } else if text_anchor == "middle" {
                halign = Center;
            } else {
                halign = Right;
            }
            geom.push(Thing::ArialText(ArialText {
                p: Point::new(x.unwrap(), y.unwrap()),
                halign: halign,
                c: Color::new(c.unwrap().0, c.unwrap().1, c.unwrap().2, o),
                t: text,
                font_size: font_size.unwrap(),
            }));
        } else {
            return None;
        }
    }
    Some(geom)
}
