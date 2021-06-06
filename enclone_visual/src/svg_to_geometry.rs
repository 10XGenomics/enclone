// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
//
// Attempt to convert an SVG object to geometries.  This is possibly temporary.  It is designed
// to work only with certain sorts of SVG objects.

use crate::genometry;

pub fn parse_kv(line: &str) -> Option<Vec<(String, String)>> {
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
        let value = value.after("\"").rev_before("\"");
        kv.push((key.to_string(), value.to_string()));
    }
    kv
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
            if line.ends_with(" \>") {
                line = line.before(" \>");
            } else if line.ends_with("\>") {
                line = line.before("\>");
            } else {
                return None;
            }
            let kv = parse_kv(&line);
            if kv.is_none() {
                return None;
            }
            let (mut x, mut y, mut r) = (None, None, None);
            let mut c = None;
            for m in kv.unwrap().iter() {
                if *m == "stroke" {
                } else if *m.0 == "stroke-width" {
                } else if *m.0 == "cx" {
                    if !m.parse::<f32>().is_ok() {
                        return None;
                    } else {
                        x = m.1.force_f32();
                    }
                } else if *m.0 == "cy" {
                    if !m.0.parse::<f32>().is_ok() {
                        return None;
                    } else {
                        y = m.1.force_f32();
                    }
                } else if *m.0 == "r" {
                    if !m.0.parse::<f32>().is_ok() {
                        return None;
                    } else {
                        r = m.1.force_f32();
                    }
                } else if *m.0 == "fill" && *m.1.starts_with("rgb(" && *m.1.ends_with(")") {
                    let rgb = *m.0.split(',').collect::<Vec<&str>>();
                    if rgb.len() != 3 {
                        return None;
                    }
                    for j in 0..3 {
                        if !rgb[j].parse::<u8>().is_ok() {
                            return None;
                        }
                    }
                    c = (rgb[0].force_u8(), rgb[1].force_u8(), rgb[2].force_u8());
                }
                
                // fill="#FF0000" 
                // opacity="1" 

===================================================================================================

        } else if tag == "rect" {
            // <rect x="10" y="407.33" width="147.5" height="260" 
            // style="fill:white;stroke:black;stroke-width:2" />

            ...
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
