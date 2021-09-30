// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use ansi_escape::*;
use io_utils::*;
use std::io::Write;
use string_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Extract the @font-face content from the current css file.

pub fn font_face_in_css() -> String {
    let f = include_str!["enclone_css_v2.css"];
    let mut x = String::new();
    let mut in_font_face = false;
    let mut count = 0;
    for line in f.lines() {
        if line.starts_with("@font-face") {
            in_font_face = true;
            count += 1;
        }
        if in_font_face {
            x += &format!("{}\n", line);
        }
        if line == "}" {
            in_font_face = false;
        }
    }
    assert_eq!(count, 2); // because there are two fonts: regular and bold
    x
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn emit_codon_color_escape(c: &[u8], log: &mut Vec<u8>) {
    let mut s = 0;
    if c == b"CTG" {
        s = 3;
    } else if c == b"AGG" {
        s = 1;
    } else if c == b"AGT" {
        s = 2;
    } else {
        for i in 0..3 {
            if c[i] == b'A' {
            } else if c[i] == b'C' {
                s += 1;
            } else if c[i] == b'G' {
                s += 2;
            } else if c[i] == b'T' {
                s += 3;
            } else {
                panic!("Illegal codon: \"{}\".", strme(c));
            }
        }
    }
    let s = s % 6;
    print_color(s, log);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn color_by_property(c: &[u8], mut log: &mut Vec<u8>) {
    for i in 0..c.len() {
        let mut color = 7;
        if c[i] == b'A'
            || c[i] == b'G'
            || c[i] == b'I'
            || c[i] == b'L'
            || c[i] == b'P'
            || c[i] == b'V'
        {
            color = 0;
        } else if c[i] == b'F' || c[i] == b'W' || c[i] == b'Y' {
            color = 1;
        } else if c[i] == b'D' || c[i] == b'E' {
            color = 2;
        } else if c[i] == b'R' || c[i] == b'H' || c[i] == b'K' {
            color = 3;
        } else if c[i] == b'S' || c[i] == b'T' {
            color = 4;
        } else if c[i] == b'C' || c[i] == b'M' {
            color = 5;
        } else if c[i] == b'N' || c[i] == b'Q' {
            color = 6;
        }
        if color < 7 {
            print_color(color, log);
        }
        fwrite!(log, "{}", c[i] as char);
        if color < 7 {
            emit_end_escape(&mut log);
        }
    }
}
