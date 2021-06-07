// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Basic geometric objects, allowing representation of simple SVG files as vectors of such objects.

pub enum Thing {
    Seg,
    ArialText,
    Circle,
    CircleWithToolTip,
    Rectangle,
}

pub struct Color {
    pub r: u8, // red
    pub g: u8, // green
    pub b: u8, // blue
    pub t: u8, // transparency (255 = not transparent at all)
}

impl Color {
    pub fn new(r: u8, g: u8, b: u8, t: u8) -> Color {
        Color {
            r: r,
            g: g,
            b: b,
            t: t,
        }
    }
}

pub struct Point {
    pub x: f32,
    pub y: f32,
}

impl Point {
    pub fn new(x: f32, y:32) -> Point {
        Point {
            x: x,
            y: y,
        }
    }
}

pub struct Circle {
    pub p: Point,
    pub r: f32,
    pub c: Color,
}

pub struct ArialText {
    pub p: Point,
    pub t: String,
    pub font_size: f32,
}

pub struct CircleWithToolTip {
    pub p: Point,
    pub r: f32,
    pub c: Color,
    pub t: String,
}

pub struct Seg {
    pub p1: Point, // start
    pub p2: Point, // stop
    pub w: f32,    // width
    pub c: Color,
}

pub struct Rectangle {
    pub p: Point,            // coordinates of upper left corner
    pub width: f32,          // width
    pub height: f32,         // height
    pub fill_color: Color,   // fill color
    pub stroke_color: Color, // stroke color
}
