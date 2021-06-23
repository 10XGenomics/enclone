// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Basic geometric objects, allowing representation of simple SVG files as vectors of such objects.

#[derive(PartialEq)]
pub enum Geometry {
    Segment(Segment),
    Text(Text),
    Circle(Circle),
    CircleWithTooltip(CircleWithTooltip),
    Rectangle(Rectangle),
    PolySegment(PolySegment),
}

#[derive(PartialEq)]
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

    pub fn from_tuple(x: (u8, u8, u8, u8)) -> Color {
        Color {
            r: x.0,
            g: x.1,
            b: x.2,
            t: x.3,
        }
    }

    pub fn from_tuple_plus(x: (u8, u8, u8), t: u8) -> Color {
        Color {
            r: x.0,
            g: x.1,
            b: x.2,
            t: t,
        }
    }
}

#[derive(Default, PartialEq)]
pub struct Point {
    pub x: f32,
    pub y: f32,
}

impl Point {
    pub fn new(x: f32, y: f32) -> Point {
        Point { x: x, y: y }
    }
}

#[derive(PartialEq)]
pub struct Circle {
    pub p: Point,
    pub r: f32,
    pub c: Color,
}

#[derive(PartialEq)]
pub enum HorizontalAlignment {
    Left,
    Center,
    Right,
}

#[derive(PartialEq)]
pub struct Text {
    pub p: Point,
    pub halign: HorizontalAlignment,
    pub t: String,
    pub font: String,
    pub font_size: f32,
    pub c: Color,
    pub rotate: [f32; 3],
}

#[derive(PartialEq)]
pub struct CircleWithTooltip {
    pub p: Point,
    pub r: f32,
    pub c: Color,
    pub t: String,
}

#[derive(PartialEq)]
pub struct Segment {
    pub p1: Point, // start
    pub p2: Point, // stop
    pub w: f32,    // width
    pub c: Color,
}

#[derive(PartialEq)]
pub struct PolySegment {
    pub p: Vec<Point>, // points
    pub w: f32,        // width
    pub c: Color,      // color
}

#[derive(PartialEq)]
pub struct Rectangle {
    pub p: Point,            // coordinates of upper left corner
    pub width: f32,          // width
    pub height: f32,         // height
    pub fill_color: Color,   // fill color
    pub stroke_width: f32,   // stroke width
    pub stroke_color: Color, // stroke color
}
