// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Basic geometric objects, allowing representation of simple SVG files as vectors of such objects.

pub enum Thing {
    Seg,
    ArialText,
    Circle,
    CircleWithToolTip,
}

pub struct Color {
    pub r: u8, // red
    pub g: u8, // green
    pub b: u8, // blue
    pub t: u8, // translucency
}

pub struct Point {
    pub x: f32,
    pub y: f32,
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
