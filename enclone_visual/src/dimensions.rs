// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::SeqCst;

pub const MAX_WIDTH: f32 = 770.0;

pub const INITIAL_WIDTH: u32 = 1100;
pub const INITIAL_HEIGHT: u32 = 1060;

pub const SVG_NULL_HEIGHT: u16 = 190;
pub const SVG_HEIGHT: u16 = 400;

pub static CURRENT_WIDTH: AtomicUsize = AtomicUsize::new(0);
pub static CURRENT_WIDTH_LAST_SEEN: AtomicUsize = AtomicUsize::new(0);

pub fn get_graphic_scale(width: f32, height: f32, empty: bool) -> f32 {
    let mut max_height = SVG_HEIGHT as f32;
    if empty {
        max_height = SVG_NULL_HEIGHT as f32;
    }
    max_height -= 5.0;
    let mut max_width = MAX_WIDTH;
    let current_width = CURRENT_WIDTH.load(SeqCst);
    if current_width > INITIAL_WIDTH as usize {
        max_width += (current_width - INITIAL_WIDTH as usize) as f32;
    }
    let scale_x = max_width / width;
    let scale_y = max_height / height;
    let mut scale = scale_y;
    if scale_x < scale_y {
        scale = scale_x;
    }
    scale
}
