// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::SeqCst;

// These are the window dimensions that are requested initially.

pub const INITIAL_WIDTH: u32 = 1100;
pub const INITIAL_HEIGHT: u32 = 1060;

// Tracking the actual screen dimensions.

pub static CURRENT_WIDTH: AtomicUsize = AtomicUsize::new(0);
pub static CURRENT_WIDTH_LAST_SEEN: AtomicUsize = AtomicUsize::new(0);
pub static CURRENT_HEIGHT: AtomicUsize = AtomicUsize::new(0);
pub static CURRENT_HEIGHT_LAST_SEEN: AtomicUsize = AtomicUsize::new(0);

// The default height of the graphic part of the screen if there is not graphic, and if there
// is a graphic.

pub const SVG_NULL_HEIGHT: u16 = 190;
pub const SVG_HEIGHT: u16 = 400;

// The default width of the graphic.

pub const MAX_WIDTH: f32 = 770.0;

// The function that determines the scale for the graphic object.

pub fn get_graphic_scale(width: f32, height: f32, empty: bool) -> f32 {
    let mut max_height = SVG_HEIGHT as f32;
    if empty {
        max_height = SVG_NULL_HEIGHT as f32;
    }
    max_height *= CURRENT_HEIGHT.load(SeqCst) as f32 / INITIAL_HEIGHT as f32;
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
