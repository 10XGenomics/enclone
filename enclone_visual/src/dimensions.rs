// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use crate::*;
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
    let current_width = CURRENT_WIDTH.load(SeqCst);
    let current_height = CURRENT_HEIGHT.load(SeqCst);
    if !GRAPHIC_MODE.load(SeqCst) {
        let mut max_height = SVG_HEIGHT as f32;
        if empty {
            max_height = SVG_NULL_HEIGHT as f32;
        }
        max_height *= current_height as f32 / INITIAL_HEIGHT as f32;
        max_height -= 5.0;
        let mut max_width = MAX_WIDTH;
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
    } else {
        const HBORDER: usize = 45;
        const VBORDER: usize = 120;
        let target_width = current_width - HBORDER;
        let target_height = current_height - VBORDER;
        let scale_x = target_width as f32 / width;
        let scale_y = target_height as f32 / height;
        let mut scale = scale_y;
        if scale_x < scale_y {
            scale = scale_x;
        }
        scale
    }
}
