// Copyright (c) 2022 10x Genomics, Inc. All rights reserved.

use crate::capture_as_file;
use crate::copy_image_to_clipboard::*;
use crate::get_window_id;
use perf_stats::*;
use std::env;
use std::fs::{remove_file, File};
use std::io::Read;
use std::thread;
use std::time::{Duration, Instant};

// Copy window image to clipboard.  If the environment variable ENCLONE_VIS_SNAPSHOT is defined,
// also save to that file.  This only works on a Mac.

pub fn snapshot(start: &Option<Instant>) {
    let mut filename = "/tmp/enclone_visual_snapshot.png".to_string();
    let mut snapshot = false;
    for (key, value) in env::vars() {
        if key == "ENCLONE_VIS_SNAPSHOT" {
            snapshot = true;
            filename = value.to_string();
        }
    }
    capture_as_file(&filename, get_window_id());
    let mut bytes = Vec::<u8>::new();
    {
        let mut f = File::open(&filename).unwrap();
        f.read_to_end(&mut bytes).unwrap();
    }
    if !snapshot {
        remove_file(&filename).unwrap();
    }
    copy_png_bytes_to_clipboard(&bytes);
    const MIN_SLEEP: f64 = 0.4;
    let used = elapsed(&start.unwrap());
    if used < MIN_SLEEP {
        let ms = ((MIN_SLEEP - used) * 1000.0).round() as u64;
        thread::sleep(Duration::from_millis(ms));
    }
}

/*

Below, there is some code that will snapshot a window under Linux.  However, it is dependent
on windows being managed by X, and most Linus distros now use Wayland, for which this won't
work.  This was an experimental main program.

extern crate libwmctl;
extern crate pretty_trace;
extern crate image;
extern crate x11;
extern crate libc;
extern crate string_utils;
extern crate arboard;
extern crate png_decoder;

// pub mod util;
pub mod xwrap;

use pretty_trace::*;

use std::fs::File;
use std::path::Path;

use std::process::Command;

use x11::xlib;

use crate::xwrap::Display;

use string_utils::strme;
use string_utils::TextUtils;

use std::io::Read;

#[cfg(target_os = "linux")]
use arboard::{Clipboard, ImageData};

#[cfg(target_os = "linux")]
pub fn copy_png_bytes_to_clipboard(bytes: &[u8]) {
    let (header, image_data) = png_decoder::decode(&bytes).unwrap();
    let (width, height) = (header.width as usize, header.height as usize);
    let mut clipboard = Clipboard::new().unwrap();
    let img_data = ImageData {
        width: width,
        height: height,
        bytes: image_data.into(),
    };
    clipboard.set_image(img_data).unwrap();
}

fn main() {
    PrettyTrace::new().on();

    // Get the enclone visual window id.

    let mut window_id = String::new();
    let new = Command::new("wmctrl")
        .arg("-l")
        .output()
        .unwrap_or_else(|_| {
            eprintln!("\nThe executable wmctrl could not be found.\n");
            eprintln!("\nYou may need to install it via:\n\
                sudo apt install wmctrl\n"
            );
            std::process::exit(1);
        });
    if new.status.code() != Some(0) {
        eprintln!("\nwmctrl -l failed\n");
        eprintln!("stderr = {}", strme(&new.stderr));
        std::process::exit(1);
    }
    let s = strme(&new.stdout);
    for line in s.lines() {
        if line.ends_with(" EncloneVisual") {
            window_id = line.before(" ").to_string();
        }
    }
    if window_id == "" {
        eprintln!("\nFailed to find enclone visual window.  The output of wmctrl -l is:\n{}", s);
        std::process::exit(1);
    }
    println!("window id = {}", window_id);

    // code taken from shotgun

    let display = match Display::open(None) {
        Some(d) => d,
        None => {
            eprintln!("Failed to open display");
            std::process::exit(1);
        }
    };

    let window = xwrap::util::parse_int::<xlib::Window>(&window_id).unwrap();

    let output_format = image::ImageOutputFormat::Png;

    let window_rect = display.get_window_rect(window);

    let sel =
        xwrap::util::Rect {
            x: 0,
            y: 0,
            w: window_rect.w,
            h: window_rect.h,
        };

    let image = match display.get_image(window, sel, xwrap::ALL_PLANES, xlib::ZPixmap) {
        Some(i) => i,
        None => {
            eprintln!("Failed to get image from X");
            std::process::exit(1);
        },
    };

    let image = match image.into_image_buffer() {
        Some(i) => image::DynamicImage::ImageRgba8(i),
        None => {
            eprintln!("Failed to convert captured framebuffer, only 24/32 \
                      bit (A)RGB8 is supported");
            std::process::exit(1);
        }
    };

    let filename = "/tmp/enclone_visual_snapshot.png";
    {
        match File::create(&Path::new(&filename)) {
            Ok(mut f) => image.write_to(&mut f, output_format).expect("Writing to file failed"),
            Err(e) => {
                eprintln!("Failed to create {}: {}", filename, e);
                std::process::exit(1);
            },
        }
    }
    let mut bytes = Vec::<u8>::new();
    {
        let mut f = File::open(&filename).unwrap();
        f.read_to_end(&mut bytes).unwrap();
    }
    std::fs::remove_file(&filename).unwrap();
    copy_png_bytes_to_clipboard(&bytes);
}

*/

// Some code that works even less.  This was supposed to be an alternative to using wmctrl.

/*

// Doesn't work:

use libwmctl::prelude::*;
let wmctl = WmCtl::connect().unwrap();
let (_, wm_name) = wmctl.winmgr().unwrap();
let win = wmctl.active_win().unwrap();
println!("X11 Information");
println!("-----------------------------------------------------------------------");
println!("Window Manager:    {}", wm_name);
println!("Composite Manager: {}", wmctl.composite_manager().unwrap());
println!("Root Window:       {}", wmctl.root());
println!("Work area:         {}x{}", wmctl.work_width(), wmctl.work_height());
println!("Screen Size:       {}x{}", wmctl.width(), wmctl.height());
println!("Desktops:          {}", wmctl.desktops().unwrap());
println!();
println!("Active Window");
println!("{:-<120}", "");

println!("{:<8} {:<3} {:<6} {:<5} {:<5} {:<4} {:<4} {:<8} {:<7} {:<18} {:<18} {}", "ID", "DSK", "PID", "X", "Y", "W", "H", "BORDERS", "TYPE", "STATE", "CLASS", "NAME");

let pid = wmctl.win_pid(win).unwrap_or(-1);
let desktop = wmctl.win_desktop(win).unwrap_or(-1);
let typ = wmctl.win_type(win).unwrap_or(WinType::Invalid);
let states = wmctl.win_state(win).unwrap_or(vec![WinState::Invalid]);
let (x, y, w, h) = wmctl.win_geometry(win).unwrap_or((0,0,0,0));
let (l, r, t, b) = wmctl.win_borders(win).unwrap_or((0, 0, 0, 0));
let class = wmctl.win_class(win).unwrap_or("".to_owned());
let name = wmctl.win_name(win).unwrap_or("".to_owned());
println!("{:<8} {:<3} {:<6} {:<5} {:<5} {:<4} {:<4} {:<8} {:<7} {:<18} {:<18} {}",
    format!("{:0>8}", win), format!("{:>2}", desktop), pid,
    format!("{:<4}", x), format!("{:<4}", y), format!("{:<4}", w), format!("{:<4}", h),
    format!("{},{},{},{}", l, r, t, b),
    typ.to_string(), format!("{:?}", states), class, name);

*/

// Cargo.toml stuff relevant to the above:

/*

[dependencies]
string_utils = { version = "0.1", git = "https://github.com/10XGenomics/rust-toolbox.git", branch = "master" }
getopts = "0.2"
libc = "0.2"
num-traits = "0.2"
libwmctl = "0.0"
pretty_trace = { version = "0.5", git = "https://github.com/10XGenomics/rust-toolbox.git", branch = "master" }
png-decoder = "0.1"

[dependencies.image]
default-features = false
version = "0.23"
features = ["png", "pnm"]

[dependencies.x11]
version = "2.18"
features = ["xlib", "xrandr"]

[target.'cfg(target_os = "linux")'.dependencies]
arboard = "2"

*/
