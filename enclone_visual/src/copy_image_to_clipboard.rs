// This code is copied with small changes from the pasteboard crate.  As implemented, it only
// works on a Mac.
//
// The same code presumably works for image types other than PNG, but has only been tested on PNG.

#[cfg(any(target_os = "macos", target_os = "ios"))]
use cocoa::{
    appkit::{NSImage, NSPasteboard},
    base::nil,
    foundation::{NSArray, NSAutoreleasePool, NSData, NSString},
};

#[cfg(any(target_os = "macos", target_os = "ios"))]
use objc::{msg_send, runtime::Object, sel, sel_impl};

#[cfg(any(target_os = "macos", target_os = "ios"))]
pub type Id = *mut Object;

#[cfg(any(target_os = "macos", target_os = "ios"))]
use libc::c_void;

#[cfg(target_os = "linux")]
use arboard::{Clipboard, ImageData};

#[cfg(target_os = "linux")]
use string_utils::*;

#[cfg(any(target_os = "macos", target_os = "ios"))]
use crate::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[cfg(any(target_os = "macos", target_os = "ios"))]
pub fn copy_png_bytes_to_clipboard(bytes: &[u8]) {
    if bytes.len() > 0 {
        unsafe {
            let pool = NSAutoreleasePool::new(nil);
            let data = NSData::dataWithBytes_length_(
                pool,
                bytes.as_ptr() as *const c_void,
                bytes.len() as u64,
            );
            let object = NSImage::initWithData_(NSImage::alloc(pool), data);
            if object != nil {
                let pasteboard = NSPasteboard::generalPasteboard(pool);
                pasteboard.clearContents();
                pasteboard.writeObjects(NSArray::arrayWithObject(pool, object));
            } else {
                xprintln!("\ncopy to pasteboard failed\n");
                std::process::exit(1);
            }
        }
    }
}

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

#[cfg(not(any(target_os = "macos", target_os = "ios", target_os = "linux")))]
pub fn copy_png_bytes_to_clipboard(_bytes: &[u8]) {}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// intended for strings:

#[cfg(any(target_os = "macos", target_os = "ios"))]
pub fn copy_bytes_to_clipboard(bytes: &[u8]) {
    if bytes.len() > 0 {
        unsafe {
            let pool = NSAutoreleasePool::new(nil);
            let data = NSData::dataWithBytes_length_(
                pool,
                bytes.as_ptr() as *const c_void,
                bytes.len() as u64,
            );
            let string: Id = msg_send![NSString::alloc(pool), initWithData:data encoding:4];
            let object = string;
            if object != nil {
                let pasteboard = NSPasteboard::generalPasteboard(pool);
                pasteboard.clearContents();
                pasteboard.writeObjects(NSArray::arrayWithObject(pool, object));
            } else {
                xprintln!("\ncopy to pasteboard failed\n");
                std::process::exit(1);
            }
        }
    }
}

#[cfg(target_os = "linux")]
pub fn copy_bytes_to_clipboard(bytes: &[u8]) {
    if bytes.len() > 0 {
        let mut clipboard = Clipboard::new().unwrap();
        let the_string = strme(&bytes);
        clipboard.set_text(the_string.into()).unwrap();
    }
}

#[cfg(not(any(target_os = "macos", target_os = "ios", target_os = "linux")))]
pub fn copy_bytes_to_clipboard(_bytes: &[u8]) {}
