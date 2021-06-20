// This code is copied with small changes from the pasteboard crate.  As implemented, it only
// works on a Mac.
//
// The same code presumably works for image types other than PNG, but has only been tested on PNG.

#[cfg(any(target_os = "macos", target_os = "ios"))]
use cocoa::{
    appkit::{NSImage, NSPasteboard},
    base::nil,
    foundation::{NSArray, NSAutoreleasePool, NSData},
};

#[cfg(any(target_os = "macos", target_os = "ios"))]
use libc::c_void;

#[cfg(any(target_os = "macos", target_os = "ios"))]
pub fn copy_png_bytes_to_mac_clipboard(bytes: &[u8]) {
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
                eprintln!("\ncopy to pasteboard failed\n");
                std::process::exit(1);
            }
        }
    }
}

#[cfg(not(any(target_os = "macos", target_os = "ios")))]
pub fn copy_png_bytes_to_mac_clipboard(_bytes: &[u8]) {}
