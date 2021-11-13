// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Takes a single argument, the name of a PNG file.  Print the pHYs chunk, which
// if present contains the physical pixel dimensions.
// See https://www.w3.org/TR/2003/REC-PNG-20031110/#11pHYs

use crc::*;
use io_utils::*;
use pretty_trace::*;
use std::env;

use std::io::{Read, Seek, SeekFrom};
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let mut f = open_for_read![&filename];
    let mut pos = 8 as u64;
    let mut found = false;
    loop {
        f.seek(SeekFrom::Start(pos)).unwrap();
        let mut x = vec![0 as u8; 4];
        let res = f.read_exact(&mut x);
        if res.is_err() {
            break;
        }
        let len = u32::from_be_bytes([x[0], x[1], x[2], x[3]]);
        pos += 4;
        f.seek(SeekFrom::Start(pos)).unwrap();
        let mut chunk_type = vec![0 as u8; 4];
        f.read_exact(&mut chunk_type).unwrap();
        pos += 4;
        println!("found chunk of type {}", strme(&chunk_type));
        if chunk_type == b"pHYs" {
            found = true;
            f.seek(SeekFrom::Start(pos)).unwrap();
            f.read_exact(&mut x).unwrap();
            let xpix = u32::from_be_bytes([x[0], x[1], x[2], x[3]]);
            f.seek(SeekFrom::Start(pos + 4)).unwrap();
            f.read_exact(&mut x).unwrap();
            let ypix = u32::from_be_bytes([x[0], x[1], x[2], x[3]]);
            f.seek(SeekFrom::Start(pos + 8)).unwrap();
            let mut unit_spec = vec![0 as u8];
            f.read_exact(&mut unit_spec).unwrap();
            println!("\npixels per x unit = {}", xpix);
            println!("pixels per y unit = {}", ypix);
            println!("unit specifier = {}\n", unit_spec[0]);
            f.seek(SeekFrom::Start(pos + 9)).unwrap();
            f.read_exact(&mut x).unwrap();
            let checksum = u32::from_be_bytes([x[0], x[1], x[2], x[3]]);
            println!("checksum in file  = {}", checksum);
            f.seek(SeekFrom::Start(pos - 4)).unwrap();
            let mut bytes = vec![0 as u8; 13];
            f.read_exact(&mut bytes).unwrap();
            let crc = Crc::<u32>::new(&CRC_32_ISO_HDLC);
            let cs = crc.checksum(&bytes);
            println!("computed checksum = {}", cs);
            // break;
        }
        pos += len as u64 + 4;
    }
    if !found {
        println!("\nThere is no pHYs block.\n");
    }
}
