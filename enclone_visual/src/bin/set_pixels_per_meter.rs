// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// set_pixels_per_meter n in out
//
// Modify or add the pHYs chunk in a PNG file, to set the pixels per meter value to n.
// See https://www.w3.org/TR/2003/REC-PNG-20031110/#11pHYs
// See also print_phys.rs.
//
// Note that Mac tools including screen capture add an iDOT chunk. That's an undocumented
// apple-ism, which appears to interfer with pHYs.  Googling will get you a partial structure
// description, obtained by reverse engineering.

use crc::*;
use io_utils::*;
use pretty_trace::*;
use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let n = args[1].force_usize() as u32;
    let (infile, outfile) = (&args[2], &args[3]);

    // Form the pHYs chunk.

    let mut bytes = Vec::<u8>::new();
    {
        let mut data = Vec::<u8>::new();
        data.append(&mut n.to_be_bytes().to_vec());
        data.append(&mut n.to_be_bytes().to_vec());
        data.append(&mut vec![1 as u8]);
        let len = 9 as u32;
        bytes.append(&mut len.to_be_bytes().to_vec());
        bytes.append(&mut b"pHYs".to_vec());
        bytes.append(&mut data);
        let crc = Crc::<u32>::new(&CRC_32_ISO_HDLC);
        let cs = crc.checksum(&bytes[4..bytes.len()]);
        bytes.append(&mut cs.to_be_bytes().to_vec());
    }

    // Copying the file, inserting a new pHYs chunk, and deleting the existing one.  We put the
    // new one right after the IHDR chunk.  According to the spec, it needs to go after IHDR and
    // before the first IDAT.

    let mut f = open_for_read![&infile];
    let mut file = open_for_write_new![&outfile];
    let mut header = vec![0 as u8; 8];
    f.read_exact(&mut header).unwrap();
    file.write_all(&header).unwrap();
    let mut first = true;
    let mut pos = 8;
    loop {
        f.seek(SeekFrom::Start(pos)).unwrap();
        let mut x = vec![0 as u8; 4];
        let res = f.read_exact(&mut x);
        if res.is_err() {
            break;
        }
        let len = u32::from_be_bytes([x[0], x[1], x[2], x[3]]);
        let total = len + 12;
        f.seek(SeekFrom::Start(pos + 4)).unwrap();
        f.read_exact(&mut x).unwrap();
        if x != b"pHYs" {
            f.seek(SeekFrom::Start(pos)).unwrap();
            let mut x = vec![0 as u8; total as usize];
            f.read_exact(&mut x).unwrap();
            file.write_all(&x).unwrap();
        }
        if first {
            file.write_all(&bytes).unwrap();
            first = false;
        }
        pos += total as u64;
    }
}
