// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// set_pixels_per_meter n in out
//
// Modify or add the pHYs chunk in a PNG file, to set the pixels per meter value to n.
// See https://www.w3.org/TR/2003/REC-PNG-20031110/#11pHYs
// See also print_phys.rs.

use crc::*;
use io_utils::*;
use pretty_trace::*;
use std::env;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
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

    // Locate the exiting pHYs chunk.

    let mut f = open_for_read![&infile];
    let mut pos = 8 as u64;
    let mut loc = None;
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
        if chunk_type == b"pHYs" {
            loc = Some(pos - 8);
        }
        pos += len as u64 + 4;
    }

    // Create the new file.

    std::fs::copy(&infile, &outfile).unwrap();
    if loc.is_none() {
        let mut file = OpenOptions::new()
            .write(true)
            .append(true)
            .open(&outfile)
            .unwrap();
        file.write_all(&bytes).unwrap();
    } else {
        let pos = loc.unwrap();
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&outfile)
            .unwrap();
        file.seek(SeekFrom::Start(pos)).unwrap();
        file.write_all(&bytes).unwrap();
    }
}
