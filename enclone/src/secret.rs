// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Get counts for secreted and membrane proteins.
//
// The basic problem with this is the relevant exon junctions are too far from the 5' end,
// so there is not enough signal to be useful.
//
// Also what is done below appears not to handle the IGHG case correctly.
//
// For IGHG, the boundaries should be
// (A)CH2-(B)CH3-CHS  [secreted]
// (A)CH2-(B)Mx [membrane].

use enclone_core::defs::*;
use std::collections::HashMap;
use std::process::Command;
use string_utils::*;
use vector_utils::*;

// copied from tenkit2/pack_dna.rs:

pub fn reverse_complement(x: &mut Vec<u8>) {
    x.reverse();
    for j in 0..x.len() {
        x[j] = match x[j] {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => x[j],
        }
    }
}

pub fn fetch_secmem(ctl: &mut EncloneControl) -> Result<(), String> {
    // Define the CH3 exon boundaries, and the sequences that could follow it, both in
    // GRCh38 or GRCm38 coordinates.

    let species = &ctl.gen_opt.species;
    let ch3;
    let fol;
    if species == "human" {
        ch3 = vec![
            ('-', "chr14:105600482-105600805"),
            ('-', "chr14:105840368-105840691"),
            ('-', "chr14:105854918-105855235"),
        ];
        fol = vec![
            ("TACCTG", "M1"),
            ("GTGAAA", "M2"),
            ("GTGAAG", "M2"),
            ("TGTGAA", "M2?"),
            ("GGCTCT", "M"),
            ("GGAAAC", "S"),
            ("TATGTA", "S"),
            ("GGCCCG", "S"),
            ("GCCCGC", "S"),
            ("GGACAG", "S"),
            ("GGGGTG", "S"),
        ];
    } else {
        ch3 = vec![
            ('-', "chr12:113414273-113414593"),
            ('-', "chr12:113271711-113272031"),
            ('-', "chr12:113421370-113421686"),
        ];
        fol = vec![
            ("GAGCTAGAC", "M1"),
            ("GAGCTGGAA", "M1"),
            ("GAGGGGGAG", "M1"),
            ("GGCATAGTC", "M1"),
            ("GGGCTAGAC", "M1"),
            ("GGGCTGCAA", "M1"),
            ("GTGAAA", "M2"),
            ("GTGAAG", "M2"),
            ("GAACGTCAA", "M"),
            ("GGAAGAGCC", "S"),
            ("GGCAGACCG", "S"),
            ("GGGCCAGTA", "S"),
            ("GGGCTAGTC", "S"),
            ("GGGTCAGTA", "S"),
            ("GTGAACACC", "S"),
            ("TGAACACCT", "S"),
            ("GAGGTGCAC", "S"),
            ("GCCAGCGCT", "S"),
            ("GGCCAGCGC", "S"),
        ];
    }

    // Traverse the datasets.

    for q in 0..ctl.origin_info.n() {
        let mut data = Vec::<(String, String, String)>::new(); // (barcode, umi, class)
        let bam = format!("{}/possorted_genome_bam.bam", ctl.origin_info.gex_path[q]);

        // Traverse the boundaries.

        for i in 0..ch3.len() {
            // Call samtools.

            let o = Command::new("samtools")
                .arg("view")
                .arg(&bam)
                .arg(&ch3[i].1)
                .output()
                .expect("failed to execute samtools");

            // Parse the output.

            let o = String::from_utf8(o.stdout).unwrap();
            for line in o.lines() {
                let fields = line.split('\t').collect::<Vec<&str>>();
                let pos = fields[3].force_usize();
                let cigar = fields[5];
                let seq = fields[9];
                let (mut barcode, mut umi) = (String::new(), String::new());
                for j in 11..fields.len() {
                    if fields[j].starts_with("CB:Z:") {
                        barcode = fields[j].after("CB:Z:").to_string();
                    } else if fields[j].starts_with("UB:Z:") {
                        umi = fields[j].after("UB:Z:").to_string();
                    }
                }
                if barcode.len() == 0 {
                    continue;
                }

                // Parse cigar string.

                let mut cg = Vec::<Vec<u8>>::new(); // pieces of cigar string
                let mut piece = Vec::<u8>::new();
                for c in cigar.chars() {
                    piece.push(c as u8);
                    if c.is_ascii_alphabetic() {
                        cg.push(piece.clone());
                        piece.clear();
                    }
                }
                if !piece.is_empty() {
                    cg.push(piece);
                }

                // Determine if the sequence is reaching off the end of the reference interval.
                // This is for the left end in the rc case, and right end otherwise.  The latter
                // case has not been tested.

                let mut ref_pos = pos;
                let mut read_pos = 1;
                let low = ch3[i].1.after(":").before("-").force_usize();
                let high = ch3[i].1.after("-").force_usize();
                let mut ext = 0;
                let mut ext_seq = Vec::<u8>::new();
                for j in 0..cg.len() {
                    let x = cg[j][cg[j].len() - 1];
                    let n = strme(&cg[j][0..cg[j].len() - 1]).force_usize();
                    if x == b'M' {
                        if ch3[i].0 == '-' {
                            if read_pos > 1 && ref_pos < high && ref_pos + n > low {
                                if read_pos + low > ref_pos + 1 {
                                    ext = read_pos + low - ref_pos - 1;
                                    ext_seq = seq.as_bytes()[0..ext].to_vec();
                                    reverse_complement(&mut ext_seq);
                                    break;
                                }
                            }
                        } else {
                            if ref_pos <= high && ref_pos + n > high {
                                ext = ref_pos + n - high;
                                ext_seq = seq.as_bytes()[seq.len() - ext..].to_vec();
                                break;
                            }
                        }
                        ref_pos += n;
                        read_pos += n;
                    } else if x == b'N' {
                        ref_pos += n;
                    } else if x == b'S' {
                        read_pos += n;
                    } else if x == b'I' {
                        read_pos += n;
                    } else if x == b'D' {
                        ref_pos += n;
                    } else {
                        return Err(format!("\nUnexpected character in cigar string.\n"));
                    }
                }

                // Check if extension long enough.

                if (species == "human" && ext < 6) || (species == "mouse" && ext < 9) {
                    continue;
                }

                // Print.

                let mut class;
                if species == "human" {
                    class = stringme(&ext_seq[0..6]);
                } else {
                    class = stringme(&ext_seq[0..9]);
                }
                for j in 0..fol.len() {
                    if strme(&ext_seq).starts_with(&fol[j].0) {
                        class = fol[j].1.to_string();
                    }
                }
                data.push((barcode, umi, class));
            }
        }

        // Fill in the map.

        let mut h = HashMap::<String, (usize, usize)>::new();
        data.sort();
        let mut i = 0;
        while i < data.len() {
            let j = next_diff1_3(&data, i as i32) as usize;
            let (mut sec, mut mem) = (0, 0);
            let mut k = i;
            while k < j {
                // let l = next_diff12_3(&data, k as i32) as usize; // crashed loader
                let mut l = k;
                while l < j {
                    if data[l].1 != data[k].1 {
                        break;
                    }
                    l += 1;
                }
                let (mut s, mut m) = (0, 0);
                for z in k..l {
                    if data[z].2.starts_with('M') {
                        m += 1;
                    } else if data[z].2.starts_with('S') {
                        s += 1;
                    }
                }
                if s > 0 && m == 0 {
                    sec += 1;
                } else if s == 0 && m > 0 {
                    mem += 1;
                }
                k = l;
            }
            h.insert(data[i].0.clone(), (sec, mem));
            i = j;
        }
        ctl.origin_info.secmem.push(h);
    }
    Ok(())
}
