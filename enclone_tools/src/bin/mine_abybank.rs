// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Extract some summary data from an abybank download.
//
// Usage:
// mine_abybank abybank/NR_LH_Combined_Martin > mine_abybank.out
// * see below for the input data abybank/NR_LH_Combined_Martin can be obtained.
//
// These are the data at http://www.abybank.org/abdb after selecting
// complete antibodies / Martin / complete dataset / NR.  The direct link for this was
// http://www.abybank.org/abdb/Data/NR_LH_Combined_Martin.tar.bz2.
// Then bunzip2 followed by tar xf - was run on these data.
//
// Distribution of pdb files when downloaded in 7/12/20.
//
// 657  human
// 522  mouse
//  15  rabbit (Oryctolagus cuniculus)
//  15  rhesus macaque
//  11  Rattus norvegicus or Rattus rattus
//   5  llama
//   4  cow
//   2  hamster (Cricetulus migratorius)
//   1  chimp
//   1  Macaca fascicularis
//   1  hare (Lepus curpaeums)
//   1  chicken.

use io_utils::*;
use itertools::Itertools;
use pretty_trace::*;
use std::env;
use std::fs::read_dir;
use std::io::BufRead;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let args: Vec<String> = env::args().collect();
    let dir = &args[1];
    let all = read_dir(&dir).unwrap();
    for f in all {
        let f = f.unwrap().path();
        let f = f.to_str().unwrap();
        let id = f.rev_after("/").before(".pdb");
        let f = open_for_read![&f];
        let mut species = String::new();
        let mut light_aa = Vec::<u8>::new();
        let mut heavy_aa = Vec::<u8>::new();
        let mut light_no = Vec::<String>::new();
        let mut heavy_no = Vec::<String>::new();
        let mut last = String::new();
        for line in f.lines() {
            let mut s = line.unwrap();
            if s.contains("SPECIES") && species.len() == 0 {
                species = s.after(":").after(": ").to_string();
            } else if s.starts_with("ATOM") {
                s = s.replace("\t", " ");
                s = s.replace("  ", " ");
                s = s.replace("  ", " ");
                s = s.replace("  ", " ");
                let fields = s.split(' ').collect::<Vec<&str>>();
                let chain = fields[4];
                if chain != "L" && chain != "H" {
                    continue;
                }
                let pos = fields[5];
                if pos.to_string() == last {
                    continue;
                }
                last = pos.to_string();
                let aa = fields[3];
                let mut a = aa.as_bytes()[0];
                if aa == "PHE" {
                    a = b'F';
                } else if aa == "TYR" {
                    a = b'Y';
                } else if aa == "GLN" {
                    a = b'Q';
                } else if aa == "ASN" {
                    a = b'N';
                } else if aa == "LYS" {
                    a = b'K';
                } else if aa == "ASP" {
                    a = b'D';
                } else if aa == "GLU" {
                    a = b'E';
                } else if aa == "TRP" {
                    a = b'W';
                } else if aa == "ARG" {
                    a = b'R';
                }
                if chain == "L" {
                    light_no.push(pos.to_string());
                    light_aa.push(a);
                } else if chain == "H" {
                    heavy_no.push(pos.to_string());
                    heavy_aa.push(a);
                }
            }
        }
        if species.contains("SYNTHETIC") || species.contains(",") || species.contains("VIRUS") {
            continue;
        }
        println!("\nid = {}", id);
        println!("species = {}", species);
        println!("light = {}", strme(&light_aa));
        println!("light_no = {}", light_no.iter().format(","));
        println!("heavy = {}", strme(&heavy_aa));
        println!("heavy_no = {}", heavy_no.iter().format(","));

        let mut fr3 = Vec::<u8>::new();
        for i in 0..light_no.len() {
            let mut s = light_no[i].clone();
            if s.contains('A') {
                s = s.before("A").to_string();
            } else if s.contains('B') {
                s = s.before("B").to_string();
            } else if s.contains('C') {
                s = s.before("C").to_string();
            } else if s.contains('D') {
                s = s.before("D").to_string();
            } else if s.contains('E') {
                s = s.before("E").to_string();
            } else if s.contains('F') {
                s = s.before("F").to_string();
            }
            let n = s.force_usize();
            if n >= 57 && n <= 88 {
                fr3.push(light_aa[i]);
            }
        }
        println!("light fr3 = {}, len = {}", strme(&fr3), fr3.len());

        let mut fr2 = Vec::<u8>::new();
        for i in 0..heavy_no.len() {
            let mut s = heavy_no[i].as_bytes();
            for j in 0..s.len() {
                if s[j] > b'9' {
                    s = &s[0..j];
                    break;
                }
            }
            let n = strme(&s).force_usize();
            if n >= 33 && n <= 49 {
                fr2.push(heavy_aa[i]);
            }
        }
        println!("heavy fr2 = {}, len = {}", strme(&fr2), fr2.len());

        let mut fr3 = Vec::<u8>::new();
        for i in 0..heavy_no.len() {
            let mut s = heavy_no[i].as_bytes();
            for j in 0..s.len() {
                if s[j] > b'9' {
                    s = &s[0..j];
                    break;
                }
            }
            let n = strme(&s).force_usize();
            // had 94???
            if n >= 57 && n <= 92 {
                fr3.push(heavy_aa[i]);
            }
        }
        println!("heavy fr3 = {}, len = {}", strme(&fr3), fr3.len());

        let mut cdr1 = Vec::<u8>::new();
        for i in 0..light_no.len() {
            let mut s = light_no[i].as_bytes();
            for j in 0..s.len() {
                if s[j] > b'9' {
                    s = &s[0..j];
                    break;
                }
            }
            let n = strme(&s).force_usize();
            if n >= 24 && n <= 34 {
                cdr1.push(light_aa[i]);
            }
        }
        println!("light cdr1 = {}, len = {}", strme(&cdr1), cdr1.len());

        let mut cdr2 = Vec::<u8>::new();
        for i in 0..light_no.len() {
            let mut s = light_no[i].as_bytes();
            for j in 0..s.len() {
                if s[j] > b'9' {
                    s = &s[0..j];
                    break;
                }
            }
            let n = strme(&s).force_usize();
            if n >= 50 && n <= 56 {
                cdr2.push(light_aa[i]);
            }
        }
        println!("light cdr2 = {}, len = {}", strme(&cdr2), cdr2.len());

        let mut cdr1 = Vec::<u8>::new();
        for i in 0..heavy_no.len() {
            let mut s = heavy_no[i].as_bytes();
            for j in 0..s.len() {
                if s[j] > b'9' {
                    s = &s[0..j];
                    break;
                }
            }
            let n = strme(&s).force_usize();
            if n >= 26 && n <= 32 {
                cdr1.push(heavy_aa[i]);
            }
        }
        println!("heavy cdr1 = {}, len = {}", strme(&cdr1), cdr1.len());

        let mut cdr2 = Vec::<u8>::new();
        for i in 0..heavy_no.len() {
            let mut s = heavy_no[i].as_bytes();
            for j in 0..s.len() {
                if s[j] > b'9' {
                    s = &s[0..j];
                    break;
                }
            }
            let n = strme(&s).force_usize();
            if n >= 52 && n <= 56 {
                cdr2.push(heavy_aa[i]);
            }
        }
        println!("heavy cdr2 = {}, len = {}", strme(&cdr2), cdr2.len());
    }
}
