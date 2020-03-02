// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

use crate::print_utils1::*;
use ansi_escape::*;
use string_utils::*;
use tables::*;

pub fn stringify(rows: Vec<Vec<&str>>) -> Vec<Vec<String>> {
    let mut r = Vec::<Vec<String>>::new();
    for i in 0..rows.len() {
        let mut x = Vec::<String>::new();
        for j in 0..rows[i].len() {
            x.push(rows[i][j].to_string());
        }
        r.push(x);
    }
    r
}

pub fn print_enclone(plain: bool) {
    if plain {
        print!("enclone");
    } else {
        let mut log = Vec::<u8>::new();
        print_color(3, &mut log);
        log.push(b'e');
        emit_end_escape(&mut log);
        print_color(1, &mut log);
        log.push(b'n');
        emit_end_escape(&mut log);
        print_color(2, &mut log);
        log.push(b'c');
        emit_end_escape(&mut log);
        print_color(0, &mut log);
        log.push(b'l');
        emit_end_escape(&mut log);
        print_color(4, &mut log);
        log.push(b'o');
        emit_end_escape(&mut log);
        print_color(5, &mut log);
        log.push(b'n');
        emit_end_escape(&mut log);
        print_color(1, &mut log);
        log.push(b'e');
        emit_end_escape(&mut log);
        print!("{}", strme(&log));
    }
}

// This encodes the color codes for each possible codon of a given amino acid
// that could be found in a BCR or TCR sequence.

pub fn colored_codon_table(plainx: bool) -> String {
    let plain = b"\
        Alanine        A  GCT GCC GCA GCG\n\
        Arginine       R  CGT CGC CGA CGG AGA AGG\n\
        Asparagine     N  AAT AAC\n\
        Aspartic Acid  D  GAT GAC\n\
        Cysteine       C  TGT TGC\n\
        Glutamine      Q  CAA CAG\n\
        Glutamic Acid  E  GAA GAG\n\
        Glycine        G  GGT GGC GGA GGG\n\
        Histidine      H  CAT CAC\n\
        Isoleucine     I  ATT ATC ATA\n\
        Leucine        L  TTA TTG CTT CTC CTA CTG\n\
        Lysine         K  AAA AAG\n\
        Methionine     M  ATG\n\
        Phenylalanine  F  TTT TTC\n\
        Proline        P  CCT CCC CCA CCG\n\
        Serine         S  TCT TCC TCA TCG AGT AGC\n\
        Threonine      T  ACT ACC ACA ACG\n\
        Tryptophan     W  TGG\n\
        Tyrosine       Y  TAT TAC\n\
        Valine         V  GTT GTC GTA GTG";
    let mut colored = Vec::<u8>::new();
    let mut p = 0;
    while p < plain.len() {
        if (plain[p] as char).is_uppercase() && (plain[p + 1] as char).is_uppercase() {
            let mut log = Vec::<u8>::new();
            if !plainx {
                emit_codon_color_escape(&plain[p..p + 3], &mut log);
            }
            for i in 0..3 {
                log.push(plain[p + i]);
            }
            if !plainx {
                emit_end_escape(&mut log);
            }
            colored.append(&mut log);
            p += 3;
        } else {
            colored.push(plain[p]);
            p += 1;
        }
    }
    stringme(&colored)
}

pub static mut PLAIN: bool = false;
pub static mut HELP_ALL: bool = false;

// Print a string, making the following conversions, the first three of which are governed
// by the state of PLAIN:
// • Change \bold{x} into a bolded string by issuing appropriate escape characters.
// • Change \red{x} into a red string by issuing appropriate escape characters.
// • Change \boldred{x} into a bold red string by issuing appropriate escape characters.
// • Fold at 99 characters.

pub fn print(x: &str) {
    print!("{}", print_to(x));
}

pub fn print_to(x: &str) -> String {
    let mut y = Vec::<char>::new();
    for c in x.chars() {
        y.push(c);
    }
    let mut s = String::new();
    let mut i = 0;
    while i < y.len() {
        if y[i..].starts_with(&['\\', 'b', 'o', 'l', 'd', '{']) {
            let mut j = i + 6;
            while j < y.len() {
                if y[j] == '}' {
                    break;
                }
                j += 1;
            }
            if j < y.len() {
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_bold_escape(&mut log);
                    }
                }
                s += &strme(&log);
                for k in i + 6..j {
                    s.push(y[k]);
                }
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_end_escape(&mut log);
                    }
                }
                s += &strme(&log);
                i = j + 1;
            }
        } else if y[i..].starts_with(&['\\', 'r', 'e', 'd', '{']) {
            let mut j = i + 5;
            while j < y.len() {
                if y[j] == '}' {
                    break;
                }
                j += 1;
            }
            if j < y.len() {
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_red_escape(&mut log);
                    }
                }
                s += &strme(&log);
                for k in i + 5..j {
                    s.push(y[k]);
                }
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_end_escape(&mut log);
                    }
                }
                s += &strme(&log);
                i = j + 1;
            } else {
                i += 1;
            }
        } else if y[i..].starts_with(&['\\', 'b', 'l', 'u', 'e', '{']) {
            let mut j = i + 6;
            while j < y.len() {
                if y[j] == '}' {
                    break;
                }
                j += 1;
            }
            if j < y.len() {
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_blue_escape(&mut log);
                    }
                }
                s += &strme(&log);
                for k in i + 6..j {
                    s.push(y[k]);
                }
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_end_escape(&mut log);
                    }
                }
                s += &strme(&log);
                i = j + 1;
            } else {
                i += 1;
            }
        } else if y[i..].starts_with(&['\\', 'g', 'r', 'e', 'e', 'n', '{']) {
            let mut j = i + 7;
            while j < y.len() {
                if y[j] == '}' {
                    break;
                }
                j += 1;
            }
            if j < y.len() {
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_green_escape(&mut log);
                    }
                }
                s += &strme(&log);
                for k in i + 7..j {
                    s.push(y[k]);
                }
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_end_escape(&mut log);
                    }
                }
                s += &strme(&log);
                i = j + 1;
            } else {
                i += 1;
            }
        } else if y[i..].starts_with(&['\\', 'b', 'o', 'l', 'd', 'r', 'e', 'd', '{']) {
            let mut j = i + 9;
            while j < y.len() {
                if y[j] == '}' {
                    break;
                }
                j += 1;
            }
            if j < y.len() {
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_bold_escape(&mut log);
                        emit_red_escape(&mut log);
                    }
                }
                s += &strme(&log);
                for k in i + 9..j {
                    s.push(y[k]);
                }
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_end_escape(&mut log);
                    }
                }
                s += &strme(&log);
                i = j + 1;
            } else {
                i += 1;
            }
        } else if y[i..].starts_with(&['\\', 'b', 'o', 'l', 'd', 'b', 'l', 'u', 'e', '{']) {
            let mut j = i + 10;
            while j < y.len() {
                if y[j] == '}' {
                    break;
                }
                j += 1;
            }
            if j < y.len() {
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_bold_escape(&mut log);
                        emit_blue_escape(&mut log);
                    }
                }
                s += &strme(&log);
                for k in i + 10..j {
                    s.push(y[k]);
                }
                let mut log = Vec::<u8>::new();
                unsafe {
                    if !PLAIN {
                        emit_end_escape(&mut log);
                    }
                }
                s += &strme(&log);
                i = j + 1;
            } else {
                i += 1;
            }
        } else {
            s.push(y[i]);
            i += 1;
        }
    }
    let mut x = Vec::<char>::new();
    for c in s.chars() {
        x.push(c);
    }
    let mut printed = 0;
    let mut escaped = false;
    let mut y = Vec::<char>::new();
    let mut i = 0;
    while i < x.len() {
        if x[i] == '' {
            escaped = true;
        }
        if escaped {
            if x[i] == 'm' {
                escaped = false;
            }
            y.push(x[i]);
            i += 1;
            continue;
        }
        if x[i] == ' ' {
            let mut j = i + 1;
            while j < x.len() {
                if x[j] == ' ' || x[j] == '\n' || x[j] == '' {
                    break;
                }
                j += 1;
            }
            if printed + j - i >= 100 - 1 {
                y.push('\n');
                printed = 0;
                i += 1;
                continue;
            }
        }
        y.push(x[i]);
        printed += 1;
        if x[i] == '\n' {
            printed = 0;
        }
        i += 1;
    }
    let mut ans = String::new();
    for i in 0..y.len() {
        ans.push(y[i]);
    }
    ans
}

pub fn print_with_box(x: &str) {
    let y = print_to(x);
    let mut rows = Vec::<Vec<String>>::new();
    let lines = y.split('\n').collect::<Vec<&str>>();
    for z in lines {
        rows.push(vec![z.to_string()]);
    }
    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows, 2, &b"l".to_vec(), false);
    println!("{}", log);
}

pub fn banner(x: &str, plain: bool) {
    let mut log = Vec::<u8>::new();
    if !plain {
        emit_blue_escape(&mut log);
    }
    print!("{}", strme(&log));
    for _ in 1..100 {
        print!("▓");
    }
    print!("{}", strme(&log));
    println!("\nenclone help {}", x);
    let mut log = Vec::<u8>::new();
    if !plain {
        emit_blue_escape(&mut log);
    }
    print!("{}", strme(&log));
    for _ in 1..100 {
        print!("▓");
    }
    let mut log = Vec::<u8>::new();
    if !plain {
        emit_end_escape(&mut log);
    }
    println!("{}", strme(&log));
}

pub fn print_tab2(rows: &Vec<Vec<String>>) {
    let mut log = String::new();
    print_tabular_vbox(&mut log, &rows, 2, &b"l|l".to_vec(), false);
    print!("{}", log);
}
