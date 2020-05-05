// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Utility for inserting html files.

use io_utils::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use string_utils::*;

pub fn insert_html(in_file: &str, out_file: &str) {
    let f = open_for_read![&in_file];
    let mut g = open_for_write_new![&out_file];
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with("#include ") {
            let h = open_for_read![&format!("../{}", s.after("#include "))];
            let mut started = false;
            for line in h.lines() {
                let t = line.unwrap();
                if t == "<body>" {
                    started = true;
                } else if t == "</body>" {
                    break;
                } else if started {
                    fwriteln!(g, "{}", t);
                }
            }
        } else {
            fwriteln!(g, "{}", s);
        }
        if s == "</head>" {
            fwriteln!(
                g,
                "
                \n
                <! ––\n
                💩 💩 💩 🔴 🔨 🔨 🔨 🔨 🔨 🔨 🔴 💩 💩 💩\n
                PUT DOWN YOUR HAMMER.
                THIS IS AN AUTO-GENERATED FILE.  PLEASE DO NOT EDIT IT.
                THANK YOU FOR YOUR COOPERATION,\n
                SINCERELY,
                THE BENEVOLENT OVERLORDS\n
                💩 💩 💩 🔴 🔨 🔨 🔨 🔨 🔨 🔨 🔴 💩 💩 💩\n
                ––>"
            );
        }
    }
}
