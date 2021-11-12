// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod export_code;
pub mod var;

use string_utils::TextUtils;
use vector_utils::sort_sync2;

pub fn sort_vars(input: &str) -> String {
    let mut preamble = String::new();
    let mut vars = Vec::<String>::new();
    let mut in_vars = false;
    let div = "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\
        ━━━━━━━━━━━━━━━━━━━━━";
    let mut groups = Vec::<String>::new();
    let mut this_group = String::new();
    for line in input.lines() {
        if !in_vars && line != div {
            preamble += &mut format!("{}\n", line);
        } else if !in_vars && line == div {
            in_vars = true;
            this_group = format!("{}\n", div);
        } else if line == div {
            groups.push(this_group.clone());
            this_group = format!("{}\n", div);
        } else {
            this_group += &mut format!("{}\n", line);
        }
    }
    for i in 0..groups.len() {
        let name = groups[i].between("name:     ", "\n");
        vars.push(name.to_string());
    }
    sort_sync2(&mut vars, &mut groups);
    let mut out = preamble;
    for i in 0..groups.len() {
        out += &mut groups[i];
    }
    out += &mut format!("{}\n", div);
    out
}

// Functions to encode and decode arithmetic operators.  Because the symbols - and / appear in gene
// names, and because these symbols are also arithmetic operators, we need a system for hiding
// them.  The system is that we encode/decode all the standard arithmetic operators + - * /
// according to the following table, but only if they appear with characters on both sides, with
// neither character being a blank.
//
// NORMAL    ENCODED
// +         ©add©
// -         ©sub©
// *         ©mul©
// /         ©div©
//
// The funny symbol © is the copyright symbol.

pub fn encode_arith(x: &str) -> String {
    let mut m = Vec::<char>::new();
    for c in x.chars() {
        m.push(c);
    }
    let mut encoded = String::new();
    for i in 0..m.len() {
        let mut saved = false;
        if i >= 1 && i < m.len() - 1 && m[i - 1] != ' ' && m[i + 1] != ' ' {
            if m[i] == '+' {
                encoded += "©add©";
                saved = true;
            } else if m[i] == '-' {
                encoded += "©sub©";
                saved = true;
            } else if m[i] == '*' {
                encoded += "©mul©";
                saved = true;
            } else if m[i] == '/' {
                encoded += "©div©";
                saved = true;
            }
        }
        if !saved {
            encoded.push(m[i]);
        }
    }
    encoded
}

pub fn decode_arith(x: &str) -> String {
    let mut x = x.replace("©add©", "+");
    x = x.replace("©sub©", "-");
    x = x.replace("©mul©", "*");
    x = x.replace("©div©", "/");
    x
}
