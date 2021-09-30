// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod var;

use string_utils::*;
use vector_utils::*;

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
