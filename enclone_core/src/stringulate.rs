// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Representation of Vec<String> objects as Strings.  We do this under that assumption that the
// strings in the vector do not contain double escapes, which seems like a reasonable assumption
// in typical circumstances.  Then the representation of v as a string is
// (double escape) v[0] (double escape) v[1] (double escape) ... (double escape).
//
// For this to be reversible, one has to know the length of v.  In fact, the way we use this
// representation is as follows.  Let have a data structure, say Widget, and we convert it to
// a Vec<String>, in such a way that the vector starts with
// "Widget", number of entries in the string, ... .
//
// Second, in practice what we want to represent are vectors of heterogeneous objects, from a
// list of types having functions to_string and from_string, using the above design.  In that
// case we do this:
// 1. Convert each of the objects to strings, except for strings, which we leave intact.
// 2. Concatenate these.

use itertools::Itertools;
use string_utils::*;

const DOUBLE: &str = "";

pub fn flatten_vec_string(v: &[String]) -> String {
    format!("{}{}{}", DOUBLE, v.iter().format(&DOUBLE), DOUBLE)
}

pub fn unflatten_string(s: &str) -> Vec<String> {
    let mut chars = Vec::new();
    for c in s.chars() {
        chars.push(c);
    }
    let mut mid = String::new();
    for i in 2..chars.len() - 2 {
        mid.push(chars[i]);
    }
    mid.split(&DOUBLE).map(str::to_owned).collect()
}

pub struct HetString {
    pub name: String,
    pub content: String,
}

pub fn unpack_to_het_string(s: &str) -> Vec<HetString> {
    let mut v = Vec::<HetString>::new();
    let fields: Vec<String> = s.split(&DOUBLE).map(str::to_owned).collect();
    let mut i = 0;
    while i < fields.len() {
        if fields[i].len() > 0 {
            v.push(HetString {
                name: "String".to_string(),
                content: fields[i].to_string(),
            });
        }
        if i + 2 < fields.len() {
            let n = fields[i + 2].force_usize();
            v.push(HetString {
                name: fields[i + 1].to_string(),
                content: flatten_vec_string(&fields[i + 1..i + 1 + n]),
            });
            i += n + 1;
        } else {
            break;
        }
    }
    v
}

// Specific implementations, should split off if significantly more are added.

pub struct FeatureBarcodeAlluvialTable {
    pub id: String,
    pub display_text: String,
    pub spreadsheet_text: String,
}

pub struct FeatureBarcodeAlluvialTableSet {
    pub s: Vec<FeatureBarcodeAlluvialTable>,
}

impl FeatureBarcodeAlluvialTableSet {
    pub fn to_string(&self) -> String {
        let mut v = Vec::<String>::new();
        v.push("FeatureBarcodeAlluvialTableSet".to_string());
        v.push((3 * self.s.len() + 2).to_string());
        for i in 0..self.s.len() {
            v.push(self.s[i].id.clone());
            v.push(self.s[i].display_text.clone());
            v.push(self.s[i].spreadsheet_text.clone());
        }
        flatten_vec_string(&v)
    }
    pub fn from_string(x: &str) -> Self {
        let v = unflatten_string(&x);
        let n = v[1].force_usize() / 3;
        let mut s = Vec::new();
        for i in 0..n {
            s.push(FeatureBarcodeAlluvialTable {
                id: v[2 + 3 * i].clone(),
                display_text: v[2 + 3 * i + 1].clone(),
                spreadsheet_text: v[2 + 3 * i + 2].clone(),
            });
        }
        FeatureBarcodeAlluvialTableSet { s: s }
    }
}
