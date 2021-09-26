// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Variable specification.  Fields are currently Strings, but could be given more structure.

use string_utils::*;

pub struct Variable {
    pub name: String,
    pub inputs: String,
    pub limits: String,
    pub class: String,
    pub level: String,
    pub val: String,
    pub doc: String,
    pub brief: String,
    pub page: String,
    pub avail: String,
    pub notes: String,
    pub code: String,
}

pub fn parse_variables(input: &str) -> Vec<Variable> {
    const FIELDS: [&str; 12] = [
        "name", "inputs", "limits", "class", "level", "val", "doc", "brief", "page", "avail",
        "notes", "code",
    ];
    const INDENT: &str = "          ";
    let mut in_vars = false;
    let div = "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\
        ━━━━━━━━━━━━━━━━━━━━━";
    let mut groups = Vec::<Vec<String>>::new();
    let mut this_group = Vec::<String>::new();
    let mut fc = 0;
    for (i, line) in input.lines().enumerate() {
        if !in_vars && line != div {
        } else if !in_vars && line == div {
            in_vars = true;
            this_group.clear();
            fc = 0;
        } else if line == div {
            if this_group.len() != FIELDS.len() {
                eprintln!("\nWrong number of fields before line {}.\n", i + 1);
                std::process::exit(1);
            }
            groups.push(this_group.clone());
            this_group.clear();
            fc = 0;
        } else {
            if fc > FIELDS.len() {
                eprintln!("\nToo many fields at line {}.\n", i + 1);
                std::process::exit(1);
            }
            let new = fc < FIELDS.len() && line.starts_with(&format!("{}:", FIELDS[fc]));
            if new {
                if line == format!("{}:", FIELDS[fc]) {
                    this_group.push(String::new());
                } else {
                    for k in FIELDS[fc].len() + 1..INDENT.len() {
                        if k >= line.as_bytes().len() || line.as_bytes()[k] != b' ' {
                            eprintln!(
                                "\nIllegal indentation or trailing blanks at line {}:\n{}\n",
                                i + 1,
                                line,
                            );
                            std::process::exit(1);
                        }
                    }
                    this_group.push(stringme(&line.as_bytes()[INDENT.len()..]));
                }
                fc += 1;
                if fc > FIELDS.len() {
                    eprintln!("\nToo many fields at line {}.\n", i + 1);
                    std::process::exit(1);
                }
            } else {
                if !line.starts_with(&INDENT) {
                    eprintln!(
                        "\nIllegal field or indentation rule violation at line {}:\n{}\n",
                        i + 1,
                        line,
                    );
                    std::process::exit(1);
                }
                let n = this_group.len();
                this_group[n - 1] += &mut format!(" {}", line.after(&INDENT));
            }
        }
    }
    let mut vars = Vec::<Variable>::new();
    for g in groups.iter() {
        vars.push(Variable {
            name: g[0].clone(),
            inputs: g[1].clone(),
            limits: g[2].clone(),
            class: g[3].clone(),
            level: g[4].clone(),
            val: g[5].clone(),
            doc: g[6].clone(),
            brief: g[7].clone(),
            page: g[8].clone(),
            avail: g[9].clone(),
            notes: g[10].clone(),
            code: g[11].clone(),
        });
    }
    vars
}
