// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Variable specification.  Fields are currently Strings, but could be given more structure.

use string_utils::{stringme, TextUtils};

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

// Parse variables, and exit if requirements are not satisfied.

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
                if !this_group[n - 1].ends_with("\n") {
                    this_group[n - 1] += "\n";
                }
                this_group[n - 1] += &mut format!(" {}", line.after(INDENT));
                if FIELDS[fc - 1] == "code" {
                    this_group[n - 1] += "\n";
                }
            }
        }
    }

    // Form variables.

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

    // Test for duplicated entries.

    for i in 1..vars.len() {
        if vars[i].name == vars[i - 1].name {
            eprintln!(
                "\nThe variable name {} appears more than once.\n",
                vars[i].name
            );
            std::process::exit(1);
        }
    }

    // Test upper-case rule.

    let classes = ["BC", "DATASET", "FEATURE", "INFO", "NAME", "REGA", "VARDEF"];
    for i in 0..vars.len() {
        let n = &vars[i].name;
        let mut chars = Vec::<char>::new();
        for c in n.chars() {
            chars.push(c);
        }
        let mut j = 0;
        while j < chars.len() {
            if chars[j] < 'A' || chars[j] > 'Z' {
                j += 1;
            } else {
                let mut k = j + 1;
                while k < chars.len() {
                    if chars[k] < 'A' || chars[k] > 'Z' {
                        break;
                    }
                    k += 1;
                }
                let mut s = String::new();
                for l in j..k {
                    s.push(chars[l]);
                }
                if !classes.contains(&s.as_str()) {
                    eprintln!("\nFound illegal class {} in variable name {}.\n", s, n);
                    std::process::exit(1);
                }
                j = k;
            }
        }
    }

    // Return.

    vars
}
