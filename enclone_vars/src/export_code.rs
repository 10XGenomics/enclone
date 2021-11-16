// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Read the vars file and export code.  This is a partial implementation.
// Output is {(filename, contents)}.
//
// This writes a temporary file.

use crate::var::parse_variables;
use io_utils::*;
use itertools::Itertools;
use std::io::{BufWriter, Write};
use std::process::Command;
use string_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Emit code that tests for a given variable, allowing for up to three bracket expressions
// in the variable.

fn emit_code_to_test_for_var<W: Write>(var: &str, f: &mut BufWriter<W>) {
    let nranges = var.matches('{').count();
    assert_eq!(nranges, var.matches('}').count());
    assert!(nranges <= 3);
    if nranges == 0 {
        fwriteln!(f, r###"}} else if var == "{}" {{"###, var);
    } else if nranges == 1 {
        let begin = var.before("{");
        let end = var.after("}").to_string();
        let low = var.after("{").before("..");
        let high = var.after("{").between("..", "}");
        let mut conditions = Vec::<String>::new();
        conditions.push(format!(r###"var.starts_with("{}")"###, begin));
        conditions.push(format!(r###"var.ends_with("{}")"###, end));
        conditions.push(format!(
            r###"var.between2("{}", "{}").parse::<i64>().is_ok()"###,
            begin, end,
        ));
        if low.len() > 0 {
            conditions.push(format!(
                r###"var.between2("{}", "{}").force_i64() >= {}"###,
                begin, end, low,
            ));
        }
        if high.len() > 0 {
            conditions.push(format!(
                r###"var.between2("{}", "{}").force_i64() <= {}"###,
                begin, end, high,
            ));
        }
        fwriteln!(f, "}} else if {} {{ ", conditions.iter().format(" && "));
        fwriteln!(
            f,
            r###"let arg1 = var.between2("{}", "{}").force_i64();"###,
            begin,
            end,
        );
    } else if nranges == 2 {
        // This code has not been exercised.
        let begin = var.before("{");
        let middle = var.between("}", "{");
        let end = var.rev_after("}").to_string();
        let low1 = var.after("{").before("..");
        let high1 = var.after("{").between("..", "}");
        let low2 = var.rev_after("{").before("..");
        let high2 = var.rev_after("{").between("..", "}");
        let mut conditions = Vec::<String>::new();
        conditions.push(format!(r###"var.starts_with("{}")"###, begin));
        conditions.push(format!(
            r###"var.after("{}").contains("{}")"###,
            begin, middle,
        ));
        conditions.push(format!(
            r###"var.after("{}").after("{}").ends_with("{}")"###,
            begin, middle, end,
        ));
        conditions.push(format!(
            r###"var.between2("{}", "{}").parse::<i64>().is_ok()"###,
            begin, middle,
        ));
        if low1.len() > 0 {
            conditions.push(format!(
                r###"var.between2("{}", "{}").force_i64() >= {}"###,
                begin, middle, low1,
            ));
        }
        if high1.len() > 0 {
            conditions.push(format!(
                r###"var.between2("{}", "{}").force_i64() <= {}"###,
                begin, middle, high1,
            ));
        }
        conditions.push(format!(
            r###"var.after("{}").between2("{}", "{}").parse::<i64>().is_ok()"###,
            begin, middle, end,
        ));
        if low2.len() > 0 {
            conditions.push(format!(
                r###"var.after("{}").between2("{}", "{}").force_i64() >= {}"###,
                begin, middle, end, low2,
            ));
        }
        if high2.len() > 0 {
            conditions.push(format!(
                r###"var.after("{}").between2("{}", "{}").force_i64() <= {}"###,
                begin, middle, end, high2,
            ));
        }
        fwriteln!(f, "}} else if {} {{ ", conditions.iter().format(" && "));
        fwriteln!(
            f,
            r###"let arg1 = var.between2("{}", "{}").force_i64();"###,
            begin,
            middle,
        );
        fwriteln!(
            f,
            r###"let arg2 = var.after("{}"),between2("{}", "{}").force_i64();"###,
            begin,
            middle,
            end,
        );
    } else {
        let begin = var.before("{");
        let mid1 = var.between("}", "{");
        let mid2 = var.after("}").between("}", "{");
        let end = var.rev_after("}").to_string();
        let low1 = var.after("{").before("..");
        let high1 = var.after("{").between("..", "}");
        let low2 = var.after("{").after("{").before("..");
        let high2 = var.after("{").after("{").between("..", "}");
        let low3 = var.rev_after("{").before("..");
        let high3 = var.rev_after("{").between("..", "}");
        let mut conditions = Vec::<String>::new();
        conditions.push(format!(r###"var.starts_with("{}")"###, begin));
        conditions.push(format!(
            r###"var.after("{}").contains("{}")"###,
            begin, mid1,
        ));
        conditions.push(format!(
            r###"var.after("{}").after("{}").contains("{}")"###,
            begin, mid1, mid2,
        ));
        conditions.push(format!(
            r###"var.after("{}").after("{}").after("{}").ends_with("{}")"###,
            begin, mid1, mid2, end,
        ));
        conditions.push(format!(
            r###"var.between("{}", "{}").parse::<i64>().is_ok()"###,
            begin, mid1,
        ));
        if low1.len() > 0 {
            conditions.push(format!(
                r###"var.between("{}", "{}").force_i64() >= {}"###,
                begin, mid1, low1,
            ));
        }
        if high1.len() > 0 {
            conditions.push(format!(
                r###"var.between("{}", "{}").force_i64() <= {}"###,
                begin, mid1, high1,
            ));
        }
        if low2.len() > 0 {
            conditions.push(format!(
                r###"var.after("{}").between("{}", "{}").force_i64() >= {}"###,
                begin, mid1, mid2, low2,
            ));
        }
        if high2.len() > 0 {
            conditions.push(format!(
                r###"var.after("{}").between("{}", "{}").force_i64() <= {}"###,
                begin, mid1, mid2, high2,
            ));
        }
        conditions.push(format!(
            r###"var.after("{}").after("{}").between("{}", "{}").parse::<i64>().is_ok()"###,
            begin, mid1, mid2, end,
        ));
        if low3.len() > 0 {
            conditions.push(format!(
                r###"var.after("{}").after("{}").between("{}", "{}").force_i64() >= {}"###,
                begin, mid1, mid2, end, low3,
            ));
        }
        if high3.len() > 0 {
            conditions.push(format!(
                r###"var.after("{}").after("{}").between("{}", "{}").force_i64() <= {}"###,
                begin, mid1, mid2, end, high3,
            ));
        }
        fwriteln!(f, "}} else if {} {{ ", conditions.iter().format(" && "));
        fwriteln!(
            f,
            r###"let arg1 = var.between("{}", "{}").force_i64();"###,
            begin,
            mid1,
        );
        fwriteln!(
            f,
            r###"let arg2 = var.after("{}").between("{}", "{}").force_i64();"###,
            begin,
            mid1,
            mid2,
        );
        fwriteln!(
            f,
            r###"let arg3 = var.after("{}").after("{}").between("{}", "{}").force_i64();"###,
            begin,
            mid1,
            mid2,
            end,
        );
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Export code.

pub fn export_code(level: usize) -> Vec<(String, String)> {
    // Do two passes, one for cvars and one for lvars.  (second pass not implemented yet)

    let mut outs = Vec::<(String, String)>::new();
    for zpass in 1..=2 {
        if zpass == 2 {
            continue;
        }

        // Define code start/stop for var_vdj.

        let cvar_vdj_start = r###"
    
            // Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
            // This file is auto-generated by the crate enclone_vars, please do not edit.
    
            use amino::*;
            use crate::print_utils1::*;
            use crate::print_utils3::*;
            use enclone_core::align_to_vdj_ref::*;
            use enclone_core::defs::*;
            use enclone_core::median::*;
            use enclone_core::opt_d::*;
            use enclone_proto::types::*;
            use itertools::Itertools;
            use stats_utils::*;
            use std::cmp::min;
            use std::collections::HashMap;
            use string_utils::*;
            use vdj_ann::refx::RefData;
            use vector_utils::*;
    
            pub fn proc_cvar_auto(
                j: usize,
                pass: usize,
                var: &String,
                ex: &ExactClonotype,
                exacts: &Vec<usize>,
                exact_clonotypes: &Vec<ExactClonotype>,
                mid: usize,
                col: usize,
                u: usize,
                rsi: &ColInfo,
                refdata: &RefData,
                dref: &Vec<DonorReferenceItem>,
                ctl: &EncloneControl,
                extra_args: &Vec<String>,
                pcols_sort: &Vec<String>,
                cx: &mut Vec<Vec<String>>,
                varmat: &Vec<Vec<Vec<u8>>>,
                out_data: &mut Vec<HashMap<String, String>>,
                stats: &mut Vec<(String, Vec<String>)>,
            ) -> Result<bool, String> {
    
                let cvars = &ctl.clono_print_opt.cvars;
    
                // Set up macros.
    
                macro_rules! speakc {
                    ($u:expr, $col:expr, $var:expr, $val:expr) => {
                        if pass == 2
                            && ((ctl.parseable_opt.pout.len() > 0
                                && (ctl.parseable_opt.pchains == "max"
                                    || col < ctl.parseable_opt.pchains.force_usize()))
                                || extra_args.len() > 0)
                        {
                            let mut v = $var.clone();
                            v = v.replace("_Σ", "_sum");
                            v = v.replace("_μ", "_mean");
            
                            // Strip escape character sequences from val.  Can happen in notes, 
                            // maybe other places.
            
                            let mut val_clean = String::new();
                            let mut chars = Vec::<char>::new();
                            let valx = format!("{}", $val);
                            for c in valx.chars() {
                                chars.push(c);
                            }
                            let mut escaped = false;
                            for l in 0..chars.len() {
                                if chars[l] == '' {
                                    escaped = true;
                                }
                                if escaped {
                                    if chars[l] == 'm' {
                                        escaped = false;
                                    }
                                    continue;
                                }
                                val_clean.push(chars[l]);
                            }
            
                            // Proceed.
            
                            let varc = format!("{}{}", v, $col + 1);
                            if pcols_sort.is_empty()
                                || bin_member(&pcols_sort, &varc)
                                || bin_member(&extra_args, &varc)
                            {
                                out_data[$u].insert(varc, val_clean);
                            }
                        }
                    };
                }
        
                // Test variable.
    
                let val =
                if false {
                    (String::new(), Vec::<String>::new())
    
            "###;

        let cvar_vdj_stop = r###"
    
                } else {
                    ("$UNDEFINED".to_string(), Vec::<String>::new())
                };
                if val.0 == "$UNDEFINED" {
                    return Ok(false);
                } else {
                    let (exact, cell) = &val;
                    let varc = format!("{}{}", var, col + 1);
                    if exact.len() > 0 {
                        if j < rsi.cvars[col].len() && cvars.contains(&var) {
                            cx[col][j] = exact.clone();
                        }
                        speakc!(u, col, var, exact);
                        if val.1.is_empty() {
                            stats.push((varc, vec![exact.to_string(); ex.ncells()]));
                        } else {
                            stats.push((varc, cell.to_vec()));
                        }
                    } else if cell.len() > 0 {
                        if pass == 2
                            && ((ctl.parseable_opt.pchains == "max"
                                || col < ctl.parseable_opt.pchains.force_usize())
                                || !extra_args.is_empty())
                        {
                            if pcols_sort.is_empty() || bin_member(pcols_sort, &varc) {
                                let vals = format!("{}", cell.iter().format(&POUT_SEP));
                                out_data[u].insert(varc, vals);
                            }
                        }
                    }
                    return Ok(true);
                }
            }
    
            "###;

        // Build cvar auto file.

        let actual_out = "enclone_print/src/proc_cvar_auto.rs".to_string();
        let mut temp_out = "enclone_exec/testx/outputs/proc_cvar_auto.rs".to_string();
        let mut vars_loc = "enclone_vars/src/vars".to_string();
        if level == 1 {
            temp_out = format!("../{}", temp_out);
            vars_loc = format!("../{}", vars_loc);
        }
        {
            let mut f = open_for_write_new![&temp_out];
            fwrite!(f, "{}", cvar_vdj_start);
            let vars = std::fs::read_to_string(&vars_loc).unwrap();
            let vars = parse_variables(&vars);
            for v in vars.iter() {
                if v.inputs == "cvar_vdj" {
                    // Parse value return lines.

                    let mut exact = "String::new()".to_string();
                    let mut cell = "Vec::new()".to_string();
                    let mut code = v.code.clone();
                    let mut lines = Vec::<String>::new();
                    for line in code.lines() {
                        lines.push(line.to_string());
                    }
                    let n = lines.len();
                    if n > 0 {
                        let mut sub = 0;
                        for i in (0..lines.len()).rev() {
                            if lines[i].contains("exact: ") {
                                exact = lines[i].after("exact: ").to_string();
                                sub += 1;
                            } else if lines[i].contains("cell: ") {
                                cell = lines[i].after("cell: ").to_string();
                                sub += 1;
                            }
                        }
                        let mut code2 = String::new();
                        for i in 0..lines.len() - sub {
                            code2 += &mut format!("{}\n", lines[i]);
                        }
                        code = code2;
                    }
                    if v.level == "cell-exact" {
                        assert!(!exact.is_empty());
                        assert!(!cell.is_empty());
                    }
                    if v.level == "cell" {
                        assert!(!cell.is_empty());
                    }

                    // Proceed.

                    // RESTRICTION 1: only allow cvar_vdj
                    let mut upper = false;
                    let var = &v.name;
                    for c in var.chars() {
                        if c.is_ascii_uppercase() {
                            upper = true;
                        }
                    }
                    // RESTRICTION 2: don't allow upper case
                    if !upper {
                        let mut passes = 1;
                        if v.level == "cell-exact" {
                            passes = 2;
                        }
                        for pass in 1..=passes {
                            let mut var = var.clone();
                            if pass == 2 {
                                var += "_cell";
                            }
                            emit_code_to_test_for_var(&var, &mut f);
                            fwriteln!(f, "{}", code);
                            if pass == 1 {
                                fwriteln!(f, "({}, {})", exact, cell);
                            } else {
                                fwriteln!(f, "let _exact = {};", exact); // to circumvent warning
                                fwriteln!(f, "(String::new(), {})", cell);
                            }
                        }
                    }
                }
            }
            fwrite!(f, "{}", cvar_vdj_stop);
        }

        // Rustfmt.

        let new = Command::new("rustfmt")
            .arg(&temp_out)
            .output()
            .unwrap_or_else(|_| panic!("{}", "failed to execute rustfmt".to_string()));
        if new.status.code() != Some(0) {
            eprintln!("\nrustfmt failed\n");
            eprintln!(
                "You can observe the problem by typing rustfmt {}.\n",
                temp_out
            );
            std::process::exit(1);
        }
        let f = std::fs::read_to_string(&temp_out).unwrap();
        outs.push((actual_out, f));
    }
    outs
}
