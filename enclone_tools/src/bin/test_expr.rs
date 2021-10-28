// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Experimental code for working with evalexpr.

use evalexpr::build_operator_tree;
use evalexpr::Function;
use evalexpr::Value;
use evalexpr::{ContextWithMutableFunctions, ContextWithMutableVariables, HashMapContext};
use string_utils::*;

// ================================================================================================

// This does not compile, so using a macro below.

/*
pub fn evalexpr_fn2(f: fn(f64, f64) -> f64) -> evalexpr::Function {
    Function::new(|t| {
        if t.is_tuple() {
            let t = t.as_tuple().unwrap();
            if t.len() == 2 {
                let x = &t[0];
                let y = &t[1];
                if x.is_number() && y.is_number() {
                    let x = x.as_number().unwrap();
                    let y = y.as_number().unwrap();
                    return Ok(Value::from(f(x, y)));
                }
            }
        }
        Ok(Value::from(""))
    })
}
*/

macro_rules! evalexpr_fn2 {
    ($f:expr) => {
        Function::new(|t| {
            if t.is_tuple() {
                let t = t.as_tuple().unwrap();
                if t.len() == 2 {
                    let x = &t[0];
                    let y = &t[1];
                    if x.is_number() && y.is_number() {
                        let x = x.as_number().unwrap();
                        let y = y.as_number().unwrap();
                        return Ok(Value::from($f(x, y)));
                    }
                }
            }
            Ok(Value::from(""))
        })
    };
}

pub fn prod(x: f64, y: f64) -> f64 {
    x * y
}

// ================================================================================================

// Given a a list of variables and values for them, define an evalexpr::HashMapContext that
// includes these variables, as well as some convenient (but arbitrarily chosen) functions.
// This can then be used to evaluate an expression.
//
// Functions take as input zero or more f64 arguments, and return an f64.  The machine allows them
// to be called on arbitrary strings, but if the strings are not all f64, then the return value is
// null.

pub fn define_evalexpr_context(vars: &Vec<String>, vals: &Vec<String>) -> evalexpr::HashMapContext {
    assert_eq!(vars.len(), vals.len());
    let mut c = HashMapContext::new();

    // Define the variable values.

    for i in 0..vars.len() {
        if vals[i].parse::<f64>().is_ok() {
            c.set_value(vars[i].clone(), evalexpr::Value::Float(vals[i].force_f64()))
                .unwrap();
        } else {
            c.set_value(vars[i].clone(), evalexpr::Value::String(vals[i].clone()))
                .unwrap();
        }
    }

    let _ = c.set_function("prod".to_string(), evalexpr_fn2![prod]);

    // Define a function.

    /*
    fn prod(x: f64, y: f64) -> f64 {
        x * y
    }
    let _ = c.set_function(
        "prod".to_string(),
        Function::new(|t| {
            if t.is_tuple() {
                let t = t.as_tuple().unwrap();
                if t.len() == 2 {
                    let x = &t[0];
                    let y = &t[1];
                    if x.is_number() && y.is_number() {
                        let x = x.as_number().unwrap();
                        let y = y.as_number().unwrap();
                        return Ok(Value::from(prod(x, y)));
                    }
                }
            }
            Ok(Value::from(""))
        }),
    );
    */

    // Define a function.

    fn square(x: f64) -> f64 {
        x * x
    }
    let _ = c.set_function(
        "square".to_string(),
        Function::new(|x| {
            if x.is_number() {
                let x = x.as_number().unwrap();
                Ok(Value::from(square(x)))
            } else {
                Ok(Value::from(""))
            }
        }),
    );

    // Done.

    c
}

// ================================================================================================

fn main() {
    let expr = "prod(x, y)";
    let compiled = build_operator_tree(&expr); // creates a Node
    if compiled.is_err() {
        println!("failed");
    } else {
        println!("ok");
    }
    let compiled = compiled.unwrap();
    let vars = vec!["x".to_string(), "y".to_string()];
    let vals = vec!["1.5".to_string(), "100".to_string()];
    let c = define_evalexpr_context(&vars, &vals);
    let res = compiled.eval_with_context(&c);
    if res.is_err() {
        println!("evaluation failed");
    }
    let res = res.unwrap();
    println!("value is {}", res);
}
