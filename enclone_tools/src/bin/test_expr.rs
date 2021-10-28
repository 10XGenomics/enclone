// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Experimental code for working with evalexpr.

use evalexpr::build_operator_tree;
use evalexpr::Function;
use evalexpr::Value;
use evalexpr::{ContextWithMutableFunctions, ContextWithMutableVariables, HashMapContext};
use statrs::distribution::ContinuousCDF;
use string_utils::*;

// ================================================================================================

// Convert a function having one of these forms:
// - fn f(x: f64) -> f64;
// - fn f(x: f64, y: f64) -> f64;
// - fn f(x: f64, y: f64, z: f64) -> f64;
// into an evalexpr::Function.
//
// This could be extended to work for zero variables or four/more.

#[macro_export]
macro_rules! evalexpr_fn1 {
    ($f:expr) => {
        Function::new(|t| {
            if t.is_number() {
                let t = t.as_number().unwrap();
                Ok(Value::from($f(t)))
            } else {
                Ok(Value::from(""))
            }
        })
    };
}

#[macro_export]
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

#[macro_export]
macro_rules! evalexpr_fn3 {
    ($f:expr) => {
        Function::new(|t| {
            if t.is_tuple() {
                let t = t.as_tuple().unwrap();
                if t.len() == 3 {
                    let x = &t[0];
                    let y = &t[1];
                    let z = &t[2];
                    if x.is_number() && y.is_number() && z.is_number() {
                        let x = x.as_number().unwrap();
                        let y = y.as_number().unwrap();
                        let z = z.as_number().unwrap();
                        return Ok(Value::from($f(x, y, z)));
                    }
                }
            }
            Ok(Value::from(""))
        })
    };
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

    // Define the beta cdf function.
    //
    // Requirements (not tested, but out of range should return *some* value):
    // - 0 <= x <= 1
    // - a > 0
    // - b > 0.
    
    fn beta_cdf(x: f64, a: f64, b: f64) -> f64 {
        let n = statrs::distribution::Beta::new(a, b).unwrap();
        n.cdf(x)
    }
    c.set_function("beta_cdf".to_string(), evalexpr_fn3![beta_cdf])
        .unwrap();

    // Done.

    c
}

// ================================================================================================

fn main() {
    // Define an expression.

    // let expr = "prod(x, y)";
    let expr = "square(x)";

    // Compile the expression.

    let compiled = build_operator_tree(&expr); // creates a Node
    if compiled.is_err() {
        println!("failed");
    } else {
        println!("ok");
    }
    let compiled = compiled.unwrap();

    // Create variables and values.

    let vars = vec!["x".to_string(), "y".to_string()];
    let vals = vec!["1.5".to_string(), "100".to_string()];

    // Set up compute context.

    let c = define_evalexpr_context(&vars, &vals);

    // Determine the value of the expression in the compute context.

    let res = compiled.eval_with_context(&c);
    if res.is_err() {
        println!("evaluation failed");
    }
    let res = res.unwrap();
    println!("value is {}", res);
}
