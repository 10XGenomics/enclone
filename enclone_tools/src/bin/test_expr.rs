// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Experimental code for working with evalexpr.

use evalexpr::build_operator_tree;
use expr_tools::*;
use itertools::Itertools;

// ================================================================================================

fn main() {
    // Define an expression.

    // let expr = "prod(x, y)";
    // let expr = "square(x)";
    let expr = "square(x.µ©¶Ƨ)";

    // Compile the expression.

    let compiled = build_operator_tree(&expr); // creates a Node
    if compiled.is_err() {
        println!("failed");
    } else {
        println!("ok");
    }
    let compiled = compiled.unwrap();
    let v = vars_of_node(&compiled);
    println!("variables in node = {}", v.iter().format(","));

    if true {
        std::process::exit(0);
    }

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
