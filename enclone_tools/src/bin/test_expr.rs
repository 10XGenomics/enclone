// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.


use evalexpr::build_operator_tree;
use evalexpr::{ContextWithMutableFunctions, ContextWithMutableVariables, HashMapContext};
use evalexpr::Function;
use evalexpr::Value;

fn main() {

    // let expr = "3 + 5";
    // let expr = "3+5";
    // let expr = "f(arglin-9,7,pot(moo)) + 5";
    // let expr = "3 > 5";
    // let expr = "3>5";
    // let expr = "f(10) - g(2)";
    // let expr = "(lambdax + 5) / 0";
    // let expr = "square(5)";
    // let expr = "square(5.0)";
    let expr = "prod(2, 3)";

    let compiled = build_operator_tree(&expr); // creates a Node
    if compiled.is_err() {
        println!("failed");
    } else {
        println!("ok");
    }
    let compiled = compiled.unwrap();

    let mut c = HashMapContext::new();
    let _ = c.set_value("lambdax".to_string(), evalexpr::Value::Float(1.234));




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
                Ok(Value::from("undefined"))
            }
        })
    );




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
            Ok(Value::from("undefined"))
        })
    );




    let res = compiled.eval_with_context(&c);
    if res.is_err() {
        println!("evaluation failed");
    }
    let res = res.unwrap();
    println!("value is {}", res);
}

/*

• use == for equality, and not =
• put string values in single quotes 
• put the entire expression in double quotes

enclone BCR=123085 BC=f KEEP_CELL_IF="nice == 'true'" 
enclone BCR=123085 BC=f KEEP_CELL_IF="nice == 'true' && rank <= 5"

*/
