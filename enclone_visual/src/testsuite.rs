// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::messages::*;
use crate::*;

const SUBMIT: fn(Result<(), std::string::String>) -> messages::Message =
    Message::SubmitButtonPressed as MsgFn;

const BACK: fn(Result<(), std::string::String>) -> messages::Message =
    Message::BackButtonPressed as MsgFn;

const FORWARD: fn(Result<(), std::string::String>) -> messages::Message =
    Message::ForwardButtonPressed as MsgFn;

#[rustfmt::skip]
pub const TESTS: [(&str, MsgFn, &str); 17] = [
    ("enclone woof", SUBMIT, ""), // this is a crash test
    ("#1", SUBMIT, "test1"),      // enclone BCR=123085 PLOT=gui MIN_CELLS=5
    ("", BACK, ""),
    ("", FORWARD, ""),
    ("10", SUBMIT, "test1b"),
    ("", BACK, ""),
    ("", FORWARD, "test1x"),
    ("enclone BCR=123085 PLOT=gui MIN_CELLS=5 G=12", SUBMIT, "test1c"),
    ("#2", SUBMIT, "test2"), // enclone BCR=123085 PLOT_BY_ISOTYPE=gui MIN_CELLS=5
    ("#3", SUBMIT, "test3"), // enclone BCR=123085 GEX=123217 PLOTXY_EXACT=HLA-A_g,CD74_g,gui
    ("", BACK, "test4"),
    ("", FORWARD, "test5"),
    ("#4", SUBMIT, "test6"), // enclone BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT NGEX
    ("200", SUBMIT, "test7"),
    ("", BACK, "test8"),
    ("", BACK, "test9"),
    ("10", SUBMIT, "test10"),
];
