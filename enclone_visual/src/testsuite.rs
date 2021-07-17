// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::messages::*;
use crate::*;

const SUBMIT: fn(Result<(), std::string::String>) -> messages::Message =
    Message::SubmitButtonPressed as MsgFn;

const BACK: fn(Result<(), std::string::String>) -> messages::Message =
    Message::BackButtonPressed as MsgFn;

const FORWARD: fn(Result<(), std::string::String>) -> messages::Message =
    Message::ForwardButtonPressed as MsgFn;

const X0: &str = "enclone woof";
const X1: &str = "enclone BCR=123085 PLOT=gui MIN_CELLS=5 G=12";
const X2: &str = "enclone BCR=123085 CHAINS=4 PLOT_BY_ISOTYPE=gui";

#[rustfmt::skip]
pub const TESTS: [(&str, MsgFn, &str); 18] = [
    (X0,    SUBMIT,  ""),        // enclone woof
    ("#1",  SUBMIT,  "test1"),   // enclone BCR=123085 PLOT=gui MIN_CELLS=5
    ("",    BACK,    ""),        // enclone woof
    ("",    FORWARD, ""),        // #1
    ("10",  SUBMIT,  "test1b"),  // 10
    ("",    BACK,    ""),        // #1
    ("",    FORWARD, "test1x"),  // 10
    (X1,    SUBMIT,  "test1c"),  // enclone BCR=123085 PLOT=gui MIN_CELLS=5 G=12
    (X2,    SUBMIT,  "test1d"),  // enclone BCR=123085 CHAINS=4 PLOT_BY_ISOTYPE=gui
    ("#2",  SUBMIT,  "test2"),   // enclone BCR=123085 PLOT_BY_ISOTYPE=gui MIN_CELLS=5
    ("#3",  SUBMIT,  "test3"),   // enclone BCR=123085 GEX=123217 PLOTXY_EXACT=HLA-A_g,CD74_g,gui
    ("",    BACK,    "test4"),   // #2
    ("",    FORWARD, "test5"),   // #3
    ("#4",  SUBMIT,  "test6"),   // enclone BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT NGEX
    ("200", SUBMIT,  "test7"),   // 200
    ("",    BACK,    "test8"),   // #4
    ("",    BACK,    "test9"),   // #3
    ("10",  SUBMIT,  "test10"),  // 10
];
