// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Automated and manual tests.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// AUTOMATED TESTS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// These are executed by bin/test_vis.

use crate::messages::*;
use crate::*;

const SUBMIT: fn(Result<(), std::string::String>) -> messages::Message =
    Message::SubmitButtonPressed as MsgFn;

const BACK: fn(Result<(), std::string::String>) -> messages::Message =
    Message::BackButtonPressed as MsgFn;

const FORWARD: fn(Result<(), std::string::String>) -> messages::Message =
    Message::ForwardButtonPressed as MsgFn;

const DEL: fn(Result<(), std::string::String>) -> messages::Message =
    Message::DelButtonPressed as MsgFn;

const X0: &str = "enclone woof";
const X1: &str = "enclone BCR=123085 PLOT=gui MIN_CELLS=5 G=12";
const X2: &str = "enclone BCR=123085 CHAINS=4 PLOT_BY_ISOTYPE=gui";

#[rustfmt::skip]
pub const TESTS: [(&str, MsgFn, &str); 19] = [
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
    ("",    DEL,     "test11"),  // #3
];

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// MANUAL TESTS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// 1. enclone VIS=b
//    (a) Type #1 and then manually scroll up/down for two minutes.  Should not slow down or crash.
//    (b) Scroll to the bottom and make sure the last clonotype is not truncated.
//    (c) Check that group number 0 and too large give errors.
//    (d) Check that group number 1 is the right clonotype.
//    (e) Try to cut and paste from clonotype tables [BROKEN].
//    (f) Test help button.
//    (g) Test that copy image button flashes and that it actually copies.
//    (h) Test cookbook commands #5 and #6.
//    (i) enclone BCR=1096354 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui
//        check tooltip functionality
//    (j) enclone BCR=123085:123089 PLOT="gui,s1->red,s2->blue" LEGEND=red,"f 085",blue,"f 089"
//        check tooltip functionality
//    (k) enclone BCR=123085 NOPRINT
//
// 2. Repeat ten times:
//    (a) type enclone VIS
//    (b) type #1
//    (c) make sure it says "thinking" briefly and produces output
//    (d) exit.
//
// 3. enclone VIS and type Help, Dismiss, Help, Dismiss, Help, Dismiss.  Do these respond instantly?
//
// 4. enclone VIS=b
//    - enclone TCR=1146724-1146731,1146751-1146758 PLOT=gui
//    - this is a large dataset, so the canvas is huge
//    - make sure scrolling is smooth and that tooltip responds essentially instantly
//
// 5. test summary
//
// 6. test that group clicks work
//
// 7. test that horizontal resizing works on enclone BCR=123085 CHAINS=4
