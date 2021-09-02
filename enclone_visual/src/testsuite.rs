// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Automated and manual tests.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// AUTOMATED TESTS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// These are executed by bin/test_vis.

use crate::messages::*;
use crate::*;
use std::env;

#[rustfmt::skip]
pub fn metatests() -> Vec<Vec<Message>> {
    let mut user = String::new();
    for (key, value) in env::vars() {
        if key == "USER" {
            user = value.to_string();
        }
    }
    vec![
        // 1 = enclone VIS
        vec![
            Message::ArchiveOpen(Ok(())),
            Message::ExpandArchiveEntry(true, 0),
            Message::SetName("expand"),   // test pushing expand button once
            Message::ExpandArchiveEntry(false, 0),
            Message::SetName("unexpand"), // test pushing expand button second time
            Message::DeleteArchiveEntry(true, 0),
            Message::SetName("delete1"),  // see if delete message shows up
            Message::ArchiveClose,
            Message::ArchiveOpen(Ok(())),
            Message::SetName("delete2"),  // see if deletion occurred 
            Message::Restore(true, 1),
            Message::SetName("restore1"), // see if restore message shows up
            Message::ArchiveClose,
            Message::SetName("restore2"), // see if restore occurs
            Message::ArchiveOpen(Ok(())),
            Message::ArchiveName("honeycomb".to_string(), 1),
            Message::ArchiveNameChange(1),
            Message::CompleteArchiveNameChange(Ok(())),
            Message::ArchiveClose,
            Message::ArchiveOpen(Ok(())),
            Message::SetName("rename"),   // see if rename occurs
        ],
        // 2 = enclone VIS=b
        vec![
            Message::ArchiveOpen(Ok(())),
            Message::ArchiveShare(true, 0),
            Message::UserName(user, 0),
            Message::UserSelected(true, 0),
            Message::DoShare(true),
            Message::CompleteDoShare(Ok(())),
            Message::ArchiveRefresh,
            Message::ArchiveRefreshComplete(Ok(())),
            Message::ExpandArchiveEntry(true, 0),
            Message::SetName("share"),    // test share
        ],
        // 3 = enclone VIS FAIL_ON_ERROR
        vec![
            Message::ArchiveOpen(Ok(())),
            Message::RestoreCookbook(true, 0),
            Message::ArchiveClose,
            Message::SubmitButtonPressed(Ok(())),
            Message::WaitCommand(Ok(())),
        ],
        // 4 = enclone VIS
        vec![
            Message::InputChanged1("enclone TCR_GEX=1175299-1175300 SUMMARY_CLEAN".to_string()),
            Message::SubmitButtonPressed(Ok(())),
            Message::WaitCommand(Ok(())),
            Message::SummaryOpen(Ok(())),
            Message::SetName("metrics"),
            Message::MetricButton(2),
            Message::MetricButton(3),
            Message::CondenseMetrics,
            Message::SetName("select_metrics"),
            Message::SummaryClose(Ok(())),
            Message::Save,
            Message::ArchiveOpen(Ok(())),
            Message::Restore(true, 0),
            Message::ArchiveClose,
            Message::SummaryOpen(Ok(())),
            Message::SetName("restore_metrics"),
        ],
    ]
}

const SUBMIT: fn(Result<(), std::string::String>) -> messages::Message =
    Message::SubmitButtonPressed as MsgFn;

const BACK: fn(Result<(), std::string::String>) -> messages::Message =
    Message::BackButtonPressed as MsgFn;

const FORWARD: fn(Result<(), std::string::String>) -> messages::Message =
    Message::ForwardButtonPressed as MsgFn;

const DEL: fn(Result<(), std::string::String>) -> messages::Message =
    Message::DelButtonPressed as MsgFn;

const HELP1: fn(Result<(), std::string::String>) -> messages::Message = Message::HelpOpen as MsgFn;
const HELP2: fn(Result<(), std::string::String>) -> messages::Message = Message::HelpClose as MsgFn;

const SUMMARY1: fn(Result<(), std::string::String>) -> messages::Message =
    Message::SummaryOpen as MsgFn;
const SUMMARY2: fn(Result<(), std::string::String>) -> messages::Message =
    Message::SummaryClose as MsgFn;

const ARCH1: fn(Result<(), std::string::String>) -> messages::Message =
    Message::ArchiveOpen as MsgFn;

const X0: &str = "enclone woof";
const X1: &str = "enclone BCR=123085 PLOT=gui MIN_CELLS=5 G=12";
const X2: &str = "enclone BCR=123085 CHAINS=4 PLOT_BY_ISOTYPE=gui";
const X3: &str = "enclone + BCR=123085 NOPRINT";
const X4: &str = "enclone BCR=123085 CHAINS=10";
const X5: &str = "enclone BCR=123085 KEEP_CLONO_IF_CELL_MAX=\"u1 >= 6000\" SEG=IGHM";
const X6: &str = "enclone BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT NGEX \
                          SIM_MAT_PLOT=gui,fb1_n,fb2_n,fb3_n,fb4_n,fb5_n SUMMARY_CLEAN";

#[rustfmt::skip]
pub const TESTS: [(&str, MsgFn, &str); 40] = [
    (X0,     SUBMIT,  ""),        // enclone woof
    ("#1",   SUBMIT,  "test1"),   // enclone BCR=123085 PLOT=gui MIN_CELLS=5
    ("#999", SUBMIT,  "test1a"),  // #999
    ("",     DEL,     ""),        // enclone BCR=123085 PLOT=gui MIN_CELLS=5
    ("",     BACK,    ""),        // enclone woof
    ("",     FORWARD, ""),        // #1
    ("10",   SUBMIT,  "test1b"),  // 10
    ("",     BACK,    ""),        // #1
    ("",     FORWARD, "test1x"),  // 10
    (X1,     SUBMIT,  "test1c"),  // enclone BCR=123085 PLOT=gui MIN_CELLS=5 G=12
    (X2,     SUBMIT,  "test1d"),  // enclone BCR=123085 CHAINS=4 PLOT_BY_ISOTYPE=gui
    ("#2",   SUBMIT,  "test2"),   // enclone BCR=123085 PLOT_BY_ISOTYPE=gui MIN_CELLS=5
    ("#3",   SUBMIT,  "test3"),   // enclone BCR=123085 GEX=123217 PLOTXY_EXACT=HLA-A_g,CD74_g,gui
    ("",     BACK,    "test4"),   // #2
    ("",     FORWARD, "test5"),   // #3
    ("#4",   SUBMIT,  "test6"),   // enclone BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT NGEX
    ("200",  SUBMIT,  "test7"),   // 200
    ("",     BACK,    "test8"),   // #4
    ("",     BACK,    "test9"),   // #3
    ("10",   SUBMIT,  "test10"),  // 10
    ("",     DEL,     "test11"),  // #3
    ("",     BACK,    ""),        // #2
    ("",     DEL,     ""),        // X2 = ...
    ("",     DEL,     ""),        // X1 = ...
    ("",     DEL,     ""),        // 10
    ("",     DEL,     ""),        // #1
    ("",     DEL,     ""),        // X0 = ...
    ("",     DEL,     "test12"),  // #3
    (X3,     SUBMIT,  ""),        // enclone + BCR=123085 NOPRINT
    ("5",    SUBMIT,  ""),        // 5
    ("",     BACK,    "test13"),  // enclone + BCR=123085 NOPRINT
    (X4,     SUBMIT,  "test14"),  // enclone BCR=123085 CHAINS=10
    (X5,     SUBMIT,  "test15"),  // enclone BCR=123085 KEEP_CLONO_IF_CELL_MAX="u1 >= 6000" SEG=IGHM
    ("#5",   SUBMIT,  "test16"),  // enclone BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT 
                                  //         NGEX LVARSP=fb1,fb1_n,fb2,fb2_n
    (X6,     SUBMIT,  "test17"),  // enclone BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT NGEX 
                                  //         SIM_MAT_PLOT=gui,fb1_n,fb2_n,fb3_n,fb4_n,fb5_n
                                  //         SUMMARY_CLEAN
    ("",     HELP1,   "test18"),  // (help)
    ("",     HELP2,   ""),        // enclone BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT NGEX
                                  //         SIM_MAT_PLOT=gui,fb1_n,fb2_n,fb3_n,fb4_n,fb5_n
                                  //         SUMMARY_CLEAN
    ("",     SUMMARY1, "test19"), // (summary)
    ("",     SUMMARY2, ""),       // enclone BCR=1145040 GEX=1142282 ALLOW_INCONSISTENT NGEX
                                  //         SIM_MAT_PLOT=gui,fb1_n,fb2_n,fb3_n,fb4_n,fb5_n
                                  //         SUMMARY_CLEAN
    ("",     ARCH1,    "test20"), // (archive)
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
//    (f) Test that copy image button flashes and that it actually copies.
//    (g) enclone BCR=1096354 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui
//        check tooltip functionality and group clicks
//    (h) enclone BCR=123085:123089 PLOT="gui,s1->red,s2->blue" LEGEND=red,"f 085",blue,"f 089"
//        check tooltip functionality and group clicks
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
// 5. test that horizontal resizing works on enclone BCR=123085 CHAINS=4
//    and test that there is no limit on horizontal resizing
