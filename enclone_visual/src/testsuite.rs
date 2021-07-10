// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::messages::*;
use crate::*;

pub const TESTS: [(&str, MsgFn, &str); 5] = [
    ("#1", Message::SubmitButtonPressed as MsgFn, "test1"),
    ("#2", Message::SubmitButtonPressed as MsgFn, "test2"),
    ("#3", Message::SubmitButtonPressed as MsgFn, "test3"),
    ("", Message::BackButtonPressed as MsgFn, "test4"),
    ("", Message::ForwardButtonPressed as MsgFn, "test5"),
];
