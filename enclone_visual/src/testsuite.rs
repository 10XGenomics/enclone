// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::messages::*;
use crate::*;

pub const TESTS: [(&str, MsgFn); 5] = [
    ("#1", Message::SubmitButtonPressed as MsgFn),
    ("#2", Message::SubmitButtonPressed as MsgFn),
    ("#3", Message::SubmitButtonPressed as MsgFn),
    ("", Message::BackButtonPressed as MsgFn),
    ("", Message::ForwardButtonPressed as MsgFn),
];
