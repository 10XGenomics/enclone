// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::*;
use crate::messages::*;

pub const TESTS: [(&str, MsgFn); 4] = [
    ("#1", Message::SubmitButtonPressed as MsgFn),
    ("#2", Message::SubmitButtonPressed as MsgFn),
    ("#3", Message::SubmitButtonPressed as MsgFn),
    ("", Message::BackButtonPressed as MsgFn),
];

