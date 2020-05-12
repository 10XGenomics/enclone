// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// The purpose of this file is to make some version information available so that it can be
// printed out at appropriate points by enclone.  This files is a slightly modified version
// of https://vallentin.dev/2019/06/06/versioning.

extern crate chrono;

use prost_build::Config;

fn main() {
    let mut config = Config::new();
    config.type_attribute(".", "#[derive(::serde::Serialize, ::serde::Deserialize)]");
    config.compile_protos(&["types.proto"], &["."]).unwrap();
}
