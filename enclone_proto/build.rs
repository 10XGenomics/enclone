// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// The purpose of this file is to auto generate `types.rs` from the `types.proto` file.

use prost_build::Config;

fn main() {
    let mut config = Config::new();
    config.type_attribute(".", "#[derive(::serde::Serialize, ::serde::Deserialize)]");
    config.compile_protos(&["types.proto"], &["."]).unwrap();
}
