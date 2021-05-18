// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // compiling protos using path on build time

    tonic_build::compile_protos("./server.proto")?;

    Ok(())
}
