// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

pub fn profiling_blacklist() -> Vec<String> {
    let blacklist = [
        "alloc",
        "build",
        "core",
        "core-arch",
        "crossbeam-deque",
        "crossbeam-epoch",
        "debruijn",
        "float-ord",
        "hashbrown",
        "hdf5-rust",
        "hdf5-types",
        "lock_api",
        "lz4",
        "ndarray",
        "parking_lot",
        "parking_lot_core",
        "rayon",
        "rayon-core",
        "regex",
        "regex-syntax",
        "rust-bio",
        "serde",
        "serde_json",
        "std",
        "superslice",
        "tokio",
        "unknown",
    ];
    let mut b = Vec::<String>::new();
    for x in blacklist.iter() {
        b.push(x.to_string());
    }
    b
}
