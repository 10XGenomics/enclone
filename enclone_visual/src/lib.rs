// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod enclone_client;
pub mod enclone_server;

pub mod proto {
    tonic::include_proto!("enclone");
}
