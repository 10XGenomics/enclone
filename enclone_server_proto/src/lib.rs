#![deny(warnings)]

pub mod proto {
    tonic::include_proto!("enclone_server");
}
