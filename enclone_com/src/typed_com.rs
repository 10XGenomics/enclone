// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Generic typed communication.
//
// A generic typed message consists of a type name plus a blob of bytes.  The name is tacked
// on with no restrictions, but one could add validity checking for some types.
//
// We include a 64-bit subchannel identifier because it seems like a port could accidentally
// receive a spurious message, and this can be used to drive that possibility to near zero.
//
// In a more sophisticated version, the message could be compressed and/or encrypted.  Encryption
// might be desirable because it is not necessarily desirable for all users on a given computer
// to have visibility of port traffic.
//
// Here is the structure of a generic typed message:
//
// NAME                        TYPE    BYTES
// subchannel identifier       u64     8
// type name length in bytes   u16     2
// type name                   [u8]    up to 256
// body length in bytes        u64     8
// body                        [u8]    up to 2^64 - 1
//
// Some hypothetical candidates for types:
//
// NAME             CONTENT
// json             any JSON file
// svg              any SVG file
// rust:<name>      any rust structure of the given name, packed using serde
// string           any UTF-8-encoded string
// colored-text     a UTF-8-encoded string representing text with coloring/bolding via ANSI escapes
// html             any html file
// request          any UTF-8-encoded string that client would sent to server
// error            any UTF-8-encoded string representing an error, sent from server to client

use std::convert::TryInto;
use std::mem::size_of_val;

type SubChannelIdentifier = u64;

pub const ENCLONE_SUBCHANNEL: SubChannelIdentifier = 123456789; // temporary

pub fn pack_object(id: SubChannelIdentifier, type_name: &str, body: &[u8]) -> Vec<u8> {
    assert!(type_name.as_bytes().len() <= u16::MAX as usize);
    let type_name_length = type_name.as_bytes().len() as u16;
    let body_length = body.len() as u64;
    let len = size_of_val(&id)
        + size_of_val(&type_name_length)
        + type_name.as_bytes().len()
        + size_of_val(&body_length)
        + body.len();
    let mut bytes = Vec::<u8>::new();
    bytes.reserve_exact(len);
    bytes.append(&mut id.to_ne_bytes().to_vec());
    bytes.append(&mut type_name_length.to_ne_bytes().to_vec());
    bytes.append(&mut type_name.as_bytes().to_vec());
    bytes.append(&mut body_length.to_ne_bytes().to_vec());
    bytes.append(&mut body.to_vec());
    assert_eq!(bytes.len(), len);
    bytes
}

pub fn unpack_message(msg: &[u8], id: &mut u64, type_name: &mut Vec<u8>, body: &mut Vec<u8>) {
    if msg.len() == 0 {
        println!("Looks like enclone_client exited, bye!\n");
        std::process::exit(0);
    }
    if msg.len() < 10 {
        eprintln!(
            "\nReceived a coded message of length {}, and thus appears to be truncated.\n",
            msg.len()
        );
        std::process::exit(1);
    }
    *id = u64::from_ne_bytes(msg[0..8].try_into().unwrap());
    let type_name_length = u16::from_ne_bytes(msg[8..10].try_into().unwrap()) as usize;
    if msg.len() < 10 + type_name_length + 8 {
        eprintln!(
            "\nReceived a coded message of length {}, which appears to be truncated.\n",
            msg.len()
        );
        std::process::exit(1);
    }
    *type_name = msg[10..10 + type_name_length].to_vec();
    let body_length = u64::from_ne_bytes(
        msg[10 + type_name_length..10 + type_name_length + 8]
            .try_into()
            .unwrap(),
    ) as usize;
    let body_length = body_length as usize;
    if msg.len() != 10 + type_name_length + 8 + body_length {
        eprintln!(
            "\nReceived a coded message of length {}, which appears to be wrong.\n",
            msg.len()
        );
        std::process::exit(1);
    }
    *body = msg[10 + type_name_length + 8..].to_vec();
}
