// Copyright (c) 2019 10X Genomics, Inc. All rights reserved.

//!
//! Code for reading and writing custom proto files
//!

use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use failure::Error;
use prost::Message;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use types::EncloneOutputs;

const BUFFER_CAPACITY: usize = 1_000_000;

pub struct ProtoWriter<W: Write> {
    encode_buffer: Vec<u8>,
    writer: W,
}

impl<W: Write> ProtoWriter<W> {
    pub fn with_writer(writer: W) -> Self {
        ProtoWriter {
            encode_buffer: Vec::with_capacity(BUFFER_CAPACITY),
            writer,
        }
    }
    pub fn encode_and_write<M>(&mut self, message: M) -> Result<usize, Error>
    where
        M: Message,
    {
        let encoded_len = message.encoded_len();
        message.encode(&mut self.encode_buffer)?;
        self.writer.write_u32::<BigEndian>(encoded_len as u32)?;
        self.writer.write_all(&self.encode_buffer)?;
        self.encode_buffer.clear();
        Ok(encoded_len)
    }
}

pub struct ProtoReader<R: Read> {
    decode_buffer: Vec<u8>,
    reader: R,
}

impl<R: Read> ProtoReader<R> {
    pub fn with_reader(reader: R) -> Self {
        ProtoReader {
            decode_buffer: Vec::with_capacity(BUFFER_CAPACITY),
            reader,
        }
    }
    pub fn read_and_decode<M>(&mut self) -> Result<M, Error>
    where
        M: Message + Default,
    {
        self.reader
            .by_ref()
            .take(4)
            .read_to_end(&mut self.decode_buffer)?;
        if self.decode_buffer.len() != 4 {
            return Err(format_err!(
                "Expected to get 4 bytes for length. Got {} bytes!",
                self.decode_buffer.len()
            ));
        }
        let decoded_len = self.decode_buffer.as_slice().read_u32::<BigEndian>()?;
        self.decode_buffer.clear();
        self.reader
            .by_ref()
            .take(decoded_len as u64)
            .read_to_end(&mut self.decode_buffer)?;
        let decoded_message = M::decode(&self.decode_buffer)?;
        self.decode_buffer.clear();
        Ok(decoded_message)
    }
}

pub fn write_proto(enclone_outputs: EncloneOutputs, path: impl AsRef<Path>) -> Result<(), Error> {
    let writer = BufWriter::new(File::create(path)?);
    let mut proto_writer = ProtoWriter::with_writer(writer);

    // Write the universal reference
    proto_writer.encode_and_write(enclone_outputs.universal_reference)?;
    // Write the donor reference
    proto_writer.encode_and_write(enclone_outputs.donor_reference)?;
    // Write the number of clonotypes. Not bothering to write this raw
    proto_writer.encode_and_write(enclone_outputs.clonotypes.len() as u32)?;
    for cl in enclone_outputs.clonotypes {
        proto_writer.encode_and_write(cl)?;
    }
    Ok(())
}

pub fn read_proto(path: impl AsRef<Path>) -> Result<EncloneOutputs, Error> {
    let reader = BufReader::new(File::open(path)?);
    let mut proto_reader = ProtoReader::with_reader(reader);

    // Read the universal reference
    let universal_reference = proto_reader.read_and_decode()?;
    // Read the donor reference
    let donor_reference = proto_reader.read_and_decode()?;
    // Number of clonotyopes
    let num_clonotypes: u32 = proto_reader.read_and_decode()?;
    let mut clonotypes = Vec::new();
    for _ in 0..num_clonotypes {
        clonotypes.push(proto_reader.read_and_decode()?);
    }

    Ok(EncloneOutputs {
        universal_reference,
        donor_reference,
        clonotypes,
    })
}
