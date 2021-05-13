// extern crate protoc_rust;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // compiling protos using path on build time
    tonic_build::compile_protos("./server.proto")?;
    Ok(())
}
