// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is a text client, which is useless except for debugging and experimentation.

use enclone_server::proto::{
    analyzer_client::AnalyzerClient,
    EncloneRequest,
};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("connecting to client"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    let mut client = AnalyzerClient::connect("http://127.0.0.1:7000").await?;

    println!("creating request"); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    let request = tonic::Request::new(EncloneRequest {
        args: "BCR=123085 MIN_CELLS=5 PLOT_BY_ISOTYPE=gui".into(),
    });

    let response = client.enclone(request).await?;

    println!("RESPONSE={:?}", response);

    Ok(())
}
