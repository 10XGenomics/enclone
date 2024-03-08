// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use serde::{Deserialize, Serialize};
use vdj_ann::annotate::ContigAnnotation;

pub mod html;
pub mod pdb;
pub mod run_test;

/// Add an optional dataset string to a ContigAnnotation.
#[derive(Serialize, Deserialize)]
pub struct AnnotationWithDataset {
    pub dataset: Option<String>,
    #[serde(flatten)]
    pub data: ContigAnnotation,
}
