// Copyright (c) 2019 10X Genomics, Inc. All rights reserved.

// Build input files for Immcantation.  This creates two files:
// filtered_contig.fasta
// filtered_contig_annotations.csv.

extern crate enclone;
use enclone::*;

use main_build_immcantation_inputs::*;

pub fn main() {
    main_build_immcantation_inputs()
}
