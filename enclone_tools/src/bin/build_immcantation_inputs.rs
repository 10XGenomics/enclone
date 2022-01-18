// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Build input files for Immcantation.  This creates two files:
// filtered_contig.fasta
// filtered_contig_annotations.csv.

use enclone_tools::main_build_immcantation_inputs::main_build_immcantation_inputs;

pub fn main() {
    main_build_immcantation_inputs()
}
