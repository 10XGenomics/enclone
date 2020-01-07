// Copyright (c) 2019 10X Genomics, Inc. All rights reserved.

#![deny(missing_docs)]
//!
//! This crate defines the data structure that would represent the clonotypes
//! computed by enclone.
//!

use serde::{Serialize, Deserialize};

/// An amino acid or nuclotide sequence represented as ASCII encoded vector of bytes
pub type Sequence = Vec<u8>;

/// Various regions within a VDJ transcript
#[derive(Clone,Serialize, Deserialize)]
pub enum Region {
    /// 5' untranslated region
    U,
    /// Variable region
    V,
    /// Diversity region
    D,
    /// Joining region
    J,
    /// Constant region
    C,
}

/// Basic cigar string operations
#[derive(Clone,Serialize, Deserialize)]
pub enum CigarOps {
    /// Match(10) represents the cigar string "10="
    Match(u16),
    /// Mismatch(2) represents the cigar string "2X"
    Mismatch(u16),
    /// Insertion(3) represents the cigar string "3I"
    Insertion(u16),
    /// Deletion(6) represents the cigar string "6D"
    Deletion(u16),
    /// SoftClip(13) represents the cigar string "13S"
    SoftClip(u16),
}

/// Representation of an alignment
#[derive(Clone,Serialize, Deserialize)]
pub struct Alignment {
    /// Start of the alignment in the reference
    pub ref_start: usize,
    /// Equivalent of a cigar string
    pub cigar: Vec<CigarOps>,
}

/// Define a chain within an exact subclonotype.
#[derive(Serialize, Deserialize)]
pub struct ExactClonotypeChain {
    /// Nucleotide sequence of the chain. This will only contain ACGT alphabets
    pub nt_sequence: Sequence,
    /// Amino acid sequence from the start codon at the beginning of the V-REGION.
    /// This can be inferred from the `nt_sequence` and `v_start`, but stored for convenience
    pub aa_sequence: Sequence,
    /// Index of the start of the V-REGION in the `nt_sequence`.
    pub v_start: usize,
    /// Index of the end of the J-REGION in the `nt_sequence` (exclusive).
    pub j_end: usize,
    /// Index of the C-REGION of this chain in the universal reference.
    /// TODO: Should we store UVDJ regions here for convenience? The reason why it is not
    /// stored at this level is because all exact subclonotypes share the same UVDJ regions.
    pub c_region_idx: Option<usize>,
    /// Index of the start of the CDR3 sequence in the `nt_sequence`. The start of the CDR3
    /// amino acid in the `aa_sequence` is `(cdr3_start - v_start)/3`
    pub cdr3_start: usize,
    /// Index of the end of the CDR3 sequence in the `nt_sequence` (exclusive). The end of the
    /// CDR3 amino acid in the `aa_sequence` is `(cdr3_end - v_start)/3`
    pub cdr3_end: usize,
    /// UMI counts of contigs associated with this exact subclonotype chain. The number of elements
    /// in this vector is equal to the number of barcodes associated with this exact subclonotype.
    pub umi_counts: Vec<u32>,
    /// Names of contigs associated with this exact subclonotype chain. The number of elements
    /// in this vector is equal to the number of barcodes associated with this exact subclonotype.
    /// The contig name would be of the form `{barcode}_contig_{id}`.
    pub contig_ids: Vec<String>,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the clonotype consensus
    /// of this chain.
    /// TODO: Do we need amino acid alignment info?
    pub clonotype_consensus_aln: Alignment,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the concatenated donor
    /// reference of this chain.
    /// TODO: Default donor reference to universal reference?
    pub donor_reference_aln: Option<Alignment>,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the concatenated universal
    /// reference of this chain.
    pub universal_reference_aln: Alignment,
}

/// Define an exact subclonotype.
///
/// All the barcodes within an exact subclonotype have the same number of productive chains,
/// the same sequence from the start of the V-REGION to the end of the J-REGION
/// as well as the same C-REGION annotation for each chain.
/// TODO: Maybe mutations outside V-J?
/// Would call this ExactSubclonotype but that name is used elsewhere and reusing it could be 
/// confusing.
#[derive(Serialize, Deserialize)]
pub struct EClonotype {
    /// The chains in an exact subclonotype. The number of elements in this vector is
    /// equal to the total number of chains in the parent clonotype. The order of chains
    /// is consistent with the order in the parent clonotype and at least one element in this
    /// vector is not a `None`.
    pub chains: Vec<Option<ExactClonotypeChain>>,
    /// List of cell barcodes in this exact subclonotype. The number of elements in this list is
    /// equal to the number of elements in the `umi_counts` and `contig_ids` vector in an
    /// `ExactClonotypeChain`. The barcodes have the appropriate gem group as the suffix.
    pub cell_barcodes: Vec<String>,
}

/// Define a clonotype chain
#[derive(Serialize, Deserialize)]
pub struct ClonotypeChain {
    /// The nuclotide sequence of this clonotype chain consensus.
    /// What we actually compute here is not the consensus across the clonotype (whose
    /// biological meaning is questionable), but rather the sequence of the first exact
    /// subclonotype which has an entry for the given chain.  Over 99% of the time, this
    /// will be the first exact subclonotype.
    pub nt_sequence: Sequence,
    /// Index of the 5' UTR region in the universal reference. The region in the universal
    /// reference is guaranteed to be `Region::U`
    pub u_idx: Option<usize>,
    /// Index of the Variable region in the universal reference. The region in the universal
    /// reference is guaranteed to be `Region::V`
    pub v_idx: usize,
    /// Index of the Diversity region in the universal reference. The region in the universal
    /// reference is guaranteed to be `Region::D`. D-REGION is not present for light/alpha
    /// chains. Even for heavy/beta chains, this might be `None` if there is ambiguity in the
    /// annotation.
    pub d_idx: Option<usize>,
    /// Index of the Joining region in the universal reference. The region in the universal
    /// reference is guaranteed to be `Region::J`
    pub j_idx: usize,
    /// Index of the Variable region in the donor reference. The region in the donor
    /// reference is guaranteed to be `Region::V` and the `universal_idx` in the donor reference
    /// item will be equal to the `v_idx`
    pub donor_v_idx: Option<usize>,
    /// Index of the Joining region in the donor reference. The region in the donor
    /// reference is guaranteed to be `Region::J` and the `universal_idx` in the donor reference
    /// item will be equal to the `j_idx`
    pub donor_j_idx: Option<usize>,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the concatenated donor
    /// reference of this chain.
    /// Concatenated donor reference =
    ///     `nt_sequence` of donor_reference[donor_v_idx] if donor_v_idx is not None +
    ///     `nt_sequence` of donor_reference[donor_j_idx] if donor_j_idx is not None
    pub donor_reference_aln: Option<Alignment>,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the concatenated universal
    /// reference of this chain.
    /// Concatenated universal reference =
    ///     `nt_sequence` of universal_reference[u_idx] if u_idx is not None +
    ///     `nt_sequence` of universal_reference[v_idx] +
    ///     `nt_sequence` of universal_reference[d_idx] if d_idx is not None+
    ///     `nt_sequence` of universal_reference[j_idx]
    pub universal_reference_aln: Alignment,
}

/// Definition of a clonotype.
///
/// A clonotype is composed of a list of exact subclonotypes
#[derive(Serialize, Deserialize)]
pub struct Clonotype {
    /// The list of chains associated with this clonotype. The ordering of the chains is important
    /// as this order is preserved in the list of chains specified under each exact subclonotype. By
    /// convention heavy chain/beta chain comes ahead of light chain/alpha chain.
    /// TODO: What is the ordering when multiple chains of same kind are present?
    pub chains: Vec<ClonotypeChain>,
    /// The list of exact subclonotypes in this clonotype ordered by the number of cell barcodes in the
    /// exact subclonotype in descending order (TODO: Verify sort order). The number of chains listed
    /// under each exact subclonotype will be equal to the number of chains in this clonotype in the
    /// same order. However, some of the exact subclonotype chains could be `None`.
    pub exact_clonotypes: Vec<EClonotype>,
    /// The total number of cell barcodes associated with this clonotype. This can be inferred by
    /// summing up the number of barcodes within each exact subclonotype, but it is stored here for
    /// convenience.
    pub frequency: u32,
}

/// A single universal reference sequence and metadata packaged in a convenient struct.
#[derive(Serialize, Deserialize)]
pub struct UniversalReferenceItem {
    /// A unique identifier for this reference sequence that traces back to the reference fasta.
    pub ref_idx: u32,
    /// The display name of this gene which will be shown in Loupe (optional allele information)
    pub display_name: String,
    /// One of the U/V/D/J/C regions
    pub region: Region,
    /// Nucleotide sequence associated with this reference item
    pub nt_sequence: Sequence,
}

/// Universal reference is a list of sequences and its associated metadata.
pub type UniversalReference = Vec<UniversalReferenceItem>;

/// A single donor reference sequence and metadata packaged in a convenient struct. In the
/// current version of enclone, the donor reference is only inferred for the V-REGIONs. But we could
/// extend it to J-REGIONs in the future.
#[derive(Clone,Serialize, Deserialize)]
pub struct DonorReferenceItem {
    /// Index of the parent sequence in the universal reference
    pub universal_idx: usize,
    /// Index of the donor associated with this reference. If there are no donors specified, this
    /// will be `0` by default.
    pub donor_idx: usize,
    /// The display name of this gene which will be shown in Loupe (optional allele information)
    /// TODO: Should this be modified to explicitly point out the donor? e.g TRAV-1 [Donor 0]?
    /// for now, like this: "TRAV-1, donor 1, alt allele 1", etc.
    pub display_name: String,
    /// Currently, the donor reference region will only be the V-REGION
    pub region: Region,
    /// The nucleotide sequence associated with this reference item
    pub nt_sequence: Sequence,
    /// Alignment of the `nt-sequence` with the nucleotide sequence of the corresponding universal
    /// reference item.
    pub universal_aln: Alignment,
}

/// Donor reference is a list of sequences and its associated metadata.
pub type DonorReference = Vec<DonorReferenceItem>;

/// Outputs from a single enclone run.
///
/// TODO: It might be easier to store the three vectors as three different files if we want to
/// stream through it.
#[derive(Serialize, Deserialize)]
pub struct EncloneOutputs {
    /// List of all clonotypes computed in this enclone run.
    pub clonotypes: Vec<Clonotype>,
    /// List of all universal reference sequences and metadata packaged in a convenient struct.
    /// UV(D)JC regions associated with a clonotype chain or an exact subclonotype chain are stored as
    /// indices into this vector
    pub universal_reference: UniversalReference,
    /// List of all donor reference sequences and metadata packaged in a convenient struct.
    /// The donor V-REGION associated with a clonotype chain is stored as an index to this vector.
    pub donor_reference: DonorReference,
}

// TODO: Donor names?
// TODO: Aggr metadata structure
