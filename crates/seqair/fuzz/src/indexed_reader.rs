//! Types for the `fuzz_reader_indexed` target — shared between the fuzz target
//! and the seed generator so seeds always match the deserialization format.

use arbitrary::Arbitrary;

/// Top-level fuzz input for the indexed reader pipeline.
#[derive(Arbitrary, Debug)]
pub struct Input {
    pub readers: ReadersInput,
    pub alignment: FastaInput,
}

/// Which alignment format to test.
#[derive(Arbitrary, Debug)]
pub enum ReadersInput {
    Bam { bam: Vec<u8>, bai: Vec<u8> },
    Sam { sam: Vec<u8>, sai: Vec<u8> },
    Cram { cram: Vec<u8>, crai: Vec<u8> },
}

/// Reference FASTA data (shared across all alignment formats).
#[derive(Arbitrary, Debug)]
pub struct FastaInput {
    pub fasta_gz: Vec<u8>,
    pub fai: String,
    pub gzi: Vec<u8>,
}
