//! Pure-Rust BAM/SAM/CRAM/FASTA reader, pileup engine, and VCF/BCF writer.
//!
//! seqair provides indexed random access to alignment files (BAM, SAM, CRAM),
//! reference sequences (FASTA), and a column-based pileup engine. It also
//! includes VCF/BCF writers with type-safe allele representation and single-pass
//! index co-production.
//!
//! # Quick start: pileup at a genomic region
//!
//! The most common workflow is opening an alignment + reference, fetching a
//! region, and iterating pileup columns:
//!
//! ```no_run
//! use seqair::reader::Readers;
//! use seqair::bam::pileup::PileupOp;
//! use seqair_types::{Pos, Zero};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Open BAM + FASTA together — auto-detects format (BAM/SAM/CRAM)
//! let mut readers = Readers::open("sample.bam".as_ref(), "reference.fa".as_ref())?;
//!
//! // Fetch a region and get a pileup iterator
//! let tid = readers.header().tid("chr19").expect("contig not found");
//! let start = Pos::<Zero>::new(6_100_000).unwrap();
//! let end = Pos::<Zero>::new(6_200_000).unwrap();
//! let mut pileup = readers.pileup(tid, start, end)?;
//!
//! // Iterate over columns (one per reference position)
//! for column in pileup.by_ref() {
//!     let pos = column.pos();
//!     let ref_base = column.reference_base();
//!     let depth = column.depth();
//!
//!     // Inspect each read's contribution at this position
//!     for aln in column.alignments() {
//!         match &aln.op {
//!             PileupOp::Match { base, qual, .. } => {
//!                 // Read shows base with quality
//!             }
//!             PileupOp::Deletion { del_len } => {
//!                 // Read has a deletion spanning this position
//!             }
//!             PileupOp::Insertion { base, qual, insert_len, .. } => {
//!                 // Read has an insertion after the previous match
//!             }
//!             PileupOp::RefSkip => {
//!                 // Spliced alignment skips this position (RNA-seq)
//!             }
//!         }
//!     }
//! }
//!
//! // Recover the RecordStore for reuse in the next region
//! readers.recover_store(&mut pileup);
//! # Ok(())
//! # }
//! ```
//!
//! # Multi-threaded pileup
//!
//! Both alignment and FASTA readers support [`fork()`](reader::Readers::fork)
//! for cheap per-thread copies that share the parsed index and header via `Arc`:
//!
//! ```no_run
//! # use seqair::reader::Readers;
//! # use seqair_types::{Pos, Zero};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let readers = Readers::open("sample.bam".as_ref(), "ref.fa".as_ref())?;
//!
//! // Each thread forks its own reader — no locking, shared index
//! std::thread::scope(|s| {
//!     for _ in 0..4 {
//!         let mut thread_reader = readers.fork().unwrap();
//!         s.spawn(move || {
//!             let pileup = thread_reader.pileup(0,
//!                 Pos::<Zero>::new(0).unwrap(),
//!                 Pos::<Zero>::new(100_000).unwrap(),
//!             ).unwrap();
//!             for column in pileup {
//!                 // process in parallel...
//!             }
//!         });
//!     }
//! });
//! # Ok(())
//! # }
//! ```
//!
//! # Writing VCF/BCF from pileup results
//!
//! After computing variant calls from pileup data, write them to VCF or BCF:
//!
//! ```
//! use seqair::vcf::*;
//! use seqair::vcf::alleles::Alleles;
//! use seqair::vcf::writer::VcfWriter;
//! use seqair::vcf::record::{VcfRecordBuilder, Genotype, SampleValue};
//! use seqair_types::{Base, Pos, One, SmolStr};
//! use std::sync::Arc;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // 1. Build a header with field definitions
//! let header = Arc::new(
//!     VcfHeader::builder()
//!         .add_contig("chr1", ContigDef { length: Some(248_956_422) })?
//!         .add_info("DP", InfoDef {
//!             number: Number::Count(1),
//!             typ: ValueType::Integer,
//!             description: SmolStr::from("Total read depth"),
//!         })?
//!         .add_format("GT", FormatDef {
//!             number: Number::Count(1),
//!             typ: ValueType::String,
//!             description: SmolStr::from("Genotype"),
//!         })?
//!         .add_format("DP", FormatDef {
//!             number: Number::Count(1),
//!             typ: ValueType::Integer,
//!             description: SmolStr::from("Sample read depth"),
//!         })?
//!         .add_sample("sample1")?
//!         .build()?
//! );
//!
//! // 2. Create a writer (VCF text to stdout, or BCF to file)
//! let mut output = Vec::new();
//! let mut writer = VcfWriter::new(&mut output, header.clone());
//! writer.write_header()?;
//!
//! // 3. Build records with type-safe alleles
//! let alleles = Alleles::snv(Base::A, Base::T)?;
//! let record = VcfRecordBuilder::new("chr1", Pos::<One>::new(12345).unwrap(), alleles)
//!     .qual(30.0)
//!     .filter_pass()
//!     .info_integer("DP", 50)
//!     .format_keys(&["GT", "DP"])
//!     .add_sample(vec![
//!         SampleValue::Genotype(Genotype::unphased(0, 1)),
//!         SampleValue::Integer(45),
//!     ])
//!     .build(&header)?;
//!
//! writer.write_record(&record)?;
//! writer.finish()?;
//!
//! // output now contains valid VCF text
//! let vcf_text = String::from_utf8(output)?;
//! assert!(vcf_text.contains("chr1\t12345\t.\tA\tT\t30.0\tPASS\tDP=50"));
//! # Ok(())
//! # }
//! ```
//!
//! # Type-safe alleles
//!
//! The [`Alleles`](vcf::Alleles) enum prevents invalid REF/ALT combinations
//! at construction time:
//!
//! ```
//! use seqair::vcf::alleles::Alleles;
//! use seqair_types::Base;
//!
//! // SNV: single ref base, single alt base
//! let snv = Alleles::snv(Base::A, Base::T).unwrap();
//! assert_eq!(snv.ref_text(), "A");
//! assert_eq!(snv.rlen(), 1);
//!
//! // Insertion: anchor + inserted bases → REF=A, ALT=ACGT
//! let ins = Alleles::insertion(Base::A, &[Base::C, Base::G, Base::T]).unwrap();
//! assert_eq!(ins.ref_text(), "A");
//!
//! // Deletion: anchor + deleted bases → REF=ACGT, ALT=A
//! let del = Alleles::deletion(Base::A, &[Base::C, Base::G, Base::T]).unwrap();
//! assert_eq!(del.ref_text(), "ACGT");
//! assert_eq!(del.rlen(), 4);
//!
//! // Compile-time safety: can't create SNV where ref == alt
//! assert!(Alleles::snv(Base::A, Base::A).is_err());
//! // Can't create empty insertion
//! assert!(Alleles::insertion(Base::A, &[]).is_err());
//! ```
//!
//! # High-performance BCF encoding
//!
//! For hot paths (millions of records), the [`encoder`](vcf::encoder) module
//! provides zero-allocation direct BCF encoding with pre-resolved typed handles:
//!
//! ```
//! use seqair::vcf::bcf_writer::BcfWriter;
//! use seqair::vcf::encoder::*;
//! use seqair::vcf::header::*;
//! use seqair::vcf::record::Genotype;
//! use seqair::vcf::alleles::Alleles;
//! use seqair_types::{Base, Pos, One, SmolStr};
//! use std::sync::Arc;
//! use std::marker::PhantomData;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let header = Arc::new(VcfHeader::builder()
//!     .add_contig("chr1", ContigDef { length: Some(1000) })?
//!     .add_info("DP", InfoDef {
//!         number: Number::Count(1), typ: ValueType::Integer,
//!         description: SmolStr::from("Depth"),
//!     })?
//!     .add_format("GT", FormatDef {
//!         number: Number::Count(1), typ: ValueType::String,
//!         description: SmolStr::from("Genotype"),
//!     })?
//!     .add_sample("sample1")?
//!     .build()?
//! );
//!
//! // Resolve handles once at setup — no string lookups per record
//! let contig = ContigHandle(0);
//! let dp = ScalarInfoHandle::<i32> {
//!     dict_idx: header.string_map().get("DP").unwrap() as u32,
//!     _marker: PhantomData,
//! };
//! let gt = GtFormatHandle {
//!     dict_idx: header.string_map().get("GT").unwrap() as u32,
//! };
//!
//! let mut buf = Vec::new();
//! let mut writer = BcfWriter::new(&mut buf, header, false);
//! writer.write_header()?;
//!
//! // Hot loop: zero allocations, pre-resolved dict indices
//! let alleles = Alleles::snv(Base::A, Base::T)?;
//! let mut enc = writer.record_encoder();
//! alleles.begin_record(&mut enc, contig, Pos::<One>::new(100).unwrap(), Some(30.0));
//! FilterHandle::PASS.encode(&mut enc);
//! dp.encode(&mut enc, 50);          // handle.encode(enc, value)
//! enc.begin_samples(1);
//! gt.encode(&mut enc, &Genotype::unphased(0, 1));
//! enc.emit()?;
//!
//! writer.finish()?;
//! assert!(!buf.is_empty()); // valid BCF bytes
//! # Ok(())
//! # }
//! ```
//!
//! # Reading FASTA reference sequences
//!
//! ```no_run
//! use seqair::fasta::IndexedFastaReader;
//! use seqair_types::{Base, Pos, Zero};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut reader = IndexedFastaReader::open("reference.fa".as_ref())?;
//!
//! // Fetch a region as raw ASCII bytes
//! let seq = reader.fetch_seq("chr1",
//!     Pos::<Zero>::new(1000).unwrap(),
//!     Pos::<Zero>::new(2000).unwrap(),
//! )?;
//!
//! // Convert to Base enum for downstream use
//! let mut bases = Base::from_ascii_vec(seq);
//! // bases is now Vec<Base> with A/C/G/T/Unknown variants
//! # Ok(())
//! # }
//! ```

pub mod bam;
pub mod cram;
pub mod fasta;
pub mod reader;
pub mod sam;
pub(crate) mod utils;
pub mod vcf;

pub use bam::{BaiError, BamError, BamHeaderError, BgzfError};
pub use cram::{CramError, CramIndexError};
pub use fasta::{FaiEntryError, FaiError, FastaError, GziError};
pub use reader::{FormatDetectionError, IndexedReader, ReaderError, Readers};
pub use sam::{SamError, SamRecordError};
pub use vcf::{AllelesError, VcfError, VcfHeaderError};
