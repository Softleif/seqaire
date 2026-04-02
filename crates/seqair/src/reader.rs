//! Format-agnostic entry point. Use [`Readers::open`] to open an alignment file (BAM/SAM/CRAM)
//! together with a FASTA reference, then call [`Readers::pileup`] to iterate columns or
//! [`Readers::fetch_base_seq`] for reference sequence. [`IndexedReader`] is the bare alignment
//! handle when no FASTA is needed. Both types are forkable for multi-threaded use.

use crate::{
    bam::{BamError, BamHeader, IndexedBamReader, pileup::PileupEngine, record_store::RecordStore},
    cram::reader::{CramError, IndexedCramReader},
    fasta::{FastaError, IndexedFastaReader},
    sam::reader::{IndexedSamReader, SamError},
};
use seqair_types::{Base, Pos, Zero};
use std::io::{Read, Seek};
use std::path::{Path, PathBuf};
use std::rc::Rc;
use tracing::instrument;

mod formats;
use formats::Format;
pub use formats::FormatDetectionError;

// r[impl io.non_exhaustive_enums]
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReaderError {
    #[error("BAM reader error")]
    Bam {
        #[from]
        source: BamError,
    },

    #[error("SAM reader error")]
    Sam {
        #[from]
        source: SamError,
    },

    #[error("CRAM reader error")]
    Cram {
        #[from]
        source: CramError,
    },

    #[error("I/O error reading {path}")]
    Open { path: PathBuf, source: std::io::Error },

    #[error(transparent)]
    Format {
        #[from]
        source: FormatDetectionError,
    },

    #[error("failed to open FASTA reference")]
    FastaOpen { source: FastaError },

    #[error("failed to fork FASTA reader")]
    FastaFork { source: FastaError },
}

/// Format-agnostic indexed alignment reader, generic over the I/O backend.
///
/// `R = File` for production (file-based I/O); `R = Cursor<Vec<u8>>` for fuzzing
/// (in-memory, no file I/O).
// r[impl unified.reader_enum]
#[non_exhaustive]
pub enum IndexedReader<R: Read + Seek = std::fs::File> {
    Bam(IndexedBamReader<R>),
    Sam(IndexedSamReader<R>),
    Cram(Box<IndexedCramReader<R>>),
}

impl<R: Read + Seek> std::fmt::Debug for IndexedReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bam(r) => r.fmt(f),
            Self::Sam(r) => r.fmt(f),
            Self::Cram(r) => r.fmt(f),
        }
    }
}

impl<R: Read + Seek> IndexedReader<R> {
    pub fn header(&self) -> &BamHeader {
        match self {
            Self::Bam(r) => r.header(),
            Self::Sam(r) => r.header(),
            Self::Cram(r) => r.header(),
        }
    }

    // r[impl unified.fetch_equivalence]
    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        match self {
            Self::Bam(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
            Self::Sam(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
            Self::Cram(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
        }
    }
}

// r[impl unified.reader_api]
impl IndexedReader<std::fs::File> {
    /// Open a BAM or bgzf-compressed SAM file, auto-detecting the format.
    ///
    /// CRAM files are detected but require a FASTA reference — use
    /// [`Readers::open`] instead. This method returns an error for CRAM
    /// with a message directing the user to provide a reference.
    // r[impl unified.readers_backward_compat]
    #[instrument(level = "debug", fields(path = %path.display()), err)]
    pub fn open(path: &Path) -> Result<Self, ReaderError> {
        match formats::detect(path)? {
            Format::Bam => {
                let reader = IndexedBamReader::open(path)?;
                Ok(IndexedReader::Bam(reader))
            }
            Format::Sam => {
                let reader = IndexedSamReader::open(path)?;
                Ok(IndexedReader::Sam(reader))
            }
            Format::Cram => Err(FormatDetectionError::CramRequiresFasta.into()),
        }
    }

    /// Open any format, with a FASTA path for CRAM support.
    fn open_with_fasta(path: &Path, fasta_path: &Path) -> Result<Self, ReaderError> {
        match formats::detect(path)? {
            Format::Bam => {
                let reader = IndexedBamReader::open(path)?;
                Ok(IndexedReader::Bam(reader))
            }
            Format::Sam => {
                let reader = IndexedSamReader::open(path)?;
                Ok(IndexedReader::Sam(reader))
            }
            Format::Cram => {
                let reader = IndexedCramReader::open(path, fasta_path)?;
                Ok(IndexedReader::Cram(Box::new(reader)))
            }
        }
    }

    // r[impl unified.fork_bam]
    // r[impl unified.fork_sam]
    // r[impl unified.fork_cram]
    pub fn fork(&self) -> Result<Self, ReaderError> {
        match self {
            Self::Bam(r) => Ok(Self::Bam(r.fork()?)),
            Self::Sam(r) => Ok(Self::Sam(r.fork()?)),
            Self::Cram(r) => Ok(Self::Cram(Box::new(r.fork()?))),
        }
    }
}

/// Alignment + reference reader bundle.
///
/// Bundles an [`IndexedReader`] (BAM/SAM/CRAM) with an [`IndexedFastaReader`]
/// so that CRAM has access to the reference it needs and all formats have
/// uniform open/fork/fetch semantics.
// r[impl unified.readers_struct]
pub struct Readers {
    alignment: IndexedReader,
    fasta: IndexedFastaReader,
    store: RecordStore,
    fasta_buf: Vec<u8>,
}

impl std::fmt::Debug for Readers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Readers").field("alignment", &self.alignment).finish()
    }
}

impl Readers {
    /// Open an alignment file (BAM/SAM/CRAM) and a FASTA reference.
    ///
    /// Auto-detects the alignment format. For CRAM, the FASTA path is passed
    /// to the CRAM reader for sequence reconstruction.
    // r[impl unified.readers_open]
    #[instrument(level = "debug", fields(alignment = %alignment_path.display(), fasta = %fasta_path.display()), err)]
    pub fn open(alignment_path: &Path, fasta_path: &Path) -> Result<Self, ReaderError> {
        let fasta = IndexedFastaReader::open(fasta_path)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        let alignment = IndexedReader::open_with_fasta(alignment_path, fasta_path)?;
        Ok(Readers { alignment, fasta, store: RecordStore::default(), fasta_buf: Vec::new() })
    }

    /// Fork both the alignment reader and the FASTA reader.
    // r[impl unified.readers_fork]
    pub fn fork(&self) -> Result<Self, ReaderError> {
        let alignment = self.alignment.fork()?;
        let fasta = self.fasta.fork().map_err(|source| ReaderError::FastaFork { source })?;
        Ok(Readers { alignment, fasta, store: RecordStore::default(), fasta_buf: Vec::new() })
    }

    // r[impl unified.readers_accessors]
    pub fn header(&self) -> &BamHeader {
        self.alignment.header()
    }

    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        self.alignment.fetch_into(tid, start, end, store)
    }

    /// Fetch records for a region and return a [`PileupEngine`] ready for iteration.
    ///
    /// Uses an internal [`RecordStore`] whose capacity is retained across calls.
    /// After iterating the engine, call [`recover_store`](Self::recover_store) to
    /// return the store for reuse. If not called, the next `pileup()` call
    /// allocates a fresh store (small perf hit, not a correctness issue).
    pub fn pileup(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
    ) -> Result<PileupEngine, ReaderError> {
        self.alignment.fetch_into(tid, start, end, &mut self.store)?;

        let store = std::mem::take(&mut self.store);
        Ok(PileupEngine::new(store, start, end))
    }

    /// Recover the [`RecordStore`] from a consumed [`PileupEngine`] for reuse.
    ///
    /// Call this after iteration is complete. The store retains its allocated
    /// capacity, avoiding ~39 MB of re-allocation on the next `pileup()` call.
    pub fn recover_store(&mut self, engine: &mut PileupEngine) {
        if let Some(store) = engine.take_store() {
            self.store = store;
        }
    }

    pub fn fasta(&self) -> &IndexedFastaReader {
        &self.fasta
    }

    pub fn fasta_mut(&mut self) -> &mut IndexedFastaReader {
        &mut self.fasta
    }

    // r[impl fasta.fetch.buffer_reuse]
    /// Fetch a reference sequence region and return it as `Rc<[Base]>`.
    ///
    /// Uses an internal buffer whose capacity is retained across calls,
    /// avoiding per-segment allocation. The returned `Rc<[Base]>` owns its
    /// own copy — the internal buffer is reused on the next call.
    pub fn fetch_base_seq(
        &mut self,
        name: &str,
        start: Pos<Zero>,
        stop: Pos<Zero>,
    ) -> Result<Rc<[Base]>, FastaError> {
        self.fasta.fetch_seq_into(name, start, stop, &mut self.fasta_buf)?;
        // Take the buffer so from_ascii_vec can reinterpret it in-place (safe),
        // then create the Rc (which copies). The buffer capacity is lost but
        // fasta_buf re-grows on the next call.
        let buf = std::mem::take(&mut self.fasta_buf);
        let bases = Base::from_ascii_vec(buf);
        Ok(Rc::from(bases))
    }

    pub fn alignment(&self) -> &IndexedReader {
        &self.alignment
    }

    pub fn alignment_mut(&mut self) -> &mut IndexedReader {
        &mut self.alignment
    }
}

/// Cursor-backed type alias for in-memory fuzzing.
#[cfg(feature = "fuzz")]
pub type CursorReader = IndexedReader<std::io::Cursor<Vec<u8>>>;

/// In-memory reader bundle for fuzzing. No file I/O — all data from byte slices.
/// Uses the same `IndexedReader` enum as production code, just with `Cursor` I/O.
#[cfg(feature = "fuzz")]
pub struct FuzzReaders {
    alignment: CursorReader,
    fasta: IndexedFastaReader<std::io::Cursor<Vec<u8>>>,
    store: RecordStore,
    fasta_buf: Vec<u8>,
}

#[cfg(feature = "fuzz")]
impl FuzzReaders {
    /// Build BAM-based readers from raw bytes: BAM + BAI + FASTA.gz + FAI + GZI.
    pub fn from_bam_bytes(
        bam_data: Vec<u8>,
        bai_data: &[u8],
        fasta_gz_data: Vec<u8>,
        fai_contents: &str,
        gzi_data: &[u8],
    ) -> Result<Self, ReaderError> {
        let bam = IndexedBamReader::from_bytes(bam_data, bai_data)?;
        let fasta = IndexedFastaReader::from_bytes(fasta_gz_data, fai_contents, gzi_data)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        Ok(FuzzReaders {
            alignment: CursorReader::Bam(bam),
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
        })
    }

    /// Build CRAM-based readers from raw bytes.
    pub fn from_cram_bytes(
        cram_data: Vec<u8>,
        crai_text: &str,
        fasta_gz_data: Vec<u8>,
        fai_contents: &str,
        gzi_data: &[u8],
    ) -> Result<Self, ReaderError> {
        let cram = IndexedCramReader::from_bytes(
            cram_data,
            crai_text.as_bytes(),
            fasta_gz_data.clone(),
            fai_contents,
            gzi_data,
        )?;
        let fasta = IndexedFastaReader::from_bytes(fasta_gz_data, fai_contents, gzi_data)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        Ok(FuzzReaders {
            alignment: CursorReader::Cram(Box::new(cram)),
            fasta,
            store: RecordStore::default(),
            fasta_buf: Vec::new(),
        })
    }

    pub fn header(&self) -> &BamHeader {
        self.alignment.header()
    }

    /// Full pileup pipeline: fetch_into → PileupEngine with reference sequence.
    pub fn pileup(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
    ) -> Result<PileupEngine, ReaderError> {
        self.alignment.fetch_into(tid, start, end, &mut self.store)?;
        let store = std::mem::take(&mut self.store);
        let mut engine = PileupEngine::new(store, start, end);

        // Try to set reference sequence
        if let Some(name) = self.alignment.header().target_name(tid) {
            let name = name.to_owned();
            if self.fasta.fetch_seq_into(&name, start, end, &mut self.fasta_buf).is_ok() {
                let buf = std::mem::take(&mut self.fasta_buf);
                let bases = Base::from_ascii_vec(buf);
                let ref_seq = crate::bam::pileup::RefSeq::new(Rc::from(bases), start);
                engine.set_reference_seq(ref_seq);
            }
        }

        Ok(engine)
    }

    pub fn recover_store(&mut self, engine: &mut PileupEngine) {
        if let Some(store) = engine.take_store() {
            self.store = store;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;
    use tempfile::NamedTempFile;

    #[test]
    fn cram_requires_fasta() {
        // A file starting with the CRAM magic bytes — IndexedReader::open (no FASTA)
        // must return CramRequiresFasta.
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"CRAM\x03\x00").expect("write");
        f.flush().expect("flush");
        let err = IndexedReader::open(f.path()).unwrap_err();
        assert!(
            matches!(err, ReaderError::Format { source: FormatDetectionError::CramRequiresFasta }),
            "unexpected error: {err}"
        );
    }
}
