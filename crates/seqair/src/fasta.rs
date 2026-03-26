//! Reference sequence access. Open a FASTA file with [`IndexedFastaReader`] to fetch
//! subsequences by name and coordinate range. Supports plain and bgzf-compressed FASTA.

mod gzi;
mod index;
mod reader;

pub use gzi::{BlockLocation, GziError, GziIndex};
pub use index::{FaiEntry, FaiEntryError, FaiError, FastaIndex};
pub use reader::{FastaError, IndexedFastaReader};
