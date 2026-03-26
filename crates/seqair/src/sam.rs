//! bgzf-compressed SAM reading. Use [`reader::IndexedSamReader`] to fetch records into a
//! [`crate::bam::RecordStore`] via a tabix index.

pub mod reader;
