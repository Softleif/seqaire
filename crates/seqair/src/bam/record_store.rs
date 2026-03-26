//! Slab-based BAM record storage with zero per-record heap allocation.
//!
//! All variable-length data is packed into contiguous byte buffers:
//! - **Name slab**: read names (qnames), accessed during dedup mate detection
//! - **Bases slab**: decoded `Base` values per record, accessed per-position in pileup
//! - **Data slab**: cigar + qual + aux per record, accessed during pileup construction

use seqair_types::Base;

use super::{cigar, record::DecodeError, seq};

/// Compact BAM record with offsets into the store's slabs.
// r[impl record_store.push_raw+2]
pub struct SlimRecord {
    pub pos: i64,
    pub end_pos: i64,
    pub flags: u16,
    pub n_cigar_ops: u16,
    pub mapq: u8,
    pub seq_len: u32,
    pub matching_bases: u32,
    pub indel_bases: u32,
    /// Offset into the name slab.
    name_off: u32,
    /// Length of the qname in the name slab.
    name_len: u16,
    /// Offset into the bases slab (seq_len Base values).
    bases_off: u32,
    /// Offset into the data slab (start of [cigar|qual|aux]).
    data_off: u32,
    /// Length of aux data in the data slab.
    aux_len: u32,
}

impl SlimRecord {
    fn cigar_off(&self) -> usize {
        self.data_off as usize
    }

    fn cigar_len(&self) -> usize {
        self.n_cigar_ops as usize * 4
    }

    fn qual_off(&self) -> usize {
        self.cigar_off() + self.cigar_len()
    }

    fn aux_off(&self) -> usize {
        self.qual_off() + self.seq_len as usize
    }
}

// r[impl record_store.push_raw+2]
// r[impl record_store.field_access]
// r[impl record_store.clear+2]
// r[impl record_store.region_scoped]
// r[impl record_store.capacity]
// r[impl record_store.no_rc]
// r[impl base_decode.slab]
/// Slab-based storage for BAM records.
///
/// Four contiguous buffers hold all data for a region:
/// - `records`: compact fixed-size structs
/// - `names`: packed qnames (for dedup)
/// - `bases`: decoded `Base` values per record (for per-position base lookup)
/// - `data`: packed [cigar|qual|aux] per record (for pileup construction)
pub struct RecordStore {
    records: Vec<SlimRecord>,
    names: Vec<u8>,
    bases: Vec<Base>,
    data: Vec<u8>,
}

impl RecordStore {
    pub fn new() -> Self {
        Self { records: Vec::new(), names: Vec::new(), bases: Vec::new(), data: Vec::new() }
    }

    // r[impl perf.arena_capacity_hint+2]
    /// Pre-allocate based on estimated compressed bytes for the region.
    /// Assumes ~3:1 compression ratio, ~400 bytes per record, ~150bp avg read.
    pub fn with_byte_hint(compressed_bytes: usize) -> Self {
        let uncompressed_est = compressed_bytes.saturating_mul(3);
        let record_count_est = (uncompressed_est / 400).max(64);
        Self {
            records: Vec::with_capacity(record_count_est),
            names: Vec::with_capacity(record_count_est * 25),
            bases: Vec::with_capacity(record_count_est * 150),
            data: Vec::with_capacity(uncompressed_est),
        }
    }

    /// Decode a raw BAM record and append to the store.
    ///
    /// Variable-length data is written directly into the slabs.
    pub fn push_raw(&mut self, raw: &[u8]) -> Result<u32, DecodeError> {
        if raw.len() < 32 {
            return Err(DecodeError::TooShort { len: raw.len() });
        }

        let idx = self.records.len() as u32;

        // Parse fixed header (same layout as BamRecord::decode)
        debug_assert!(raw.len() >= 32, "raw record too short for fixed fields: {}", raw.len());
        let tid = i32::from_le_bytes(read4(raw, 0));
        let pos = i64::from(i32::from_le_bytes(read4(raw, 4)));
        #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
        let name_len = raw[8] as usize;
        #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
        let mapq = raw[9];
        let n_cigar_ops = u16::from_le_bytes(read2(raw, 12));
        let flags = u16::from_le_bytes(read2(raw, 14));
        let seq_len = u32::from_le_bytes(read4(raw, 16));
        let _ = tid; // stored for filtering, but SlimRecord doesn't need it (filtered before push)

        let cigar_bytes = usize::from(n_cigar_ops) * 4;
        let seq_bytes = (seq_len as usize).div_ceil(2);

        let var_start = 32 + name_len;
        let cigar_end = var_start + cigar_bytes;
        let seq_end = cigar_end + seq_bytes;
        let qual_end = seq_end + seq_len as usize;

        if raw.len() < qual_end {
            return Err(DecodeError::TooShort { len: raw.len() });
        }

        // All slice bounds ≤ qual_end ≤ raw.len()
        debug_assert!(qual_end <= raw.len(), "qual_end overrun: {qual_end} > {}", raw.len());
        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let qname_raw = &raw[32..var_start];
        let qname_actual_len = qname_raw.iter().position(|&b| b == 0).unwrap_or(qname_raw.len());

        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let cigar_slice = &raw[var_start..cigar_end];
        let end_pos = super::record::compute_end_pos(pos, cigar_slice);
        let (matching_bases, indel_bases) = cigar::calc_matches_indels(cigar_slice);

        // --- Write into name slab ---
        let name_off = self.names.len() as u32;
        #[allow(clippy::indexing_slicing, reason = "qname_actual_len ≤ qname_raw.len()")]
        self.names.extend_from_slice(&qname_raw[..qname_actual_len]);

        // --- Write into bases slab (4-bit → Base via SIMD) ---
        let bases_off = self.bases.len() as u32;
        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        let packed_seq = &raw[cigar_end..seq_end];
        let decoded_bases = seq::decode_bases(packed_seq, seq_len as usize);
        self.bases.extend_from_slice(&decoded_bases);

        // --- Write into data slab [cigar|qual|aux] ---
        let data_off = self.data.len() as u32;

        // Cigar bytes (raw packed u32s)
        self.data.extend_from_slice(cigar_slice);

        // Quality scores
        #[allow(clippy::indexing_slicing, reason = "all bounds ≤ qual_end ≤ raw.len()")]
        self.data.extend_from_slice(&raw[seq_end..qual_end]);

        // Aux data (everything after qual)
        #[allow(clippy::indexing_slicing, reason = "qual_end ≤ raw.len()")]
        let aux_slice = &raw[qual_end..];
        self.data.extend_from_slice(aux_slice);

        self.records.push(SlimRecord {
            pos,
            end_pos,
            flags,
            n_cigar_ops,
            mapq,
            seq_len,
            matching_bases,
            indel_bases,
            name_off,
            name_len: qname_actual_len as u16,
            bases_off,
            data_off,
            aux_len: aux_slice.len() as u32,
        });

        Ok(idx)
    }

    // r[impl unified.record_store_push]
    // r[impl record_store.push_fields]
    // r[impl unified.push_fields_equivalence]
    /// Append a record from pre-parsed fields (for SAM/CRAM readers).
    ///
    /// Writes directly into the slabs without going through BAM binary encoding.
    /// CIGAR must be in BAM packed u32 format (`len << 4 | op`).
    #[allow(clippy::too_many_arguments)]
    pub fn push_fields(
        &mut self,
        pos: i64,
        end_pos: i64,
        flags: u16,
        mapq: u8,
        matching_bases: u32,
        indel_bases: u32,
        qname: &[u8],
        cigar_packed: &[u8],
        bases: &[Base],
        qual: &[u8],
        aux: &[u8],
    ) -> u32 {
        let idx = self.records.len() as u32;
        let n_cigar_ops = (cigar_packed.len() / 4) as u16;
        let seq_len = bases.len() as u32;

        // Name slab
        let name_off = self.names.len() as u32;
        self.names.extend_from_slice(qname);

        // Bases slab
        let bases_off = self.bases.len() as u32;
        self.bases.extend_from_slice(bases);

        // Data slab [cigar|qual|aux]
        let data_off = self.data.len() as u32;
        self.data.extend_from_slice(cigar_packed);
        self.data.extend_from_slice(qual);
        self.data.extend_from_slice(aux);

        self.records.push(SlimRecord {
            pos,
            end_pos,
            flags,
            n_cigar_ops,
            mapq,
            seq_len,
            matching_bases,
            indel_bases,
            name_off,
            name_len: qname.len() as u16,
            bases_off,
            data_off,
            aux_len: aux.len() as u32,
        });

        idx
    }

    // --- Accessors ---

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    #[allow(clippy::indexing_slicing, reason = "idx is always a valid index returned by push_raw")]
    pub fn record(&self, idx: u32) -> &SlimRecord {
        debug_assert!(
            (idx as usize) < self.records.len(),
            "record idx {idx} out of bounds (len={})",
            self.records.len()
        );
        &self.records[idx as usize]
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn qname(&self, idx: u32) -> &[u8] {
        let rec = self.record(idx);
        let start = rec.name_off as usize;
        let end = start + rec.name_len as usize;
        debug_assert!(end <= self.names.len(), "qname slab overrun: {end} > {}", self.names.len());
        &self.names[start..end]
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn cigar(&self, idx: u32) -> &[u8] {
        let rec = self.record(idx);
        let end = rec.cigar_off() + rec.cigar_len();
        debug_assert!(end <= self.data.len(), "cigar slab overrun: {end} > {}", self.data.len());
        &self.data[rec.cigar_off()..end]
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn seq(&self, idx: u32) -> &[Base] {
        let rec = self.record(idx);
        let start = rec.bases_off as usize;
        let end = start + rec.seq_len as usize;
        debug_assert!(end <= self.bases.len(), "bases slab overrun: {end} > {}", self.bases.len());
        &self.bases[start..end]
    }

    pub fn seq_at(&self, idx: u32, pos: usize) -> Base {
        let rec = self.record(idx);
        self.bases.get(rec.bases_off as usize + pos).copied().unwrap_or(Base::Unknown)
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn qual(&self, idx: u32) -> &[u8] {
        let rec = self.record(idx);
        let end = rec.qual_off() + rec.seq_len as usize;
        debug_assert!(end <= self.data.len(), "qual slab overrun: {end} > {}", self.data.len());
        &self.data[rec.qual_off()..end]
    }

    #[allow(clippy::indexing_slicing, reason = "offsets written by push_raw; within slab bounds")]
    pub fn aux(&self, idx: u32) -> &[u8] {
        let rec = self.record(idx);
        let end = rec.aux_off() + rec.aux_len as usize;
        debug_assert!(end <= self.data.len(), "aux slab overrun: {end} > {}", self.data.len());
        &self.data[rec.aux_off()..end]
    }

    pub fn clear(&mut self) {
        self.records.clear();
        self.names.clear();
        self.bases.clear();
        self.data.clear();
    }

    pub fn records_capacity(&self) -> usize {
        self.records.capacity()
    }

    pub fn names_capacity(&self) -> usize {
        self.names.capacity()
    }

    pub fn bases_capacity(&self) -> usize {
        self.bases.capacity()
    }

    pub fn data_capacity(&self) -> usize {
        self.data.capacity()
    }
}

impl Default for RecordStore {
    fn default() -> Self {
        Self::new()
    }
}

// offset + 3 < buf.len() ensured by the caller's bounds checks (raw.len() >= 32).
#[allow(clippy::indexing_slicing, reason = "offset + 3 < buf.len() ensured by caller")]
fn read4(buf: &[u8], offset: usize) -> [u8; 4] {
    debug_assert!(
        offset + 3 < buf.len(),
        "read4 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [buf[offset], buf[offset + 1], buf[offset + 2], buf[offset + 3]]
}

#[allow(clippy::indexing_slicing, reason = "offset + 1 < buf.len() ensured by caller")]
fn read2(buf: &[u8], offset: usize) -> [u8; 2] {
    debug_assert!(
        offset + 1 < buf.len(),
        "read2 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [buf[offset], buf[offset + 1]]
}
