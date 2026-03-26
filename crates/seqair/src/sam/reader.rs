//! Open and query bgzf-compressed SAM files. [`IndexedSamReader`] uses a tabix index to fetch
//! records for a region into a [`RecordStore`], with the same API as [`crate::bam::reader::IndexedBamReader`].

use crate::bam::{
    BaiError, BamHeader, BamHeaderError, BamIndex, cigar, flags::FLAG_UNMAPPED,
    record::compute_end_pos, record_store::RecordStore, region_buf::RegionBuf,
};
use seqair_types::Base;
use std::{
    fs::File,
    path::{Path, PathBuf},
    sync::Arc,
};
use tracing::instrument;

#[derive(Debug, thiserror::Error)]
pub enum SamRecordError {
    #[error("expected >= 11 TAB-separated fields, got {found}")]
    TooFewFields { found: usize },

    #[error("invalid FLAG field: {value:?}")]
    InvalidFlag { value: Box<[u8]> },

    #[error("RNAME is not valid UTF-8: {value:?}")]
    InvalidRname { value: Box<[u8]> },

    #[error("invalid POS field: {value:?}")]
    InvalidPos { value: Box<[u8]> },

    #[error("invalid MAPQ field: {value:?}")]
    InvalidMapq { value: Box<[u8]> },

    #[error("invalid CIGAR operation length: {value:?}")]
    InvalidCigarLength { value: Box<[u8]> },

    #[error("unknown CIGAR operation: {op}")]
    UnknownCigarOp { op: char },

    #[error("SAM header is not valid UTF-8")]
    HeaderNotUtf8 { source: std::string::FromUtf8Error },
}

#[derive(Debug, thiserror::Error)]
pub enum SamError {
    #[error("I/O error opening {path}")]
    Open { path: PathBuf, source: std::io::Error },

    #[error("SAM header error")]
    Header {
        #[from]
        source: BamHeaderError,
    },

    #[error("tabix index error")]
    Index {
        #[from]
        source: BaiError,
    },

    #[error("tabix index not found for {sam_path} (tried .tbi and .bai)")]
    IndexNotFound { sam_path: PathBuf },

    #[error(
        "found CSI index at {path} but CSI indexes are not yet supported. \
         Re-index with `tabix -p sam` to create a .tbi index instead.",
        path = path.display()
    )]
    CsiNotSupported { path: PathBuf },

    #[error(transparent)]
    MalformedRecord {
        #[from]
        source: SamRecordError,
    },

    #[error("BGZF error")]
    Bgzf {
        #[from]
        source: crate::bam::BgzfError,
    },

    #[error(
        "plain (uncompressed) SAM cannot be indexed. \
         Compress with `bgzip {path}` then index with `tabix -p sam {path}.gz`.",
        path = path.display()
    )]
    UncompressedSam { path: PathBuf },
}

pub struct SamShared {
    index: BamIndex,
    header: BamHeader,
    sam_path: PathBuf,
}

pub struct IndexedSamReader {
    bulk_reader: File,
    shared: Arc<SamShared>,
}

impl std::fmt::Debug for IndexedSamReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedSamReader").field("sam_path", &self.shared.sam_path).finish()
    }
}

impl IndexedSamReader {
    // r[impl sam.reader.open]
    #[instrument(level = "debug", fields(path = %path.display()), err)]
    pub fn open(path: &Path) -> Result<Self, SamError> {
        // r[impl sam.reader.reject_uncompressed]
        {
            let mut f = File::open(path)
                .map_err(|source| SamError::Open { path: path.to_path_buf(), source })?;
            let mut magic = [0u8; 1];
            use std::io::Read;
            if f.read_exact(&mut magic).is_ok() && magic[0] == b'@' {
                return Err(SamError::UncompressedSam { path: path.to_path_buf() });
            }
        }

        // Parse SAM header from BGZF stream
        let mut bgzf = crate::bam::bgzf::BgzfReader::open(path)?;
        let header_text = read_sam_header(&mut bgzf)?;
        let header = BamHeader::from_sam_text(&header_text)?;
        // r[impl unified.sort_order]
        header.validate_sort_order()?;

        // Find tabix index
        let index = find_tabix_index(path)?;

        let bulk_file = File::open(path)
            .map_err(|source| SamError::Open { path: path.to_path_buf(), source })?;

        Ok(IndexedSamReader {
            bulk_reader: bulk_file,
            shared: Arc::new(SamShared { index, header, sam_path: path.to_path_buf() }),
        })
    }

    pub fn fork(&self) -> Result<Self, SamError> {
        let bulk_file = File::open(&self.shared.sam_path)
            .map_err(|source| SamError::Open { path: self.shared.sam_path.clone(), source })?;

        Ok(IndexedSamReader { bulk_reader: bulk_file, shared: Arc::clone(&self.shared) })
    }

    pub fn shared(&self) -> &Arc<SamShared> {
        &self.shared
    }

    pub fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    // r[impl sam.reader.fetch_into]
    // r[impl sam.reader.sorted_order]
    #[instrument(level = "debug", skip(self, store), fields(tid, start, end), err)]
    // r[impl sam.perf.bulk_read]
    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: u64,
        end: u64,
        store: &mut RecordStore,
    ) -> Result<usize, SamError> {
        store.clear();

        let chunks = self.shared.index.query(tid, start, end);
        if chunks.is_empty() {
            return Ok(0);
        }

        let mut region = RegionBuf::load(&mut self.bulk_reader, &chunks)?;

        let start_i64 = start as i64;
        let end_i64 = end as i64;
        let tid_i32 = tid as i32;

        // Buffer for accumulating lines that span BGZF block boundaries
        let mut line_buf = Vec::with_capacity(1024);
        // Scratch buffers reused across records
        let mut cigar_buf = Vec::with_capacity(256);
        let mut bases_buf = Vec::with_capacity(256);
        let mut qual_buf = Vec::with_capacity(256);
        let mut aux_buf = Vec::with_capacity(256);

        for chunk in &chunks {
            region.seek_virtual(chunk.begin)?;
            line_buf.clear();

            // TODO(perf): r[sam.perf.text_parsing] — use memchr for TAB/newline scanning instead of byte-at-a-time
            loop {
                let current_voff = region.virtual_offset();
                if current_voff >= chunk.end {
                    break;
                }

                // Read one byte at a time, accumulating into line_buf until newline
                let byte = match region.read_byte() {
                    Ok(b) => b,
                    Err(_) => break,
                };

                if byte == b'\n' {
                    // Strip \r for Windows line endings
                    if line_buf.last() == Some(&b'\r') {
                        line_buf.pop();
                    }

                    // r[impl sam.edge.empty_lines]
                    if line_buf.is_empty() || line_buf.first() == Some(&b'@') {
                        line_buf.clear();
                        continue;
                    }

                    // Parse the SAM line
                    if let Some(()) = parse_sam_line(
                        &line_buf,
                        &self.shared.header,
                        tid_i32,
                        start_i64,
                        end_i64,
                        store,
                        &mut cigar_buf,
                        &mut bases_buf,
                        &mut qual_buf,
                        &mut aux_buf,
                    )? {
                        // record was added to store
                    }

                    line_buf.clear();
                } else {
                    line_buf.push(byte);
                }
            }

            // Handle partial line at end of chunk (line without trailing newline)
            if !line_buf.is_empty() && line_buf.first() != Some(&b'@') {
                if line_buf.last() == Some(&b'\r') {
                    line_buf.pop();
                }
                if !line_buf.is_empty() {
                    parse_sam_line(
                        &line_buf,
                        &self.shared.header,
                        tid_i32,
                        start_i64,
                        end_i64,
                        store,
                        &mut cigar_buf,
                        &mut bases_buf,
                        &mut qual_buf,
                        &mut aux_buf,
                    )?;
                }
                line_buf.clear();
            }
        }

        Ok(store.len())
    }
}

// r[impl sam.record.parse]
#[allow(clippy::too_many_arguments)]
fn parse_sam_line(
    line: &[u8],
    header: &BamHeader,
    tid_filter: i32,
    start: i64,
    end: i64,
    store: &mut RecordStore,
    cigar_buf: &mut Vec<u8>,
    bases_buf: &mut Vec<Base>,
    qual_buf: &mut Vec<u8>,
    aux_buf: &mut Vec<u8>,
) -> Result<Option<()>, SamError> {
    let fields: Vec<&[u8]> = line.splitn(12, |&b| b == b'\t').collect();
    if fields.len() < 11 {
        return Err(SamRecordError::TooFewFields { found: fields.len() }.into());
    }

    // Field 2: FLAG
    let flag_field = fields.get(1).copied().unwrap_or(b"0");
    let flags = parse_u16(flag_field)
        .ok_or_else(|| SamRecordError::InvalidFlag { value: flag_field.into() })?;

    // Skip unmapped
    if flags & FLAG_UNMAPPED != 0 {
        return Ok(None);
    }

    // r[impl sam.edge.rname_star]
    let rname = fields.get(2).copied().unwrap_or(b"*");
    if rname == b"*" {
        return Ok(None);
    }

    let rname_str = std::str::from_utf8(rname)
        .map_err(|_| SamRecordError::InvalidRname { value: rname.into() })?;
    let rec_tid = match header.tid(rname_str) {
        Some(t) => t as i32,
        None => return Ok(None), // unknown contig, skip
    };

    if rec_tid != tid_filter {
        return Ok(None);
    }

    // r[impl sam.record.coordinate_conversion]
    let pos_field = fields.get(3).copied().unwrap_or(b"0");
    let pos_1based = parse_i64(pos_field)
        .ok_or_else(|| SamRecordError::InvalidPos { value: pos_field.into() })?;
    let pos = pos_1based - 1;

    // Field 5: MAPQ
    let mapq_field = fields.get(4).copied().unwrap_or(b"0");
    let mapq = parse_u8(mapq_field)
        .ok_or_else(|| SamRecordError::InvalidMapq { value: mapq_field.into() })?;

    // Field 6: CIGAR
    let cigar_str = fields.get(5).copied().unwrap_or(b"*");
    cigar_buf.clear();
    let cigar_available = parse_cigar(cigar_str, cigar_buf)?;

    // Compute end_pos from packed CIGAR
    let end_pos = if cigar_available { compute_end_pos(pos, cigar_buf) } else { pos };

    // r[impl sam.reader.overlap_filter]
    if pos > end || end_pos < start {
        return Ok(None);
    }

    // Field 1: QNAME
    let qname = fields.first().copied().unwrap_or(b"*");

    // r[impl sam.record.seq_decode]
    // r[impl sam.edge.missing_seq]
    let seq_field = fields.get(9).copied().unwrap_or(b"*");
    if seq_field == b"*" {
        bases_buf.clear();
    } else {
        *bases_buf = Base::from_ascii_vec(seq_field.to_vec());
    }

    // r[impl sam.record.qual_decode]
    // r[impl sam.edge.missing_qual]
    let qual_field = fields.get(10).copied().unwrap_or(b"*");
    qual_buf.clear();
    if qual_field == b"*" {
        qual_buf.resize(bases_buf.len(), 0xFF);
    } else {
        qual_buf.reserve(qual_field.len());
        for &b in qual_field {
            qual_buf.push(b.wrapping_sub(33));
        }
    }

    // Optional fields (12+): aux tags
    aux_buf.clear();
    if fields.len() > 11 {
        // The 12th element contains everything after the 11th TAB (because of splitn(12,...))
        let aux_text = fields.get(11).copied().unwrap_or(b"");
        parse_aux_tags(aux_text, aux_buf)?;
    }

    let (matching_bases, indel_bases) = cigar::calc_matches_indels(cigar_buf);

    store.push_fields(
        pos,
        end_pos,
        flags,
        mapq,
        matching_bases,
        indel_bases,
        qname,
        cigar_buf,
        bases_buf,
        qual_buf,
        aux_buf,
    );

    Ok(Some(()))
}

// r[impl sam.record.cigar_parse]
// r[impl sam.edge.long_cigar]
fn parse_cigar(cigar_str: &[u8], buf: &mut Vec<u8>) -> Result<bool, SamError> {
    if cigar_str == b"*" {
        return Ok(false);
    }

    let mut num_start = 0;
    for (i, &b) in cigar_str.iter().enumerate() {
        if b.is_ascii_digit() {
            continue;
        }

        let len_bytes = cigar_str.get(num_start..i).unwrap_or(b"");
        let len = parse_u32(len_bytes)
            .ok_or_else(|| SamRecordError::InvalidCigarLength { value: len_bytes.into() })?;

        let op = match b {
            b'M' => 0u32,
            b'I' => 1,
            b'D' => 2,
            b'N' => 3,
            b'S' => 4,
            b'H' => 5,
            b'P' => 6,
            b'=' => 7,
            b'X' => 8,
            _ => {
                return Err(SamRecordError::UnknownCigarOp { op: b as char }.into());
            }
        };

        let packed = (len << 4) | op;
        buf.extend_from_slice(&packed.to_le_bytes());
        num_start = i + 1;
    }

    Ok(true)
}

// r[impl sam.record.aux_tags]
fn parse_aux_tags(text: &[u8], buf: &mut Vec<u8>) -> Result<(), SamError> {
    for field in text.split(|&b| b == b'\t') {
        if field.len() < 5 {
            continue; // need at least "XX:T:V"
        }
        // field[0..2] = tag name, field[2] = ':', field[3] = type, field[4] = ':'
        if field.get(2) != Some(&b':') || field.get(4) != Some(&b':') {
            continue;
        }
        let tag = field.get(..2).unwrap_or(b"??");
        let type_char = field.get(3).copied().unwrap_or(b'Z');
        let value = field.get(5..).unwrap_or(b"");

        buf.extend_from_slice(tag);
        match type_char {
            b'A' => {
                buf.push(b'A');
                buf.push(value.first().copied().unwrap_or(b' '));
            }
            b'i' => {
                let val = parse_i64(value).unwrap_or(0);
                serialize_bam_int(buf, val);
            }
            b'f' => {
                buf.push(b'f');
                let f: f32 =
                    std::str::from_utf8(value).ok().and_then(|s| s.parse().ok()).unwrap_or(0.0);
                buf.extend_from_slice(&f.to_le_bytes());
            }
            b'Z' => {
                buf.push(b'Z');
                buf.extend_from_slice(value);
                buf.push(0);
            }
            b'H' => {
                buf.push(b'H');
                buf.extend_from_slice(value);
                buf.push(0);
            }
            b'B' => {
                buf.push(b'B');
                // B:type,v1,v2,...
                let subtype = value.first().copied().unwrap_or(b'C');
                buf.push(subtype);
                let values_str = value.get(2..).unwrap_or(b"");
                let values: Vec<&[u8]> = values_str.split(|&b| b == b',').collect();
                buf.extend_from_slice(&(values.len() as u32).to_le_bytes());
                for v in &values {
                    match subtype {
                        b'c' | b'C' => {
                            let n = parse_i64(v).unwrap_or(0) as u8;
                            buf.push(n);
                        }
                        b's' | b'S' => {
                            let n = parse_i64(v).unwrap_or(0) as u16;
                            buf.extend_from_slice(&n.to_le_bytes());
                        }
                        b'i' | b'I' => {
                            let n = parse_i64(v).unwrap_or(0) as u32;
                            buf.extend_from_slice(&n.to_le_bytes());
                        }
                        b'f' => {
                            let f: f32 = std::str::from_utf8(v)
                                .ok()
                                .and_then(|s| s.parse().ok())
                                .unwrap_or(0.0);
                            buf.extend_from_slice(&f.to_le_bytes());
                        }
                        _ => {}
                    }
                }
            }
            _ => {
                // Unknown type — store as Z-string
                buf.push(b'Z');
                buf.extend_from_slice(value);
                buf.push(0);
            }
        }
    }
    Ok(())
}

/// Serialize an integer value using the smallest BAM integer type.
fn serialize_bam_int(buf: &mut Vec<u8>, val: i64) {
    if (-128..=127).contains(&val) {
        buf.push(b'c');
        buf.push(val as u8);
    } else if (0..=255).contains(&val) {
        buf.push(b'C');
        buf.push(val as u8);
    } else if (-32768..=32767).contains(&val) {
        buf.push(b's');
        buf.extend_from_slice(&(val as i16).to_le_bytes());
    } else if (0..=65535).contains(&val) {
        buf.push(b'S');
        buf.extend_from_slice(&(val as u16).to_le_bytes());
    } else if (-2_147_483_648..=2_147_483_647).contains(&val) {
        buf.push(b'i');
        buf.extend_from_slice(&(val as i32).to_le_bytes());
    } else {
        buf.push(b'I');
        buf.extend_from_slice(&(val as u32).to_le_bytes());
    }
}

// r[impl sam.header.parse]
// r[impl sam.header.sq_required]
// r[impl sam.header.preserve_full_text]
// r[impl sam.edge.line_spanning_blocks]
fn read_sam_header(bgzf: &mut crate::bam::bgzf::BgzfReader) -> Result<String, SamError> {
    let mut header = Vec::with_capacity(4096);
    let mut line = Vec::with_capacity(256);

    while let Ok(b) = bgzf.read_byte() {
        if b == b'\n' {
            if line.first() == Some(&b'@') {
                header.extend_from_slice(&line);
                header.push(b'\n');
                line.clear();
            } else {
                // First non-header line — stop
                break;
            }
        } else {
            line.push(b);
        }
    }

    String::from_utf8(header).map_err(|source| SamRecordError::HeaderNotUtf8 { source }.into())
}

// r[impl sam.index.tabix]
// r[impl unified.detect_index]
// r[impl unified.detect_error]
// TODO: r[sam.index.locate+2] — CSI not yet supported, only TBI and BAI
fn find_tabix_index(sam_path: &Path) -> Result<BamIndex, SamError> {
    // Try .gz.tbi first
    let tbi_path = sam_path.with_extension("gz.tbi");
    if tbi_path.exists() {
        return BamIndex::from_tabix_path(&tbi_path).map_err(SamError::from);
    }

    // Try appending .tbi
    let mut path_with_tbi = sam_path.to_path_buf();
    let mut name = path_with_tbi.file_name().unwrap_or_default().to_os_string();
    name.push(".tbi");
    path_with_tbi.set_file_name(name);
    if path_with_tbi.exists() {
        return BamIndex::from_tabix_path(&path_with_tbi).map_err(SamError::from);
    }

    // Try .gz.bai (BAI format, not tabix)
    let bai_path = sam_path.with_extension("gz.bai");
    if bai_path.exists() {
        return BamIndex::from_path(&bai_path).map_err(SamError::from);
    }

    // Try appending .bai
    let mut path_with_bai = sam_path.to_path_buf();
    let mut name = path_with_bai.file_name().unwrap_or_default().to_os_string();
    name.push(".bai");
    path_with_bai.set_file_name(name);
    if path_with_bai.exists() {
        return BamIndex::from_path(&path_with_bai).map_err(SamError::from);
    }

    // Check if a CSI index exists — not yet supported
    let csi_path = sam_path.with_extension("gz.csi");
    if csi_path.exists() {
        return Err(SamError::CsiNotSupported { path: csi_path });
    }

    let mut path_with_csi = sam_path.to_path_buf();
    let mut name = path_with_csi.file_name().unwrap_or_default().to_os_string();
    name.push(".csi");
    path_with_csi.set_file_name(name);
    if path_with_csi.exists() {
        return Err(SamError::CsiNotSupported { path: path_with_csi });
    }

    Err(SamError::IndexNotFound { sam_path: sam_path.to_path_buf() })
}

// --- Integer parsing from byte slices (no String allocation) ---

fn parse_u8(bytes: &[u8]) -> Option<u8> {
    let mut val = 0u16;
    for &b in bytes {
        if !b.is_ascii_digit() {
            return None;
        }
        val = val.checked_mul(10)?.checked_add(u16::from(b - b'0'))?;
    }
    u8::try_from(val).ok()
}

fn parse_u16(bytes: &[u8]) -> Option<u16> {
    let mut val = 0u32;
    for &b in bytes {
        if !b.is_ascii_digit() {
            return None;
        }
        val = val.checked_mul(10)?.checked_add(u32::from(b - b'0'))?;
    }
    u16::try_from(val).ok()
}

fn parse_u32(bytes: &[u8]) -> Option<u32> {
    let mut val = 0u64;
    for &b in bytes {
        if !b.is_ascii_digit() {
            return None;
        }
        val = val.checked_mul(10)?.checked_add(u64::from(b - b'0'))?;
    }
    u32::try_from(val).ok()
}

fn parse_i64(bytes: &[u8]) -> Option<i64> {
    if bytes.is_empty() {
        return None;
    }
    let (negative, digits) =
        if bytes.first() == Some(&b'-') { (true, bytes.get(1..)?) } else { (false, bytes) };
    let mut val = 0i64;
    for &b in digits {
        if !b.is_ascii_digit() {
            return None;
        }
        val = val.checked_mul(10)?.checked_add(i64::from(b - b'0'))?;
    }
    if negative { Some(-val) } else { Some(val) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::record_store::RecordStore;
    use seqair_types::Base;
    use std::io::Write;

    fn make_header() -> BamHeader {
        BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:1000\n").expect("failed to build test header")
    }

    fn make_store_and_bufs() -> (RecordStore, Vec<u8>, Vec<Base>, Vec<u8>, Vec<u8>) {
        (RecordStore::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new())
    }

    fn call_parse(line: &[u8], header: &BamHeader) -> Result<Option<()>, SamError> {
        let (mut store, mut cigar_buf, mut bases_buf, mut qual_buf, mut aux_buf) =
            make_store_and_bufs();
        parse_sam_line(
            line,
            header,
            0,
            0,
            10000,
            &mut store,
            &mut cigar_buf,
            &mut bases_buf,
            &mut qual_buf,
            &mut aux_buf,
        )
    }

    #[test]
    fn error_too_few_fields() {
        let header = make_header();
        let line = b"read1\t0\tchr1";
        let err = call_parse(line, &header).expect_err("expected TooFewFields error");
        assert!(
            matches!(err, SamError::MalformedRecord { source: SamRecordError::TooFewFields { found } } if found < 11),
            "expected TooFewFields, got: {err}"
        );
    }

    #[test]
    fn error_invalid_flag() {
        let header = make_header();
        let line = b"r\tNOTANUM\tchr1\t100\t60\t10M\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected InvalidFlag error");
        assert!(
            matches!(err, SamError::MalformedRecord { source: SamRecordError::InvalidFlag { .. } }),
            "expected InvalidFlag, got: {err}"
        );
    }

    #[test]
    fn error_invalid_rname_not_utf8() {
        let header = make_header();
        // Build a line with invalid UTF-8 bytes in the RNAME field (field[2])
        let mut line = b"r\t0\t".to_vec();
        line.extend_from_slice(b"\xFF\xFE"); // invalid UTF-8
        line.extend_from_slice(b"\t100\t60\t10M\t*\t0\t0\tACGT\t~~~~");
        let err = call_parse(&line, &header).expect_err("expected InvalidRname error");
        assert!(
            matches!(
                err,
                SamError::MalformedRecord { source: SamRecordError::InvalidRname { .. } }
            ),
            "expected InvalidRname, got: {err}"
        );
    }

    #[test]
    fn error_invalid_pos() {
        let header = make_header();
        let line = b"r\t0\tchr1\tNOTANUM\t60\t10M\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected InvalidPos error");
        assert!(
            matches!(err, SamError::MalformedRecord { source: SamRecordError::InvalidPos { .. } }),
            "expected InvalidPos, got: {err}"
        );
    }

    #[test]
    fn error_invalid_mapq() {
        let header = make_header();
        // 999 overflows u8 (max 255)
        let line = b"r\t0\tchr1\t100\t999\t10M\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected InvalidMapq error");
        assert!(
            matches!(err, SamError::MalformedRecord { source: SamRecordError::InvalidMapq { .. } }),
            "expected InvalidMapq, got: {err}"
        );
    }

    #[test]
    fn error_invalid_cigar_length() {
        let header = make_header();
        // Length 4294967296 exceeds u32::MAX, causing parse_u32 to return None
        let line = b"r\t0\tchr1\t100\t60\t4294967296M\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected InvalidCigarLength error");
        assert!(
            matches!(
                err,
                SamError::MalformedRecord { source: SamRecordError::InvalidCigarLength { .. } }
            ),
            "expected InvalidCigarLength, got: {err}"
        );
    }

    #[test]
    fn error_unknown_cigar_op() {
        let header = make_header();
        // "10Z" — Z is not a valid CIGAR operation
        let line = b"r\t0\tchr1\t100\t60\t10Z\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected UnknownCigarOp error");
        assert!(
            matches!(
                err,
                SamError::MalformedRecord { source: SamRecordError::UnknownCigarOp { op: 'Z' } }
            ),
            "expected UnknownCigarOp {{ op: 'Z' }}, got: {err}"
        );
    }

    // TODO: HeaderNotUtf8 requires crafting a BGZF stream with invalid UTF-8 header bytes —
    // it is raised by read_sam_header() which reads from a BgzfReader, not from parse_sam_line.

    // r[verify sam.reader.reject_uncompressed]
    #[test]
    fn reject_uncompressed_sam() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let sam_path = dir.path().join("test.sam");
        {
            let mut f = File::create(&sam_path).expect("failed to create temp file");
            f.write_all(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n").expect("failed to write");
        }

        let err = IndexedSamReader::open(&sam_path).expect_err("should reject uncompressed SAM");
        assert!(
            matches!(err, SamError::UncompressedSam { .. }),
            "expected UncompressedSam, got: {err}"
        );
        let msg = err.to_string();
        assert!(msg.contains("bgzip"), "error message should mention bgzip: {msg}");
        assert!(msg.contains("test.sam"), "error message should contain file name: {msg}");
    }

    // r[verify sam.index.csi+2]
    #[test]
    fn csi_index_returns_clear_error() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let sam_path = dir.path().join("test.sam.gz");

        // Create a fake .csi file next to the SAM
        let csi_path = dir.path().join("test.sam.gz.csi");
        File::create(&csi_path).expect("failed to create .csi file");

        let err = find_tabix_index(&sam_path).expect_err("should reject CSI index");
        assert!(
            matches!(err, SamError::CsiNotSupported { .. }),
            "expected CsiNotSupported, got: {err}"
        );
        let msg = err.to_string();
        assert!(msg.contains("CSI"), "error message should mention CSI: {msg}");
        assert!(msg.contains("tabix"), "error message should mention tabix: {msg}");
    }
}
