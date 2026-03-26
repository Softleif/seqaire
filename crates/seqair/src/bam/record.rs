//! Decode a BAM record from raw bytes into a [`BamRecord`] with owned variable-length fields.
//! Used transiently during region loading; the pileup engine works from [`crate::bam::RecordStore`] instead.

use super::{
    flags::{
        BamFlags, FLAG_FIRST_IN_TEMPLATE, FLAG_REVERSE, FLAG_SECOND_IN_TEMPLATE, FLAG_UNMAPPED,
    },
    seq,
};

/// A decoded BAM record with owned variable-length data.
///
/// Wrapped in `Rc` for cheap sharing between the arena and pileup columns.
/// All accessor methods live directly on the struct — no separate `RecordRef`.
// r[impl bam.record.fields]
// r[impl bam.record.decode]
#[derive(Debug, Clone)]
pub struct BamRecord {
    pub pos: i64,
    pub end_pos: i64,
    pub tid: i32,
    pub seq_len: u32,
    pub flags: u16,
    pub n_cigar_ops: u16,
    pub mapq: u8,
    // r[impl perf.precompute_matches_indels]
    pub matching_bases: u32,
    pub indel_bases: u32,
    pub qname: Box<[u8]>,
    pub cigar: Box<[u8]>,
    pub seq: Box<[u8]>,
    pub qual: Box<[u8]>,
    pub aux: Box<[u8]>,
}

impl BamRecord {
    /// Decode from raw BAM bytes (after the 4-byte block_size prefix).
    pub fn decode(raw: &[u8]) -> Result<Self, DecodeError> {
        if raw.len() < 32 {
            return Err(DecodeError::TooShort { len: raw.len() });
        }

        debug_assert!(raw.len() >= 32, "raw record too short for fixed fields: {}", raw.len());
        let tid = i32::from_le_bytes(read4(raw, 0));
        let pos = i64::from(i32::from_le_bytes(read4(raw, 4)));
        // raw.len() >= 32, so raw[8] and raw[9] are in bounds
        #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
        let name_len = raw[8] as usize;
        #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
        let mapq = raw[9];
        let n_cigar_ops = u16::from_le_bytes(read2(raw, 12));
        let flags = u16::from_le_bytes(read2(raw, 14));
        let seq_len = u32::from_le_bytes(read4(raw, 16));
        let seq_len_usize = seq_len as usize;

        let cigar_bytes = usize::from(n_cigar_ops) * 4;
        let seq_bytes = seq_len_usize.div_ceil(2);

        let var_start = 32 + name_len;
        let cigar_end = var_start + cigar_bytes;
        let seq_end = cigar_end + seq_bytes;
        let qual_end = seq_end + seq_len_usize;

        if raw.len() < qual_end {
            return Err(DecodeError::TooShort { len: raw.len() });
        }

        // All slice bounds (32..var_start, var_start..cigar_end, etc.) are ≤ qual_end ≤ raw.len()
        debug_assert!(qual_end <= raw.len(), "qual_end overrun: {qual_end} > {}", raw.len());
        #[allow(
            clippy::indexing_slicing,
            reason = "all bounds ≤ qual_end ≤ raw.len() checked above"
        )]
        let qname_raw = &raw[32..var_start];
        let qname_actual_len = qname_raw.iter().position(|&b| b == 0).unwrap_or(qname_raw.len());

        #[allow(
            clippy::indexing_slicing,
            reason = "all bounds ≤ qual_end ≤ raw.len() checked above"
        )]
        let cigar_slice = &raw[var_start..cigar_end];
        let end_pos = compute_end_pos(pos, cigar_slice);
        let (matching_bases, indel_bases) = super::cigar::calc_matches_indels(cigar_slice);

        #[allow(
            clippy::indexing_slicing,
            reason = "all bounds ≤ qual_end ≤ raw.len() checked above"
        )]
        Ok(BamRecord {
            pos,
            end_pos,
            tid,
            seq_len,
            flags,
            n_cigar_ops,
            mapq,
            matching_bases,
            indel_bases,
            qname: qname_raw[..qname_actual_len].into(),
            cigar: cigar_slice.into(),
            // r[impl bam.record.seq_4bit]
            // r[impl bam.record.seq_at_simd+2]
            seq: seq::decode_seq(&raw[cigar_end..seq_end], seq_len_usize).into_boxed_slice(),
            qual: raw[seq_end..qual_end].into(),
            aux: raw[qual_end..].into(),
        })
    }

    // --- accessors ---

    pub fn bam_flags(&self) -> BamFlags {
        BamFlags::new(self.flags)
    }

    // r[impl bam.record.flag_reverse]
    pub fn is_reverse(&self) -> bool {
        self.flags & FLAG_REVERSE != 0
    }

    // r[impl bam.record.flag_first]
    pub fn is_first_in_template(&self) -> bool {
        self.flags & FLAG_FIRST_IN_TEMPLATE != 0
    }

    // r[impl bam.record.flag_second]
    pub fn is_second_in_template(&self) -> bool {
        self.flags & FLAG_SECOND_IN_TEMPLATE != 0
    }

    // r[impl bam.record.flag_unmapped]
    pub fn is_unmapped(&self) -> bool {
        self.flags & FLAG_UNMAPPED != 0
    }

    // r[impl bam.record.seq_at]
    pub fn seq_at(&self, pos: usize) -> u8 {
        self.seq.get(pos).copied().unwrap_or(b'N')
    }

    // r[impl bam.record.aux_parse]
    pub fn aux(&self, tag: &[u8; 2]) -> Option<AuxValue<'_>> {
        find_aux_tag(&self.aux, tag)
    }
}

/// Compute end_pos from raw BAM record bytes (before full decode).
pub fn compute_end_pos_from_raw(raw: &[u8]) -> Option<i64> {
    if raw.len() < 32 {
        return None;
    }
    debug_assert!(raw.len() >= 32, "raw record too short for fixed fields: {}", raw.len());
    let pos = i64::from(i32::from_le_bytes(read4(raw, 4)));
    // raw.len() >= 32, so raw[8] is in bounds
    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
    let name_len = raw[8] as usize;
    let n_cigar_ops = u16::from_le_bytes(read2(raw, 12)) as usize;
    let cigar_start = 32 + name_len;
    let cigar_end = cigar_start + n_cigar_ops * 4;
    if raw.len() < cigar_end {
        return None;
    }
    // cigar_end ≤ raw.len() checked above
    debug_assert!(cigar_end <= raw.len(), "cigar overrun: {cigar_end} > {}", raw.len());
    #[allow(clippy::indexing_slicing, reason = "cigar_end ≤ raw.len() checked above")]
    Some(compute_end_pos(pos, &raw[cigar_start..cigar_end]))
}

// r[impl bam.record.end_pos]
// r[impl bam.record.zero_refspan]
pub(crate) fn compute_end_pos(pos: i64, cigar_bytes: &[u8]) -> i64 {
    use super::cigar::{CIGAR_D, CIGAR_EQ, CIGAR_M, CIGAR_N, CIGAR_X};

    let mut ref_len: i64 = 0;
    let n_ops = cigar_bytes.len() / 4;
    for i in 0..n_ops {
        let op = u32::from_le_bytes(read4(cigar_bytes, i * 4));
        let op_len = i64::from(op >> 4);
        let op_type = (op & 0xF) as u8;
        match op_type {
            CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X => ref_len += op_len,
            _ => {}
        }
    }
    if ref_len == 0 { pos } else { pos + ref_len - 1 }
}

// Callers only invoke read2/read4 after validating raw.len() >= 32 or equivalent.
#[allow(
    clippy::indexing_slicing,
    reason = "offset + 1 < raw.len() ensured by caller's length check"
)]
fn read2(buf: &[u8], offset: usize) -> [u8; 2] {
    debug_assert!(
        offset + 1 < buf.len(),
        "read2 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [buf[offset], buf[offset + 1]]
}

#[allow(
    clippy::indexing_slicing,
    reason = "offset + 3 < raw.len() ensured by caller's length check"
)]
fn read4(buf: &[u8], offset: usize) -> [u8; 4] {
    debug_assert!(
        offset + 3 < buf.len(),
        "read4 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [buf[offset], buf[offset + 1], buf[offset + 2], buf[offset + 3]]
}

#[derive(Debug, thiserror::Error)]
pub enum DecodeError {
    #[error("BAM record too short: {len} bytes")]
    TooShort { len: usize },
}

// r[impl bam.record.raw_aux]
/// Parsed auxiliary tag value.
#[derive(Debug, Clone, PartialEq)]
pub enum AuxValue<'a> {
    Char(u8),
    I8(i8),
    U8(u8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    Float(f32),
    Double(f64),
    String(&'a [u8]),
    Hex(&'a [u8]),
    ArrayI8(&'a [u8]),
    ArrayU8(&'a [u8]),
    ArrayI16(&'a [u8]),
    ArrayU16(&'a [u8]),
    ArrayI32(&'a [u8]),
    ArrayU32(&'a [u8]),
    ArrayFloat(&'a [u8]),
}

// aux[0..2] valid since while loop checks aux.len() >= 3; aux[3..] valid for same reason.
#[allow(clippy::indexing_slicing, reason = "aux.len() >= 3 enforced by while condition")]
pub fn find_aux_tag<'a>(mut aux: &'a [u8], tag: &[u8; 2]) -> Option<AuxValue<'a>> {
    while aux.len() >= 3 {
        debug_assert!(aux.len() >= 3, "aux too short for tag+type: {}", aux.len());
        let t = [aux[0], aux[1]];
        let type_byte = aux[2];
        aux = &aux[3..];
        let (value, consumed) = parse_aux_value(type_byte, aux)?;
        if t == *tag {
            return Some(value);
        }
        aux = aux.get(consumed..)?;
    }
    None
}

// Indexing is guarded by explicit length checks within each match arm.
#[allow(clippy::indexing_slicing, reason = "each arm validates data.len() before indexing")]
fn parse_aux_value(type_byte: u8, data: &[u8]) -> Option<(AuxValue<'_>, usize)> {
    match type_byte {
        b'A' => Some((AuxValue::Char(*data.first()?), 1)),
        b'c' => Some((AuxValue::I8(*data.first()? as i8), 1)),
        b'C' => Some((AuxValue::U8(*data.first()?), 1)),
        b's' => {
            let v = i16::from_le_bytes([*data.first()?, *data.get(1)?]);
            Some((AuxValue::I16(v), 2))
        }
        b'S' => {
            let v = u16::from_le_bytes([*data.first()?, *data.get(1)?]);
            Some((AuxValue::U16(v), 2))
        }
        b'i' => {
            let v = i32::from_le_bytes(read4_checked(data)?);
            Some((AuxValue::I32(v), 4))
        }
        b'I' => {
            let v = u32::from_le_bytes(read4_checked(data)?);
            Some((AuxValue::U32(v), 4))
        }
        b'f' => {
            let v = f32::from_le_bytes(read4_checked(data)?);
            Some((AuxValue::Float(v), 4))
        }
        b'd' => {
            if data.len() < 8 {
                return None;
            }
            debug_assert!(data.len() >= 8, "double aux data too short: {}", data.len());
            let mut bytes = [0u8; 8];
            bytes.copy_from_slice(&data[..8]);
            Some((AuxValue::Double(f64::from_le_bytes(bytes)), 8))
        }
        b'Z' => {
            let end = data.iter().position(|&b| b == 0)?;
            Some((AuxValue::String(&data[..end]), end + 1))
        }
        b'H' => {
            let end = data.iter().position(|&b| b == 0)?;
            Some((AuxValue::Hex(&data[..end]), end + 1))
        }
        b'B' => {
            if data.len() < 5 {
                return None;
            }
            debug_assert!(data.len() >= 5, "array aux data too short: {}", data.len());
            let sub_type = data[0];
            let count = u32::from_le_bytes(read4_checked(&data[1..])?) as usize;
            let elem_size = match sub_type {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            let total = 5 + count * elem_size;
            if data.len() < total {
                return None;
            }
            debug_assert!(total <= data.len(), "array data overrun: {total} > {}", data.len());
            let array_data = &data[5..total];
            let value = match sub_type {
                b'c' => AuxValue::ArrayI8(array_data),
                b'C' => AuxValue::ArrayU8(array_data),
                b's' => AuxValue::ArrayI16(array_data),
                b'S' => AuxValue::ArrayU16(array_data),
                b'i' => AuxValue::ArrayI32(array_data),
                b'I' => AuxValue::ArrayU32(array_data),
                b'f' => AuxValue::ArrayFloat(array_data),
                _ => return None,
            };
            Some((value, total))
        }
        _ => None,
    }
}

fn read4_checked(data: &[u8]) -> Option<[u8; 4]> {
    Some([*data.first()?, *data.get(1)?, *data.get(2)?, *data.get(3)?])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_end_pos() {
        let op = 50u32 << 4;
        assert_eq!(compute_end_pos(100, &op.to_le_bytes()), 149);
    }

    #[test]
    fn test_decode_record() {
        let mut raw = [0u8; 64];
        raw[0..4].copy_from_slice(&0i32.to_le_bytes());
        raw[4..8].copy_from_slice(&100i32.to_le_bytes());
        raw[8] = 5;
        raw[9] = 60;
        raw[10..12].copy_from_slice(&0u16.to_le_bytes());
        raw[12..14].copy_from_slice(&1u16.to_le_bytes());
        raw[14..16].copy_from_slice(&99u16.to_le_bytes());
        raw[16..20].copy_from_slice(&4u32.to_le_bytes());
        raw[20..32].fill(0);
        raw[32..37].copy_from_slice(b"read\0");
        raw[37..41].copy_from_slice(&(4u32 << 4).to_le_bytes());
        raw[41] = 0x12;
        raw[42] = 0x48;
        raw[43..47].copy_from_slice(&[30, 30, 30, 30]);

        let rec = BamRecord::decode(&raw[..47]).unwrap();
        assert_eq!(rec.pos, 100);
        assert_eq!(rec.end_pos, 103);
        assert_eq!(rec.flags, 99);
        assert_eq!(rec.mapq, 60);
        assert_eq!(&*rec.qname, b"read");
        assert_eq!(rec.seq_at(0), b'A');
        assert_eq!(rec.seq_at(3), b'T');
        assert_eq!(rec.matching_bases, 4);
    }
}
