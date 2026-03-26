//! Shared test helpers for building synthetic BAM records.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
#![allow(dead_code)]

pub fn cigar_op(len: u32, op: u8) -> u32 {
    (len << 4) | u32::from(op)
}

pub fn cigar_bytes(ops: &[u32]) -> Vec<u8> {
    ops.iter().flat_map(|op| op.to_le_bytes()).collect()
}

/// Build a synthetic BAM record with a simple N×M CIGAR.
pub fn make_record(tid: i32, pos: i32, flags: u16, mapq: u8, seq_len: u32) -> Vec<u8> {
    make_record_with_cigar(tid, pos, flags, mapq, &[cigar_op(seq_len, 0)], seq_len)
}

/// Build a synthetic BAM record with a custom CIGAR.
pub fn make_record_with_cigar(
    tid: i32,
    pos: i32,
    flags: u16,
    mapq: u8,
    cigar_ops: &[u32],
    seq_len: u32,
) -> Vec<u8> {
    let name = b"read\0";
    let name_len = name.len();
    let n_cigar_ops = u16::try_from(cigar_ops.len()).unwrap();
    let cigar_bytes_len = cigar_ops.len() * 4;
    let seq_bytes = (seq_len as usize).div_ceil(2);

    let total = 32 + name_len + cigar_bytes_len + seq_bytes + seq_len as usize;
    let mut raw = vec![0u8; total];

    raw[0..4].copy_from_slice(&tid.to_le_bytes());
    raw[4..8].copy_from_slice(&pos.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = mapq;
    raw[10..12].copy_from_slice(&0u16.to_le_bytes());
    raw[12..14].copy_from_slice(&n_cigar_ops.to_le_bytes());
    raw[14..16].copy_from_slice(&flags.to_le_bytes());
    raw[16..20].copy_from_slice(&seq_len.to_le_bytes());
    raw[20..32].fill(0);
    raw[32..32 + name_len].copy_from_slice(name);

    let cigar_start = 32 + name_len;
    for (i, op) in cigar_ops.iter().enumerate() {
        raw[cigar_start + i * 4..cigar_start + (i + 1) * 4].copy_from_slice(&op.to_le_bytes());
    }

    let seq_start = cigar_start + cigar_bytes_len;
    for i in 0..seq_bytes {
        raw[seq_start + i] = 0x11; // A,A
    }

    let qual_start = seq_start + seq_bytes;
    for i in 0..seq_len as usize {
        raw[qual_start + i] = 30;
    }

    raw
}

/// Build a BAM record with a custom qname and packed sequence bytes.
/// `seq_packed` is the 4-bit packed sequence. `seq_len` is the number of bases.
pub fn make_named_record(
    qname: &[u8],
    tid: i32,
    pos: i32,
    flags: u16,
    mapq: u8,
    seq_len: u32,
    seq_packed: &[u8],
) -> Vec<u8> {
    let mut name_with_nul = qname.to_vec();
    name_with_nul.push(0);
    // Pad to 4-byte alignment as BAM spec requires
    while !name_with_nul.len().is_multiple_of(4) {
        name_with_nul.push(0);
    }
    let name_len = name_with_nul.len();
    let n_cigar_ops = 1u16;
    let cigar = cigar_op(seq_len, 0); // simple NM
    let seq_bytes = (seq_len as usize).div_ceil(2);

    let total = 32 + name_len + 4 + seq_bytes + seq_len as usize;
    let mut raw = vec![0u8; total];

    raw[0..4].copy_from_slice(&tid.to_le_bytes());
    raw[4..8].copy_from_slice(&pos.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = mapq;
    raw[10..12].copy_from_slice(&0u16.to_le_bytes());
    raw[12..14].copy_from_slice(&n_cigar_ops.to_le_bytes());
    raw[14..16].copy_from_slice(&flags.to_le_bytes());
    raw[16..20].copy_from_slice(&seq_len.to_le_bytes());
    raw[20..32].fill(0);
    raw[32..32 + name_len].copy_from_slice(&name_with_nul);

    let cigar_start = 32 + name_len;
    raw[cigar_start..cigar_start + 4].copy_from_slice(&cigar.to_le_bytes());

    let seq_start = cigar_start + 4;
    let copy_len = seq_packed.len().min(seq_bytes);
    raw[seq_start..seq_start + copy_len].copy_from_slice(&seq_packed[..copy_len]);

    let qual_start = seq_start + seq_bytes;
    for i in 0..seq_len as usize {
        raw[qual_start + i] = 30;
    }

    raw
}

/// BAM 4-bit encoding for common bases
pub const BASE_A: u8 = 1;
pub const BASE_C: u8 = 2;
pub const BASE_G: u8 = 4;
pub const BASE_T: u8 = 8;

/// Pack two 4-bit bases into one byte
pub fn pack_bases(hi: u8, lo: u8) -> u8 {
    (hi << 4) | lo
}
