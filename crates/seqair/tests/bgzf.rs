//! Tests for BGZF virtual offset representation and ordering.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use proptest::prelude::*;
use seqair::bam::bgzf::VirtualOffset;

// r[verify bgzf.virtual_offset]
proptest! {
    #[test]
    fn virtual_offset_roundtrips(
        block_offset in 0u64..=(1u64 << 48) - 1,
        within_block in 0u16..=u16::MAX,
    ) {
        let vo = VirtualOffset::new(block_offset, within_block);
        prop_assert_eq!(vo.block_offset(), block_offset);
        prop_assert_eq!(vo.within_block(), within_block);
    }
}

// r[verify bgzf.virtual_offset]
proptest! {
    #[test]
    fn ordering_matches_file_order(b1 in 0u64..1_000_000, b2 in 0u64..1_000_000, w1 in 0u16..=1000, w2 in 0u16..=1000) {
        let vo1 = VirtualOffset::new(b1, w1);
        let vo2 = VirtualOffset::new(b2, w2);
        if b1 != b2 {
            prop_assert_eq!(vo1.cmp(&vo2), b1.cmp(&b2));
        } else {
            prop_assert_eq!(vo1.cmp(&vo2), w1.cmp(&w2));
        }
    }
}

// r[verify io.errors]
#[test]
fn error_types_implement_std_error_and_capture_context() {
    use std::error::Error;

    let err = seqair::BamHeaderError::InvalidMagic;
    let _: &dyn Error = &err;
    assert!(format!("{err}").contains("BAM"));

    let err = seqair::BamError::ContigNotFound { name: seqair_types::SmolStr::new("chr99") };
    assert!(format!("{err}").contains("chr99"));

    let err = seqair::BamError::RegionOutOfBounds {
        contig: seqair_types::SmolStr::new("chr1"),
        start: 100,
        end: 200,
        contig_len: 50,
    };
    let msg = format!("{err}");
    assert!(msg.contains("chr1") && msg.contains("100") && msg.contains("50"));
}

// r[verify bgzf.crc32]
// r[verify bgzf.fast_header]
// r[verify bgzf.resize_uninit]
// r[verify bgzf.block_offset_tracking]
#[test]
fn read_bam_verifies_crc32_and_tracks_offsets() {
    use seqair::bam::reader::IndexedBamReader;
    use seqair::bam::record_store::RecordStore;
    use std::path::Path;

    let bam_path = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
    let mut reader = IndexedBamReader::open(bam_path).expect("open BAM");
    let tid = reader.header().tid("chr19").expect("tid");

    let mut store = RecordStore::new();
    let count = reader
        .fetch_into(tid, 6_105_700, 6_105_800, &mut store)
        .expect("CRC32 verification should pass on valid BAM data");

    // If CRC32 verification failed, fetch_into would have returned an error.
    // If block_offset tracking broke, seeking via the index would fail.
    // If resize_uninit corrupted data, record decoding would fail.
    // If fast_header parsing failed, the block wouldn't decompress.
    assert!(count > 0, "should have read some records");
}

// r[verify bgzf.crc32]
#[test]
fn checksum_mismatch_error_formats_correctly() {
    let err = seqair::BgzfError::ChecksumMismatch { expected: 0xDEADBEEF, found: 0xCAFEBABE };
    let msg = format!("{err}");
    assert!(msg.contains("0xdeadbeef"), "should contain expected CRC: {msg}");
    assert!(msg.contains("0xcafebabe"), "should contain found CRC: {msg}");
}

// r[verify io.named_constants]
#[test]
fn bam_flag_constants_match_sam_spec() {
    use seqair::bam::flags::*;
    assert_eq!(FLAG_PAIRED, 0x1);
    assert_eq!(FLAG_PROPER_PAIR, 0x2);
    assert_eq!(FLAG_UNMAPPED, 0x4);
    assert_eq!(FLAG_MATE_UNMAPPED, 0x8);
    assert_eq!(FLAG_REVERSE, 0x10);
    assert_eq!(FLAG_MATE_REVERSE, 0x20);
    assert_eq!(FLAG_FIRST_IN_TEMPLATE, 0x40);
    assert_eq!(FLAG_SECOND_IN_TEMPLATE, 0x80);
    assert_eq!(FLAG_SECONDARY, 0x100);
    assert_eq!(FLAG_FAILED_QC, 0x200);
    assert_eq!(FLAG_DUPLICATE, 0x400);
    assert_eq!(FLAG_SUPPLEMENTARY, 0x800);
}

// r[verify io.typed_flags]
#[test]
fn bam_flags_newtype_predicate_methods() {
    use seqair::bam::flags::BamFlags;

    // 0x53 = paired(0x1) + proper_pair(0x2) + reverse(0x10) + first_in_template(0x40)
    let flags = BamFlags::new(0x53);
    assert!(flags.is_paired());
    assert!(flags.is_proper_pair());
    assert!(flags.is_reverse());
    assert!(flags.is_first_in_template());
    assert!(!flags.is_unmapped());
    assert!(!flags.is_secondary());
    assert!(!flags.is_supplementary());

    let unmapped = BamFlags::new(0x4);
    assert!(unmapped.is_unmapped());
    assert!(!unmapped.is_reverse());

    assert_eq!(flags.raw(), 0x53);
}
