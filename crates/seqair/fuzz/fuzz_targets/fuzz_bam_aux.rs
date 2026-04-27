#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::aux::{self, AuxValue};
use seqair::bam::aux_data::AuxData;

/// Round-trip an `AuxValue` through `AuxData` and verify it comes back intact.
///
/// Integer types may select a wider BAM encoding (e.g., i8 42 → U8 42).
/// We accept `as_i64()` equivalence for integers. Floats compare bit patterns
/// since NaN ≠ NaN under `PartialEq`.
fn assert_value_roundtrip(original: &AuxValue<'_>, roundtripped: Option<AuxValue<'_>>) {
    let rt = roundtripped.expect("tag should survive round-trip");
    match (original, &rt) {
        (AuxValue::Float(a), AuxValue::Float(b)) => {
            assert_eq!(a.to_bits(), b.to_bits(), "float bit-pattern mismatch");
        }
        (AuxValue::Double(a), AuxValue::Double(b)) => {
            assert_eq!(a.to_bits(), b.to_bits(), "double bit-pattern mismatch");
        }
        (AuxValue::I8(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::U8(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::I16(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::U16(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::I32(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::U32(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        _ => assert_eq!(original, &rt, "non-integer round-trip mismatch"),
    }
}

/// Set a single tag value into `AuxData`.
fn set_one(aux: &mut AuxData, tag: [u8; 2], value: &AuxValue<'_>) {
    match value {
        AuxValue::Char(v) => aux.set_char(tag, *v),
        AuxValue::String(s) => aux.set_string(tag, s),
        AuxValue::Hex(h) => aux.set_string(tag, h),
        AuxValue::I8(v) => drop(aux.set_int(tag, i64::from(*v))),
        AuxValue::U8(v) => drop(aux.set_int(tag, i64::from(*v))),
        AuxValue::I16(v) => drop(aux.set_int(tag, i64::from(*v))),
        AuxValue::U16(v) => drop(aux.set_int(tag, i64::from(*v))),
        AuxValue::I32(v) => drop(aux.set_int(tag, i64::from(*v))),
        AuxValue::U32(v) => drop(aux.set_int(tag, i64::from(*v))),
        AuxValue::Float(v) => aux.set_float(tag, *v),
        AuxValue::Double(v) => aux.set_double(tag, *v),
        AuxValue::ArrayI8(a) => drop(aux.set_array_i16(tag, a)),
        AuxValue::ArrayU8(a) => drop(aux.set_array_u8(tag, a)),
        AuxValue::ArrayI16(a) => drop(aux.set_array_i16(tag, a)),
        AuxValue::ArrayU16(a) => drop(aux.set_array_u16(tag, a)),
        AuxValue::ArrayI32(a) => drop(aux.set_array_i32(tag, a)),
        AuxValue::ArrayU32(a) => drop(aux.set_array_u32(tag, a)),
        AuxValue::ArrayFloat(a) => drop(aux.set_array_f32(tag, a)),
    }
}

fuzz_target!(|data: &[u8]| {
    // ── 1. Parse raw bytes (no panics) ──
    let tags: Vec<_> = aux::iter_tags(data).collect();

    // Also exercise find_tag for a few fixed tags
    let _ = aux::find_tag(data, *b"RG");
    let _ = aux::find_tag(data, *b"NM");
    let _ = aux::find_tag(data, *b"ZZ");

    // ── 2. Deduplicate: AuxData overwrites on set, so only the LAST
    //    occurrence of each tag name survives the round-trip.
    //    Skip Hex (no dedicated setter) and ArrayI8 (alignment mismatch
    //    with set_array_i16).
    let mut last_value: rustc_hash::FxHashMap<[u8; 2], &AuxValue<'_>> =
        rustc_hash::FxHashMap::default();
    for (tag, value) in &tags {
        if matches!(value, AuxValue::Hex(_) | AuxValue::ArrayI8(_)) {
            continue;
        }
        last_value.insert(*tag, value);
    }

    // ── 3. Build AuxData from deduplicated tags (last-wins) ──
    let mut aux = AuxData::new();
    let mut stored: Vec<[u8; 2]> = Vec::new();
    for (tag, value) in &last_value {
        set_one(&mut aux, *tag, value);
        stored.push(*tag);
    }

    // ── 4. Verify round-trip for every stored tag ──
    for tag in &stored {
        let expected = last_value.get(tag).unwrap();
        assert_value_roundtrip(expected, aux.get(*tag));
    }

    // ── 5. Sequential removal → empty ──
    for tag in &stored {
        aux.remove(*tag);
    }
    assert!(aux.is_empty(), "all tags removed but aux not empty");
});
