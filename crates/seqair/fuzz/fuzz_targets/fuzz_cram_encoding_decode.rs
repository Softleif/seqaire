#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::cram::{
    bitstream::BitReader,
    encoding::{ExternalCursor, HuffmanTable},
};

#[derive(Arbitrary, Debug)]
struct EncodingInput {
    alphabet: Vec<i32>,
    bit_lengths: Vec<u32>,
    core_data: Vec<u8>,
    external_data: Vec<u8>,
    external_ops: Vec<ExternalOp>,
}

#[derive(Arbitrary, Debug)]
enum ExternalOp {
    ReadByte,
    ReadItf8,
    ReadBytesUntil(u8),
    ReadBytes(u8),
    Remaining,
}

fuzz_target!(|input: EncodingInput| {
    // Fuzz HuffmanTable construction and decoding
    if let Ok(table) = HuffmanTable::new(&input.alphabet, &input.bit_lengths) {
        let mut reader = BitReader::new(&input.core_data);
        let _ = table.decode(&mut reader);
        // Try decoding multiple symbols
        let mut reader2 = BitReader::new(&input.core_data);
        for _ in 0..8 {
            if table.decode(&mut reader2).is_none() {
                break;
            }
        }
    }

    // Fuzz ExternalCursor operations
    let mut cursor = ExternalCursor::new(input.external_data);
    let mut buf = Vec::new();
    for op in &input.external_ops {
        match op {
            ExternalOp::ReadByte => {
                let _ = cursor.read_byte();
            }
            ExternalOp::ReadItf8 => {
                let _ = cursor.read_itf8();
            }
            ExternalOp::ReadBytesUntil(stop) => {
                buf.clear();
                let _ = cursor.read_bytes_until_into(*stop, &mut buf);
            }
            ExternalOp::ReadBytes(n) => {
                // Clamp to avoid huge allocations
                let n = (*n as usize) % 64;
                buf.clear();
                let _ = cursor.read_bytes_into(n, &mut buf);
            }
            ExternalOp::Remaining => {
                let _ = cursor.remaining();
            }
        }
    }
});
