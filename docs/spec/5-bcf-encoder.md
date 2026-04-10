# BCF Encoding Primitives

> **Sources:** [BCF2] — typed value encoding, record layout, field-major FORMAT encoding, GT encoding. See [BCF Writer](./5-bcf-writer.md) for the underlying BCF format rules. These primitives are used internally by the unified [`RecordEncoder`](./5-record-encoder.md) for BCF output. Also see [References](./99-references.md).

## BcfValue trait

r[bcf_encoder.bcf_value]
The `BcfValue` trait MUST define: `bcf_type_code() -> u8` (BCF type code for this Rust type), `encode_bcf(&self, buf: &mut Vec<u8>)` (write value bytes), `encode_missing(buf: &mut Vec<u8>)` (write missing sentinel), `encode_end_of_vector(buf: &mut Vec<u8>)` (write EOV sentinel).

r[bcf_encoder.bcf_value_int]
Scalar `i32` values MUST select the smallest BCF integer type that fits, matching `r[bcf_writer.smallest_int_type]`. For arrays, the type is determined by scanning all concrete (non-missing) values first. Missing values within integer arrays MUST use the per-type sentinel (int8=0x80, int16=0x8000, int32=0x80000000) matching the selected type, not a fixed i32::MIN.

r[bcf_encoder.bcf_value_float]
`f32` values MUST be encoded as IEEE 754 single-precision LE bytes. Missing sentinel `0x7F800001` and EOV sentinel `0x7F800002` MUST be written as raw bytes (never through float arithmetic) per `r[bcf_writer.missing_sentinels]`.

## BCF Record Encoding

r[bcf_encoder.begin_record]
Beginning a BCF record MUST write the 24-byte fixed header (with placeholder n_info/n_fmt), ID as `.`, and REF/ALT allele strings using zero-alloc `write_ref_into`/`write_alts_into`. It MUST set `n_allele` and `n_alt` for downstream validation.

r[bcf_encoder.checked_casts]
All conversions from `usize`/`u32` to narrower integer types (`i32`, `u16`, `u8`) for user-controlled values (contig tid, position, allele count, sample count) MUST use checked conversion (`try_from`) and return a typed error on overflow rather than silently truncating via `as`.

r[bcf_encoder.emit]
BCF `emit()` MUST patch the n_info, n_fmt, n_sample fields in the 24-byte fixed header, then flush the record to BGZF (with `flush_if_needed`), then push to the index builder if present.

r[bcf_encoder.info_counting]
Each INFO field encoded MUST increment the encoder's `n_info` counter. Each FORMAT field encoded MUST increment `n_fmt`.

r[bcf_encoder.format_field_major]
FORMAT fields MUST be encoded in field-major order per `r[bcf_writer.indiv_field_major]`. For single-sample records (the common case), each format encode call writes the key + type descriptor + 1 value.
