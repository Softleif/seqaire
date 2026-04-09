# Record Encoder

> **Sources:** The `RecordEncoder` trait unifies the [BCF Direct Encoder](./5-bcf-encoder.md) and the [VCF Text Writer](./5-vcf-writer.md) behind a single format-agnostic API. Field definitions provide a single source of truth for header construction, key resolution, and documentation generation. See also [VCF Header](./5-vcf-header.md) for string dictionary assignment.

## Motivation

The `VcfRecord` intermediate allocates per-record (string keys, SmallVec values, enum boxing). For high-throughput callers (millions of records), a direct encoder avoids this entirely. The existing `BcfRecordEncoder` does this for BCF but not for VCF text. The `RecordEncoder` trait provides format-agnostic encoding so the same calling code works for VCF, VCF.gz, and BCF output without allocating per-record.

## Field Definitions

r[record_encoder.field_def]
Field definitions MUST carry all metadata needed for header construction, key resolution, and documentation: name (`&'static str`), Number, ValueType, and description (`&'static str`). Field definitions MUST be parameterized by a value type marker to enable type-safe key resolution.

> r[record_encoder.field_def_types]
> Three definition types MUST be provided:
>
> - `InfoFieldDef<V>` — INFO field definition, produces `InfoKey<V>` on registration
> - `FormatFieldDef<V>` — FORMAT field definition, produces `FormatKey<V>` on registration
> - `FilterFieldDef` — FILTER definition, produces `FilterId` on registration

r[record_encoder.field_def_const]
Field definitions MUST be constructible as `const` values so that downstream crates can define them at compile time.

r[record_encoder.field_description]
A `FieldDescription` trait MUST provide read access to name, number, value type, and description. Both `InfoFieldDef<V>` and `FormatFieldDef<V>` MUST implement this trait (via a type-erased form) to enable documentation generation from heterogeneous collections of field definitions.

## Field Registration

r[record_encoder.register]
`VcfHeaderBuilder` MUST provide `register_info`, `register_format`, and `register_filter` methods that accept a field definition, add the corresponding header entry, and return a resolved typed key. This combines header construction and key resolution in a single step, eliminating the possibility of declaring a header field without resolving its key or vice versa.

r[record_encoder.register_contig]
`VcfHeaderBuilder` MUST provide a `register_contig` method that adds a contig to the header and returns a `ContigId` carrying both the integer tid and the contig name.

## Field Identifiers and Keys

r[record_encoder.field_id]
`FieldId` MUST carry both a BCF dictionary index (`u32`) and a string name (`SmolStr`). The BCF encoder uses the dictionary index; the VCF text encoder uses the name. Both are resolved at registration time, not per-record.

r[record_encoder.contig_id]
`ContigId` MUST carry both the integer tid (`u32`) and the contig name (`SmolStr`).

r[record_encoder.filter_id]
`FilterId` MUST carry both the BCF dictionary index and the filter name. A `FilterId::PASS` constant MUST be provided with dictionary index 0 and name `"PASS"`.

### Typed Keys

r[record_encoder.typed_keys]
Typed keys MUST encode the VCF value type at the Rust type level via uninhabited marker enums. This prevents passing a float value where an integer key is expected, or using an INFO key in a FORMAT context.

> r[record_encoder.key_types]
> The following key types MUST be provided:
>
> - `InfoKey<Scalar<i32>>` — single integer INFO value
> - `InfoKey<Scalar<f32>>` — single float INFO value
> - `InfoKey<Arr<i32>>` — integer array INFO value
> - `InfoKey<Arr<f32>>` — float array INFO value
> - `InfoKey<Flag>` — flag INFO value (no data)
> - `InfoKey<Str>` — string INFO value
> - `InfoKey<OptArr<i32>>` — integer array with optional (missing) elements
> - `FormatKey<Gt>` — genotype FORMAT value
> - `FormatKey<Scalar<i32>>` — single integer FORMAT value
> - `FormatKey<Scalar<f32>>` — single float FORMAT value

r[record_encoder.key_encode]
Each typed key MUST provide an `encode` method that accepts `&mut impl RecordEncoder` and the appropriate value type, delegating to the corresponding `RecordEncoder` trait method. The method signature MUST vary by key kind:

- `InfoKey<Scalar<i32>>::encode(&self, enc, value: i32)`
- `InfoKey<Scalar<f32>>::encode(&self, enc, value: f32)`
- `InfoKey<Arr<T>>::encode(&self, enc, values: &[T])`
- `InfoKey<Flag>::encode(&self, enc)` — no value argument
- `InfoKey<Str>::encode(&self, enc, value: &str)`
- `InfoKey<OptArr<i32>>::encode(&self, enc, values: &[Option<i32>])`
- `FormatKey<Gt>::encode(&self, enc, gt: &Genotype)`
- `FormatKey<Scalar<T>>::encode(&self, enc, value: T)`

## RecordEncoder Trait

r[record_encoder.trait]
The `RecordEncoder` trait MUST be object-safe (all methods use `&mut self` and concrete parameter types, no generics) to enable both monomorphized and dynamic dispatch.

r[record_encoder.begin]
`begin()` MUST accept a `ContigId`, 1-based position, `Alleles`, and optional quality. It MUST clear all per-record state and write/buffer the fixed fields (CHROM, POS, ID, REF, ALT, QUAL). Overflow checks on position and allele counts MUST use checked conversions and return typed errors.

r[record_encoder.filters]
Exactly one filter method MUST be called per record: `filter_pass()` for PASS, or `filter_fail(&[&FilterId])` for one or more failed filters.

r[record_encoder.info_methods]
INFO methods (`info_int`, `info_float`, `info_ints`, `info_floats`, `info_flag`, `info_string`, `info_int_opts`) MUST accept a `&FieldId` and the appropriate value. These methods MUST be infallible — they write to in-memory buffers which cannot fail.

r[record_encoder.format_methods]
`begin_samples(n)` MUST be called before any FORMAT method. FORMAT methods (`format_gt`, `format_int`, `format_float`) MUST accept a `&FieldId` and the appropriate value. These methods MUST be infallible.

r[record_encoder.emit]
`emit()` MUST finalize the record and write it to the output. For BCF, this patches the fixed header and flushes to BGZF. For VCF text, this writes accumulated FORMAT fields and flushes the line buffer. `emit()` is the only encoding method that performs I/O and MAY return an error.

r[record_encoder.state_queries]
`n_allele()` and `n_alt()` MUST return the number of alleles and alternate alleles for the current record, set during `begin()`.

## Custom Type Encoding

r[record_encoder.encode_info_trait]
An `EncodeInfo` trait MUST be provided with an associated `Key` type and an `encode_info(&self, enc: &mut dyn RecordEncoder, key: &Self::Key)` method. This allows domain types to encapsulate their VCF encoding logic, including multi-field expansion (e.g., a strand-specific type producing separate OT and OB fields).

r[record_encoder.encode_format_trait]
An `EncodeFormat` trait MUST be provided with the same pattern as `EncodeInfo`. Implementations that produce no output for certain values (e.g., unknown methylation status) simply do not call any encoder methods, and the FORMAT key does not appear.

r[record_encoder.encode_dyn]
`EncodeInfo` and `EncodeFormat` MUST use `&mut dyn RecordEncoder` (not `&mut impl RecordEncoder`) so that implementations do not need to be generic over the encoder type.

## VCF Text Record Encoder

r[record_encoder.vcf_text_encoder]
`VcfRecordEncoder` MUST implement `RecordEncoder` for VCF text output. It MUST borrow reusable buffers from the `VcfWriter` to avoid per-record allocation.

r[record_encoder.vcf_text_begin]
`begin()` MUST clear the line buffer and FORMAT accumulators, then write the fixed columns (CHROM through QUAL) as tab-separated text per `r[vcf_writer.tab_delimited]`.

r[record_encoder.vcf_text_info]
INFO methods MUST append to the line buffer using semicolon separators per `r[vcf_writer.info_serialization]`. The first INFO field has no leading semicolon. If no INFO fields are written, `emit()` MUST write `.`.

r[record_encoder.vcf_text_format_accumulation]
FORMAT methods MUST accumulate key names and formatted values into separate reusable buffers. Keys accumulate in order of first `format_*` call. Values are formatted into a per-sample buffer with `:` separators.

r[record_encoder.vcf_text_emit]
`emit()` MUST write the FILTER and INFO columns already in the line buffer, then the FORMAT key column and sample value column(s) from the accumulators, then a newline. The complete line MUST be flushed to the output (plain write or BGZF with index co-production).

r[record_encoder.vcf_text_output]
A `VcfOutput` trait MUST abstract over plain and BGZF output modes, providing `write_line()` and `push_index()` methods. This mirrors the `BgzfWrite` trait used by the BCF encoder.

r[record_encoder.vcf_text_buffer_reuse]
All buffers (line buffer, FORMAT key accumulator, FORMAT value buffer) MUST be reused across records per `r[vcf_writer.buffer_reuse]`. After warmup, encoding MUST NOT allocate.

## BCF Record Encoder

r[record_encoder.bcf_impl]
`BcfRecordEncoder` MUST implement `RecordEncoder` by delegating to its existing buffer operations. INFO methods write to `shared_buf`, FORMAT methods write to `indiv_buf`, following all rules from the [BCF Direct Encoder](./5-bcf-encoder.md) spec.

r[record_encoder.bcf_backwards_compat]
The existing handle-based API (`ScalarInfoHandle<T>::encode`, etc.) MUST remain available alongside the `RecordEncoder` trait implementation for backward compatibility.

## Equivalence

r[record_encoder.vcf_bcf_equivalence]
For the same logical record, encoding via `BcfRecordEncoder` and `VcfRecordEncoder` MUST produce output that, when read back, yields identical field values (within floating-point formatting precision).

r[record_encoder.record_path_equivalence]
Encoding via `RecordEncoder` trait methods MUST produce output identical to encoding via the `VcfRecord` path (`VcfRecordBuilder` → `VcfWriter::write_record` / `BcfWriter::write_vcf_record`) for the same logical record.
