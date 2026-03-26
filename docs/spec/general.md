# Seqair — General Design

These are cross-cutting design rules that apply to all modules in the Seqair crate. They establish conventions for error handling, type safety, and platform portability.

> **Sources:** This file contains seqair-specific design rules with no direct upstream spec counterpart. Named constants (FLAG bits, CIGAR op codes) are derived from [SAM1] §1.4 and §4.2. See [references.md](references.md).

## Error handling

Seqair aims to provide helpful and complete error messages at all layers of the system.

r[io.errors]
Error types MUST implement `std::error::Error` using `thiserror` and capture necessary context and provide causes.

r[io.errors.typed_variants]
Every failure mode MUST be a distinct enum variant with typed fields — never a catch-all `String` or `reason: String`. Callers MUST be able to pattern-match on failure modes without inspecting display text. When a single enum variant is constructed at multiple call sites with different format strings, split it into separate variants. Paths, byte arrays, and counts MUST be stored in typed fields (e.g. `PathBuf`, `[u8; 4]`, `usize`), not interpolated into strings.

## Type safety

BAM files are full of magic numbers: flag bits (0x4 = unmapped, 0x10 = reverse strand), CIGAR operation codes (0 = M, 1 = I, ...), format constants. Using raw literals is error-prone and hard to read.

r[io.named_constants]
Magic numbers (bit flags, format constants, operation codes) MUST be defined as named constants, not used as raw literals. BAM flags MUST use named constants (e.g., `FLAG_UNMAPPED` instead of `0x4`). CIGAR operation types MUST use named constants (e.g., `CIGAR_M` instead of `0`).

r[io.typed_flags]
BAM flags MUST have a `BamFlags` newtype wrapping `u16` with named predicate methods (e.g., `is_unmapped()`, `is_reverse()`). The public API (`RecordRef`) MUST expose `BamFlags` rather than raw `u16` where appropriate.

r[io.typed_cigar_ops]
CIGAR operations MUST have a `CigarOpType` enum with variants for all 9 SAM-spec operations (M, I, D, N, S, H, P, =, X). The enum MUST provide `consumes_ref()` and `consumes_query()` methods. Invalid operation codes MUST be represented as `None` via `from_bam(u8) -> Option<CigarOpType>`.

## Platform portability

r[io.platform_optimizations]
Every platform-specific optimization (SIMD, architecture-specific intrinsics) MUST have a scalar fallback that produces identical results on all platforms. The optimized and fallback paths MUST be tested with property-based tests that verify equivalence for arbitrary inputs.
