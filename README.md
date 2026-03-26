# seqair

Pure-Rust BAM/SAM/CRAM/FASTA reader with a pileup engine. I/O backend for [rastair](https://github.com/bsblabludwig/rastair).

This repo has two crates:

- [`seqair`](crates/seqair): Indexed readers and pileup engine
- [`seqair-types`](crates/seqair-types): Core types: `Base`, `Strand`, `Phred`, `Probability`, `RegionString`

## Highlights

- Spec-driven development using [tracey](https://tracey.bearcove.eu/)
- Slab-based `RecordStore` — zero per-record heap allocation
- `RegionBuf` bulk I/O — one large read per region, then decompress from memory (good for NFS/Lustre)
- CRAM v3.0/v3.1 (rANS, tok3, bzip2, lzma, multi-ref slices, embedded refs, MD5 verification)
- Overlapping read-pair dedup built into the pileup engine
- `fork()` for cheap per-thread readers sharing index and header via `Arc`
- `ChunkCache` for distant BAI bins — loads wide-spanning chunks once per chromosome

## Testing

```sh
cargo test
```

Comparison tests in `crates/seqair/tests/` validate against `rust-htslib` (and sometimes `noodles`).
Property tests use `proptest`.
