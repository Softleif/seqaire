# htslib test data

Test files copied from the htslib test suite for cross-validation.

- **Source**: https://github.com/samtools/htslib
- **Commit**: `ca58677156180196aa74091f39619389d81cbf32`
- **License**: MIT/Expat (see htslib LICENSE file)
- **Date copied**: 2026-04-15

## Contents

| Directory  | Files | Description                                    |
|------------|------:|------------------------------------------------|
| `sam/`     |    45 | SAM edge cases (unmapped, clipping, padding, aux tags, etc.) |
| `bam/`     |     6 | BAM files (colons in names, range queries, missing SQ) |
| `cram/`    |     6 | Java-encoded CRAMs, range query CRAM + indexes |
| `fasta/`   |    12 | Reference sequences + FAI indexes              |
| `vcf/`     |    12 | VCF/BCF edge cases (missing values, headers, VCF 4.4) |
| `bgzf/`    |     2 | BGZF compressed test file + GZI index          |
| `index/`   |     7 | BAI, CSI, TBI index files                      |
| `mpileup/` |    26 | Pileup SAM inputs + expected outputs           |
| `longrefs/`|     6 | Long reference (>2GB position) test data       |
| `tlen/`    |    63 | Template length CRAM/SAM pairs                 |
