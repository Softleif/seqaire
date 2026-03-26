#!/usr/bin/env python3
"""Analyze seqair profile JSON output.

Usage:
    RASTAIR_PROFILE_JSON=/tmp/profile.jsonl rastair call ...
    python3 tools/analyze_profile.py /tmp/profile.jsonl
"""

import json
import sys
from collections import defaultdict


def percentile(values: list[float], p: float) -> float:
    if not values:
        return 0.0
    k = (len(values) - 1) * p
    f = int(k)
    c = f + 1
    if c >= len(values):
        return values[-1]
    return values[f] + (k - f) * (values[c] - values[f])


def fmt_bytes(n: float) -> str:
    if n >= 1024 * 1024 * 1024:
        return f"{n / (1024 * 1024 * 1024):.2f} GB"
    if n >= 1024 * 1024:
        return f"{n / (1024 * 1024):.1f} MB"
    if n >= 1024:
        return f"{n / 1024:.1f} KB"
    return f"{n:.0f} B"


def fmt_num(n: float) -> str:
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f}M"
    if n >= 1_000:
        return f"{n / 1_000:.1f}K"
    return f"{n:.0f}"


def fmt_pct(n: float) -> str:
    return f"{n:.1f}%"


def distribution(values: list[float], unit_fn=str) -> str:
    if not values:
        return "  (no data)"
    values.sort()
    lines = []
    lines.append(f"  n={len(values)}")
    lines.append(
        f"  min={unit_fn(values[0])}  p25={unit_fn(percentile(values, 0.25))}  "
        f"median={unit_fn(percentile(values, 0.5))}  p75={unit_fn(percentile(values, 0.75))}  "
        f"p95={unit_fn(percentile(values, 0.95))}  max={unit_fn(values[-1])}"
    )
    lines.append(
        f"  total={unit_fn(sum(values))}  mean={unit_fn(sum(values) / len(values))}"
    )
    return "\n".join(lines)


def analyze(events: list[dict]) -> None:
    by_type: dict[str, list[dict]] = defaultdict(list)
    for e in events:
        by_type[e["fields"]["message"]].append(e["fields"])

    n_regions = len(by_type.get("region_buf::load", []))
    print(f"=== seqair Profile ({n_regions} regions) ===\n")

    # --- Bin 0 cache ---
    cache_loads = by_type.get("chunk_cache::load", [])
    if cache_loads:
        print("--- Chunk Cache (distant bins, levels 0-2) ---")
        print(f"  Loaded {len(cache_loads)} times (once per thread per tid).")
        records_per = [f["cached_records"] for f in cache_loads]
        if records_per:
            print(f"  Records cached: {records_per[0]} per chromosome")
        print()

    # --- RegionBuf Load ---
    # Exclude bin 0 cache loads (they're small and one-time)
    all_loads = by_type.get("region_buf::load", [])
    loads = [l for l in all_loads if l["total_read_bytes"] > 200_000]
    if not loads:
        loads = all_loads
    if loads:
        print("--- RegionBuf I/O ---")
        print(f"  How data is read from the BAM file into memory.\n")

        print("  Merged ranges per region (each = 1 seek + read):")
        print(distribution([f["merged_ranges"] for f in loads], fmt_num))
        print()

        print("  Input index chunks per region:")
        print(distribution([f["input_chunks"] for f in loads], fmt_num))
        print()

        print("  Total bytes read per region:")
        print(distribution([f["total_read_bytes"] for f in loads], fmt_bytes))
        print()

        print("  File span (distance between first and last range):")
        spans = [f.get("file_span", 0) for f in loads]
        print(distribution(spans, fmt_bytes))
        print()

        # Check if there's a consistent distant offset
        last_offsets = [f.get("last_offset", 0) for f in loads]
        if last_offsets and len(set(last_offsets)) <= 3:
            print(
                f"  ** All regions reach the same distant offset: {last_offsets[0]:,}"
            )
            first_offsets = [f.get("first_offset", 0) for f in loads]
            print(
                f"     Main data cluster: {min(first_offsets):,} - {max(first_offsets):,}"
            )
            print(f"     This is likely BAI bin 0 (whole-chromosome catch-all).")
            print()

        print("  Load time per region (us):")
        print(distribution([f["elapsed_us"] for f in loads], fmt_num))
        print()

        print("  Slowest single range read (us):")
        print(distribution([f["max_range_us"] for f in loads], fmt_num))
        print()

    # --- Fetch record filtering ---
    fetches = by_type.get("fetch_into", [])
    if fetches:
        print("--- Record Filtering (fetch_into) ---")
        print(f"  Records decoded from BAM vs accepted into RecordStore.\n")

        total_accepted = sum(f["accepted"] for f in fetches)
        total_tid = sum(f["skipped_tid"] for f in fetches)
        total_unmapped = sum(f["skipped_unmapped"] for f in fetches)
        total_range = sum(f["skipped_out_of_range"] for f in fetches)
        total_decoded = total_accepted + total_tid + total_unmapped + total_range

        print(f"  Total decoded:  {fmt_num(total_decoded)}")
        print(
            f"  Accepted:       {fmt_num(total_accepted)} ({total_accepted / total_decoded * 100:.1f}%)"
        )
        if total_tid > 0:
            print(
                f"  Skipped (tid):  {fmt_num(total_tid)} ({total_tid / total_decoded * 100:.1f}%)"
            )
        if total_unmapped > 0:
            print(
                f"  Skipped (unmap):{fmt_num(total_unmapped)} ({total_unmapped / total_decoded * 100:.1f}%)"
            )
        if total_range > 0:
            print(
                f"  Skipped (range):{fmt_num(total_range)} ({total_range / total_decoded * 100:.1f}%)"
            )

        total_bin0 = sum(f.get("cache_injected", 0) for f in fetches)
        if total_bin0 > 0:
            print(
                f"  From bin0 cache:{fmt_num(total_bin0)} ({total_bin0 / total_decoded * 100:.1f}%)"
            )
        print()

        waste_pct = [
            (
                f["skipped_out_of_range"]
                / (
                    f["accepted"]
                    + f["skipped_tid"]
                    + f["skipped_unmapped"]
                    + f["skipped_out_of_range"]
                )
            )
            * 100
            for f in fetches
        ]
        print("  Out-of-range waste per region:")
        print(distribution(waste_pct, fmt_pct))
        print()

    # --- RegionBuf Decompression ---
    summaries = by_type.get("region_buf summary", [])
    if summaries:
        print("--- BGZF Decompression ---")

        print("  Blocks decompressed per region:")
        print(distribution([f["blocks"] for f in summaries], fmt_num))
        print()

        ratios = []
        for f in summaries:
            if f["compressed_bytes"] > 0:
                ratios.append(f["decompressed_bytes"] / f["compressed_bytes"])
        if ratios:
            ratios.sort()
            print(
                f"  Compression ratio: "
                f"min={ratios[0]:.2f}x  median={percentile(ratios, 0.5):.2f}x  "
                f"max={ratios[-1]:.2f}x"
            )
            print()

        print("  Max gap between disjoint ranges:")
        print(distribution([f["max_gap_bytes"] for f in summaries], fmt_bytes))
        print()

    # --- RecordStore ---
    stores = by_type.get("record_store", [])
    if stores:
        print("--- RecordStore (per region) ---")
        print(f"  Slab-based storage for decoded BAM records.\n")

        print("  Records stored:")
        print(distribution([f["records"] for f in stores], fmt_num))
        print()

        print("  Memory breakdown (capacity = allocated, including slack):")
        for slab_name, slab_key, elem_size in [
            ("records", "records_cap", 48),  # SlimRecord is ~48 bytes
            ("names", "names_bytes", 1),
            ("bases", "bases_bytes", 1),
            ("data", "data_bytes", 1),
        ]:
            values = [f[slab_key] * elem_size for f in stores]
            print(f"    {slab_name}: {distribution(values, fmt_bytes).strip()}")

        total_per_region = []
        for f in stores:
            total = (
                f["records_cap"] * 48
                + f["names_bytes"]
                + f["bases_bytes"]
                + f["data_bytes"]
            )
            total_per_region.append(total)
        print(f"\n  Total memory per region:")
        print(distribution(total_per_region, fmt_bytes))
        print()

    # --- PileupEngine ---
    pileups = by_type.get("pileup_engine", [])
    if pileups:
        print("--- PileupEngine (per region) ---")

        print("  Columns produced:")
        print(distribution([f["columns"] for f in pileups], fmt_num))
        print()

        print("  Max active depth (reads overlapping a single position):")
        print(distribution([f["max_depth"] for f in pileups], fmt_num))
        print()

        print("  Active vec capacity (peak allocation):")
        print(distribution([f["active_cap"] for f in pileups], fmt_num))
        print()

        dedup_count = sum(1 for f in pileups if f["dedup"])
        print(f"  Dedup enabled: {dedup_count}/{len(pileups)} regions")
        if dedup_count > 0:
            print("  Dedup scratch capacities:")
            print(
                f"    remove_buf: "
                + distribution(
                    [f["dedup_remove_cap"] for f in pileups if f["dedup"]], fmt_num
                ).strip()
            )
            print(
                f"    seen_map:   "
                + distribution(
                    [f["dedup_seen_cap"] for f in pileups if f["dedup"]], fmt_num
                ).strip()
            )
        print()

    # --- Cross-cutting: memory per thread ---
    if stores and summaries:
        print("--- Per-Region Total Memory ---")
        print(f"  Includes compressed I/O buffer + RecordStore slabs.\n")

        combined = []
        for i in range(min(len(stores), len(summaries))):
            s = stores[i]
            r = summaries[i]
            io_buf = r["compressed_bytes"]
            store_mem = (
                s["records_cap"] * 48
                + s["names_bytes"]
                + s["bases_bytes"]
                + s["data_bytes"]
            )
            combined.append(io_buf + store_mem)
        print(distribution(combined, fmt_bytes))
        print()


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <profile.jsonl>", file=sys.stderr)
        sys.exit(1)

    path = sys.argv[1]
    with open(path) as f:
        events = [json.loads(line) for line in f if line.strip()]

    if not events:
        print(f"No events found in {path}", file=sys.stderr)
        sys.exit(1)

    analyze(events)


if __name__ == "__main__":
    main()
