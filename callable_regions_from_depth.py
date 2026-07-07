#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import statistics
import sys
from collections import Counter
from typing import Iterable, Iterator, Tuple


def iter_depth_lines(handle: Iterable[str]) -> Iterator[Tuple[str, int, int]]:
    for raw in handle:
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split("\t")
        if len(parts) < 3:
            raise ValueError(f"Malformed depth line: {raw!r}")
        chrom = parts[0]
        pos = int(parts[1])
        depth = int(parts[2])
        yield chrom, pos, depth


def compute_center(method: str) -> None:
    depth_hist: Counter[int] = Counter()
    total_bp = 0
    total_depth = 0
    nonzero_bp = 0
    nonzero_depth = 0

    for _chrom, _pos, depth in iter_depth_lines(sys.stdin):
        total_bp += 1
        total_depth += depth
        depth_hist[depth] += 1
        if depth > 0:
            nonzero_bp += 1
            nonzero_depth += depth

    if total_bp == 0:
        raise SystemExit("No depth records received on stdin")

    if method == "nonzero_median":
        if nonzero_bp == 0:
            center = 0
        else:
            target = (nonzero_bp + 1) // 2
            seen = 0
            center = 0
            for depth in sorted(d for d in depth_hist if d > 0):
                seen += depth_hist[depth]
                if seen >= target:
                    center = depth
                    break
    elif method == "nonzero_mean":
        center = int(round(nonzero_depth / nonzero_bp)) if nonzero_bp else 0
    elif method == "all_mean":
        center = int(round(total_depth / total_bp)) if total_bp else 0
    else:
        raise SystemExit(f"Unsupported center method: {method}")

    mean_all = total_depth / total_bp
    mean_nonzero = (nonzero_depth / nonzero_bp) if nonzero_bp else 0.0
    zero_bp = total_bp - nonzero_bp

    rows = [
        ("depth_center_method", method),
        ("depth_center_x", center),
        ("total_bp", total_bp),
        ("nonzero_bp", nonzero_bp),
        ("zero_bp", zero_bp),
        ("mean_all_depth_x", f"{mean_all:.6f}"),
        ("mean_nonzero_depth_x", f"{mean_nonzero:.6f}"),
    ]
    print("metric\tvalue")
    for key, value in rows:
        print(f"{key}\t{value}")


def emit_bed(low: int, high: int) -> None:
    current_chrom = None
    start = None
    prev_pos = None

    def flush(chrom: str | None, start0: int | None, prev1: int | None) -> None:
        if chrom is None or start0 is None or prev1 is None:
            return
        sys.stdout.write(f"{chrom}\t{start0}\t{prev1}\n")

    for chrom, pos, depth in iter_depth_lines(sys.stdin):
        in_range = low <= depth <= high
        if in_range:
            start0 = pos - 1
            if current_chrom == chrom and start is not None and prev_pos is not None and pos == prev_pos + 1:
                prev_pos = pos
            else:
                flush(current_chrom, start, prev_pos)
                current_chrom = chrom
                start = start0
                prev_pos = pos
        else:
            flush(current_chrom, start, prev_pos)
            current_chrom = None
            start = None
            prev_pos = None

    flush(current_chrom, start, prev_pos)


def main() -> None:
    ap = argparse.ArgumentParser(description="Derive callable-region summaries from samtools depth output")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_center = sub.add_parser("center", help="estimate a depth center from samtools depth stdin")
    ap_center.add_argument("--method", default="nonzero_median", choices=["nonzero_median", "nonzero_mean", "all_mean"])

    ap_bed = sub.add_parser("bed", help="emit callable BED from samtools depth stdin")
    ap_bed.add_argument("--low", type=int, required=True)
    ap_bed.add_argument("--high", type=int, required=True)

    args = ap.parse_args()
    if args.cmd == "center":
        compute_center(args.method)
    elif args.cmd == "bed":
        emit_bed(args.low, args.high)
    else:
        raise SystemExit(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    main()
