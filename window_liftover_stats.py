#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(
        description="Classify fixed windows from a query assembly against a target assembly using a PAF file."
    )
    ap.add_argument("--windows", required=True, help="BED file of fixed query windows")
    ap.add_argument("--paf", required=True, help="PAF alignment file; query must match the BED assembly")
    ap.add_argument("--pair-label", required=True, help="Output label for the comparison")
    ap.add_argument("--window-size", required=True, type=int, help="Nominal window size used to make the BED")
    ap.add_argument(
        "--min-overlap-frac",
        type=float,
        default=0.5,
        help="Minimum fraction of the nominal window size required to count as a hit",
    )
    ap.add_argument("--summary-out", required=True, help="TSV file for one-line summary output")
    ap.add_argument("--per-window-out", required=True, help="TSV file for per-window classifications")
    return ap.parse_args()


def load_windows(path):
    by_chrom = defaultdict(list)
    ordered = []
    with open(path) as fh:
        for line in fh:
            chrom, start, end = line.rstrip().split("\t")[:3]
            start = int(start)
            end = int(end)
            rec = (chrom, start, end)
            by_chrom[chrom].append(rec)
            ordered.append(rec)
    return by_chrom, ordered


def main():
    args = parse_args()
    min_overlap = int(args.window_size * args.min_overlap_frac)
    windows_by_chrom, ordered_windows = load_windows(args.windows)
    hits = defaultdict(set)

    with open(args.paf) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip().split("\t")
            qname = f[0]
            qstart = int(f[2])
            qend = int(f[3])
            tname = f[5]
            tstart = int(f[7])
            tend = int(f[8])
            if qname not in windows_by_chrom:
                continue
            first_idx = qstart // args.window_size
            last_idx = max((qend - 1) // args.window_size, first_idx)
            chrom_windows = windows_by_chrom[qname]
            if first_idx >= len(chrom_windows):
                continue
            last_idx = min(last_idx, len(chrom_windows) - 1)
            target_interval = f"{tname}:{tstart}-{tend}"
            for idx in range(first_idx, last_idx + 1):
                chrom, wstart, wend = chrom_windows[idx]
                overlap = min(wend, qend) - max(wstart, qstart)
                if overlap >= min_overlap:
                    hits[(chrom, wstart, wend)].add(target_interval)

    total = len(ordered_windows)
    one_to_one = 0
    one_to_many = 0
    unmapped = 0
    with open(args.per_window_out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["query_chrom", "window_start", "window_end", "status", "n_target_intervals", "target_intervals"])
        for rec in ordered_windows:
            targets = sorted(hits.get(rec, set()))
            if not targets:
                status = "unmapped"
                unmapped += 1
            elif len(targets) == 1:
                status = "one_to_one"
                one_to_one += 1
            else:
                status = "one_to_many"
                one_to_many += 1
            w.writerow([rec[0], rec[1], rec[2], status, len(targets), ",".join(targets)])

    with open(args.summary_out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(
            [
                "projection_pair",
                "window_size_bp",
                "total_windows",
                "one_to_one_windows",
                "one_to_many_windows",
                "unmapped_windows",
                "one_to_one_pct",
                "one_to_many_pct",
                "unmapped_pct",
            ]
        )
        w.writerow(
            [
                args.pair_label,
                args.window_size,
                total,
                one_to_one,
                one_to_many,
                unmapped,
                f"{(100 * one_to_one / total) if total else 0:.4f}",
                f"{(100 * one_to_many / total) if total else 0:.4f}",
                f"{(100 * unmapped / total) if total else 0:.4f}",
            ]
        )


if __name__ == "__main__":
    main()
