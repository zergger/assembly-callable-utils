#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(
        description="Lift BED intervals from AGP component coordinates to object coordinates."
    )
    ap.add_argument("--agp", required=True, help="AGP file describing object/component placement")
    ap.add_argument("--bed", required=True, help="BED file in component coordinates")
    ap.add_argument("--output", required=True, help="Lifted BED file in object coordinates")
    return ap.parse_args()


def load_agp(path):
    placements = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if f[4] != "W":
                continue
            obj = f[0]
            obj_beg = int(f[1])
            obj_end = int(f[2])
            comp = f[5]
            comp_beg = int(f[6])
            comp_end = int(f[7])
            orient = f[8]
            placements[comp].append((obj, obj_beg, obj_end, comp_beg, comp_end, orient))
    return placements


def main():
    args = parse_args()
    placements = load_agp(args.agp)
    with open(args.output, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        with open(args.bed) as fh:
            for line in fh:
                chrom, start, end, *rest = line.rstrip().split("\t")
                start = int(start)
                end = int(end)
                for obj, obj_beg, obj_end, comp_beg, comp_end, orient in placements.get(chrom, []):
                    comp_zero_start = comp_beg - 1
                    comp_zero_end = comp_end
                    inter_start = max(start, comp_zero_start)
                    inter_end = min(end, comp_zero_end)
                    if inter_start >= inter_end:
                        continue
                    if orient == "+":
                        obj_start = (obj_beg - 1) + (inter_start - comp_zero_start)
                        obj_end_excl = obj_start + (inter_end - inter_start)
                    else:
                        obj_start = obj_end - (inter_end - comp_zero_start)
                        obj_end_excl = obj_end - (inter_start - comp_zero_start)
                    row = [obj, obj_start, obj_end_excl]
                    if rest:
                        row.extend(rest)
                    w.writerow(row)


if __name__ == "__main__":
    main()
