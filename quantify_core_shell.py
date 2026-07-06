#!/usr/bin/env python3
import argparse
import csv
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


Interval = Tuple[int, int]
ChromIntervals = Dict[str, List[Interval]]


@dataclass
class WindowClassification:
    primary: ChromIntervals
    core: ChromIntervals
    shell: ChromIntervals
    primary_window_n: int
    core_window_n: int
    shell_window_n: int
    primary_bp: int
    core_bp: int
    shell_bp: int


class WindowIndex:
    def __init__(self, intervals_by_chrom: ChromIntervals):
        self.intervals_by_chrom = intervals_by_chrom
        self.ends_by_chrom = {
            chrom: [end for _, end in intervals]
            for chrom, intervals in intervals_by_chrom.items()
        }

    def overlaps(self, chrom: str, start: int, end: int) -> List[Interval]:
        intervals = self.intervals_by_chrom.get(chrom, [])
        if not intervals:
            return []
        ends = self.ends_by_chrom[chrom]
        idx = bisect_right(ends, start)
        out: List[Interval] = []
        while idx < len(intervals) and intervals[idx][0] < end:
            ov_start = max(start, intervals[idx][0])
            ov_end = min(end, intervals[idx][1])
            if ov_start < ov_end:
                out.append((ov_start, ov_end))
            idx += 1
        return out


class MergedBedWriter:
    def __init__(self, path: Path):
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = self.path.open("w")
        self.last_chrom: Optional[str] = None
        self.last_start: Optional[int] = None
        self.last_end: Optional[int] = None

    def write(self, chrom: str, start: int, end: int) -> None:
        if start >= end:
            return
        if self.last_chrom == chrom and self.last_end is not None and start <= self.last_end:
            if end > self.last_end:
                self.last_end = end
            return
        self.flush()
        self.last_chrom = chrom
        self.last_start = start
        self.last_end = end

    def flush(self) -> None:
        if self.last_chrom is None:
            return
        self.fh.write(f"{self.last_chrom}\t{self.last_start}\t{self.last_end}\n")
        self.last_chrom = None
        self.last_start = None
        self.last_end = None

    def close(self) -> None:
        self.flush()
        self.fh.close()


@dataclass
class ReferenceSpec:
    projection: str
    display_name: str
    window_class_tsv: Path
    callable_bed: Path
    watch_chroms: Sequence[str]


@dataclass
class CallableSummary:
    callable_core_bp: int = 0
    callable_shell_bp: int = 0
    callable_outside_primary_bp: int = 0
    watch_callable_shell_bp: int = 0


FMT_DIGITS = 4


def fmt_pct(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{100.0 * numerator / denominator:.{FMT_DIGITS}f}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Explicitly quantify reference-specific coordinate core, reference-sensitive shell, "
            "and callable overlap from dual-reference window classifications and callable BEDs."
        )
    )
    parser.add_argument(
        "--callable-table",
        required=True,
        help="TSV with total_bp/callable_bp keyed by projection (e.g. results/ragtag_compare/tableS_callable.tsv)",
    )
    parser.add_argument(
        "--reference",
        action="append",
        required=True,
        help=(
            "Reference spec: projection::display_name::window_classification_primary.tsv::callable_scaffold.bed"
            "::watch_chr1,watch_chr2"
        ),
    )
    parser.add_argument(
        "--bed-root",
        required=True,
        help="Directory for explicit BED outputs (primary/core/shell/callable overlaps)",
    )
    parser.add_argument(
        "--summary-output",
        required=True,
        help="TSV summary output",
    )
    return parser.parse_args()


def parse_reference_specs(items: Sequence[str]) -> List[ReferenceSpec]:
    specs: List[ReferenceSpec] = []
    for item in items:
        parts = item.split("::")
        if len(parts) < 4 or len(parts) > 5:
            raise SystemExit(f"Invalid --reference spec: {item}")
        projection, display_name, window_class_tsv, callable_bed = parts[:4]
        watch_raw = parts[4] if len(parts) == 5 else ""
        watch_chroms = [x for x in watch_raw.split(",") if x]
        specs.append(
            ReferenceSpec(
                projection=projection,
                display_name=display_name,
                window_class_tsv=Path(window_class_tsv),
                callable_bed=Path(callable_bed),
                watch_chroms=watch_chroms,
            )
        )
    return specs


def load_callable_totals(path: Path) -> Dict[str, Dict[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            projection = row["projection"]
            out[projection] = {
                "total_bp": int(float(row["total_bp"])),
                "callable_bp": int(float(row["callable_bp"])),
            }
    return out


def load_window_classification(path: Path) -> WindowClassification:
    primary: ChromIntervals = {}
    core: ChromIntervals = {}
    shell: ChromIntervals = {}
    primary_window_n = 0
    core_window_n = 0
    shell_window_n = 0
    primary_bp = 0
    core_bp = 0
    shell_bp = 0

    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row["query_chrom"]
            start = int(row["window_start"])
            end = int(row["window_end"])
            status = row["status"]
            primary.setdefault(chrom, []).append((start, end))
            primary_window_n += 1
            primary_bp += end - start
            if status == "one_to_one":
                core.setdefault(chrom, []).append((start, end))
                core_window_n += 1
                core_bp += end - start
            elif status in {"one_to_many", "unmapped"}:
                shell.setdefault(chrom, []).append((start, end))
                shell_window_n += 1
                shell_bp += end - start
            else:
                raise SystemExit(f"Unexpected window status in {path}: {status}")

    return WindowClassification(
        primary=primary,
        core=core,
        shell=shell,
        primary_window_n=primary_window_n,
        core_window_n=core_window_n,
        shell_window_n=shell_window_n,
        primary_bp=primary_bp,
        core_bp=core_bp,
        shell_bp=shell_bp,
    )


def write_interval_bed(intervals_by_chrom: ChromIntervals, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as out:
        for chrom in sorted(intervals_by_chrom):
            for start, end in intervals_by_chrom[chrom]:
                out.write(f"{chrom}\t{start}\t{end}\n")


def sum_intervals(intervals: Iterable[Interval]) -> int:
    return sum(end - start for start, end in intervals)


def subtract_segments(start: int, end: int, overlaps: Sequence[Interval]) -> List[Interval]:
    if not overlaps:
        return [(start, end)]
    out: List[Interval] = []
    cursor = start
    for ov_start, ov_end in overlaps:
        if cursor < ov_start:
            out.append((cursor, ov_start))
        if ov_end > cursor:
            cursor = ov_end
    if cursor < end:
        out.append((cursor, end))
    return out


def scan_callable(
    callable_bed: Path,
    primary_idx: WindowIndex,
    core_idx: WindowIndex,
    shell_idx: WindowIndex,
    watch_chroms: Sequence[str],
    callable_core_path: Path,
    callable_shell_path: Path,
    callable_outside_primary_path: Path,
    watch_shell_path: Optional[Path] = None,
) -> CallableSummary:
    watch_set = set(watch_chroms)
    core_writer = MergedBedWriter(callable_core_path)
    shell_writer = MergedBedWriter(callable_shell_path)
    outside_writer = MergedBedWriter(callable_outside_primary_path)
    watch_writer = MergedBedWriter(watch_shell_path) if watch_shell_path is not None else None
    summary = CallableSummary()

    with callable_bed.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])

            primary_ovs = primary_idx.overlaps(chrom, start, end)
            core_ovs = core_idx.overlaps(chrom, start, end)
            shell_ovs = shell_idx.overlaps(chrom, start, end)

            for ov_start, ov_end in core_ovs:
                core_writer.write(chrom, ov_start, ov_end)
                summary.callable_core_bp += ov_end - ov_start

            for ov_start, ov_end in shell_ovs:
                shell_writer.write(chrom, ov_start, ov_end)
                summary.callable_shell_bp += ov_end - ov_start
                if chrom in watch_set:
                    summary.watch_callable_shell_bp += ov_end - ov_start
                    if watch_writer is not None:
                        watch_writer.write(chrom, ov_start, ov_end)

            for ov_start, ov_end in subtract_segments(start, end, primary_ovs):
                outside_writer.write(chrom, ov_start, ov_end)
                summary.callable_outside_primary_bp += ov_end - ov_start

    core_writer.close()
    shell_writer.close()
    outside_writer.close()
    if watch_writer is not None:
        watch_writer.close()
    return summary


def build_summary_row(
    spec: ReferenceSpec,
    windows: WindowClassification,
    callable_totals: Dict[str, int],
    callable_summary: CallableSummary,
    watch_shell_bp: int,
) -> Dict[str, str]:
    total_bp = callable_totals["total_bp"]
    callable_total_bp = callable_totals["callable_bp"]
    excluded_nonprimary_bp = total_bp - windows.primary_bp
    callable_primary_bp = callable_summary.callable_core_bp + callable_summary.callable_shell_bp
    return {
        "projection": spec.projection,
        "reference": spec.display_name,
        "total_bp": str(total_bp),
        "primary_window_n": str(windows.primary_window_n),
        "primary_window_bp": str(windows.primary_bp),
        "primary_window_pct_of_total": fmt_pct(windows.primary_bp, total_bp),
        "coordinate_core_window_n": str(windows.core_window_n),
        "coordinate_core_bp": str(windows.core_bp),
        "coordinate_core_pct_of_primary": fmt_pct(windows.core_bp, windows.primary_bp),
        "coordinate_shell_window_n": str(windows.shell_window_n),
        "coordinate_shell_bp": str(windows.shell_bp),
        "coordinate_shell_pct_of_primary": fmt_pct(windows.shell_bp, windows.primary_bp),
        "excluded_nonprimary_bp": str(excluded_nonprimary_bp),
        "excluded_nonprimary_pct_of_total": fmt_pct(excluded_nonprimary_bp, total_bp),
        "callable_total_bp": str(callable_total_bp),
        "callable_total_pct_of_total": fmt_pct(callable_total_bp, total_bp),
        "callable_primary_bp": str(callable_primary_bp),
        "callable_primary_pct_of_callable": fmt_pct(callable_primary_bp, callable_total_bp),
        "callable_core_bp": str(callable_summary.callable_core_bp),
        "callable_core_pct_of_total": fmt_pct(callable_summary.callable_core_bp, total_bp),
        "callable_core_pct_of_callable": fmt_pct(callable_summary.callable_core_bp, callable_total_bp),
        "callable_shell_bp": str(callable_summary.callable_shell_bp),
        "callable_shell_pct_of_total": fmt_pct(callable_summary.callable_shell_bp, total_bp),
        "callable_shell_pct_of_callable": fmt_pct(callable_summary.callable_shell_bp, callable_total_bp),
        "callable_outside_primary_bp": str(callable_summary.callable_outside_primary_bp),
        "callable_outside_primary_pct_of_total": fmt_pct(callable_summary.callable_outside_primary_bp, total_bp),
        "callable_outside_primary_pct_of_callable": fmt_pct(callable_summary.callable_outside_primary_bp, callable_total_bp),
        "watch_shell_bp": str(watch_shell_bp),
        "watch_shell_pct_of_shell": fmt_pct(watch_shell_bp, windows.shell_bp),
        "watch_callable_shell_bp": str(callable_summary.watch_callable_shell_bp),
        "watch_callable_shell_pct_of_callable_shell": fmt_pct(callable_summary.watch_callable_shell_bp, callable_summary.callable_shell_bp),
        "note": (
            "coordinate core = one-to-one pseudomolecule windows; "
            "coordinate shell = one-to-many + unmapped pseudomolecule windows; "
            "callable overlaps are reported in scaffold coordinates"
        ),
    }


def main() -> None:
    args = parse_args()
    callable_table = load_callable_totals(Path(args.callable_table))
    specs = parse_reference_specs(args.reference)
    bed_root = Path(args.bed_root)
    bed_root.mkdir(parents=True, exist_ok=True)
    summary_rows: List[Dict[str, str]] = []

    fieldnames = [
        "projection",
        "reference",
        "total_bp",
        "primary_window_n",
        "primary_window_bp",
        "primary_window_pct_of_total",
        "coordinate_core_window_n",
        "coordinate_core_bp",
        "coordinate_core_pct_of_primary",
        "coordinate_shell_window_n",
        "coordinate_shell_bp",
        "coordinate_shell_pct_of_primary",
        "excluded_nonprimary_bp",
        "excluded_nonprimary_pct_of_total",
        "callable_total_bp",
        "callable_total_pct_of_total",
        "callable_primary_bp",
        "callable_primary_pct_of_callable",
        "callable_core_bp",
        "callable_core_pct_of_total",
        "callable_core_pct_of_callable",
        "callable_shell_bp",
        "callable_shell_pct_of_total",
        "callable_shell_pct_of_callable",
        "callable_outside_primary_bp",
        "callable_outside_primary_pct_of_total",
        "callable_outside_primary_pct_of_callable",
        "watch_shell_bp",
        "watch_shell_pct_of_shell",
        "watch_callable_shell_bp",
        "watch_callable_shell_pct_of_callable_shell",
        "note",
    ]

    for spec in specs:
        if spec.projection not in callable_table:
            raise SystemExit(
                f"Projection {spec.projection} not found in callable table {args.callable_table}"
            )
        if not spec.window_class_tsv.is_file():
            raise SystemExit(f"Missing window classification TSV: {spec.window_class_tsv}")
        if not spec.callable_bed.is_file():
            raise SystemExit(f"Missing callable BED: {spec.callable_bed}")

        ref_dir = bed_root / spec.display_name
        ref_dir.mkdir(parents=True, exist_ok=True)

        windows = load_window_classification(spec.window_class_tsv)
        write_interval_bed(windows.primary, ref_dir / f"{spec.display_name}.primary_windows.bed")
        write_interval_bed(windows.core, ref_dir / f"{spec.display_name}.coordinate_core_windows.bed")
        write_interval_bed(windows.shell, ref_dir / f"{spec.display_name}.coordinate_shell_windows.bed")

        watch_shell_bp = 0
        if spec.watch_chroms:
            watch_shell: ChromIntervals = {}
            for chrom in spec.watch_chroms:
                intervals = windows.shell.get(chrom, [])
                if intervals:
                    watch_shell[chrom] = intervals
                    watch_shell_bp += sum_intervals(intervals)
            if watch_shell:
                write_interval_bed(watch_shell, ref_dir / f"{spec.display_name}.watch_shell_windows.bed")

        callable_summary = scan_callable(
            callable_bed=spec.callable_bed,
            primary_idx=WindowIndex(windows.primary),
            core_idx=WindowIndex(windows.core),
            shell_idx=WindowIndex(windows.shell),
            watch_chroms=spec.watch_chroms,
            callable_core_path=ref_dir / f"{spec.display_name}.callable_core.bed",
            callable_shell_path=ref_dir / f"{spec.display_name}.callable_shell.bed",
            callable_outside_primary_path=ref_dir / f"{spec.display_name}.callable_outside_primary.bed",
            watch_shell_path=(ref_dir / f"{spec.display_name}.watch_callable_shell.bed") if spec.watch_chroms else None,
        )

        summary_rows.append(
            build_summary_row(
                spec=spec,
                windows=windows,
                callable_totals=callable_table[spec.projection],
                callable_summary=callable_summary,
                watch_shell_bp=watch_shell_bp,
            )
        )

    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)


if __name__ == "__main__":
    main()
