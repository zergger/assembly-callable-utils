#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULT_CONFIG="$ROOT_DIR/configs/callable_assembly.env"

usage() {
  cat <<EOF
Usage: bash scripts/run_callable_assembly.sh [config.env]

Materialize a callable-constrained assembly view from a scaffold FASTA plus
an existing callable BED or a BAM-derived callable-region calculation.

Environment/config variables:
  SAMPLE_TAG            output tag; default from FASTA basename
  MODE                  scaffold_bam|corrected_query_agp (default: scaffold_bam)
  FASTA                 required scaffold FASTA
  BAM                   BAM used to derive callable regions (required if CALLABLE_BED unset)
  AGP                   required in corrected_query_agp mode; lifts query BED to scaffold coordinates
  RAGTAG_CORRECT_LOG    optional explicit log path for auto-parsing corrected-query depth center
  CALLABLE_BED          precomputed callable BED (optional; bypasses BAM derivation)
  OUTDIR                output directory (default: results/callable_assembly/<sample>)
  MQ_CUTOFF             samtools depth -Q cutoff (default: 30)
  DEPTH_CENTER_METHOD   nonzero_median|nonzero_mean|all_mean (default: nonzero_median)
  DEPTH_CENTER          explicit center override (optional)
  DEPTH_LOW             explicit low cutoff override (optional)
  DEPTH_HIGH            explicit high cutoff override (optional)
  DEPTH_LOW_MULT        multiplier for low cutoff if DEPTH_LOW unset (default: 0.5)
  DEPTH_HIGH_MULT       multiplier for high cutoff if DEPTH_HIGH unset (default: 2.0)
  MASK_CHAR             masking character for noncallable sequence (default: N)
  OUTPUT_CALLABLE_ONLY  1 to also emit callable-only FASTA intervals (default: 0)
  FORCE                 1 to overwrite existing outputs (default: 0)
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

CONFIG_FILE="${1:-$DEFAULT_CONFIG}"
[[ -f "$CONFIG_FILE" ]] || { echo "[callable-assembly] missing config: $CONFIG_FILE" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG_FILE"

PYTHON_BIN="${PYTHON_BIN:-python3}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
FORCE="${FORCE:-0}"
MODE="${MODE:-scaffold_bam}"
MQ_CUTOFF="${MQ_CUTOFF:-30}"
DEPTH_CENTER_METHOD="${DEPTH_CENTER_METHOD:-nonzero_median}"
DEPTH_LOW_MULT="${DEPTH_LOW_MULT:-0.5}"
DEPTH_HIGH_MULT="${DEPTH_HIGH_MULT:-2.0}"
MASK_CHAR="${MASK_CHAR:-N}"
OUTPUT_CALLABLE_ONLY="${OUTPUT_CALLABLE_ONLY:-0}"

log() {
  printf '[callable-assembly] %s\n' "$*" >&2
}

resolve_path() {
  "$PYTHON_BIN" - "$1" "$ROOT_DIR" <<'PY'
import os, sys
path, root = sys.argv[1:3]
if not path:
    print("")
elif os.path.isabs(path):
    print(os.path.realpath(path))
else:
    print(os.path.realpath(os.path.join(root, path)))
PY
}

require_file() {
  [[ -f "$1" ]] || { log "missing file: $1"; exit 1; }
}

should_refresh() {
  local out="$1"
  [[ "$FORCE" == "1" || ! -s "$out" ]]
}

FASTA="$(resolve_path "${FASTA:-}")"
require_file "$FASTA"
FAI="$FASTA.fai"
if [[ ! -s "$FAI" ]]; then
  log "index FASTA: $FASTA"
  "$SAMTOOLS" faidx "$FASTA"
fi

SAMPLE_TAG="${SAMPLE_TAG:-$(basename "$FASTA")}"
SAMPLE_TAG="${SAMPLE_TAG%%.*}"
OUTDIR="$(resolve_path "${OUTDIR:-results/callable_assembly/$SAMPLE_TAG}")"
mkdir -p "$OUTDIR" "$OUTDIR/stats" "$OUTDIR/bed" "$OUTDIR/fasta"

CALLABLE_BED_RAW="${CALLABLE_BED:-}"
BAM_RAW="${BAM:-}"
AGP_RAW="${AGP:-}"
RAGTAG_CORRECT_LOG_RAW="${RAGTAG_CORRECT_LOG:-}"
CALLABLE_BED=""
BAM=""
AGP=""
RAGTAG_CORRECT_LOG=""
if [[ -n "$CALLABLE_BED_RAW" ]]; then
  CALLABLE_BED="$(resolve_path "$CALLABLE_BED_RAW")"
  require_file "$CALLABLE_BED"
fi
if [[ -n "$BAM_RAW" ]]; then
  BAM="$(resolve_path "$BAM_RAW")"
  require_file "$BAM"
fi
if [[ -n "$AGP_RAW" ]]; then
  AGP="$(resolve_path "$AGP_RAW")"
  require_file "$AGP"
fi
if [[ -n "$RAGTAG_CORRECT_LOG_RAW" ]]; then
  RAGTAG_CORRECT_LOG="$(resolve_path "$RAGTAG_CORRECT_LOG_RAW")"
  require_file "$RAGTAG_CORRECT_LOG"
fi
if [[ -z "$CALLABLE_BED" && -z "$BAM" ]]; then
  log "set CALLABLE_BED or BAM"
  exit 1
fi
if [[ -z "$CALLABLE_BED" && "$MODE" == "corrected_query_agp" && -z "$AGP" ]]; then
  log "AGP is required in corrected_query_agp mode"
  exit 1
fi
if [[ "$MODE" != "scaffold_bam" && "$MODE" != "corrected_query_agp" ]]; then
  log "unsupported MODE: $MODE"
  exit 1
fi

CENTER_TSV="$OUTDIR/stats/${SAMPLE_TAG}.depth_center.tsv"
THRESH_TSV="$OUTDIR/stats/${SAMPLE_TAG}.callable_thresholds.tsv"
RAW_BED="$OUTDIR/bed/${SAMPLE_TAG}.callable.raw.bed"
LIFTED_RAW_BED="$OUTDIR/bed/${SAMPLE_TAG}.callable.lifted.raw.bed"
CALLABLE_OUT="$OUTDIR/bed/${SAMPLE_TAG}.callable.bed"
NONCALLABLE_OUT="$OUTDIR/bed/${SAMPLE_TAG}.noncallable.bed"
MASKED_FASTA="$OUTDIR/fasta/${SAMPLE_TAG}.callable_masked.fa"
CALLABLE_ONLY_FASTA="$OUTDIR/fasta/${SAMPLE_TAG}.callable_only.fa"
SUMMARY_TSV="$OUTDIR/stats/${SAMPLE_TAG}.callable_assembly_summary.tsv"

infer_ragtag_correct_log() {
  local bam_path="$1"
  local workdir
  workdir="$(dirname "$(dirname "$bam_path")")"
  printf '%s\n' "$workdir/logs/ragtag_correct.log"
}

extract_depth_center_from_log() {
  local log_path="$1"
  "$PYTHON_BIN" - "$log_path" <<'PY'
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
pattern = re.compile(r'global median read coverage is\s+([0-9]+(?:\.[0-9]+)?)X', re.IGNORECASE)
for line in path.read_text(encoding='utf-8', errors='ignore').splitlines():
    m = pattern.search(line)
    if m:
        value = m.group(1)
        if '.' in value:
            print(int(float(value)))
        else:
            print(int(value))
        raise SystemExit(0)
raise SystemExit(1)
PY
}

if [[ -n "$CALLABLE_BED" ]]; then
  log "reuse provided callable BED: $CALLABLE_BED"
  if should_refresh "$CALLABLE_OUT"; then
    sort -k1,1 -k2,2n "$CALLABLE_BED" | "$BEDTOOLS" merge -i - > "$CALLABLE_OUT"
  fi
  printf 'metric\tvalue\nsource\tprovided_callable_bed\n' > "$CENTER_TSV"
  {
    printf 'metric\tvalue\n'
    printf 'source\tprovided_callable_bed\n'
    printf 'mode\t%s\n' "$MODE"
    printf 'mq_cutoff\t%s\n' "${MQ_CUTOFF}"
    printf 'depth_center_method\t%s\n' "provided"
    printf 'depth_center_x\t%s\n' "${DEPTH_CENTER:-}"
    printf 'depth_low\t%s\n' "${DEPTH_LOW:-}"
    printf 'depth_high\t%s\n' "${DEPTH_HIGH:-}"
  } > "$THRESH_TSV"
else
  [[ -n "$BAM" ]] || { log "BAM resolution failure"; exit 1; }

  if [[ -z "${DEPTH_CENTER:-}" ]]; then
    if [[ "$MODE" == "corrected_query_agp" ]]; then
      if [[ -z "$RAGTAG_CORRECT_LOG" ]]; then
        inferred_log="$(infer_ragtag_correct_log "$BAM")"
        if [[ -f "$inferred_log" ]]; then
          RAGTAG_CORRECT_LOG="$inferred_log"
        fi
      fi
      if [[ -n "$RAGTAG_CORRECT_LOG" ]]; then
        log "parse depth center from RagTag correct log: $RAGTAG_CORRECT_LOG"
        if ! DEPTH_CENTER="$(extract_depth_center_from_log "$RAGTAG_CORRECT_LOG")"; then
          log "failed to parse DEPTH_CENTER from $RAGTAG_CORRECT_LOG"
          exit 1
        fi
        {
          printf 'metric\tvalue\n'
          printf 'source\tragtag_correct_log\n'
          printf 'mode\t%s\n' "$MODE"
          printf 'depth_center_method\t%s\n' "ragtag_correct_global_median"
          printf 'depth_center_x\t%s\n' "$DEPTH_CENTER"
          printf 'ragtag_correct_log\t%s\n' "$RAGTAG_CORRECT_LOG"
        } > "$CENTER_TSV"
      else
        log "MODE=corrected_query_agp requires DEPTH_CENTER or a parsable ragtag_correct.log"
        exit 1
      fi
    else
      log "estimate depth center from BAM"
      "$SAMTOOLS" depth -aa -Q "$MQ_CUTOFF" -d 0 "$BAM" \
        | "$PYTHON_BIN" "$ROOT_DIR/scripts/callable_regions_from_depth.py" center --method "$DEPTH_CENTER_METHOD" \
        > "$CENTER_TSV"
      DEPTH_CENTER="$(awk -F '\t' '$1=="depth_center_x" {print $2}' "$CENTER_TSV")"
    fi
  else
    {
      printf 'metric\tvalue\n'
      printf 'mode\t%s\n' "$MODE"
      printf 'source\texplicit\n'
      printf 'depth_center_method\t%s\n' "explicit"
      printf 'depth_center_x\t%s\n' "$DEPTH_CENTER"
    } > "$CENTER_TSV"
  fi

  if [[ -z "${DEPTH_LOW:-}" || -z "${DEPTH_HIGH:-}" ]]; then
    eval "$($PYTHON_BIN - <<'PY' "$DEPTH_CENTER" "$DEPTH_LOW_MULT" "$DEPTH_HIGH_MULT"
import math, sys
center = float(sys.argv[1])
low_mult = float(sys.argv[2])
high_mult = float(sys.argv[3])
print(f'DERIVED_LOW={math.floor(center * low_mult)}')
print(f'DERIVED_HIGH={math.floor(center * high_mult)}')
PY
)"
    DEPTH_LOW="${DEPTH_LOW:-$DERIVED_LOW}"
    DEPTH_HIGH="${DEPTH_HIGH:-$DERIVED_HIGH}"
  fi

  {
    printf 'metric\tvalue\n'
    printf 'source\tbam_depth\n'
    printf 'mode\t%s\n' "$MODE"
    printf 'mq_cutoff\t%s\n' "$MQ_CUTOFF"
    printf 'depth_center_method\t%s\n' "$DEPTH_CENTER_METHOD"
    printf 'depth_center_x\t%s\n' "$DEPTH_CENTER"
    printf 'depth_low\t%s\n' "$DEPTH_LOW"
    printf 'depth_high\t%s\n' "$DEPTH_HIGH"
  } > "$THRESH_TSV"

  if should_refresh "$CALLABLE_OUT"; then
    log "derive callable BED from BAM"
    "$SAMTOOLS" depth -aa -Q "$MQ_CUTOFF" -d 0 "$BAM" \
      | "$PYTHON_BIN" "$ROOT_DIR/scripts/callable_regions_from_depth.py" bed --low "$DEPTH_LOW" --high "$DEPTH_HIGH" \
      > "$RAW_BED"
    if [[ "$MODE" == "corrected_query_agp" ]]; then
      log "lift callable BED from corrected-query to scaffold coordinates"
      "$PYTHON_BIN" "$ROOT_DIR/scripts/agp_liftover_bed.py" \
        --agp "$AGP" \
        --bed "$RAW_BED" \
        --output "$LIFTED_RAW_BED"
      sort -k1,1 -k2,2n "$LIFTED_RAW_BED" | "$BEDTOOLS" merge -i - > "$CALLABLE_OUT"
    else
      sort -k1,1 -k2,2n "$RAW_BED" | "$BEDTOOLS" merge -i - > "$CALLABLE_OUT"
    fi
  fi
fi

if should_refresh "$NONCALLABLE_OUT"; then
  log "derive noncallable complement"
  "$BEDTOOLS" complement -i "$CALLABLE_OUT" -g "$FAI" > "$NONCALLABLE_OUT"
fi

if should_refresh "$MASKED_FASTA"; then
  log "mask noncallable sequence"
  "$BEDTOOLS" maskfasta -fi "$FASTA" -bed "$NONCALLABLE_OUT" -fo "$MASKED_FASTA" -mc "$MASK_CHAR"
  "$SAMTOOLS" faidx "$MASKED_FASTA"
fi

if [[ "$OUTPUT_CALLABLE_ONLY" == "1" ]]; then
  if should_refresh "$CALLABLE_ONLY_FASTA"; then
    log "emit callable-only FASTA"
    "$BEDTOOLS" getfasta -fi "$FASTA" -bed "$CALLABLE_OUT" -fo "$CALLABLE_ONLY_FASTA" -name
  fi
fi

log "write summary"
"$PYTHON_BIN" - <<'PY' "$FAI" "$CALLABLE_OUT" "$NONCALLABLE_OUT" "$THRESH_TSV" "$SUMMARY_TSV" "$SAMPLE_TAG" "$FASTA" "$BAM" "$CALLABLE_BED" "$MODE" "$AGP"
import csv, sys
from pathlib import Path
fai_path, callable_bed, noncallable_bed, thresh_tsv, out_tsv, sample_tag, fasta, bam, callable_src, mode, agp = sys.argv[1:]

def bed_bp(path: str) -> int:
    total = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, start, end, *_ = line.rstrip('\n').split('\t')
            total += int(end) - int(start)
    return total

total_bp = 0
with open(fai_path) as fh:
    for line in fh:
        if not line.strip():
            continue
        total_bp += int(line.split('\t')[1])
callable_bp = bed_bp(callable_bed)
noncallable_bp = bed_bp(noncallable_bed)
metrics = {}
with open(thresh_tsv) as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        metrics[row['metric']] = row['value']
rows = [{
    'sample': sample_tag,
    'fasta': fasta,
    'bam': bam,
    'callable_source': callable_src or 'bam_depth',
    'mode': mode,
    'agp': agp,
    'mq_cutoff': metrics.get('mq_cutoff', ''),
    'depth_center_method': metrics.get('depth_center_method', ''),
    'depth_center_x': metrics.get('depth_center_x', ''),
    'depth_low': metrics.get('depth_low', ''),
    'depth_high': metrics.get('depth_high', ''),
    'total_bp': total_bp,
    'callable_bp': callable_bp,
    'noncallable_bp': noncallable_bp,
    'callable_pct': f"{(100.0 * callable_bp / total_bp) if total_bp else 0.0:.4f}",
    'noncallable_pct': f"{(100.0 * noncallable_bp / total_bp) if total_bp else 0.0:.4f}",
}]
with open(out_tsv, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
PY

log "done: $SUMMARY_TSV"
log "outputs: $CALLABLE_OUT $NONCALLABLE_OUT $MASKED_FASTA"
