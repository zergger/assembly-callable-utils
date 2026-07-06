#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
umask 0022

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
DEFAULT_CONFIG="$ROOT_DIR/configs/map_for_callable.env"

usage() {
  cat <<USAGE
Usage: bash scripts/map_for_callable.sh [config.env]

Map paired short reads back to one assembly FASTA and produce the sorted/indexed
BAM expected by run_callable_assembly.sh.

Required inputs:
  FASTA                 target assembly FASTA; BAM coordinates will match this FASTA
  READ1                 read 1 FASTQ/FASTA, optionally gzip-compressed
  READ2                 read 2 FASTQ/FASTA, optionally gzip-compressed

Common variables:
  SAMPLE_TAG            output label [FASTA basename]
  OUTDIR                output root [results/callable_mapping/<SAMPLE_TAG>]
  BAM_OUT               explicit output BAM path [OUTDIR/map/<SAMPLE_TAG>.sorted.bam]
  THREADS               mapping/sort/stats threads [10]
  ALIGNER               auto|bwa-mem2|bwa [auto]
  FORCE                 recompute BAM/stats [0]

Tool variables:
  SAMTOOLS              samtools executable [samtools]
  BWA                   bwa executable [bwa]
  BWA_MEM2              bwa-mem2 executable [bwa-mem2]
  PYTHON_BIN            python executable for depth summaries [python3]

Optional summaries:
  RUN_FLAGSTAT          write stats/flagstat.txt [1]
  RUN_SAMTOOLS_STATS    write stats/samtools.stats.txt [1]
  RUN_DEPTH_SUMMARY     write stats/depth_summary.tsv [0]
  SAVE_DEPTH_TSV        keep full depth TSV if RUN_DEPTH_SUMMARY=1 [0]
  DEPTH_MQ              optional mapping-quality cutoff for depth summary [0]

Outputs:
  OUTDIR/map/<SAMPLE_TAG>.sorted.bam
  OUTDIR/map/<SAMPLE_TAG>.sorted.bam.bai
  OUTDIR/stats/flagstat.txt
  OUTDIR/stats/samtools.stats.txt
  OUTDIR/map_for_callable.outputs.env
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

CONFIG_FILE="${1:-}"
if [[ -n "$CONFIG_FILE" ]]; then
  if [[ ! -f "$CONFIG_FILE" ]]; then
    if [[ -f "$ROOT_DIR/$CONFIG_FILE" ]]; then
      CONFIG_FILE="$ROOT_DIR/$CONFIG_FILE"
    else
      echo "ERROR: config file not found: $CONFIG_FILE" >&2
      exit 1
    fi
  fi
  # shellcheck disable=SC1090
  set -a
  source "$CONFIG_FILE"
  set +a
elif [[ -f "$DEFAULT_CONFIG" ]]; then
  # shellcheck disable=SC1090
  set -a
  source "$DEFAULT_CONFIG"
  set +a
fi

log() {
  printf '[%(%F %T)T] %s\n' -1 "$*" >&2
}

resolve_path() {
  local p="$1"
  if [[ -z "$p" ]]; then
    printf '\n'
  elif [[ "$p" = /* ]]; then
    printf '%s\n' "$p"
  else
    printf '%s/%s\n' "$ROOT_DIR" "$p"
  fi
}

require_file() {
  local p="$1"
  [[ -f "$p" ]] || { echo "ERROR: required file not found: $p" >&2; exit 1; }
}

has_exec() {
  local c="$1"
  if [[ "$c" == */* ]]; then
    [[ -x "$c" ]]
  else
    command -v "$c" >/dev/null 2>&1
  fi
}

require_exec() {
  local c="$1"
  has_exec "$c" || { echo "ERROR: required executable missing from PATH or not executable: $c" >&2; exit 1; }
}

exists_nonempty() {
  [[ -s "$1" ]]
}

samtools_threads_arg() {
  local subcmd="$1"
  local threads="$2"
  if [[ "$threads" -gt 1 ]] && { "$SAMTOOLS" "$subcmd" 2>&1 || true; } | grep -q -- ' -@'; then
    printf -- '-@ %s' "$threads"
  fi
}

choose_aligner() {
  case "$ALIGNER" in
    auto)
      if has_exec "$BWA_MEM2"; then
        ALIGNER_NAME=bwa-mem2
        ALIGNER_BIN=$BWA_MEM2
      elif has_exec "$BWA"; then
        ALIGNER_NAME=bwa
        ALIGNER_BIN=$BWA
      else
        echo "ERROR: neither bwa-mem2 nor bwa is executable" >&2
        exit 1
      fi
      ;;
    bwa-mem2)
      require_exec "$BWA_MEM2"
      ALIGNER_NAME=bwa-mem2
      ALIGNER_BIN=$BWA_MEM2
      ;;
    bwa)
      require_exec "$BWA"
      ALIGNER_NAME=bwa
      ALIGNER_BIN=$BWA
      ;;
    *)
      echo "ERROR: ALIGNER must be one of: auto, bwa-mem2, bwa" >&2
      exit 1
      ;;
  esac
}

ensure_index() {
  local fa="$1"
  local need=0
  if [[ "$ALIGNER_NAME" == "bwa-mem2" ]]; then
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
      if [[ ! -f "${fa}${suffix}" ]]; then
        need=1
        break
      fi
    done
  else
    [[ -f "${fa}.bwt" ]] || need=1
  fi
  if [[ "$need" == "1" ]]; then
    log "Index FASTA with $ALIGNER_NAME: $fa"
    "$ALIGNER_BIN" index "$fa"
  fi
}

write_depth_summary() {
  local bam="$1"
  local out_tsv="$2"
  local depth_tsv="$3"
  local mq_arg=()
  local thread_arg=()
  local depth_cmd

  if [[ "$DEPTH_MQ" -gt 0 ]]; then
    mq_arg=(-Q "$DEPTH_MQ")
  fi
  read -r -a thread_arg <<<"$(samtools_threads_arg depth "$THREADS")"

  mkdir -p "$(dirname "$out_tsv")"
  if [[ "$SAVE_DEPTH_TSV" == "1" ]]; then
    if [[ "$FORCE" == "1" || ! -s "$depth_tsv" ]]; then
      log "Write full depth TSV: $depth_tsv"
      "$SAMTOOLS" depth -aa "${mq_arg[@]}" "${thread_arg[@]}" "$bam" > "$depth_tsv"
    fi
    depth_cmd=(cat "$depth_tsv")
  else
    depth_cmd=("$SAMTOOLS" depth -aa "${mq_arg[@]}" "${thread_arg[@]}" "$bam")
  fi

  log "Write depth summary: $out_tsv"
  "${depth_cmd[@]}" | "$PYTHON_BIN" -c '
import collections
import csv
import sys

sample, out_path = sys.argv[1:3]
hist = collections.Counter()
total_positions = 0
total_depth = 0
covered_positions = 0
for line in sys.stdin:
    if not line.strip():
        continue
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 3:
        continue
    depth = int(float(parts[2]))
    hist[depth] += 1
    total_positions += 1
    total_depth += depth
    if depth > 0:
        covered_positions += 1

def median_from_hist(counter, include_zero=True):
    if include_zero:
        n = sum(counter.values())
        items = sorted(counter.items())
    else:
        n = sum(v for k, v in counter.items() if k > 0)
        items = sorted((k, v) for k, v in counter.items() if k > 0)
    if n == 0:
        return "NA"
    mid1 = (n - 1) // 2
    mid2 = n // 2
    seen = 0
    vals = []
    for depth, count in items:
        nxt = seen + count
        if seen <= mid1 < nxt:
            vals.append(depth)
        if seen <= mid2 < nxt:
            vals.append(depth)
            break
        seen = nxt
    return sum(vals) / len(vals)

def mode_from_hist(counter, include_zero=True):
    items = counter.items() if include_zero else ((k, v) for k, v in counter.items() if k > 0)
    best = None
    for depth, count in items:
        if best is None or count > best[1] or (count == best[1] and depth < best[0]):
            best = (depth, count)
    return "NA" if best is None else best[0]

mean_depth = (total_depth / total_positions) if total_positions else 0
mean_nonzero = (total_depth / covered_positions) if covered_positions else 0
row = {
    "sample": sample,
    "total_positions": total_positions,
    "covered_positions": covered_positions,
    "zero_positions": total_positions - covered_positions,
    "mean_depth": f"{mean_depth:.6f}",
    "mean_depth_nonzero": f"{mean_nonzero:.6f}",
    "median_depth_all": median_from_hist(hist, True),
    "median_depth_nonzero": median_from_hist(hist, False),
    "modal_depth_all": mode_from_hist(hist, True),
    "modal_depth_nonzero": mode_from_hist(hist, False),
}
with open(out_path, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=list(row), delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
' "$SAMPLE_TAG" "$out_tsv"
}

THREADS=${THREADS:-10}
ALIGNER=${ALIGNER:-auto}
FORCE=${FORCE:-0}
RUN_FLAGSTAT=${RUN_FLAGSTAT:-1}
RUN_SAMTOOLS_STATS=${RUN_SAMTOOLS_STATS:-1}
RUN_DEPTH_SUMMARY=${RUN_DEPTH_SUMMARY:-0}
SAVE_DEPTH_TSV=${SAVE_DEPTH_TSV:-0}
DEPTH_MQ=${DEPTH_MQ:-0}
SAMTOOLS=${SAMTOOLS:-samtools}
BWA=${BWA:-bwa}
BWA_MEM2=${BWA_MEM2:-bwa-mem2}
PYTHON_BIN=${PYTHON_BIN:-python3}

FASTA=$(resolve_path "${FASTA:-}")
READ1=$(resolve_path "${READ1:-}")
READ2=$(resolve_path "${READ2:-}")
[[ -n "$FASTA" ]] || { echo "ERROR: FASTA is required" >&2; exit 1; }
[[ -n "$READ1" ]] || { echo "ERROR: READ1 is required" >&2; exit 1; }
[[ -n "$READ2" ]] || { echo "ERROR: READ2 is required" >&2; exit 1; }
require_file "$FASTA"
require_file "$READ1"
require_file "$READ2"
require_exec "$SAMTOOLS"
require_exec "$PYTHON_BIN"

SAMPLE_TAG=${SAMPLE_TAG:-$(basename "$FASTA")}
SAMPLE_TAG=${SAMPLE_TAG%%.*}
OUTDIR=$(resolve_path "${OUTDIR:-results/callable_mapping/$SAMPLE_TAG}")
MAP_DIR=${MAP_DIR:-$OUTDIR/map}
STATS_DIR=${STATS_DIR:-$OUTDIR/stats}
BAM_OUT=$(resolve_path "${BAM_OUT:-$MAP_DIR/${SAMPLE_TAG}.sorted.bam}")
OUTPUT_ENV=${OUTPUT_ENV:-$OUTDIR/map_for_callable.outputs.env}
FLAGSTAT_OUT=${FLAGSTAT_OUT:-$STATS_DIR/flagstat.txt}
SAMSTATS_OUT=${SAMSTATS_OUT:-$STATS_DIR/samtools.stats.txt}
DEPTH_TSV=${DEPTH_TSV:-$STATS_DIR/depth.tsv}
DEPTH_SUMMARY=${DEPTH_SUMMARY:-$STATS_DIR/depth_summary.tsv}

mkdir -p "$(dirname "$BAM_OUT")" "$STATS_DIR"
choose_aligner
ensure_index "$FASTA"
[[ -s "$FASTA.fai" ]] || "$SAMTOOLS" faidx "$FASTA"

log "Mapping reads for callable BAM"
log "Sample tag:  $SAMPLE_TAG"
log "FASTA:       $FASTA"
log "READ1:       $READ1"
log "READ2:       $READ2"
log "Aligner:     $ALIGNER_NAME"
log "Threads:     $THREADS"
log "BAM output:  $BAM_OUT"

sort_thread_arg=()
read -r -a sort_thread_arg <<<"$(samtools_threads_arg sort "$THREADS")"
if [[ "$FORCE" == "1" || ! -s "$BAM_OUT" ]]; then
  tmp_bam="${BAM_OUT}.tmp.$$"
  rm -f "$tmp_bam" "$tmp_bam.bai" "$tmp_bam.csi"
  log "Run alignment and sort BAM"
  "$ALIGNER_BIN" mem -t "$THREADS" "$FASTA" "$READ1" "$READ2" \
    | "$SAMTOOLS" sort "${sort_thread_arg[@]}" -o "$tmp_bam" -
  index_thread_arg=()
  read -r -a index_thread_arg <<<"$(samtools_threads_arg index "$THREADS")"
  "$SAMTOOLS" index "${index_thread_arg[@]}" "$tmp_bam"
  mv -f "$tmp_bam" "$BAM_OUT"
  if [[ -s "$tmp_bam.bai" ]]; then
    mv -f "$tmp_bam.bai" "$BAM_OUT.bai"
  elif [[ -s "$tmp_bam.csi" ]]; then
    mv -f "$tmp_bam.csi" "$BAM_OUT.csi"
  elif [[ -s "${tmp_bam}.bai" ]]; then
    mv -f "${tmp_bam}.bai" "$BAM_OUT.bai"
  elif [[ -s "${tmp_bam}.csi" ]]; then
    mv -f "${tmp_bam}.csi" "$BAM_OUT.csi"
  fi
else
  log "BAM exists; skipping mapping: $BAM_OUT"
fi

if [[ ! -s "$BAM_OUT.bai" && ! -s "$BAM_OUT.csi" ]]; then
  log "BAM exists but index is missing; indexing: $BAM_OUT"
  index_thread_arg=()
  read -r -a index_thread_arg <<<"$(samtools_threads_arg index "$THREADS")"
  "$SAMTOOLS" index "${index_thread_arg[@]}" "$BAM_OUT"
fi

if [[ ! -s "$BAM_OUT" ]]; then
  echo "ERROR: BAM not found after mapping: $BAM_OUT" >&2
  exit 1
fi
if [[ ! -s "$BAM_OUT.bai" && ! -s "$BAM_OUT.csi" ]]; then
  echo "ERROR: BAM index not found after mapping: $BAM_OUT.bai or $BAM_OUT.csi" >&2
  exit 1
fi

if [[ "$RUN_FLAGSTAT" == "1" && ( "$FORCE" == "1" || ! -s "$FLAGSTAT_OUT" ) ]]; then
  log "Write flagstat: $FLAGSTAT_OUT"
  "$SAMTOOLS" flagstat "$BAM_OUT" > "$FLAGSTAT_OUT"
fi

if [[ "$RUN_SAMTOOLS_STATS" == "1" && ( "$FORCE" == "1" || ! -s "$SAMSTATS_OUT" ) ]]; then
  log "Write samtools stats: $SAMSTATS_OUT"
  stats_thread_arg=()
  read -r -a stats_thread_arg <<<"$(samtools_threads_arg stats "$THREADS")"
  "$SAMTOOLS" stats "${stats_thread_arg[@]}" "$BAM_OUT" > "$SAMSTATS_OUT"
fi

if [[ "$RUN_DEPTH_SUMMARY" == "1" && ( "$FORCE" == "1" || ! -s "$DEPTH_SUMMARY" ) ]]; then
  write_depth_summary "$BAM_OUT" "$DEPTH_SUMMARY" "$DEPTH_TSV"
fi

cat > "$OUTPUT_ENV" <<OUT
SAMPLE_TAG="$SAMPLE_TAG"
FASTA="$FASTA"
READ1="$READ1"
READ2="$READ2"
OUTDIR="$OUTDIR"
BAM="$BAM_OUT"
FLAGSTAT="$FLAGSTAT_OUT"
SAMTOOLS_STATS="$SAMSTATS_OUT"
DEPTH_SUMMARY="$DEPTH_SUMMARY"
ALIGNER="$ALIGNER_NAME"
THREADS="$THREADS"
OUT

log "Done. Key outputs:"
log "  BAM:       $BAM_OUT"
log "  Index:     $BAM_OUT.bai"
log "  Flagstat:  $FLAGSTAT_OUT"
log "  Stats:     $SAMSTATS_OUT"
log "  Output env: $OUTPUT_ENV"
