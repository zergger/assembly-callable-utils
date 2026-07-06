#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
umask 0022

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [[ -d "$SCRIPT_DIR/../configs" ]]; then
  ROOT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
else
  ROOT_DIR="$SCRIPT_DIR"
fi
DEFAULT_CONFIG="$ROOT_DIR/configs/ab_window_core_shell.env"
[[ -f "$DEFAULT_CONFIG" ]] || DEFAULT_CONFIG="$SCRIPT_DIR/ab_window_core_shell.env"

usage() {
  cat <<USAGE
Usage: bash scripts/ab_window_core_shell.sh [config.env]

Run A/B assembly-window classification and quantify coordinate core/shell
against branch-specific callable BED files.

Required inputs:
  REF_A_FASTA              branch A scaffold FASTA
  REF_B_FASTA              branch B scaffold FASTA
  REF_A_CALLABLE_BED       callable BED in branch A FASTA coordinates
  REF_B_CALLABLE_BED       callable BED in branch B FASTA coordinates

Common variables:
  REF_A_PROJECTION         callable-table key for branch A [refA]
  REF_B_PROJECTION         callable-table key for branch B [refB]
  REF_A_NAME               display name for branch A [reference_A]
  REF_B_NAME               display name for branch B [reference_B]
  OUTDIR                   output root [results/ab_window_core_shell]
  COMPARE_OUT              table output directory [results/ragtag_compare]
  WINDOW_SIZE              fixed window size [50000]
  WINDOW_MIN_OVERLAP_FRAC  minimum nominal-window overlap for a PAF hit [0.5]
  PRIMARY_REGEX_A          regex selecting primary A windows [.*]
  PRIMARY_REGEX_B          regex selecting primary B windows [.*]
  THREADS                  default threads [10]
  FORCE                    recompute outputs [0]

Tool variables:
  PYTHON_BIN               python executable [python3]
  SAMTOOLS                 samtools executable [samtools]
  MINIMAP2                 minimap2 executable [minimap2]
  MINIMAP2_PRESET          minimap2 preset [asm5]
  MINIMAP2_THREADS         minimap2 threads [THREADS]

Outputs:
  OUTDIR/stats/A_to_B_window_classification_primary.tsv
  OUTDIR/stats/B_to_A_window_classification_primary.tsv
  COMPARE_OUT/tableS_callable.tsv
  COMPARE_OUT/tableS_window_primary.tsv
  COMPARE_OUT/tableS_core_shell.tsv
  OUTDIR/core_shell/summary/core_shell_summary.tsv
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

die() {
  echo "ERROR: $*" >&2
  exit 1
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
  [[ -f "$p" ]] || die "required file not found: $p"
}

require_exec() {
  local c="$1"
  if [[ "$c" == */* ]]; then
    [[ -x "$c" ]] || die "required executable not found or not executable: $c"
  else
    command -v "$c" >/dev/null 2>&1 || die "required executable missing from PATH: $c"
  fi
}

find_helper() {
  local name="$1"
  if [[ -f "$SCRIPT_DIR/$name" ]]; then
    printf '%s/%s\n' "$SCRIPT_DIR" "$name"
  elif [[ -f "$ROOT_DIR/scripts/$name" ]]; then
    printf '%s/scripts/%s\n' "$ROOT_DIR" "$name"
  else
    die "helper script not found: $name"
  fi
}

THREADS="${THREADS:-10}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SAMTOOLS="${SAMTOOLS:-samtools}"
MINIMAP2="${MINIMAP2:-minimap2}"
MINIMAP2_PRESET="${MINIMAP2_PRESET:-asm5}"
MINIMAP2_THREADS="${MINIMAP2_THREADS:-$THREADS}"
WINDOW_SIZE="${WINDOW_SIZE:-50000}"
WINDOW_MIN_OVERLAP_FRAC="${WINDOW_MIN_OVERLAP_FRAC:-0.5}"
PRIMARY_REGEX_A="${PRIMARY_REGEX_A:-.*}"
PRIMARY_REGEX_B="${PRIMARY_REGEX_B:-.*}"
FORCE="${FORCE:-0}"

REF_A_PROJECTION="${REF_A_PROJECTION:-refA}"
REF_B_PROJECTION="${REF_B_PROJECTION:-refB}"
REF_A_NAME="${REF_A_NAME:-reference_A}"
REF_B_NAME="${REF_B_NAME:-reference_B}"
REF_A_WATCH_CHROMS="${REF_A_WATCH_CHROMS:-}"
REF_B_WATCH_CHROMS="${REF_B_WATCH_CHROMS:-}"

OUTDIR="$(resolve_path "${OUTDIR:-results/ab_window_core_shell}")"
COMPARE_OUT="$(resolve_path "${COMPARE_OUT:-results/ragtag_compare}")"
CALLABLE_TABLE="$(resolve_path "${CALLABLE_TABLE:-$COMPARE_OUT/tableS_callable.tsv}")"
CORE_SHELL_OUTDIR="$(resolve_path "${CORE_SHELL_OUTDIR:-$OUTDIR/core_shell}")"

REF_A_FASTA="$(resolve_path "${REF_A_FASTA:?REF_A_FASTA must be set}")"
REF_B_FASTA="$(resolve_path "${REF_B_FASTA:?REF_B_FASTA must be set}")"
REF_A_CALLABLE_BED="$(resolve_path "${REF_A_CALLABLE_BED:?REF_A_CALLABLE_BED must be set}")"
REF_B_CALLABLE_BED="$(resolve_path "${REF_B_CALLABLE_BED:?REF_B_CALLABLE_BED must be set}")"

WINDOW_HELPER="$(find_helper window_liftover_stats.py)"
QUANTIFY_HELPER="$(find_helper quantify_core_shell.py)"

require_file "$REF_A_FASTA"
require_file "$REF_B_FASTA"
require_file "$REF_A_CALLABLE_BED"
require_file "$REF_B_CALLABLE_BED"
require_file "$WINDOW_HELPER"
require_file "$QUANTIFY_HELPER"
require_exec "$PYTHON_BIN"
require_exec "$SAMTOOLS"
require_exec "$MINIMAP2"

STATS_DIR="$OUTDIR/stats"
ASM2ASM_DIR="$OUTDIR/asm2asm"
A_WINDOWS_DIR="$OUTDIR/A/windows"
B_WINDOWS_DIR="$OUTDIR/B/windows"
CORE_SUMMARY_DIR="$CORE_SHELL_OUTDIR/summary"
CORE_BED_ROOT="$CORE_SHELL_OUTDIR/bed"
mkdir -p "$STATS_DIR" "$ASM2ASM_DIR" "$A_WINDOWS_DIR" "$B_WINDOWS_DIR" "$COMPARE_OUT" "$CORE_SUMMARY_DIR" "$CORE_BED_ROOT"

A_FAI="$REF_A_FASTA.fai"
B_FAI="$REF_B_FASTA.fai"

ensure_fai() {
  local fasta="$1"
  local fai="$2"
  if [[ "$FORCE" == "1" || ! -s "$fai" ]]; then
    log "samtools faidx: $fasta"
    "$SAMTOOLS" faidx "$fasta"
  fi
}

make_windows() {
  local fai="$1"
  local out="$2"
  local tmp="$out.tmp.$$"
  if [[ "$FORCE" == "1" || ! -s "$out" ]]; then
    log "make windows: $out"
    awk -v w="$WINDOW_SIZE" 'BEGIN{OFS="\t"} {for (s=0; s<$2; s+=w) {e=s+w; if (e>$2) e=$2; print $1,s,e}}' "$fai" > "$tmp"
    mv "$tmp" "$out"
  fi
}

filter_primary_windows() {
  local all_windows="$1"
  local regex="$2"
  local out="$3"
  local tmp="$out.tmp.$$"
  if [[ "$FORCE" == "1" || ! -s "$out" ]]; then
    log "filter primary windows: $out"
    awk -v pat="$regex" '$1 ~ pat' "$all_windows" > "$tmp"
    mv "$tmp" "$out"
  fi
  [[ -s "$out" ]] || die "primary window BED is empty: $out; check PRIMARY_REGEX"
}

run_minimap_pair() {
  local label="$1"
  local target="$2"
  local query="$3"
  local out="$4"
  local done="$out.done"
  if [[ "$FORCE" == "1" || ! -f "$done" ]]; then
    log "minimap2 $label"
    rm -f "$out" "$done"
    /usr/bin/time -f "${label}\t%E\t%M" -o "$out.time" \
      "$MINIMAP2" -t "$MINIMAP2_THREADS" -x "$MINIMAP2_PRESET" "$target" "$query" > "$out"
    touch "$done"
  fi
}

run_window_stats() {
  local label="$1"
  local windows="$2"
  local paf="$3"
  local summary="$4"
  local per_window="$5"
  if [[ "$FORCE" == "1" || ! -s "$summary" || ! -s "$per_window" ]]; then
    log "window classification: $label"
    "$PYTHON_BIN" "$WINDOW_HELPER" \
      --windows "$windows" \
      --paf "$paf" \
      --pair-label "$label" \
      --window-size "$WINDOW_SIZE" \
      --min-overlap-frac "$WINDOW_MIN_OVERLAP_FRAC" \
      --summary-out "$summary" \
      --per-window-out "$per_window"
  fi
}

sum_fai_bp() {
  awk '{s += $2} END{printf "%.0f", s + 0}' "$1"
}

sum_bed_bp() {
  awk 'BEGIN{s=0} !/^#/ && NF>=3 {s += $3 - $2} END{printf "%.0f", s + 0}' "$1"
}

pct() {
  awk -v n="$1" -v d="$2" 'BEGIN{if (d>0) printf "%.4f", 100*n/d; else printf "0.0000"}'
}

write_callable_table() {
  local out="$1"
  local a_total a_callable b_total b_callable
  a_total=$(sum_fai_bp "$A_FAI")
  b_total=$(sum_fai_bp "$B_FAI")
  a_callable=$(sum_bed_bp "$REF_A_CALLABLE_BED")
  b_callable=$(sum_bed_bp "$REF_B_CALLABLE_BED")

  log "write callable table: $out"
  {
    printf 'projection\tmq_cutoff\tdepth_low\tdepth_high\ttotal_bp\tcallable_bp\tcallable_pct\tnote\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$REF_A_PROJECTION" "${REF_A_MQ_CUTOFF:-${MQ_CUTOFF:-}}" "${REF_A_DEPTH_LOW:-${DEPTH_LOW:-}}" "${REF_A_DEPTH_HIGH:-${DEPTH_HIGH:-}}" \
      "$a_total" "$a_callable" "$(pct "$a_callable" "$a_total")" "Generated from REF_A_CALLABLE_BED"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$REF_B_PROJECTION" "${REF_B_MQ_CUTOFF:-${MQ_CUTOFF:-}}" "${REF_B_DEPTH_LOW:-${DEPTH_LOW:-}}" "${REF_B_DEPTH_HIGH:-${DEPTH_HIGH:-}}" \
      "$b_total" "$b_callable" "$(pct "$b_callable" "$b_total")" "Generated from REF_B_CALLABLE_BED"
  } > "$out"
}

combine_two_tables() {
  local first="$1"
  local second="$2"
  local out="$3"
  if [[ "$FORCE" == "1" || ! -s "$out" ]]; then
    {
      head -n 1 "$first"
      tail -n +2 "$first"
      tail -n +2 "$second"
    } > "$out"
  fi
}

ensure_fai "$REF_A_FASTA" "$A_FAI"
ensure_fai "$REF_B_FASTA" "$B_FAI"

A_WINDOWS="$A_WINDOWS_DIR/A.${WINDOW_SIZE}.bed"
B_WINDOWS="$B_WINDOWS_DIR/B.${WINDOW_SIZE}.bed"
A_PRIMARY_WINDOWS="$A_WINDOWS_DIR/A.${WINDOW_SIZE}.primary.bed"
B_PRIMARY_WINDOWS="$B_WINDOWS_DIR/B.${WINDOW_SIZE}.primary.bed"
make_windows "$A_FAI" "$A_WINDOWS"
make_windows "$B_FAI" "$B_WINDOWS"
filter_primary_windows "$A_WINDOWS" "$PRIMARY_REGEX_A" "$A_PRIMARY_WINDOWS"
filter_primary_windows "$B_WINDOWS" "$PRIMARY_REGEX_B" "$B_PRIMARY_WINDOWS"

A_TO_B_PAF="$ASM2ASM_DIR/A_to_B.paf"
B_TO_A_PAF="$ASM2ASM_DIR/B_to_A.paf"
run_minimap_pair "A_to_B" "$REF_B_FASTA" "$REF_A_FASTA" "$A_TO_B_PAF"
run_minimap_pair "B_to_A" "$REF_A_FASTA" "$REF_B_FASTA" "$B_TO_A_PAF"

run_window_stats "A_to_B" "$A_WINDOWS" "$A_TO_B_PAF" \
  "$STATS_DIR/A_to_B_window_liftover_stats.tsv" \
  "$STATS_DIR/A_to_B_window_classification.tsv"
run_window_stats "B_to_A" "$B_WINDOWS" "$B_TO_A_PAF" \
  "$STATS_DIR/B_to_A_window_liftover_stats.tsv" \
  "$STATS_DIR/B_to_A_window_classification.tsv"
combine_two_tables "$STATS_DIR/A_to_B_window_liftover_stats.tsv" "$STATS_DIR/B_to_A_window_liftover_stats.tsv" \
  "$STATS_DIR/window_liftover_stats.tsv"
cp "$STATS_DIR/window_liftover_stats.tsv" "$COMPARE_OUT/tableS_window_all.tsv"

run_window_stats "A_to_B_primary" "$A_PRIMARY_WINDOWS" "$A_TO_B_PAF" \
  "$STATS_DIR/A_to_B_window_liftover_primary.tsv" \
  "$STATS_DIR/A_to_B_window_classification_primary.tsv"
run_window_stats "B_to_A_primary" "$B_PRIMARY_WINDOWS" "$B_TO_A_PAF" \
  "$STATS_DIR/B_to_A_window_liftover_primary.tsv" \
  "$STATS_DIR/B_to_A_window_classification_primary.tsv"
combine_two_tables "$STATS_DIR/A_to_B_window_liftover_primary.tsv" "$STATS_DIR/B_to_A_window_liftover_primary.tsv" \
  "$STATS_DIR/window_liftover_stats_primary.tsv"
cp "$STATS_DIR/window_liftover_stats_primary.tsv" "$COMPARE_OUT/tableS_window_primary.tsv"

write_callable_table "$CALLABLE_TABLE"

CORE_SUMMARY="$CORE_SUMMARY_DIR/core_shell_summary.tsv"
if [[ "$FORCE" == "1" || ! -s "$CORE_SUMMARY" ]]; then
  log "core/shell quantification"
  "$PYTHON_BIN" "$QUANTIFY_HELPER" \
    --callable-table "$CALLABLE_TABLE" \
    --bed-root "$CORE_BED_ROOT" \
    --summary-output "$CORE_SUMMARY" \
    --reference "$REF_A_PROJECTION::$REF_A_NAME::$STATS_DIR/A_to_B_window_classification_primary.tsv::$REF_A_CALLABLE_BED::$REF_A_WATCH_CHROMS" \
    --reference "$REF_B_PROJECTION::$REF_B_NAME::$STATS_DIR/B_to_A_window_classification_primary.tsv::$REF_B_CALLABLE_BED::$REF_B_WATCH_CHROMS"
fi

cp "$CORE_SUMMARY" "$COMPARE_OUT/tableS_core_shell.tsv"
cp "$CORE_SUMMARY" "$STATS_DIR/core_shell_summary.tsv"

cat > "$OUTDIR/ab_window_core_shell.outputs.env" <<EOF
OUTDIR="$OUTDIR"
COMPARE_OUT="$COMPARE_OUT"
CALLABLE_TABLE="$CALLABLE_TABLE"
TABLES_WINDOW_ALL="$COMPARE_OUT/tableS_window_all.tsv"
TABLES_WINDOW_PRIMARY="$COMPARE_OUT/tableS_window_primary.tsv"
TABLES_CORE_SHELL="$COMPARE_OUT/tableS_core_shell.tsv"
A_TO_B_WINDOW_CLASS="$STATS_DIR/A_to_B_window_classification_primary.tsv"
B_TO_A_WINDOW_CLASS="$STATS_DIR/B_to_A_window_classification_primary.tsv"
CORE_SHELL_SUMMARY="$CORE_SUMMARY"
EOF

log "Done: $COMPARE_OUT/tableS_window_primary.tsv"
log "Done: $COMPARE_OUT/tableS_callable.tsv"
log "Done: $COMPARE_OUT/tableS_core_shell.tsv"
