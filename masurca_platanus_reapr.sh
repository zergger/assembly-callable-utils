#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
umask 0022

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [[ -f "$SCRIPT_DIR/README.md" && -f "$SCRIPT_DIR/masurca_platanus_reapr.sh" ]]; then
  ROOT_DIR=$SCRIPT_DIR
else
  ROOT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
fi
DEFAULT_CONFIG="$ROOT_DIR/masurca_platanus_reapr.env"
export PATH="$SCRIPT_DIR:$PATH"

kill_descendants() {
  local parent=$1
  local sig=${2:-TERM}
  local child
  while read -r child; do
    [[ -n "$child" ]] || continue
    kill_descendants "$child" "$sig"
    kill "-$sig" "$child" 2>/dev/null || true
  done < <(ps -o pid= --ppid "$parent" 2>/dev/null)
}

release_global_lock() {
  if [[ -n "${GLOBAL_LOCK_FD:-}" ]]; then
    flock -u "$GLOBAL_LOCK_FD" 2>/dev/null || true
    exec {GLOBAL_LOCK_FD}>&- 2>/dev/null || true
    unset GLOBAL_LOCK_FD
    log "Released global lock: $GLOBAL_LOCK_FILE"
  fi
}

cleanup() {
  local sig=${1:-INT}
  set +e +u
  release_global_lock
  echo "Interrupted (${sig}); cleaning up child processes..." >&2
  kill_descendants "$$" TERM
  sleep 1
  kill_descendants "$$" KILL
  wait 2>/dev/null || true
  exit 130
}

trap 'cleanup INT' INT
trap 'cleanup TERM' TERM

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'USAGE'
Usage: masurca_platanus_reapr.sh [config.env]

One-sample short-read assembly refinement pipeline:
  1) MaSuRCA draft assembly
  2) Platanus donor draft assembly
  3) REAPR breakpoint detection and breaking
  4) RagTag patch using the Platanus donor assembly
  5) NextPolish Illumina polishing by default

This script stops at a polished, unanchored assembly. Run ragtag_project.sh
separately for reference-guided RagTag correct/scaffold.

Common variables:
  SAMPLE_ID=sample1
  WORKDIR=./work_masurca_platanus_reapr
  R1=/path/to/reads_R1.fastq.gz
  R2=/path/to/reads_R2.fastq.gz
  THREADS=10
  SEED=42
  RUN_DRAFT_ASSEMBLY=1
  RUN_REAPR_POLISH=1

Draft assembly variables:
  MASURCA_CFG=./masurca.cfg
  MASURCA_BIN=masurca
  MASURCA_THREADS=$THREADS
  PLATANUS_BIN=platanus
  PLATANUS_PREFIX=platanus_alt
  PLATANUS_THREADS=$THREADS
  PLATANUS_MEM_GB=64
  DOWNSAMPLE_FRAC=1.0 (gzip reads are staged as plain FASTQ for Platanus only)
  PLATANUS_DECOMPRESSOR=auto (prefer pigz, fall back to gzip)
  PLATANUS_DECOMPRESS_THREADS=$PLATANUS_THREADS

REAPR/RagTag/NextPolish variables:
  PRIMARY_ASM=$WORKDIR/masurca/CA/primary.genome.scf.fasta
  DONOR_ASM=$WORKDIR/platanus/platanus_alt_contig.fa
  REAPR_ALIGNER=bwa-mem2 (bwa-mem2|smalt|minibwa)
  POLISHER=nextpolish (nextpolish|polca|none)
  NEXTPOLISH_OUT=$WORKDIR/nextpolish/01_rundir/genome.nextpolish.fasta

Resume behavior:
  Each stage is skipped if its primary output exists and is non-empty.
  Logs are written under WORKDIR/logs/.
USAGE
  exit 0
fi

if [[ $# -gt 1 ]]; then
  echo "ERROR: expected at most one config file argument" >&2
  exit 1
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
  set -a
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
  set +a
elif [[ -f "$DEFAULT_CONFIG" ]]; then
  set -a
  # shellcheck disable=SC1090
  source "$DEFAULT_CONFIG"
  set +a
fi

log() { printf '[%(%F %T)T] %s\n' -1 "$*" >&2; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

exists_nonempty() {
  local p=$1
  [[ -s "$p" ]]
}

require_file() {
  local p=$1
  if [[ -z "$p" ]]; then
    die "required file path is empty. Set it in the environment or config file."
  fi
  [[ -f "$p" ]] || die "required file not found: $p"
}

require_cmd() {
  local c=$1
  if [[ "$c" == */* ]]; then
    [[ -x "$c" ]] || die "required executable not found or not executable: $c"
  else
    command -v "$c" >/dev/null 2>&1 || die "required command not found on PATH: $c"
  fi
}

run_with_log() {
  local logfile=$1
  shift
  mkdir -p "$(dirname "$logfile")"
  local cmdline="$*"
  local had_errexit=0
  [[ $- == *e* ]] && had_errexit=1
  local t0=$SECONDS

  log "Running: $cmdline"
  set +e
  {
    printf '[%(%F %T)T] CMD: %s\n' -1 "$cmdline"
    "$@"
  } >>"$logfile" 2>&1
  local rc=$?
  (( had_errexit )) && set -e

  local dt=$((SECONDS - t0))
  if [[ $rc -eq 0 ]]; then
    log "Finished (${dt}s): $cmdline"
  else
    log "ERROR (exit=$rc, ${dt}s): $cmdline ; see $logfile"
  fi
  return $rc
}

acquire_global_lock() {
  if ! command -v flock >/dev/null 2>&1; then
    die "flock is required for pipeline locking but was not found"
  fi
  mkdir -p "$(dirname "$GLOBAL_LOCK_FILE")"
  exec {GLOBAL_LOCK_FD}> "$GLOBAL_LOCK_FILE"
  if [[ "$GLOBAL_LOCK_WAIT_SEC" -gt 0 ]]; then
    flock -w "$GLOBAL_LOCK_WAIT_SEC" "$GLOBAL_LOCK_FD" || {
      die "lock busy ($GLOBAL_LOCK_FILE), waited ${GLOBAL_LOCK_WAIT_SEC}s"
    }
  else
    flock -n "$GLOBAL_LOCK_FD" || {
      die "lock busy ($GLOBAL_LOCK_FILE). Another run may be active."
    }
  fi
  log "Acquired global lock: $GLOBAL_LOCK_FILE"
}

bam_ready() {
  local bam=$1
  [[ -s "$bam" ]] || return 1
  [[ -s "${bam}.bai" || -s "${bam}.csi" ]]
}

cleanup_incomplete_bam() {
  local bam=$1
  if [[ -s "$bam" ]] && ! bam_ready "$bam"; then
    log "Found incomplete BAM; removing: $bam"
    rm -f "$bam" "${bam}.bai" "${bam}.csi"
  fi
}

samtools_supports_threads() {
  local subcmd=$1
  local out
  out=$(samtools "$subcmd" --help 2>&1 || true)
  grep -q -- " -@" <<<"$out"
}

samtools_sort_supports_output() {
  local out
  out=$(samtools sort --help 2>&1 || true)
  grep -q -- " -o" <<<"$out"
}

samtools_supports_markdup() {
  local out
  out=$(samtools 2>&1 || true)
  grep -qE '(^|[[:space:]])markdup([[:space:]]|$)' <<<"$out"
}

samtools_fixmate_supports_m() {
  local out
  out=$(samtools fixmate 2>&1 || true)
  grep -q -- " -m" <<<"$out"
}

check_nextpolish_prereqs() {
  if ! samtools_supports_markdup; then
    echo "ERROR: NextPolish requires samtools 'markdup', but current samtools does not provide it." >&2
    samtools 2>&1 | head -n 3 >&2 || true
    exit 1
  fi
  if ! samtools_fixmate_supports_m; then
    echo "ERROR: NextPolish requires 'samtools fixmate -m', but current samtools does not support '-m'." >&2
    samtools 2>&1 | head -n 3 >&2 || true
    exit 1
  fi
}

samtools_sort_index() {
  local in_bam=$1
  local out_bam=$2
  local style=$3
  local sort_threads=()
  local index_threads=()
  if bam_ready "$out_bam"; then
    return 0
  fi
  cleanup_incomplete_bam "$out_bam"
  local tmp_out="${out_bam}.tmp.$$.bam"
  local sort_tmp_prefix="${tmp_out}.sort"
  if [[ "${SAMTOOLS_THREADS:-1}" -gt 1 ]] && samtools_supports_threads sort; then
    sort_threads=(-@ "$SAMTOOLS_THREADS")
  fi
  if [[ "${SAMTOOLS_THREADS:-1}" -gt 1 ]] && samtools_supports_threads index; then
    index_threads=(-@ "$SAMTOOLS_THREADS")
  fi
  if [[ "$style" == "auto" ]]; then
    if samtools_sort_supports_output; then
      style="new"
    else
      style="old"
    fi
  fi
  local sort_rc=0
  if [[ "$style" == "new" ]]; then
    samtools sort "${sort_threads[@]}" -T "$sort_tmp_prefix" -o "$tmp_out" "$in_bam" || sort_rc=$?
  else
    local prefix=${tmp_out%.bam}
    samtools sort "${sort_threads[@]}" "$in_bam" "$prefix" || sort_rc=$?
    tmp_out="${prefix}.bam"
  fi
  if [[ "$sort_rc" -ne 0 ]]; then
    log "samtools sort failed (exit=$sort_rc): $in_bam"
    rm -f "$tmp_out" "${tmp_out}.bai" "${tmp_out}.csi" "${sort_tmp_prefix}".*.bam
    return "$sort_rc"
  fi

  local index_rc=0
  samtools index "${index_threads[@]}" "$tmp_out" || index_rc=$?
  if [[ "$index_rc" -ne 0 ]]; then
    log "samtools index failed (exit=$index_rc): $tmp_out"
    rm -f "$tmp_out" "${tmp_out}.bai" "${tmp_out}.csi" "${sort_tmp_prefix}".*.bam
    return "$index_rc"
  fi
  if [[ -s "${tmp_out}.bai" ]]; then
    mv -f "$tmp_out" "$out_bam"
    mv -f "${tmp_out}.bai" "${out_bam}.bai"
  elif [[ -s "${tmp_out}.csi" ]]; then
    mv -f "$tmp_out" "$out_bam"
    mv -f "${tmp_out}.csi" "${out_bam}.csi"
  else
    log "samtools index did not produce .bai/.csi for $tmp_out"
    rm -f "$tmp_out" "${tmp_out}.bai" "${tmp_out}.csi" "${sort_tmp_prefix}".*.bam
    return 1
  fi
}

bwa_mem2_index_ready() {
  local prefix=$1
  local suffix
  for suffix in .0123 .amb .ann .pac; do
    [[ -s "${prefix}${suffix}" ]] || return 1
  done
  [[ -s "${prefix}.bwt.2bit" || -s "${prefix}.bwt.2bit.64" ]]
}

require_reapr_aligner() {
  case "$REAPR_ALIGNER" in
    smalt) require_cmd smalt ;;
    bwa-mem2) require_cmd bwa-mem2 ;;
    minibwa) require_cmd "$MINIBWA_BIN" ;;
    *) die "REAPR_ALIGNER must be smalt, bwa-mem2, or minibwa; got: $REAPR_ALIGNER" ;;
  esac
}

link_or_downsample_reads() {
  local frac=$1
  local tag=$2
  local out_r1="$READS_DIR/${SAMPLE_ID}_${tag}_R1.ds.fq"
  local out_r2="$READS_DIR/${SAMPLE_ID}_${tag}_R2.ds.fq"
  if awk "BEGIN{exit !($frac < 1.0)}"; then
    require_cmd "$SEQTK_BIN"
    if exists_nonempty "$out_r1" && exists_nonempty "$out_r2"; then
      log "Downsampled reads already exist for $tag; skipping"
    else
      run_with_log "$LOGDIR/${tag}_downsample.log" bash -lc \
        "$SEQTK_BIN sample -s${SEED} '$R1' '$frac' > '$out_r1'"
      run_with_log "$LOGDIR/${tag}_downsample.log" bash -lc \
        "$SEQTK_BIN sample -s${SEED} '$R2' '$frac' > '$out_r2'"
    fi
    printf '%s\t%s\n' "$out_r1" "$out_r2"
  else
    local suffix=""
    case "$R1" in
      *.fq.gz) suffix=".fq.gz" ;;
      *.fastq.gz) suffix=".fastq.gz" ;;
      *.fq) suffix=".fq" ;;
      *.fastq) suffix=".fastq" ;;
      *) suffix="" ;;
    esac
    local link_r1="$READS_DIR/${SAMPLE_ID}_${tag}_R1${suffix}"
    local link_r2="$READS_DIR/${SAMPLE_ID}_${tag}_R2${suffix}"
    log "DOWNSAMPLE_FRAC=${frac}; linking original reads for $tag"
    ln -sf "$R1" "$link_r1"
    ln -sf "$R2" "$link_r2"
    printf '%s\t%s\n' "$link_r1" "$link_r2"
  fi
}

resolve_platanus_decompressor() {
  local requested=$PLATANUS_DECOMPRESSOR
  local selected

  if [[ "$requested" == "auto" ]]; then
    if command -v pigz >/dev/null 2>&1; then
      selected=pigz
    else
      selected=gzip
    fi
    require_cmd "$selected"
    printf '%s\n' "$selected"
    return 0
  fi

  case "${requested##*/}" in
    pigz|gzip)
      require_cmd "$requested"
      printf '%s\n' "$requested"
      ;;
    *)
      log "PLATANUS_DECOMPRESSOR must be auto, pigz, gzip, or a path ending in pigz/gzip; got: $requested"
      return 1
      ;;
  esac
}

decompress_gzip_atomic() {
  local input=$1
  local output=$2
  local decompressor=$3
  local threads=$4
  local tmp="${output}.tmp.$$"

  rm -f "${output}.tmp."*
  case "${decompressor##*/}" in
    pigz)
      if ! "$decompressor" -dc -p "$threads" -- "$input" > "$tmp"; then
        rm -f "$tmp"
        return 1
      fi
      ;;
    gzip)
      if ! "$decompressor" -cd -- "$input" > "$tmp"; then
        rm -f "$tmp"
        return 1
      fi
      ;;
    *)
      rm -f "$tmp"
      return 1
      ;;
  esac
  if [[ ! -s "$tmp" ]]; then
    rm -f "$tmp"
    return 1
  fi
  mv -f "$tmp" "$output"
}

stage_platanus_read() {
  local input=$1
  local mate=$2
  local decompressor=$3

  case "$input" in
    *.gz)
      local output="$READS_DIR/${SAMPLE_ID}_platanus_${mate}.fq"
      if exists_nonempty "$output"; then
        log "Uncompressed Platanus read exists; skipping: $output"
      else
        log "Platanus cannot read gzip input; staging uncompressed $mate"
        if ! run_with_log "$LOGDIR/platanus_reads.log" \
          decompress_gzip_atomic "$input" "$output" "$decompressor" "$PLATANUS_DECOMPRESS_THREADS"; then
          log "Failed to stage uncompressed Platanus read: $input"
          return 1
        fi
      fi
      printf '%s\n' "$output"
      ;;
    *)
      printf '%s\n' "$input"
      ;;
  esac
}

prepare_platanus_reads() {
  local input_r1=$1
  local input_r2=$2
  local output_r1
  local output_r2
  local decompressor=""

  if [[ "$input_r1" == *.gz || "$input_r2" == *.gz ]]; then
    if ! [[ "$PLATANUS_DECOMPRESS_THREADS" =~ ^[1-9][0-9]*$ ]]; then
      log "PLATANUS_DECOMPRESS_THREADS must be a positive integer; got: $PLATANUS_DECOMPRESS_THREADS"
      return 1
    fi
    decompressor=$(resolve_platanus_decompressor) || return 1
    log "Platanus decompressor: $decompressor (threads=$PLATANUS_DECOMPRESS_THREADS)"
  fi

  output_r1=$(stage_platanus_read "$input_r1" R1 "$decompressor") || return 1
  output_r2=$(stage_platanus_read "$input_r2" R2 "$decompressor") || return 1
  printf '%s\t%s\n' "$output_r1" "$output_r2"
}

cleanup_platanus_full_read_staging() {
  local frac=$1
  local staged

  if ! awk "BEGIN{exit !($frac == 1.0)}"; then
    return 0
  fi
  for staged in \
    "$READS_DIR/${SAMPLE_ID}_platanus_R1.fq" \
    "$READS_DIR/${SAMPLE_ID}_platanus_R2.fq"; do
    if [[ -e "$staged" || -L "$staged" ]]; then
      rm -f "$staged"
      log "Removed Platanus full-read staging FASTQ: $staged"
    fi
  done
}

find_masurca_asm() {
  local cand
  for cand in \
    "$MASURCA_RUN_DIR/genome.scf.fasta" \
    "$MASURCA_RUN_DIR/CA/primary.genome.scf.fasta" \
    "$MASURCA_RUN_DIR/CA/primary.genome.ctg.fasta" \
    "$MASURCA_RUN_DIR/CA/primary.genome.contigs.fasta"; do
    if exists_nonempty "$cand"; then
      printf '%s\n' "$cand"
      return 0
    fi
  done
  return 1
}

THREADS=${THREADS:-10}
SAMPLE_ID=${SAMPLE_ID:-sample1}
WORKDIR=${WORKDIR:-$ROOT_DIR/work_masurca_platanus_reapr}
LOGDIR=${LOGDIR:-$WORKDIR/logs}
R1=${R1:-}
R2=${R2:-}
SEED=${SEED:-42}
RUN_DRAFT_ASSEMBLY=${RUN_DRAFT_ASSEMBLY:-1}
RUN_REAPR_POLISH=${RUN_REAPR_POLISH:-1}

MASURCA_THREADS=${MASURCA_THREADS:-$THREADS}
MASURCA_CFG=${MASURCA_CFG:-$ROOT_DIR/masurca.cfg}
MASURCA_BIN=${MASURCA_BIN:-masurca}
PLATANUS_BIN=${PLATANUS_BIN:-platanus}
PLATANUS_PREFIX=${PLATANUS_PREFIX:-platanus_alt}
PLATANUS_THREADS=${PLATANUS_THREADS:-$THREADS}
PLATANUS_MEM_GB=${PLATANUS_MEM_GB:-64}
DOWNSAMPLE_FRAC=${DOWNSAMPLE_FRAC:-1.0}
SEQTK_BIN=${SEQTK_BIN:-seqtk}
PLATANUS_DECOMPRESSOR=${PLATANUS_DECOMPRESSOR:-auto}
PLATANUS_DECOMPRESS_THREADS=${PLATANUS_DECOMPRESS_THREADS:-$PLATANUS_THREADS}

REAPR_THREADS=${REAPR_THREADS:-$THREADS}
POLCA_THREADS=${POLCA_THREADS:-$THREADS}
RAGTAG_THREADS=${RAGTAG_THREADS:-$THREADS}
SAMTOOLS_THREADS=${SAMTOOLS_THREADS:-1}
RUN_REAPR_PERFECTMAP=${RUN_REAPR_PERFECTMAP:-0}
REAPR_INSERT_SIZE=${REAPR_INSERT_SIZE:-220}
REAPR_SMALT_K=${REAPR_SMALT_K:-13}
REAPR_SMALT_S=${REAPR_SMALT_S:-6}
REAPR_PERFECT_PREFIX=${REAPR_PERFECT_PREFIX:-perfect}
REAPR_BREAK_ENABLE=${REAPR_BREAK_ENABLE:-1}
REAPR_BREAK_B=${REAPR_BREAK_B:-0}
REAPR_BREAK_E=${REAPR_BREAK_E:-0.5}
REAPR_SCORE_THREADS=${REAPR_SCORE_THREADS:-$REAPR_THREADS}
REAPR_FILTER_ENABLE=${REAPR_FILTER_ENABLE:-0}
REAPR_FILTER_F=${REAPR_FILTER_F:-2}
REAPR_FILTER_FF=${REAPR_FILTER_FF:-904}
REAPR_FILTER_MAPQ=${REAPR_FILTER_MAPQ:-0}
SAMTOOLS_SORT_STYLE=${SAMTOOLS_SORT_STYLE:-auto}
REAPR_ALIGNER=${REAPR_ALIGNER:-bwa-mem2}
REAPR_BAM=${REAPR_BAM:-}
BWA_MEM2_ARGS=${BWA_MEM2_ARGS:--SP}
MINIBWA_BIN=${MINIBWA_BIN:-minibwa}
MINIBWA_ARGS=${MINIBWA_ARGS:-}

RAGTAG_PATCH_NUCMER_PARAMS=${RAGTAG_PATCH_NUCMER_PARAMS:---maxmatch -l 100 -c 500}
GLOBAL_LOCK_WAIT_SEC=${GLOBAL_LOCK_WAIT_SEC:-${PATCH_LOCK_WAIT_SEC:-0}}
GLOBAL_LOCK_FILE=${GLOBAL_LOCK_FILE:-${PATCH_LOCK_FILE:-$WORKDIR/.pipeline.lock}}

REAPR_BIN=${REAPR_BIN:-reapr}
RAGTAG_BIN=${RAGTAG_BIN:-ragtag.py}
POLCA_BIN=${POLCA_BIN:-polca.sh}
PY3=${PY3:-python3}

POLCA_MEM_PER_THREAD=${POLCA_MEM_PER_THREAD:-1G}
POLISHER=${POLISHER:-nextpolish}
NEXTPOLISH_BIN=${NEXTPOLISH_BIN:-nextPolish}
NEXTPOLISH_CFG=${NEXTPOLISH_CFG:-}
NEXTPOLISH_CMD=${NEXTPOLISH_CMD:-}
NEXTPOLISH_PARALLEL_JOBS=${NEXTPOLISH_PARALLEL_JOBS:-1}
NEXTPOLISH_THREADS=${NEXTPOLISH_THREADS:-$THREADS}
NEXTPOLISH_TASK=${NEXTPOLISH_TASK:-best}
NEXTPOLISH_SGS_OPTIONS=${NEXTPOLISH_SGS_OPTIONS:--max_depth 100 -bwa}

READS_DIR="$WORKDIR/reads"
MASURCA_RUN_DIR="$WORKDIR/masurca"
PLATANUS_OUT_DIR="$WORKDIR/platanus"
PLATANUS_TMP=${PLATANUS_TMP:-$PLATANUS_OUT_DIR/tmp}
PLATANUS_PREFIX_PATH="$PLATANUS_OUT_DIR/$PLATANUS_PREFIX"
PLATANUS_CONTIG="$PLATANUS_PREFIX_PATH"_contig.fa

REAPR_DIR=$WORKDIR/reapr
REAPR_OUT=$REAPR_DIR/reapr_out
PATCH_DIR=$WORKDIR/patch
POLCA_DIR=$WORKDIR/polca
NEXTPOLISH_DIR=$WORKDIR/nextpolish
NEXTPOLISH_WORKDIR=${NEXTPOLISH_WORKDIR:-$NEXTPOLISH_DIR/01_rundir}
NEXTPOLISH_SGS_FOFN=${NEXTPOLISH_SGS_FOFN:-$NEXTPOLISH_DIR/sgs.fofn}
NEXTPOLISH_INPUT=$NEXTPOLISH_DIR/nextpolish_input.fa
NEXTPOLISH_OUT=${NEXTPOLISH_OUT:-$NEXTPOLISH_WORKDIR/genome.nextpolish.fasta}

PRIMARY_ASM=${PRIMARY_ASM:-$MASURCA_RUN_DIR/CA/primary.genome.scf.fasta}
DONOR_ASM=${DONOR_ASM:-$PLATANUS_CONTIG}

reapr_checked_prefix=$REAPR_DIR/primary.checked
REAPR_CHECKED=${REAPR_CHECKED:-${reapr_checked_prefix}.fa}
REAPR_MAP_BAM=$REAPR_DIR/reapr.reads.sorted.bam
REAPR_EXT_BAM=$REAPR_DIR/reapr.external.bam
REAPR_MAP_FILTERED=$REAPR_DIR/reapr.reads.filtered.bam
REAPR_BROKEN=$REAPR_OUT/04.break.broken_assembly.fa
PATCH_FASTA=$PATCH_DIR/ragtag.patch.fasta
POLCA_INPUT=$POLCA_DIR/polca_input.fa
POLCA_RAW_OUT=${POLCA_INPUT}.Polca.fa
POLCA_FINAL=$POLCA_DIR/polca.polished.fa

mkdir -p "$WORKDIR" "$LOGDIR" "$READS_DIR" "$MASURCA_RUN_DIR" "$PLATANUS_OUT_DIR" \
  "$REAPR_DIR" "$PATCH_DIR" "$POLCA_DIR" "$NEXTPOLISH_DIR"
acquire_global_lock

require_file "$R1"
require_file "$R2"
prepared_r1=""
prepared_r2=""
platanus_r1=""
platanus_r2=""

if [[ "$RUN_DRAFT_ASSEMBLY" == "1" ]]; then
  read -r prepared_r1 prepared_r2 < <(link_or_downsample_reads "$DOWNSAMPLE_FRAC" draft)

  if MASURCA_ASM=$(find_masurca_asm); then
    log "MaSuRCA assembly exists; skipping assemble.sh"
  else
    if [[ ! -s "$MASURCA_RUN_DIR/assemble.sh" ]]; then
      require_file "$MASURCA_CFG"
      require_cmd "$MASURCA_BIN"
      MASURCA_CFG_COPY="$MASURCA_RUN_DIR/masurca.cfg"
      if [[ ! -s "$MASURCA_CFG_COPY" ]]; then
        log "Copying MaSuRCA config to run dir"
        cp "$MASURCA_CFG" "$MASURCA_CFG_COPY"
        if grep -qE '^NUM_THREADS' "$MASURCA_CFG_COPY"; then
          sed -i -E "s/^NUM_THREADS.*/NUM_THREADS = ${MASURCA_THREADS}/" "$MASURCA_CFG_COPY"
        else
          printf 'NUM_THREADS = %s\n' "$MASURCA_THREADS" >> "$MASURCA_CFG_COPY"
        fi
      fi
      run_with_log "$LOGDIR/masurca_config.log" bash -lc \
        "cd '$MASURCA_RUN_DIR' && '$MASURCA_BIN' '$MASURCA_CFG_COPY'"
    fi
    require_file "$MASURCA_RUN_DIR/assemble.sh"
    run_with_log "$LOGDIR/masurca_assemble.log" bash -lc \
      "cd '$MASURCA_RUN_DIR' && bash assemble.sh"
    MASURCA_ASM=$(find_masurca_asm) || die "MaSuRCA assembly not found after assemble.sh"
  fi
  PRIMARY_ASM=${PRIMARY_ASM:-$MASURCA_ASM}
  if [[ ! -s "$PRIMARY_ASM" && -s "$MASURCA_ASM" ]]; then
    PRIMARY_ASM="$MASURCA_ASM"
  fi
  log "MaSuRCA primary assembly: $PRIMARY_ASM"

  if exists_nonempty "$PLATANUS_CONTIG"; then
    log "Platanus contig exists; skipping"
  else
    require_cmd "$PLATANUS_BIN"
    mkdir -p "$PLATANUS_TMP"
    if ! read -r platanus_r1 platanus_r2 < <(prepare_platanus_reads "$prepared_r1" "$prepared_r2"); then
      die "failed to prepare plain FASTQ inputs for Platanus"
    fi
    log "Platanus reads: $platanus_r1 , $platanus_r2"
    run_with_log "$LOGDIR/platanus_assemble.log" \
        "$PLATANUS_BIN" assemble \
        -o "$PLATANUS_PREFIX_PATH" \
        -f "$platanus_r1" "$platanus_r2" \
        -t "$PLATANUS_THREADS" \
        -m "$PLATANUS_MEM_GB" \
        -tmp "$PLATANUS_TMP"
  fi
  exists_nonempty "$PLATANUS_CONTIG" || die "Platanus contig missing or empty after assembly: $PLATANUS_CONTIG"
  cleanup_platanus_full_read_staging "$DOWNSAMPLE_FRAC"
  DONOR_ASM=${DONOR_ASM:-$PLATANUS_CONTIG}
  log "Platanus donor assembly: $DONOR_ASM"
else
  log "Draft assembly stage disabled (RUN_DRAFT_ASSEMBLY=0)"
fi

if [[ "$RUN_REAPR_POLISH" != "1" ]]; then
  log "REAPR/RagTag/NextPolish stage disabled (RUN_REAPR_POLISH=0)"
  release_global_lock
  exit 0
fi

require_file "$PRIMARY_ASM"
require_file "$DONOR_ASM"
require_cmd "$REAPR_BIN"
if [[ -n "$REAPR_BAM" ]]; then
  require_file "$REAPR_BAM"
else
  require_reapr_aligner
fi
require_cmd samtools
require_cmd "$RAGTAG_BIN"
require_cmd minimap2
require_cmd nucmer
require_cmd bwa
if [[ "$POLISHER" == "polca" ]]; then
  require_cmd "$POLCA_BIN"
elif [[ "$POLISHER" == "nextpolish" ]]; then
  if [[ -z "$NEXTPOLISH_CMD" ]]; then
    require_cmd "$NEXTPOLISH_BIN"
  fi
  check_nextpolish_prereqs
elif [[ "$POLISHER" != "none" ]]; then
  die "POLISHER must be polca, nextpolish, or none; got: $POLISHER"
fi

log "Pipeline: MaSuRCA -> Platanus -> REAPR -> RagTag patch -> polishing"
log "Primary assembly: $PRIMARY_ASM"
log "Donor assembly:   $DONOR_ASM"
log "Reads:             $R1 , $R2"
log "Workdir:           $WORKDIR"
if [[ -z "$prepared_r1" || -z "$prepared_r2" ]]; then
  read -r prepared_r1 prepared_r2 < <(link_or_downsample_reads "$DOWNSAMPLE_FRAC" draft)
fi
log "REAPR reads: $prepared_r1 , $prepared_r2"

if exists_nonempty "$REAPR_CHECKED"; then
  log "REAPR facheck output exists; skipping facheck"
else
  run_with_log "$LOGDIR/reapr_facheck.log" \
    "$REAPR_BIN" facheck "$PRIMARY_ASM" "$reapr_checked_prefix"
fi

if [[ -n "$REAPR_BAM" ]]; then
  if bam_ready "$REAPR_EXT_BAM"; then
    log "REAPR external BAM already staged; skipping copy/sort"
  else
    cleanup_incomplete_bam "$REAPR_EXT_BAM"
    run_with_log "$LOGDIR/reapr_external_bam.log" \
      samtools_sort_index "$REAPR_BAM" "$REAPR_EXT_BAM" "$SAMTOOLS_SORT_STYLE"
  fi
  REAPR_MAP_BAM="$REAPR_EXT_BAM"
else
  if bam_ready "$REAPR_MAP_BAM"; then
    log "REAPR mapping BAM exists; skipping mapping"
  else
    cleanup_incomplete_bam "$REAPR_MAP_BAM"
    if [[ "$REAPR_ALIGNER" == "smalt" ]]; then
      tmp_bam=$REAPR_DIR/reapr.smaltmap.tmp.bam
      if ! exists_nonempty "$tmp_bam"; then
        run_with_log "$LOGDIR/reapr_smaltmap.log" \
          "$REAPR_BIN" smaltmap -n "$REAPR_THREADS" -k "$REAPR_SMALT_K" -s "$REAPR_SMALT_S" \
          "$REAPR_CHECKED" "$prepared_r1" "$prepared_r2" "$tmp_bam"
      fi
      run_with_log "$LOGDIR/reapr_smaltmap.log" \
        samtools_sort_index "$tmp_bam" "$REAPR_MAP_BAM" "$SAMTOOLS_SORT_STYLE"
      rm -f "$tmp_bam"
    elif [[ "$REAPR_ALIGNER" == "bwa-mem2" ]]; then
      bwa2_prefix=$REAPR_DIR/reapr_bwa2
      if ! bwa_mem2_index_ready "$bwa2_prefix"; then
        run_with_log "$LOGDIR/reapr_bwa2_index.log" \
          bwa-mem2 index -p "$bwa2_prefix" "$REAPR_CHECKED"
      fi
      tmp_bam=$REAPR_DIR/reapr.paired.tmp.bam
      if ! exists_nonempty "$tmp_bam"; then
        run_with_log "$LOGDIR/reapr_bwa2_map.log" bash -lc \
          "bwa-mem2 mem -t $REAPR_THREADS -a $BWA_MEM2_ARGS '$bwa2_prefix' '$prepared_r1' '$prepared_r2' | samtools view -b - > '$tmp_bam'"
      fi
      run_with_log "$LOGDIR/reapr_bwa2_sort.log" \
        samtools_sort_index "$tmp_bam" "$REAPR_MAP_BAM" "$SAMTOOLS_SORT_STYLE"
      rm -f "$tmp_bam"
    elif [[ "$REAPR_ALIGNER" == "minibwa" ]]; then
      tmp_bam=$REAPR_DIR/reapr.minibwa.tmp.bam
      if ! exists_nonempty "$tmp_bam"; then
        if [[ ! -s "${REAPR_CHECKED}.mbw" || ! -s "${REAPR_CHECKED}.l2b" ]]; then
          run_with_log "$LOGDIR/reapr_minibwa_index.log" bash -lc \
            "$MINIBWA_BIN index -t $REAPR_THREADS '$REAPR_CHECKED'"
        fi
        run_with_log "$LOGDIR/reapr_minibwa_map.log" bash -lc \
          "$MINIBWA_BIN map -t $REAPR_THREADS $MINIBWA_ARGS '$REAPR_CHECKED' '$prepared_r1' '$prepared_r2' | samtools view -b - > '$tmp_bam'"
      fi
      run_with_log "$LOGDIR/reapr_minibwa_sort.log" \
        samtools_sort_index "$tmp_bam" "$REAPR_MAP_BAM" "$SAMTOOLS_SORT_STYLE"
      rm -f "$tmp_bam"
    fi
  fi
fi

reapr_bam_for_pipeline=$REAPR_MAP_BAM
if [[ "$REAPR_FILTER_ENABLE" == "1" ]]; then
  if bam_ready "$REAPR_MAP_FILTERED"; then
    log "REAPR filtered BAM exists; skipping filter"
  else
    cleanup_incomplete_bam "$REAPR_MAP_FILTERED"
    view_args=(-b)
    if [[ "${SAMTOOLS_THREADS:-1}" -gt 1 ]] && samtools_supports_threads view; then
      view_args+=(-@ "$SAMTOOLS_THREADS")
    fi
    [[ -n "${REAPR_FILTER_F:-}" ]] && view_args+=(-f "$REAPR_FILTER_F")
    [[ -n "${REAPR_FILTER_FF:-}" ]] && view_args+=(-F "$REAPR_FILTER_FF")
    [[ "${REAPR_FILTER_MAPQ}" != "0" ]] && view_args+=(-q "$REAPR_FILTER_MAPQ")
    tmp_bam=$REAPR_DIR/reapr.reads.filtered.tmp.bam
    run_with_log "$LOGDIR/reapr_filter.log" bash -lc \
      "samtools view ${view_args[*]} '$REAPR_MAP_BAM' > '$tmp_bam'"
    run_with_log "$LOGDIR/reapr_filter.log" \
      samtools_sort_index "$tmp_bam" "$REAPR_MAP_FILTERED" "$SAMTOOLS_SORT_STYLE"
    rm -f "$tmp_bam"
  fi
  reapr_bam_for_pipeline=$REAPR_MAP_FILTERED
fi

reapr_perfect_arg=()
if [[ "$RUN_REAPR_PERFECTMAP" == "1" ]]; then
  perfect_prefix=$REAPR_DIR/$REAPR_PERFECT_PREFIX
  perfect_plot=${perfect_prefix}.plot
  if exists_nonempty "$perfect_plot"; then
    log "REAPR perfectmap output exists; skipping perfectmap"
  else
    run_with_log "$LOGDIR/reapr_perfectmap.log" \
      "$REAPR_BIN" perfectmap "$REAPR_CHECKED" "$prepared_r1" "$prepared_r2" "$REAPR_INSERT_SIZE" "$perfect_prefix"
  fi
  reapr_perfect_arg=("$perfect_prefix")
fi

if exists_nonempty "$REAPR_BROKEN"; then
  log "REAPR broken assembly exists; skipping pipeline"
else
  cmd=("$REAPR_BIN" pipeline "$REAPR_CHECKED" "$reapr_bam_for_pipeline" "$REAPR_OUT")
  [[ -n "${REAPR_SCORE_THREADS:-}" ]] && cmd+=(-score "t=${REAPR_SCORE_THREADS}")
  [[ ${#reapr_perfect_arg[@]} -gt 0 ]] && cmd+=("${reapr_perfect_arg[@]}")
  if [[ "$REAPR_BREAK_ENABLE" == "1" ]]; then
    break_parts=()
    [[ "${REAPR_BREAK_B:-0}" == "1" ]] && break_parts+=("b=1")
    [[ -n "${REAPR_BREAK_E:-}" ]] && break_parts+=("e=${REAPR_BREAK_E}")
    [[ ${#break_parts[@]} -gt 0 ]] && cmd+=(-break "${break_parts[*]}")
  fi
  run_with_log "$LOGDIR/reapr_pipeline.log" "${cmd[@]}"
fi
require_file "$REAPR_BROKEN"

if exists_nonempty "$PATCH_FASTA"; then
  log "RagTag patch output exists; skipping patch"
else
  ragtag_patch_nucmer_params="$RAGTAG_PATCH_NUCMER_PARAMS"
  if [[ ! "$ragtag_patch_nucmer_params" =~ (^|[[:space:]])-t([[:space:]]|$) ]]; then
    ragtag_patch_nucmer_params="${ragtag_patch_nucmer_params} -t ${RAGTAG_THREADS}"
  fi
  patch_cmd=("$RAGTAG_BIN" patch -u -o "$PATCH_DIR" --aligner nucmer --nucmer-params "$ragtag_patch_nucmer_params")
  patch_cmd+=("$REAPR_BROKEN" "$DONOR_ASM")
  run_with_log "$LOGDIR/ragtag_patch.log" "${patch_cmd[@]}"
fi
require_file "$PATCH_FASTA"

polished_assembly=$PATCH_FASTA
if [[ "$POLISHER" == "polca" ]]; then
  if exists_nonempty "$POLCA_FINAL"; then
    log "POLCA polished assembly exists; skipping POLCA"
  else
    ln -sf "$PATCH_FASTA" "$POLCA_INPUT"
    reads_arg="${R1} ${R2}"
    run_with_log "$LOGDIR/polca.log" \
      "$POLCA_BIN" -a "$POLCA_INPUT" -r "$reads_arg" -t "$POLCA_THREADS" -m "$POLCA_MEM_PER_THREAD"
    [[ -s "$POLCA_RAW_OUT" ]] && cp -f "$POLCA_RAW_OUT" "$POLCA_FINAL"
  fi
  require_file "$POLCA_FINAL"
  polished_assembly=$POLCA_FINAL
elif [[ "$POLISHER" == "nextpolish" ]]; then
  if exists_nonempty "$NEXTPOLISH_OUT"; then
    log "NextPolish output exists; skipping NextPolish"
  else
    ln -sf "$PATCH_FASTA" "$NEXTPOLISH_INPUT"
    [[ -z "$NEXTPOLISH_CFG" ]] && NEXTPOLISH_CFG=$NEXTPOLISH_DIR/run.cfg
    mkdir -p "$(dirname "$NEXTPOLISH_CFG")" "$NEXTPOLISH_WORKDIR"
    [[ -s "$NEXTPOLISH_SGS_FOFN" ]] || printf '%s\n%s\n' "$R1" "$R2" > "$NEXTPOLISH_SGS_FOFN"
    if [[ ! -s "$NEXTPOLISH_CFG" ]]; then
      cat > "$NEXTPOLISH_CFG" <<EOF
[General]
job_type = local
job_prefix = nextPolish
task = $NEXTPOLISH_TASK
rewrite = yes
rerun = 3
parallel_jobs = $NEXTPOLISH_PARALLEL_JOBS
multithread_jobs = $NEXTPOLISH_THREADS
genome = $NEXTPOLISH_INPUT
genome_size = auto
workdir = $NEXTPOLISH_WORKDIR
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = $NEXTPOLISH_SGS_FOFN
sgs_options = $NEXTPOLISH_SGS_OPTIONS
EOF
    fi
    [[ -z "$NEXTPOLISH_CMD" ]] && NEXTPOLISH_CMD="$NEXTPOLISH_BIN $NEXTPOLISH_CFG"
    run_with_log "$LOGDIR/nextpolish.log" bash -lc "$NEXTPOLISH_CMD"
  fi
  require_file "$NEXTPOLISH_OUT"
  polished_assembly=$NEXTPOLISH_OUT
else
  log "Polishing disabled (POLISHER=$POLISHER); using patch output"
fi

log "Done. Key outputs:"
log "  MaSuRCA primary:       $PRIMARY_ASM"
log "  Platanus donor:        $DONOR_ASM"
log "  REAPR broken assembly: $REAPR_BROKEN"
log "  RagTag patch assembly: $PATCH_FASTA"
log "  Final polished assembly: $polished_assembly"
release_global_lock
