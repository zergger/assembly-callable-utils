#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
umask 0022

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# Prefer local wrapper scripts in this directory (e.g. ragtag_rename.py override).
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

cleanup() {
  local sig=${1:-INT}
  set +e +u
  if [[ -n "${GLOBAL_LOCK_FD:-}" ]]; then
    flock -u "$GLOBAL_LOCK_FD" 2>/dev/null || true
    exec {GLOBAL_LOCK_FD}>&- 2>/dev/null || true
  fi
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
Usage: reapr_patch_nextpolish_ragtag.sh

This pipeline intentionally follows a conservative,
read-backed sequence:
  1) REAPR breakpoint detection and breaking
  2) RagTag patch using a donor assembly
  3) NextPolish Illumina polishing (default; POLCA optional)
  4) RagTag correct with read validation (optional; requires reference)
  5) RagTag reference-guided scaffolding (optional; requires reference)

Key environment variables (defaults in brackets):
  THREADS=10
  REAPR_THREADS=10
  POLCA_THREADS=10
  RAGTAG_THREADS=10
  SAMTOOLS_THREADS=1

Inputs:
  PRIMARY_ASM=/home/assem/SJF/work_masurca_platanus_ragtag/masurca/CA/primary.genome.scf.fasta
  DONOR_ASM=/home/assem/SJF/work_masurca_platanus_ragtag/platanus/platanus_alt_contig.fa
  REFERENCE_ASM= (optional; enables RagTag scaffolding)
  R1=/home/fqdata/clean/SJF/SJF_R1.fastq.gz
  R2=/home/fqdata/clean/SJF/SJF_R2.fastq.gz

Workdir and toggles:
  WORKDIR=/home/assem/SJF/work_reapr_patch_nextpolish_ragtag
  DOWNSAMPLE_FRAC=1.0 (set <1.0 to downsample reads with seqtk)
  SEED=42
  SAMPLE_ID=SJF
  RUN_SORTNR=0
  SORTNR_BASE=$WORKDIR/dedup/${SAMPLE_ID}_reapr_patch_nextpolish_ragtag

RagTag correct options (effective only when REFERENCE_ASM is set):
  RUN_RAGTAG_CORRECT=1
  RAGTAG_CORRECT_VALIDATE=1
  RAGTAG_READTYPE=sr
  RAGTAG_READ_ALIGNER=minimap2

RagTag patch options:
  RAGTAG_PATCH_NUCMER_PARAMS="--maxmatch -l 100 -c 500"
  RAGTAG_FILTER_PROCS=10 (ragtag_patch.py filter workers; default follows THREADS)
  (script invokes `ragtag patch -u` to avoid AGP ID collisions)
  (script will auto-append "-t RAGTAG_THREADS" if -t is not present)
  GLOBAL_LOCK_FILE=$WORKDIR/.pipeline.lock
  GLOBAL_LOCK_WAIT_SEC=0 (0 means fail-fast if lock is busy)
  (backward compatible: PATCH_LOCK_FILE/PATCH_LOCK_WAIT_SEC are still accepted)

Polishing options:
  POLISHER=nextpolish (nextpolish|polca|none)
  POLCA_MEM_PER_THREAD=1G
  NEXTPOLISH_BIN=nextpolish
  NEXTPOLISH_CFG= (auto-generated if empty)
  NEXTPOLISH_CMD= (optional; overrides NEXTPOLISH_CFG)
  NEXTPOLISH_OUT= (expected polished FASTA output)
  NEXTPOLISH_WORKDIR= (default: workdir/nextpolish/01_rundir)
  NEXTPOLISH_SGS_FOFN= (default: workdir/nextpolish/sgs.fofn)
  NEXTPOLISH_PARALLEL_JOBS=1
  NEXTPOLISH_THREADS=10
  NEXTPOLISH_TASK=best
  NEXTPOLISH_SGS_OPTIONS="-max_depth 100 -bwa"

REAPR options:
  RUN_REAPR_PERFECTMAP=0 (set 1 to run perfectmap using the same reads)
  REAPR_INSERT_SIZE=220 (used only if RUN_REAPR_PERFECTMAP=1)
  REAPR_SMALT_K=13
  REAPR_SMALT_S=6
  REAPR_BREAK_ENABLE=1
  REAPR_BREAK_B=0
  REAPR_BREAK_E=0.5
  REAPR_SCORE_THREADS=10
  REAPR_FILTER_ENABLE=0
  REAPR_FILTER_F=2
  REAPR_FILTER_FF=904
  REAPR_FILTER_MAPQ=0
  SAMTOOLS_SORT_STYLE=auto (auto|new|old)
  REAPR_ALIGNER=smalt (smalt|bwa-mem2|strobealign)
  REAPR_BAM= (optional; precomputed paired BAM)
  BWA_MEM2_ARGS=-SP
  STROBEALIGN_ARGS="--mcs=always -R 0"

Tool overrides:
  PY3=python3
  DIPLOIDOCUS_PY=/app/diploidocus/code/diploidocus.py
  DIPLOIDOCUS_INTERACTIVE=-1
  DIPLOIDOCUS_FORKS=$THREADS
  REAPR_BIN=reapr
  RAGTAG_BIN=ragtag.py
  RAGTAG_PATCH_BIN=$SCRIPT_DIR/ragtag_patch.py
  POLCA_BIN=polca.sh
  SEQTK_BIN=seqtk
  NOTE: scripts/ is prepended to PATH, so scripts/ragtag_rename.py is preferred.

Resume behavior:
  Each step is skipped if its primary output FASTA/BAM already exists
  and is non-empty. Logs are written under WORKDIR/logs/.
USAGE
  exit 0
fi

THREADS=${THREADS:-10}
REAPR_THREADS=${REAPR_THREADS:-$THREADS}
POLCA_THREADS=${POLCA_THREADS:-$THREADS}
RAGTAG_THREADS=${RAGTAG_THREADS:-$THREADS}
SAMTOOLS_THREADS=${SAMTOOLS_THREADS:-1}

PRIMARY_ASM=${PRIMARY_ASM:-/home/assem/SJF/work_masurca_platanus_ragtag/masurca/CA/primary.genome.scf.fasta}
DONOR_ASM=${DONOR_ASM:-/home/assem/SJF/work_masurca_platanus_ragtag/platanus/platanus_alt_contig.fa}
REFERENCE_ASM=${REFERENCE_ASM:-}
R1=${R1:-/home/fqdata/clean/SJF/SJF_R1.fastq.gz}
R2=${R2:-/home/fqdata/clean/SJF/SJF_R2.fastq.gz}

WORKDIR=${WORKDIR:-/home/assem/SJF/work_reapr_patch_nextpolish_ragtag}
LOGDIR=${LOGDIR:-$WORKDIR/logs}
READS_DIR=$WORKDIR/reads
REAPR_DIR=$WORKDIR/reapr
REAPR_OUT=$REAPR_DIR/reapr_out
PATCH_DIR=$WORKDIR/patch
POLCA_DIR=$WORKDIR/polca
NEXTPOLISH_DIR=$WORKDIR/nextpolish
SCAFFOLD_DIR=$WORKDIR/scaffold
CORRECT_DIR=$WORKDIR/correct
DEDUP_DIR=$WORKDIR/dedup

DOWNSAMPLE_FRAC=${DOWNSAMPLE_FRAC:-1.0}
SEED=${SEED:-42}
SAMPLE_ID=${SAMPLE_ID:-SJF}
RUN_SORTNR=${RUN_SORTNR:-0}
SORTNR_BASE=${SORTNR_BASE:-$DEDUP_DIR/${SAMPLE_ID}_reapr_patch_nextpolish_ragtag}

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
REAPR_ALIGNER=${REAPR_ALIGNER:-smalt}
REAPR_BAM=${REAPR_BAM:-}
BWA_MEM2_ARGS=${BWA_MEM2_ARGS:--SP}
STROBEALIGN_BIN=${STROBEALIGN_BIN:-strobealign}
STROBEALIGN_ARGS=${STROBEALIGN_ARGS:---mcs=always -R 0}

RUN_RAGTAG_CORRECT=${RUN_RAGTAG_CORRECT:-1}
RAGTAG_CORRECT_VALIDATE=${RAGTAG_CORRECT_VALIDATE:-1}
RAGTAG_READTYPE=${RAGTAG_READTYPE:-sr}
RAGTAG_READ_ALIGNER=${RAGTAG_READ_ALIGNER:-minimap2}
RAGTAG_PATCH_NUCMER_PARAMS=${RAGTAG_PATCH_NUCMER_PARAMS:---maxmatch -l 100 -c 500}
RAGTAG_FILTER_PROCS=${RAGTAG_FILTER_PROCS:-$THREADS}
GLOBAL_LOCK_WAIT_SEC=${GLOBAL_LOCK_WAIT_SEC:-${PATCH_LOCK_WAIT_SEC:-0}}
GLOBAL_LOCK_FILE=${GLOBAL_LOCK_FILE:-${PATCH_LOCK_FILE:-$WORKDIR/.pipeline.lock}}

REAPR_BIN=${REAPR_BIN:-reapr}
RAGTAG_BIN=${RAGTAG_BIN:-ragtag.py}
RAGTAG_PATCH_BIN=${RAGTAG_PATCH_BIN:-$SCRIPT_DIR/ragtag_patch.py}
POLCA_BIN=${POLCA_BIN:-polca.sh}
SEQTK_BIN=${SEQTK_BIN:-seqtk}
PY3=${PY3:-python3}
DIPLOIDOCUS_PY=${DIPLOIDOCUS_PY:-/app/diploidocus/code/diploidocus.py}
DIPLOIDOCUS_INTERACTIVE=${DIPLOIDOCUS_INTERACTIVE:--1}
DIPLOIDOCUS_FORKS=${DIPLOIDOCUS_FORKS:-$THREADS}

POLCA_MEM_PER_THREAD=${POLCA_MEM_PER_THREAD:-1G}
POLISHER=${POLISHER:-nextpolish}
NEXTPOLISH_BIN=${NEXTPOLISH_BIN:-nextPolish}
NEXTPOLISH_CFG=${NEXTPOLISH_CFG:-}
NEXTPOLISH_CMD=${NEXTPOLISH_CMD:-}
NEXTPOLISH_WORKDIR=${NEXTPOLISH_WORKDIR:-$NEXTPOLISH_DIR/01_rundir}
NEXTPOLISH_SGS_FOFN=${NEXTPOLISH_SGS_FOFN:-$NEXTPOLISH_DIR/sgs.fofn}
NEXTPOLISH_PARALLEL_JOBS=${NEXTPOLISH_PARALLEL_JOBS:-1}
NEXTPOLISH_THREADS=${NEXTPOLISH_THREADS:-$THREADS}
NEXTPOLISH_TASK=${NEXTPOLISH_TASK:-best}
NEXTPOLISH_SGS_OPTIONS=${NEXTPOLISH_SGS_OPTIONS:--max_depth 100 -bwa}
NEXTPOLISH_INPUT=$NEXTPOLISH_DIR/nextpolish_input.fa
NEXTPOLISH_OUT=${NEXTPOLISH_OUT:-$NEXTPOLISH_WORKDIR/genome.nextpolish.fasta}

reapr_checked_prefix=$REAPR_DIR/primary.checked
REAPR_CHECKED=${REAPR_CHECKED:-${reapr_checked_prefix}.fa}
REAPR_MAP_BAM=$REAPR_DIR/reapr.smaltmap.bam
REAPR_EXT_BAM=$REAPR_DIR/reapr.external.bam
REAPR_MAP_FILTERED=$REAPR_DIR/reapr.smaltmap.filtered.bam
REAPR_BROKEN=$REAPR_OUT/04.break.broken_assembly.fa
REAPR_ERRORS_GFF=$REAPR_OUT/03.score.errors.gff.gz

PATCH_FASTA=$PATCH_DIR/ragtag.patch.fasta
POLCA_INPUT=$POLCA_DIR/polca_input.fa
POLCA_RAW_OUT=${POLCA_INPUT}.Polca.fa
POLCA_FINAL=$POLCA_DIR/polca.polished.fa
CORRECT_FASTA=$CORRECT_DIR/ragtag.correct.fasta
SCAFFOLD_FASTA=$SCAFFOLD_DIR/ragtag.scaffold.fasta
SORTNR_OUT=${SORTNR_BASE}.nr.fasta

log() {
  printf '[%(%F %T)T] %s\n' -1 "$*" >&2
}

exists_nonempty() {
  local p=$1
  [[ -s "$p" ]]
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

require_file() {
  local p=$1
  if [[ ! -f "$p" ]]; then
    echo "ERROR: required file not found: $p" >&2
    exit 1
  fi
}

require_cmd() {
  local c=$1
  if ! command -v "$c" >/dev/null 2>&1; then
    echo "ERROR: required command not found on PATH: $c" >&2
    exit 1
  fi
}

ragtag_patch_supports_filter_procs() {
  "$RAGTAG_PATCH_BIN" -h 2>&1 | grep -q -- "--filter-procs"
}

samtools_supports_threads() {
  local subcmd=$1
  samtools "$subcmd" 2>&1 | grep -q -- " -@"
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
    echo "ERROR: Please use a newer samtools (recommended >=1.10, ideally 1.17/1.18)." >&2
    echo "ERROR: current samtools path: $(command -v samtools || echo 'not found')" >&2
    samtools 2>&1 | head -n 3 >&2 || true
    exit 1
  fi
  if ! samtools_fixmate_supports_m; then
    echo "ERROR: NextPolish requires 'samtools fixmate -m', but current samtools does not support '-m'." >&2
    echo "ERROR: Please use a newer samtools (recommended >=1.10, ideally 1.17/1.18)." >&2
    echo "ERROR: current samtools path: $(command -v samtools || echo 'not found')" >&2
    samtools 2>&1 | head -n 3 >&2 || true
    exit 1
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
    echo "ERROR: 'flock' is required for pipeline locking but not found on PATH." >&2
    exit 1
  fi
  mkdir -p "$(dirname "$GLOBAL_LOCK_FILE")"
  exec {GLOBAL_LOCK_FD}> "$GLOBAL_LOCK_FILE"
  if [[ "$GLOBAL_LOCK_WAIT_SEC" -gt 0 ]]; then
    flock -w "$GLOBAL_LOCK_WAIT_SEC" "$GLOBAL_LOCK_FD" || {
      echo "ERROR: pipeline lock busy ($GLOBAL_LOCK_FILE), waited ${GLOBAL_LOCK_WAIT_SEC}s." >&2
      exit 1
    }
  else
    flock -n "$GLOBAL_LOCK_FD" || {
      echo "ERROR: pipeline lock busy ($GLOBAL_LOCK_FILE). Another pipeline run may be running." >&2
      exit 1
    }
  fi
  log "Acquired global pipeline lock: $GLOBAL_LOCK_FILE"
}

release_global_lock() {
  if [[ -n "${GLOBAL_LOCK_FD:-}" ]]; then
    flock -u "$GLOBAL_LOCK_FD" 2>/dev/null || true
    exec {GLOBAL_LOCK_FD}>&- 2>/dev/null || true
    unset GLOBAL_LOCK_FD
    log "Released global pipeline lock: $GLOBAL_LOCK_FILE"
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
  local tmp_out="${out_bam}.tmp.$$"
  if [[ "${SAMTOOLS_THREADS:-1}" -gt 1 ]] && samtools_supports_threads sort; then
    sort_threads=(-@ "$SAMTOOLS_THREADS")
  fi
  if [[ "${SAMTOOLS_THREADS:-1}" -gt 1 ]] && samtools_supports_threads index; then
    index_threads=(-@ "$SAMTOOLS_THREADS")
  fi
  if [[ "$style" == "auto" ]]; then
    if samtools sort 2>&1 | grep -q -- " -o"; then
      style="new"
    else
      style="old"
    fi
  fi
  if [[ "$style" == "new" ]]; then
    samtools sort "${sort_threads[@]}" -o "$tmp_out" "$in_bam"
  else
    local prefix=${tmp_out%.bam}
    samtools sort "${sort_threads[@]}" "$in_bam" "$prefix"
    tmp_out="${prefix}.bam"
  fi
  samtools index "${index_threads[@]}" "$tmp_out"
  if [[ -s "${tmp_out}.bai" ]]; then
    mv -f "$tmp_out" "$out_bam"
    mv -f "${tmp_out}.bai" "${out_bam}.bai"
  elif [[ -s "${tmp_out}.csi" ]]; then
    mv -f "$tmp_out" "$out_bam"
    mv -f "${tmp_out}.csi" "${out_bam}.csi"
  else
    echo "ERROR: samtools index did not produce .bai/.csi for $tmp_out" >&2
    exit 1
  fi
}

require_reapr_aligner() {
  case "$REAPR_ALIGNER" in
    smalt)
      require_cmd smalt
      ;;
    bwa-mem2)
      require_cmd bwa-mem2
      ;;
    strobealign)
      require_cmd "$STROBEALIGN_BIN"
      ;;
    *)
      echo "ERROR: REAPR_ALIGNER must be smalt, bwa-mem2, or strobealign; got: $REAPR_ALIGNER" >&2
      exit 1
      ;;
  esac
}

mkdir -p "$WORKDIR" "$LOGDIR" "$READS_DIR" "$REAPR_DIR" "$PATCH_DIR" "$POLCA_DIR" "$NEXTPOLISH_DIR" "$CORRECT_DIR" "$SCAFFOLD_DIR" "$DEDUP_DIR"
acquire_global_lock

require_file "$PRIMARY_ASM"
require_file "$DONOR_ASM"
require_file "$R1"
require_file "$R2"
if [[ -n "$REFERENCE_ASM" ]]; then
  require_file "$REFERENCE_ASM"
fi

require_cmd "$REAPR_BIN"
if [[ -n "$REAPR_BAM" ]]; then
  require_file "$REAPR_BAM"
else
  require_reapr_aligner
fi
require_cmd samtools
require_cmd "$RAGTAG_BIN"
require_cmd "$RAGTAG_PATCH_BIN"
require_cmd minimap2
require_cmd nucmer
if [[ "$POLISHER" == "polca" ]]; then
  require_cmd "$POLCA_BIN"
elif [[ "$POLISHER" == "nextpolish" ]]; then
  if [[ -z "$NEXTPOLISH_CMD" ]]; then
    require_cmd "$NEXTPOLISH_BIN"
  fi
  check_nextpolish_prereqs
elif [[ "$POLISHER" == "none" ]]; then
  :
else
  echo "ERROR: POLISHER must be polca, nextpolish, or none; got: $POLISHER" >&2
  exit 1
fi
require_cmd bwa
if [[ "$RUN_SORTNR" == "1" ]]; then
  require_cmd "$PY3"
  require_file "$DIPLOIDOCUS_PY"
fi

log "Pipeline: REAPR -> RagTag patch -> polishing (NextPolish default; POLCA optional) -> RagTag correct -> RagTag scaffold"
log "Primary assembly: $PRIMARY_ASM"
log "Donor assembly:   $DONOR_ASM"
if [[ -n "$REFERENCE_ASM" ]]; then
  log "Reference:         $REFERENCE_ASM"
else
  log "Reference:         (none; scaffolding will be skipped)"
fi
log "Reads:             $R1 , $R2"
log "Workdir:           $WORKDIR"
if [[ -z "$REAPR_BAM" && "$REAPR_ALIGNER" != "smalt" ]]; then
  log "WARNING: REAPR_ALIGNER=$REAPR_ALIGNER is not equivalent to SMALT -x; paired scoring/MAPQ may differ from SMALT results."
fi

current_r1=$R1
current_r2=$R2

if awk "BEGIN{exit !($DOWNSAMPLE_FRAC < 1.0)}"; then
  require_cmd "$SEQTK_BIN"
  ds_r1=$READS_DIR/${SAMPLE_ID}_R1.ds.fq
  ds_r2=$READS_DIR/${SAMPLE_ID}_R2.ds.fq
  if exists_nonempty "$ds_r1" && exists_nonempty "$ds_r2"; then
    log "Downsampled reads already exist; skipping downsampling"
  else
    run_with_log "$LOGDIR/downsample.log" bash -lc \
      "$SEQTK_BIN sample -s${SEED} $R1 $DOWNSAMPLE_FRAC > $ds_r1"
    run_with_log "$LOGDIR/downsample.log" bash -lc \
      "$SEQTK_BIN sample -s${SEED} $R2 $DOWNSAMPLE_FRAC > $ds_r2"
  fi
  current_r1=$ds_r1
  current_r2=$ds_r2
else
  suffix=""
  case "$R1" in
    *.fq.gz) suffix=".fq.gz" ;;
    *.fastq.gz) suffix=".fastq.gz" ;;
    *.fq) suffix=".fq" ;;
    *.fastq) suffix=".fastq" ;;
    *) suffix="" ;;
  esac
  ln -sf "$R1" "$READS_DIR/${SAMPLE_ID}_R1${suffix}"
  ln -sf "$R2" "$READS_DIR/${SAMPLE_ID}_R2${suffix}"
  current_r1=$READS_DIR/${SAMPLE_ID}_R1${suffix}
  current_r2=$READS_DIR/${SAMPLE_ID}_R2${suffix}
fi

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
    log "Using external paired BAM for REAPR: $REAPR_BAM"
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
      if exists_nonempty "$tmp_bam"; then
        log "REAPR smaltmap temp BAM exists; skipping mapping"
      else
        run_with_log "$LOGDIR/reapr_smaltmap.log" \
          "$REAPR_BIN" smaltmap -n "$REAPR_THREADS" -k "$REAPR_SMALT_K" -s "$REAPR_SMALT_S" \
          "$REAPR_CHECKED" "$current_r1" "$current_r2" "$tmp_bam"
      fi
      run_with_log "$LOGDIR/reapr_smaltmap.log" \
        samtools_sort_index "$tmp_bam" "$REAPR_MAP_BAM" "$SAMTOOLS_SORT_STYLE"
      rm -f "$tmp_bam"
    elif [[ "$REAPR_ALIGNER" == "bwa-mem2" ]]; then
      bwa2_prefix=$REAPR_DIR/reapr_bwa2
      if [[ ! -s "${bwa2_prefix}.bwt.2bit" ]]; then
        run_with_log "$LOGDIR/reapr_bwa2_index.log" \
          bwa-mem2 index -p "$bwa2_prefix" "$REAPR_CHECKED"
      fi
      tmp_bam=$REAPR_DIR/reapr.paired.tmp.bam
      if exists_nonempty "$tmp_bam"; then
        log "REAPR bwa-mem2 temp BAM exists; skipping mapping"
      else
        run_with_log "$LOGDIR/reapr_bwa2_map.log" bash -lc \
          "bwa-mem2 mem -t $REAPR_THREADS -a $BWA_MEM2_ARGS '$bwa2_prefix' '$current_r1' '$current_r2' | samtools view -b - > '$tmp_bam'"
      fi
      run_with_log "$LOGDIR/reapr_bwa2_sort.log" \
        samtools_sort_index "$tmp_bam" "$REAPR_MAP_BAM" "$SAMTOOLS_SORT_STYLE"
      rm -f "$tmp_bam"
    else
      tmp_bam=$REAPR_DIR/reapr.paired.tmp.bam
      if exists_nonempty "$tmp_bam"; then
        log "REAPR strobealign temp BAM exists; skipping mapping"
      else
        run_with_log "$LOGDIR/reapr_strobealign_map.log" bash -lc \
          "$STROBEALIGN_BIN -t $REAPR_THREADS $STROBEALIGN_ARGS '$REAPR_CHECKED' '$current_r1' '$current_r2' | samtools view -b - > '$tmp_bam'"
      fi
      run_with_log "$LOGDIR/reapr_strobealign_sort.log" \
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
    if [[ "${REAPR_FILTER_MAPQ}" != "0" ]]; then
      view_args+=(-q "$REAPR_FILTER_MAPQ")
    fi
    tmp_bam=$REAPR_DIR/reapr.smaltmap.filtered.tmp.bam
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
      "$REAPR_BIN" perfectmap "$REAPR_CHECKED" "$current_r1" "$current_r2" "$REAPR_INSERT_SIZE" "$perfect_prefix"
  fi
  reapr_perfect_arg=("$perfect_prefix")
fi

if exists_nonempty "$REAPR_BROKEN"; then
  log "REAPR broken assembly exists; skipping pipeline"
else
  cmd=("$REAPR_BIN" pipeline "$REAPR_CHECKED" "$reapr_bam_for_pipeline" "$REAPR_OUT")
  if [[ -n "${REAPR_SCORE_THREADS:-}" ]]; then
    cmd+=(-score "t=${REAPR_SCORE_THREADS}")
  fi
  if [[ ${#reapr_perfect_arg[@]} -gt 0 ]]; then
    cmd+=("${reapr_perfect_arg[@]}")
  fi
  if [[ "$REAPR_BREAK_ENABLE" == "1" ]]; then
    break_parts=()
    # REAPR pipeline expects "-break b=1" only for enabling -b.
    # Passing b=0 becomes "reapr break -b 0" and breaks argument parsing.
    if [[ "${REAPR_BREAK_B:-0}" == "1" ]]; then
      break_parts+=("b=1")
    fi
    [[ -n "${REAPR_BREAK_E:-}" ]] && break_parts+=("e=${REAPR_BREAK_E}")
    if [[ ${#break_parts[@]} -gt 0 ]]; then
      break_arg="${break_parts[*]}"
      cmd+=(-break "$break_arg")
    fi
  fi
  run_with_log "$LOGDIR/reapr_pipeline.log" "${cmd[@]}"
fi

if ! exists_nonempty "$REAPR_BROKEN"; then
  echo "ERROR: expected REAPR broken assembly not found: $REAPR_BROKEN" >&2
  exit 1
fi

if exists_nonempty "$PATCH_FASTA"; then
  log "RagTag patch output exists; skipping patch"
else
  ragtag_patch_nucmer_params="$RAGTAG_PATCH_NUCMER_PARAMS"
  if [[ ! "$ragtag_patch_nucmer_params" =~ (^|[[:space:]])-t([[:space:]]|$) ]]; then
    ragtag_patch_nucmer_params="${ragtag_patch_nucmer_params} -t ${RAGTAG_THREADS}"
  fi
  patch_cmd=("$RAGTAG_PATCH_BIN" -u -o "$PATCH_DIR" --aligner nucmer --nucmer-params "$ragtag_patch_nucmer_params")
  if ragtag_patch_supports_filter_procs; then
    patch_cmd+=(--filter-procs "$RAGTAG_FILTER_PROCS")
  else
    log "WARNING: current RagTag patch does not support --filter-procs; running without it"
  fi
  patch_cmd+=("$REAPR_BROKEN" "$DONOR_ASM")
  run_with_log "$LOGDIR/ragtag_patch.log" "${patch_cmd[@]}"
fi

if ! exists_nonempty "$PATCH_FASTA"; then
  echo "ERROR: expected RagTag patch output not found: $PATCH_FASTA" >&2
  exit 1
fi

polished_assembly=$PATCH_FASTA
if [[ "$POLISHER" == "polca" ]]; then
  if exists_nonempty "$POLCA_FINAL"; then
    log "POLCA polished assembly exists; skipping POLCA"
  else
    ln -sf "$PATCH_FASTA" "$POLCA_INPUT"
    reads_arg="${R1} ${R2}"
    run_with_log "$LOGDIR/polca.log" \
      "$POLCA_BIN" -a "$POLCA_INPUT" -r "$reads_arg" -t "$POLCA_THREADS" -m "$POLCA_MEM_PER_THREAD"
    if exists_nonempty "$POLCA_RAW_OUT"; then
      cp -f "$POLCA_RAW_OUT" "$POLCA_FINAL"
    fi
  fi
  if ! exists_nonempty "$POLCA_FINAL"; then
    echo "ERROR: expected POLCA output not found: $POLCA_FINAL" >&2
    exit 1
  fi
  polished_assembly=$POLCA_FINAL
elif [[ "$POLISHER" == "nextpolish" ]]; then
  if exists_nonempty "$NEXTPOLISH_OUT"; then
    log "NextPolish output exists; skipping NextPolish"
  else
    ln -sf "$PATCH_FASTA" "$NEXTPOLISH_INPUT"
    if [[ -z "$NEXTPOLISH_CFG" ]]; then
      NEXTPOLISH_CFG=$NEXTPOLISH_DIR/run.cfg
    fi
    mkdir -p "$(dirname "$NEXTPOLISH_CFG")" "$NEXTPOLISH_WORKDIR"
    if [[ ! -s "$NEXTPOLISH_SGS_FOFN" ]]; then
      printf '%s\n%s\n' "$R1" "$R2" > "$NEXTPOLISH_SGS_FOFN"
    fi
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
    if [[ -z "$NEXTPOLISH_CMD" ]]; then
      NEXTPOLISH_CMD="$NEXTPOLISH_BIN $NEXTPOLISH_CFG"
    fi
    run_with_log "$LOGDIR/nextpolish.log" bash -lc "$NEXTPOLISH_CMD"
  fi
  if ! exists_nonempty "$NEXTPOLISH_OUT"; then
    echo "ERROR: expected NextPolish output not found: $NEXTPOLISH_OUT" >&2
    exit 1
  fi
  polished_assembly=$NEXTPOLISH_OUT
else
  log "Polishing disabled (POLISHER=$POLISHER); using patch output"
fi

corrected_assembly=$polished_assembly
if [[ -n "$REFERENCE_ASM" && "$RUN_RAGTAG_CORRECT" == "1" ]]; then
  correct_validate_arg=()
  if [[ "$RAGTAG_CORRECT_VALIDATE" == "1" ]]; then
    if [[ -n "$current_r2" && "$current_r2" != "$current_r1" ]]; then
      correct_validate_fofn=$CORRECT_DIR/ragtag_correct.reads.fofn
      printf '%s\n%s\n' "$current_r1" "$current_r2" > "$correct_validate_fofn"
      correct_validate_arg=(-F "$correct_validate_fofn" -T "$RAGTAG_READTYPE")
    else
      correct_validate_arg=(-R "$current_r1" -T "$RAGTAG_READTYPE")
    fi
    if [[ -n "$RAGTAG_READ_ALIGNER" ]]; then
      correct_validate_arg+=(--read-aligner "$RAGTAG_READ_ALIGNER")
    fi
  fi

  if exists_nonempty "$CORRECT_FASTA"; then
    log "RagTag correct output exists; skipping correct"
  else
    run_with_log "$LOGDIR/ragtag_correct.log" \
      "$RAGTAG_BIN" correct -o "$CORRECT_DIR" -t "$RAGTAG_THREADS" \
      "${correct_validate_arg[@]}" \
      "$REFERENCE_ASM" "$polished_assembly"
  fi

  if exists_nonempty "$CORRECT_FASTA"; then
    corrected_assembly=$CORRECT_FASTA
  else
    echo "ERROR: expected RagTag correct output not found: $CORRECT_FASTA" >&2
    exit 1
  fi
elif [[ -n "$REFERENCE_ASM" && "$RUN_RAGTAG_CORRECT" != "1" ]]; then
  log "RagTag correct disabled (RUN_RAGTAG_CORRECT=$RUN_RAGTAG_CORRECT); scaffolding will use polished output"
fi

final_assembly=$corrected_assembly
if [[ -n "$REFERENCE_ASM" ]]; then
  if exists_nonempty "$SCAFFOLD_FASTA"; then
    log "RagTag scaffold output exists; skipping scaffold"
  else
    run_with_log "$LOGDIR/ragtag_scaffold.log" \
      "$RAGTAG_BIN" scaffold -o "$SCAFFOLD_DIR" -t "$RAGTAG_THREADS" "$REFERENCE_ASM" "$corrected_assembly"
  fi
  if exists_nonempty "$SCAFFOLD_FASTA"; then
    final_assembly=$SCAFFOLD_FASTA
  else
    echo "WARNING: RagTag scaffold did not produce $SCAFFOLD_FASTA; using polished output" >&2
  fi
fi

if [[ "$RUN_SORTNR" == "1" ]]; then
  if exists_nonempty "$SORTNR_OUT"; then
    log "Diploidocus sortnr output exists; skipping sortnr"
  else
    log "Running diploidocus sortnr"
    run_with_log "$LOGDIR/diploidocus_sortnr.log" \
      "$PY3" "$DIPLOIDOCUS_PY" runmode=sortnr seqin="$final_assembly" basefile="$SORTNR_BASE" "i=$DIPLOIDOCUS_INTERACTIVE" "forks=$DIPLOIDOCUS_FORKS"
  fi
  if ! exists_nonempty "$SORTNR_OUT"; then
    echo "ERROR: expected Diploidocus sortnr output not found: $SORTNR_OUT" >&2
    exit 1
  fi
  final_assembly=$SORTNR_OUT
fi

log "Done. Key outputs:"
log "  REAPR broken assembly: $REAPR_BROKEN"
log "  RagTag patch assembly: $PATCH_FASTA"
if [[ "$POLISHER" == "polca" ]]; then
  log "  POLCA polished:        $POLCA_FINAL"
elif [[ "$POLISHER" == "nextpolish" ]]; then
  log "  NextPolish polished:   $NEXTPOLISH_OUT"
else
  log "  Polishing:             (disabled)"
fi
if [[ -n "$REFERENCE_ASM" ]]; then
  log "  RagTag correct:        $CORRECT_FASTA"
  log "  RagTag scaffold:       $SCAFFOLD_FASTA"
fi
if [[ "$RUN_SORTNR" == "1" ]]; then
  log "  Diploidocus sortnr:    $SORTNR_OUT"
fi
log "  Final assembly:        $final_assembly"
release_global_lock
