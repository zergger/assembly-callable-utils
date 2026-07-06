#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
umask 0022

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
DEFAULT_CONFIG="$ROOT_DIR/configs/ragtag_project.env"
export PATH="$SCRIPT_DIR:$PATH"

usage() {
  cat <<USAGE
Usage: bash scripts/ragtag_project.sh [config.env]

Project an already polished assembly onto one reference with RagTag.
This is the clean post-polish projection stage for paper-style A/B branches.
It does not run assembly refinement, patching, polishing, sortnr, or callable masking.

Required inputs:
  POLISHED_ASM          polished/unanchored assembly FASTA
  REFERENCE_ASM         reference FASTA used by RagTag correct/scaffold

Read-backed correction inputs when RAGTAG_CORRECT_VALIDATE=1:
  R1                    read 1 FASTQ/FASTA, optionally gzip-compressed
  R2                    read 2 FASTQ/FASTA, optionally gzip-compressed

Common variables:
  SAMPLE_TAG            output label [basename of WORKDIR]
  WORKDIR               output workdir [results/ragtag_project/<SAMPLE_TAG>]
  THREADS               default thread count [10]
  RAGTAG_THREADS        RagTag thread count [THREADS]
  RUN_RAGTAG_CORRECT    run ragtag correct before scaffold [1]
  RAGTAG_CORRECT_VALIDATE  pass reads to ragtag correct [1]
  RUN_RAGTAG_SCAFFOLD   run ragtag scaffold [1]
  FORCE                 remove and recompute correct/scaffold dirs [0]

RagTag options:
  RAGTAG_READTYPE       read type for RagTag correct validation [sr]
  RAGTAG_READ_ALIGNER   read aligner passed to RagTag correct [minimap2]
  RAGTAG_BIN            RagTag executable [ragtag.py]

Tool/log options:
  SAMTOOLS              samtools executable [samtools]
  LOGDIR                log directory [WORKDIR/logs]
  GLOBAL_LOCK_FILE      lock file [WORKDIR/.ragtag_project.lock]
  GLOBAL_LOCK_WAIT_SEC  seconds to wait for lock; 0 means fail-fast [0]

Primary outputs:
  WORKDIR/correct/ragtag.correct.fasta
  WORKDIR/correct/ragtag.correct.reads.s.bam  (when read validation is enabled)
  WORKDIR/scaffold/ragtag.scaffold.fasta
  WORKDIR/scaffold/ragtag.scaffold.agp
  WORKDIR/ragtag_project.outputs.env
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

require_cmd() {
  local c="$1"
  has_exec "$c" || { echo "ERROR: required command not found or not executable: $c" >&2; exit 1; }
}

exists_nonempty() {
  [[ -s "$1" ]]
}

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
    echo "ERROR: flock is required for locking but was not found." >&2
    exit 1
  fi
  mkdir -p "$(dirname "$GLOBAL_LOCK_FILE")"
  exec {GLOBAL_LOCK_FD}> "$GLOBAL_LOCK_FILE"
  if [[ "$GLOBAL_LOCK_WAIT_SEC" -gt 0 ]]; then
    flock -w "$GLOBAL_LOCK_WAIT_SEC" "$GLOBAL_LOCK_FD" || {
      echo "ERROR: lock busy ($GLOBAL_LOCK_FILE), waited ${GLOBAL_LOCK_WAIT_SEC}s." >&2
      exit 1
    }
  else
    flock -n "$GLOBAL_LOCK_FD" || {
      echo "ERROR: lock busy ($GLOBAL_LOCK_FILE). Another run may be active." >&2
      exit 1
    }
  fi
  log "Acquired global lock: $GLOBAL_LOCK_FILE"
}

safe_remove_dir_for_force() {
  local d="$1"
  [[ -n "$d" && "$d" != "/" ]] || { echo "ERROR: unsafe empty/root directory removal request" >&2; exit 1; }
  [[ "$d" == "$WORKDIR"/* ]] || { echo "ERROR: refusing to remove directory outside WORKDIR: $d" >&2; exit 1; }
  rm -rf "$d"
}

THREADS=${THREADS:-10}
RAGTAG_THREADS=${RAGTAG_THREADS:-$THREADS}
FORCE=${FORCE:-0}
RUN_RAGTAG_CORRECT=${RUN_RAGTAG_CORRECT:-1}
RAGTAG_CORRECT_VALIDATE=${RAGTAG_CORRECT_VALIDATE:-1}
RUN_RAGTAG_SCAFFOLD=${RUN_RAGTAG_SCAFFOLD:-1}
RAGTAG_READTYPE=${RAGTAG_READTYPE:-sr}
RAGTAG_READ_ALIGNER=${RAGTAG_READ_ALIGNER:-minimap2}
RAGTAG_BIN=${RAGTAG_BIN:-ragtag.py}
SAMTOOLS=${SAMTOOLS:-samtools}

POLISHED_ASM=$(resolve_path "${POLISHED_ASM:-}")
REFERENCE_ASM=$(resolve_path "${REFERENCE_ASM:-}")
R1=$(resolve_path "${R1:-}")
R2=$(resolve_path "${R2:-}")

SAMPLE_TAG=${SAMPLE_TAG:-}
if [[ -z "$SAMPLE_TAG" ]]; then
  if [[ -n "${WORKDIR:-}" ]]; then
    SAMPLE_TAG=$(basename "$WORKDIR")
  else
    SAMPLE_TAG=ragtag_project
  fi
fi
WORKDIR=$(resolve_path "${WORKDIR:-results/ragtag_project/$SAMPLE_TAG}")
LOGDIR=$(resolve_path "${LOGDIR:-$WORKDIR/logs}")
CORRECT_DIR=${CORRECT_DIR:-$WORKDIR/correct}
SCAFFOLD_DIR=${SCAFFOLD_DIR:-$WORKDIR/scaffold}
GLOBAL_LOCK_WAIT_SEC=${GLOBAL_LOCK_WAIT_SEC:-0}
GLOBAL_LOCK_FILE=${GLOBAL_LOCK_FILE:-$WORKDIR/.ragtag_project.lock}
OUTPUT_ENV=${OUTPUT_ENV:-$WORKDIR/ragtag_project.outputs.env}

[[ -n "$POLISHED_ASM" ]] || { echo "ERROR: POLISHED_ASM is required" >&2; exit 1; }
[[ -n "$REFERENCE_ASM" ]] || { echo "ERROR: REFERENCE_ASM is required" >&2; exit 1; }
require_file "$POLISHED_ASM"
require_file "$REFERENCE_ASM"
require_cmd "$RAGTAG_BIN"
require_cmd "$SAMTOOLS"

if [[ "$RAGTAG_CORRECT_VALIDATE" == "1" ]]; then
  [[ -n "$R1" ]] || { echo "ERROR: R1 is required when RAGTAG_CORRECT_VALIDATE=1" >&2; exit 1; }
  require_file "$R1"
  if [[ -n "$R2" ]]; then
    require_file "$R2"
  fi
fi

CORRECT_FASTA=$CORRECT_DIR/ragtag.correct.fasta
CORRECT_AGP=$CORRECT_DIR/ragtag.correct.agp
CORRECT_BAM=$CORRECT_DIR/ragtag.correct.reads.s.bam
SCAFFOLD_FASTA=$SCAFFOLD_DIR/ragtag.scaffold.fasta
SCAFFOLD_AGP=$SCAFFOLD_DIR/ragtag.scaffold.agp

mkdir -p "$WORKDIR" "$LOGDIR"
acquire_global_lock

log "Pipeline: polished assembly -> RagTag correct -> RagTag scaffold"
log "Sample tag:       $SAMPLE_TAG"
log "Polished asm:     $POLISHED_ASM"
log "Reference asm:    $REFERENCE_ASM"
log "Workdir:          $WORKDIR"
log "Threads:          $RAGTAG_THREADS"

if [[ ! -s "$REFERENCE_ASM.fai" ]]; then
  run_with_log "$LOGDIR/faidx_reference.log" "$SAMTOOLS" faidx "$REFERENCE_ASM"
fi
if [[ ! -s "$POLISHED_ASM.fai" ]]; then
  run_with_log "$LOGDIR/faidx_polished.log" "$SAMTOOLS" faidx "$POLISHED_ASM"
fi

corrected_assembly=$POLISHED_ASM
if [[ "$RUN_RAGTAG_CORRECT" == "1" ]]; then
  if [[ "$FORCE" == "1" && -d "$CORRECT_DIR" ]]; then
    log "FORCE=1; removing correct dir: $CORRECT_DIR"
    safe_remove_dir_for_force "$CORRECT_DIR"
  fi
  mkdir -p "$CORRECT_DIR"

  correct_validate_arg=()
  if [[ "$RAGTAG_CORRECT_VALIDATE" == "1" ]]; then
    if [[ -n "$R2" && "$R2" != "$R1" ]]; then
      correct_validate_fofn=$CORRECT_DIR/ragtag_correct.reads.fofn
      printf '%s\n%s\n' "$R1" "$R2" > "$correct_validate_fofn"
      correct_validate_arg=(-F "$correct_validate_fofn" -T "$RAGTAG_READTYPE")
    else
      correct_validate_arg=(-R "$R1" -T "$RAGTAG_READTYPE")
    fi
    if [[ -n "$RAGTAG_READ_ALIGNER" ]]; then
      correct_validate_arg+=(--read-aligner "$RAGTAG_READ_ALIGNER")
    fi
  fi

  if exists_nonempty "$CORRECT_FASTA"; then
    log "RagTag correct output exists; skipping correct: $CORRECT_FASTA"
  else
    run_with_log "$LOGDIR/ragtag_correct.log" \
      "$RAGTAG_BIN" correct -o "$CORRECT_DIR" -t "$RAGTAG_THREADS" \
      "${correct_validate_arg[@]}" \
      "$REFERENCE_ASM" "$POLISHED_ASM"
  fi
  if ! exists_nonempty "$CORRECT_FASTA"; then
    echo "ERROR: expected RagTag correct output not found: $CORRECT_FASTA" >&2
    exit 1
  fi
  corrected_assembly=$CORRECT_FASTA
else
  log "RagTag correct disabled; scaffolding polished assembly directly"
fi

final_assembly=$corrected_assembly
if [[ "$RUN_RAGTAG_SCAFFOLD" == "1" ]]; then
  if [[ "$FORCE" == "1" && -d "$SCAFFOLD_DIR" ]]; then
    log "FORCE=1; removing scaffold dir: $SCAFFOLD_DIR"
    safe_remove_dir_for_force "$SCAFFOLD_DIR"
  fi
  mkdir -p "$SCAFFOLD_DIR"

  if exists_nonempty "$SCAFFOLD_FASTA"; then
    log "RagTag scaffold output exists; skipping scaffold: $SCAFFOLD_FASTA"
  else
    run_with_log "$LOGDIR/ragtag_scaffold.log" \
      "$RAGTAG_BIN" scaffold -o "$SCAFFOLD_DIR" -t "$RAGTAG_THREADS" \
      "$REFERENCE_ASM" "$corrected_assembly"
  fi
  if ! exists_nonempty "$SCAFFOLD_FASTA"; then
    echo "ERROR: expected RagTag scaffold output not found: $SCAFFOLD_FASTA" >&2
    exit 1
  fi
  final_assembly=$SCAFFOLD_FASTA
else
  log "RagTag scaffold disabled; final assembly is: $final_assembly"
fi

cat > "$OUTPUT_ENV" <<OUT
SAMPLE_TAG="$SAMPLE_TAG"
WORKDIR="$WORKDIR"
POLISHED_ASM="$POLISHED_ASM"
REFERENCE_ASM="$REFERENCE_ASM"
CORRECT_FASTA="$CORRECT_FASTA"
CORRECT_AGP="$CORRECT_AGP"
CORRECT_BAM="$CORRECT_BAM"
SCAFFOLD_FASTA="$SCAFFOLD_FASTA"
SCAFFOLD_AGP="$SCAFFOLD_AGP"
RAGTAG_CORRECT_LOG="$LOGDIR/ragtag_correct.log"
RAGTAG_SCAFFOLD_LOG="$LOGDIR/ragtag_scaffold.log"
FINAL_ASSEMBLY="$final_assembly"
OUT

log "Done. Key outputs:"
log "  Correct FASTA:   $CORRECT_FASTA"
log "  Correct BAM:     $CORRECT_BAM"
log "  Scaffold FASTA:  $SCAFFOLD_FASTA"
log "  Scaffold AGP:    $SCAFFOLD_AGP"
log "  Output env:      $OUTPUT_ENV"
log "  Final assembly:  $final_assembly"
release_global_lock
