#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
PIPELINE="$ROOT_DIR/map_for_callable.sh"
TMPDIR_TEST=$(mktemp -d)
trap 'rm -rf "$TMPDIR_TEST"' EXIT

mkdir -p "$TMPDIR_TEST/bin" "$TMPDIR_TEST/input"

cat >"$TMPDIR_TEST/bin/samtools" <<'FAKE'
#!/usr/bin/env bash
set -euo pipefail

subcmd=${1:?missing subcommand}
shift
{
  printf '%s' "$subcmd"
  printf ' %s' "$@"
  printf '\n'
} >>"$FAKE_SAMTOOLS_LOG"

case "$subcmd" in
  sort)
    if [[ "${1:-}" == "--help" ]]; then
      printf 'Usage: samtools sort [options]\n  -@ INT\n'
      exit 1
    fi
    out=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) out=$2; shift 2 ;;
        -@) shift 2 ;;
        *) shift ;;
      esac
    done
    [[ -n "$out" ]]
    cat >"$FAKE_SORT_STDIN"
    printf 'fake bam\n' >"$out"
    ;;
  index)
    if [[ "${1:-}" == "--help" ]]; then
      printf 'Usage: samtools index [options]\n  -@ INT\n'
      exit 1
    fi
    bam="${*: -1}"
    printf 'fake index\n' >"${bam}.bai"
    ;;
  flagstat)
    printf '2 + 0 in total\n'
    ;;
  stats)
    if [[ "${1:-}" == "--help" ]]; then
      printf 'Usage: samtools stats [options]\n  -@ INT\n'
      exit 1
    fi
    printf 'SN\traw total sequences:\t2\n'
    ;;
  faidx)
    printf 'chr1\t4\t6\t4\t5\n' >"${1}.fai"
    ;;
  *)
    exit 64
    ;;
esac
FAKE

cat >"$TMPDIR_TEST/bin/bwa-mem2" <<'FAKE'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
  mem)
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    printf '@SQ\tSN:chr1\tLN:4\n'
    ;;
  index)
    exit 65
    ;;
  *)
    exit 64
    ;;
esac
FAKE
chmod 0755 "$TMPDIR_TEST/bin/samtools" "$TMPDIR_TEST/bin/bwa-mem2"

printf '>chr1\nACGT\n' >"$TMPDIR_TEST/input/assembly.fa"
printf '@r1\nA\n+\nI\n' >"$TMPDIR_TEST/input/R1.fq"
printf '@r1\nT\n+\nI\n' >"$TMPDIR_TEST/input/R2.fq"
for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
  printf 'index\n' >"$TMPDIR_TEST/input/assembly.fa${suffix}"
done

cat >"$TMPDIR_TEST/config.env" <<EOF
SAMPLE_TAG="probe_test"
FASTA="$TMPDIR_TEST/input/assembly.fa"
READ1="$TMPDIR_TEST/input/R1.fq"
READ2="$TMPDIR_TEST/input/R2.fq"
OUTDIR="$TMPDIR_TEST/output"
THREADS=3
ALIGNER="bwa-mem2"
RUN_FLAGSTAT=1
RUN_SAMTOOLS_STATS=1
RUN_DEPTH_SUMMARY=0
SAVE_DEPTH_TSV=0
SAMTOOLS="$TMPDIR_TEST/bin/samtools"
BWA_MEM2="$TMPDIR_TEST/bin/bwa-mem2"
FORCE=0
EOF

export FAKE_SAMTOOLS_LOG="$TMPDIR_TEST/samtools.log"
export FAKE_SORT_STDIN="$TMPDIR_TEST/sort.stdin"
bash "$PIPELINE" "$TMPDIR_TEST/config.env"

grep -Fx 'sort --help' "$FAKE_SAMTOOLS_LOG"
grep -Fx 'index --help' "$FAKE_SAMTOOLS_LOG"
grep -Fx 'stats --help' "$FAKE_SAMTOOLS_LOG"
if grep -Eq '^(sort|index|stats)$' "$FAKE_SAMTOOLS_LOG"; then
  echo "FAIL: samtools capability probe invoked a subcommand without --help" >&2
  cat "$FAKE_SAMTOOLS_LOG" >&2
  exit 1
fi
test -s "$TMPDIR_TEST/output/map/probe_test.sorted.bam"
test -s "$TMPDIR_TEST/output/map/probe_test.sorted.bam.bai"
test -s "$TMPDIR_TEST/output/stats/flagstat.txt"
test -s "$TMPDIR_TEST/output/stats/samtools.stats.txt"

echo "PASS: map_for_callable samtools probes use --help and the smoke workflow completes"
