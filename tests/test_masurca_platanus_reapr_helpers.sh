#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
PIPELINE="$ROOT_DIR/masurca_platanus_reapr.sh"
TMPDIR_TEST=$(mktemp -d)
trap 'rm -rf "$TMPDIR_TEST"' EXIT

mkdir -p "$TMPDIR_TEST/bin"
cat >"$TMPDIR_TEST/bin/samtools" <<'FAKE'
#!/usr/bin/env bash
set -euo pipefail

subcmd=${1:?missing subcommand}
shift
case "$subcmd" in
  sort)
    if [[ "${1:-}" == "--help" ]]; then
      printf 'Usage: samtools sort [options] [in.bam]\n  -o FILE\n  -@ INT\n'
      exit 1
    fi
    printf 'sort\t%s\n' "$*" >>"$FAKE_SAMTOOLS_LOG"
    if [[ "${FAKE_SORT_FAIL:-0}" -ne 0 ]]; then
      exit "$FAKE_SORT_FAIL"
    fi
    out=""
    input=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) out=$2; shift 2 ;;
        -T|-@) shift 2 ;;
        *) input=$1; shift ;;
      esac
    done
    [[ -n "$out" && -n "$input" ]]
    cp "$input" "$out"
    ;;
  index)
    if [[ "${1:-}" == "--help" ]]; then
      printf 'Usage: samtools index [options] in.bam\n  -@ INT\n'
      exit 1
    fi
    printf 'index\t%s\n' "$*" >>"$FAKE_SAMTOOLS_LOG"
    bam=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -@) shift 2 ;;
        *) bam=$1; shift ;;
      esac
    done
    [[ -n "$bam" ]]
    printf 'fake index\n' >"${bam}.bai"
    ;;
  *)
    exit 64
    ;;
esac
FAKE
chmod 0755 "$TMPDIR_TEST/bin/samtools"

helpers=$(awk '
  /^log\(\)/ { emit=1 }
  /^require_reapr_aligner\(\)/ { emit=0 }
  emit { print }
' "$PIPELINE")
eval "$helpers"

export PATH="$TMPDIR_TEST/bin:$PATH"
export FAKE_SAMTOOLS_LOG="$TMPDIR_TEST/samtools.log"
export SAMTOOLS_THREADS=3

bwa2_prefix="$TMPDIR_TEST/reapr_bwa2"
for suffix in .0123 .amb .ann .pac .bwt.2bit.64; do
  printf 'index data\n' >"${bwa2_prefix}${suffix}"
done
bwa_mem2_index_ready "$bwa2_prefix"
find "${bwa2_prefix}.pac" -maxdepth 0 -type f -delete
if bwa_mem2_index_ready "$bwa2_prefix"; then
  echo "FAIL: incomplete bwa-mem2 index accepted" >&2
  exit 1
fi
printf 'index data\n' >"${bwa2_prefix}.pac"
find "${bwa2_prefix}.bwt.2bit.64" -maxdepth 0 -type f -delete
printf 'index data\n' >"${bwa2_prefix}.bwt.2bit"
bwa_mem2_index_ready "$bwa2_prefix"

printf 'fake bam\n' >"$TMPDIR_TEST/input.bam"
samtools_sort_index "$TMPDIR_TEST/input.bam" "$TMPDIR_TEST/output.bam" auto
cmp "$TMPDIR_TEST/input.bam" "$TMPDIR_TEST/output.bam"
test -s "$TMPDIR_TEST/output.bam.bai"
grep -q $'^sort\t-@ 3 .* -o ' "$FAKE_SAMTOOLS_LOG"
grep -q $'^index\t-@ 3 ' "$FAKE_SAMTOOLS_LOG"

: >"$FAKE_SAMTOOLS_LOG"
set +e
FAKE_SORT_FAIL=42 samtools_sort_index \
  "$TMPDIR_TEST/input.bam" "$TMPDIR_TEST/failed.bam" auto
rc=$?
set -e
test "$rc" -eq 42
test ! -e "$TMPDIR_TEST/failed.bam"
test ! -e "$TMPDIR_TEST/failed.bam.bai"
if grep -q '^index' "$FAKE_SAMTOOLS_LOG"; then
  echo "FAIL: index ran after sort failure" >&2
  exit 1
fi

echo "PASS: samtools sort/index handling and bwa-mem2 index detection"
