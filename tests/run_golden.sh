#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools on the files in data/.
# Usage: ./tests/run_golden.sh
#
# This is the single entry point CI runs. The cases themselves live one file per
# subcommand -- golden_sort.sh, golden_merge.sh, golden_intersect.sh,
# golden_subtract.sh -- so that work on separate subcommands does not collide.
# Each of those is also runnable on its own while you are working on one.
#
# Exits non-zero if any suite fails. Set MYTOOLS to test a different build.
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
export MYTOOLS=${MYTOOLS:-$ROOT/mytools}

if ! command -v bedtools >/dev/null 2>&1; then
  echo "bedtools not found -- it is the oracle, so these tests cannot run without it" >&2
  exit 1
fi

suites=(golden_sort.sh golden_merge.sh golden_intersect.sh golden_subtract.sh)

total_pass=0
total_fail=0
failed_suites=()

for s in "${suites[@]}"; do
  echo "=== $s ==="
  out=$(bash "$HERE/$s")
  rc=$?
  echo "$out"

  # Each suite's last line is "<n> passed, <n> failed".
  summary=$(printf '%s\n' "$out" | tail -1)
  p=$(sed -n 's/^\([0-9]\+\) passed.*/\1/p' <<<"$summary")
  f=$(sed -n 's/.*, \([0-9]\+\) failed$/\1/p' <<<"$summary")
  total_pass=$((total_pass + ${p:-0}))
  total_fail=$((total_fail + ${f:-0}))

  # A suite that dies before printing its summary still has to count as failure.
  if [[ $rc -ne 0 || -z $p ]]; then
    failed_suites+=("$s")
    [[ -z $p ]] && total_fail=$((total_fail + 1))
  fi
  echo
done

echo "========================================"
echo "TOTAL: $total_pass passed, $total_fail failed"
if (( ${#failed_suites[@]} )); then
  echo "failing suites: ${failed_suites[*]}"
fi
(( total_fail == 0 && ${#failed_suites[@]} == 0 ))
