#!/usr/bin/env bash
# Unit tests: every tests/test_*.R file. No bedtools required, runs in seconds.
# Usage: ./tests/run_unit.sh
#
# These pin the behaviour that has no oracle -- our exit codes, merge on
# unsorted input (SPEC.md section 8) -- and the edge cases we had to reason
# about, so a failure names the function rather than just the disagreement.
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

total_pass=0
total_fail=0
failed=()

for t in "$HERE"/test_*.R; do
  name=$(basename "$t")
  out=$(Rscript "$t" 2>&1)
  rc=$?
  summary=$(printf '%s\n' "$out" | tail -1)
  p=$(sed -n 's/^\([0-9]\+\) passed.*/\1/p' <<<"$summary")
  f=$(sed -n 's/.*, \([0-9]\+\) failed$/\1/p' <<<"$summary")

  if [[ $rc -ne 0 || -z $p ]]; then
    printf '%-24s FAIL  %s\n' "$name" "$summary"
    printf '%s\n' "$out" | grep '^FAIL' | sed 's/^/    /'
    failed+=("$name")
    [[ -z $p ]] && total_fail=$((total_fail + 1))
  else
    printf '%-24s ok    %s\n' "$name" "$summary"
  fi
  total_pass=$((total_pass + ${p:-0}))
  total_fail=$((total_fail + ${f:-0}))
done

echo "========================================"
echo "TOTAL: $total_pass passed, $total_fail failed"
if (( ${#failed[@]} )); then
  echo "failing files: ${failed[*]}"
fi
(( total_fail == 0 && ${#failed[@]} == 0 ))
