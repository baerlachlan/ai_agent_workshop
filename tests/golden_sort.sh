#!/usr/bin/env bash
# Golden tests for `mytools sort`: diff it against real bedtools.
# Usage: ./tests/golden_sort.sh  (or via ./tests/run_golden.sh, which runs all)
#
# Only cases where bedtools and mytools are meant to agree byte-for-byte live
# here (SPEC.md section 8). Deliberate deviations -- our exit codes for usage
# errors, merge on unsorted input -- have no oracle and belong in the unit
# tests next door.
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
MYTOOLS=${MYTOOLS:-$ROOT/mytools}   # override to test a different build
DATA=$ROOT/data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

report() {
  local name=$1 got_rc=$2 want_rc=$3
  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", diffs them
check() {
  local name=$1; shift; shift        # drop the literal --
  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  report "$name" $got_rc $want_rc
}

# check_stdin <name> <file> -- <args...>
#   same, but <file> is fed to stdin of each -- both get their own copy, so
#   the first run cannot eat the input out from under the second.
check_stdin() {
  local name=$1 infile=$2; shift 2; shift
  "$MYTOOLS" "$@" < "$infile" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" < "$infile" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  report "$name" $got_rc $want_rc
}

# --- sort (issue #5) -------------------------------------------------------
# a.bed and b.bed are deliberately unsorted; row order is the answer, so the
# outputs are diffed as-is and never re-sorted first.

check "sort a.bed"                    -- sort -i "$DATA/a.bed"
check "sort b.bed"                    -- sort -i "$DATA/b.bed"
check "sort genes.bed"                -- sort -i "$DATA/genes.bed"
check_stdin "sort a.bed from stdin" "$DATA/a.bed" -- sort -i -


echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
