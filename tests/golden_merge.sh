#!/usr/bin/env bash
# Golden tests for `mytools merge`: diff it against real bedtools.
# Usage: ./tests/golden_merge.sh
#
# One file per subcommand so that parallel work does not collide;
# tests/run_golden.sh (issue #3) is expected to call this.
#
# Every case here runs on PRE-SORTED input. bedtools refuses unsorted input to
# merge and exits 1 where we sort and exit 0 (SPEC.md section 8), so the
# unsorted case has no oracle and lives in tests/test_merge.R instead.
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
MYTOOLS=${MYTOOLS:-$HERE/../mytools}   # override to test a different build
DATA=$HERE/../data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", diffs stdout and exit codes
check() {
  local name=$1; shift; shift        # drop the literal --
  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" > "$tmp/want" 2>/dev/null
  local want_rc=$?

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

# a.bed and b.bed are deliberately unsorted -- sort them with bedtools first,
# so the fixture the oracle and we both read is in the order bedtools demands.
bedtools sort -i "$DATA/a.bed" > "$tmp/a.sorted"
bedtools sort -i "$DATA/b.bed" > "$tmp/b.sorted"
bedtools sort -i "$DATA/genes.bed" > "$tmp/genes.sorted"
bedtools sort -i "$DATA/hg002.highconf.bed" > "$tmp/highconf.sorted"

# The two cases SPEC.md section 8 names. a.sorted carries the zero-length
# features, so `merge` alone pins the leftward expansion (chr1 499 600) and
# `-d 10` pins that a lone zero-length feature out of reach stays as it is.
check "merge a.bed"                  -- merge        -i "$tmp/a.sorted"
check "merge -d 10 a.bed"            -- merge -d 10  -i "$tmp/a.sorted"

# -d 100 pulls chr2 0 0 into its neighbour, and the expansion is not clamped:
# bedtools prints a start of -1. Worth a case of its own.
check "merge -d 100 a.bed"           -- merge -d 100 -i "$tmp/a.sorted"
check "merge -d 1000 a.bed"          -- merge -d 1000 -i "$tmp/a.sorted"

# Negative -d demands that many bases of overlap, and it is the one setting
# where a nested feature opens a cluster ending before its container does
# (a06 chr1 320 350 inside a05 chr1 300 400).
check "merge -d -1 a.bed"            -- merge -d -1   -i "$tmp/a.sorted"
check "merge -d -100 a.bed"          -- merge -d -100 -i "$tmp/a.sorted"

# b.bed has its own zero-length features; genes.bed and hg002.highconf.bed are
# real coordinates in the millions, with no edge cases arranged on purpose.
check "merge b.bed"                  -- merge        -i "$tmp/b.sorted"
check "merge -d 25 b.bed"            -- merge -d 25  -i "$tmp/b.sorted"
check "merge genes.bed"              -- merge        -i "$tmp/genes.sorted"
check "merge -d 5000 genes.bed"      -- merge -d 5000 -i "$tmp/genes.sorted"
check "merge highconf.bed"           -- merge        -i "$tmp/highconf.sorted"
check "merge -d 500 highconf.bed"    -- merge -d 500 -i "$tmp/highconf.sorted"

# Empty output is a real answer: print nothing, exit 0 (SPEC.md section 5).
: > "$tmp/empty.bed"
check "merge an empty file"          -- merge -i "$tmp/empty.bed"
printf '# a comment\ntrack name=x\n\n' > "$tmp/skipped.bed"
check "merge a file of only comments" -- merge -i "$tmp/skipped.bed"

# stdin: `-i -` means the two runs take identical arguments, so each needs its
# own copy of the input -- the first must not eat it (SPEC.md section 2).
name="merge -i - (stdin)"
"$MYTOOLS" merge -i - < "$tmp/a.sorted" > "$tmp/got" 2>/dev/null
got_rc=$?
bedtools   merge -i - < "$tmp/a.sorted" > "$tmp/want" 2>/dev/null
want_rc=$?
if [[ $got_rc -eq $want_rc ]] && diff -q "$tmp/want" "$tmp/got" >/dev/null; then
  echo "ok   $name"; (( pass++ ))
else
  echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
  diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
  (( fail++ ))
fi

# Usage errors are ours, not bedtools' -- it exits 1 and prints its whole help
# text for a bad -d. SPEC.md section 7: exit 2, message on stderr, stdout clean.
usage_error() {
  local name=$1; shift
  "$MYTOOLS" "$@" > "$tmp/got" 2>"$tmp/got.err"
  local rc=$?
  if [[ $rc -eq 2 && -s "$tmp/got.err" && ! -s "$tmp/got" ]]; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name (exit $rc, stderr $(wc -c < "$tmp/got.err") bytes, stdout $(wc -c < "$tmp/got") bytes)"
    (( fail++ ))
  fi
}

usage_error "merge without -i"        merge
usage_error "merge unknown flag"      merge -q -i "$tmp/a.sorted"
usage_error "merge -d without value"  merge -i "$tmp/a.sorted" -d
usage_error "merge -d not an integer" merge -d abc -i "$tmp/a.sorted"
usage_error "merge -i missing file"   merge -i "$DATA/nope.bed"

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
