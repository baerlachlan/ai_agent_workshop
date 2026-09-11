#!/usr/bin/env bash
# Golden tests for `mytools intersect`: diff it against real bedtools.
# Usage: ./tests/golden_intersect.sh
#
# One file per subcommand so that parallel work does not collide; tests/run_golden.sh
# (issue #3) is expected to call this.
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

# The four output modes (SPEC.md s5).
check "intersect"                    -- intersect     -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -u"                 -- intersect -u  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -v"                 -- intersect -v  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -wa"                -- intersect -wa -a "$DATA/a.bed" -b "$DATA/b.bed"

# Roles swapped, so the zero-length features sit in -b as well as in -a. bedtools
# refuses a zero-length -b feature at position 0 (a12 is chr2 0 0) and exits 1 with
# empty stdout, because widening it to index it gives a start of -1. We match that,
# so this case is a real comparison and not just an empty diff -- it checks the exit
# code as much as the output.
check "intersect b-as-a"             -- intersect     -a "$DATA/b.bed" -b "$DATA/a.bed"
check "intersect b-as-a -v"          -- intersect -v  -a "$DATA/b.bed" -b "$DATA/a.bed"

# Real coordinates, in the millions: hits then span several bin levels, and bedtools
# reports the hits for one -a feature finest-bin-first rather than in -b file order.
# a.bed and b.bed cannot catch that -- every coordinate in them lands in one bin.
check "intersect genes x highconf"   -- intersect -a "$DATA/genes.bed" -b "$DATA/hg002.highconf.bed"
check "intersect highconf x genes"   -- intersect -a "$DATA/hg002.highconf.bed" -b "$DATA/genes.bed"

# stdin: -a is "-", so the two runs do not take identical arguments and check() cannot
# express it (SPEC.md s2, and bedtools reads -a the same way).
name="intersect -a - (stdin)"
"$MYTOOLS" intersect -a - -b "$DATA/b.bed" < "$DATA/a.bed" > "$tmp/got" 2>/dev/null
got_rc=$?
bedtools intersect -a - -b "$DATA/b.bed" < "$DATA/a.bed" > "$tmp/want" 2>/dev/null
want_rc=$?
if [[ $got_rc -eq $want_rc ]] && diff -q "$tmp/want" "$tmp/got" >/dev/null; then
  echo "ok   $name"; (( pass++ ))
else
  echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
  diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
  (( fail++ ))
fi

# Usage errors are ours, not bedtools' -- its wording and exit codes for these are its
# own business (SPEC.md s7: exit 2, message on stderr, stdout stays clean).
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

usage_error "intersect -u -v is exit 2"  intersect -u -v -a "$DATA/a.bed" -b "$DATA/b.bed"
usage_error "intersect -u -wa is exit 2" intersect -u -wa -a "$DATA/a.bed" -b "$DATA/b.bed"
usage_error "intersect unknown flag"     intersect -q -a "$DATA/a.bed" -b "$DATA/b.bed"
usage_error "intersect without -b"       intersect -a "$DATA/a.bed"
usage_error "intersect without -a"       intersect -b "$DATA/b.bed"
usage_error "intersect two stdins"       intersect -a - -b -

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
