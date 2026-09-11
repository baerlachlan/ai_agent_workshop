#!/usr/bin/env bash
# Golden tests for "mytools subtract": diff it against real bedtools.
# Usage: ./tests/golden_subtract.sh
#
# Kept in its own file so that the four subcommands can be built in parallel
# without four agents editing one script. tests/run_golden.sh (issue #3) should
# call this, or fold these cases in once it exists.
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
MYTOOLS=${MYTOOLS:-$ROOT/mytools}     # override to test a different build
DATA=$ROOT/data

if ! command -v bedtools >/dev/null; then
  echo "SKIP subtract golden tests: bedtools is not installed" >&2
  exit 0
fi

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", diffs stdout and exit code
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

# The SPEC.md section 8 case. a.bed and b.bed are deliberately unsorted and
# carry every zero-length and bookended edge case; subtract is the subcommand
# they hit hardest.
check "subtract a.bed - b.bed" -- subtract -a "$DATA/a.bed" -b "$DATA/b.bed"

# a.bed as -b is the refusal case: widening a12 (chr2 0 0) gives it a start of
# -1, which bedtools' bin index rejects, so it prints nothing and exits 1. Both
# sides of this check must agree on that -- check() compares the exit code as
# well as stdout, and stderr is not part of the contract (tests/README.md).
check "subtract b.bed - a.bed (both refuse: a12 at 0)" \
  -- subtract -a "$DATA/b.bed" -b "$DATA/a.bed"

# The same pair with a12 dropped, so that -a still carries b.bed's zero-length
# features (b02, b07) through a run that produces real output.
grep -v $'^chr2\t0\t0\t' "$DATA/a.bed" > "$tmp/a-no-zero-at-0.bed"
check "subtract b.bed - a.bed (minus a12)" \
  -- subtract -a "$DATA/b.bed" -b "$tmp/a-no-zero-at-0.bed"

# Real data, no arranged edge cases: 25 MANE genes against GIAB's
# high-confidence regions, which have real gaps.
check "subtract genes.bed - hg002.highconf.bed" \
  -- subtract -a "$DATA/genes.bed" -b "$DATA/hg002.highconf.bed"
check "subtract hg002.highconf.bed - genes.bed" \
  -- subtract -a "$DATA/hg002.highconf.bed" -b "$DATA/genes.bed"

# -a from stdin. check() cannot do this: both processes would race for the
# same stdin, so run them one at a time.
"$MYTOOLS" subtract -a - -b "$DATA/b.bed" < "$DATA/a.bed" > "$tmp/got" 2>/dev/null
got_rc=$?
bedtools subtract -a - -b "$DATA/b.bed" < "$DATA/a.bed" > "$tmp/want" 2>/dev/null
want_rc=$?
if [[ $got_rc -eq $want_rc ]] && diff -q "$tmp/want" "$tmp/got" >/dev/null; then
  echo "ok   subtract -a - (stdin)"; (( pass++ ))
else
  echo "FAIL subtract -a - (stdin) (exit $got_rc, bedtools gave $want_rc)"
  diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
  (( fail++ ))
fi

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
