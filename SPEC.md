# SPEC.md — mytools

A small reimplementation of a subset of bedtools, in R. Real `bedtools` is the oracle:
if our output differs from it on the same input, we are wrong — except where §8 records
a deliberate deviation.

## 1. Scope

v1 ships four subcommands:

| Subcommand  | Flags in v1         |
|-------------|---------------------|
| `sort`      | (none)              |
| `merge`     | `-d <int>`          |
| `intersect` | `-u`, `-v`, `-wa`   |
| `subtract`  | (none)              |

Explicitly **not** in v1: `closest`, `-s`/`-S` strand awareness, `-f` minimum overlap,
BED12, GFF3/VCF input, compressed (`.gz`) input, `-header`. Four subcommands finished
beats seven half-built ones; say no to all of the above.

## 2. Invocation

    mytools sort      -i <file|->
    mytools merge     -i <file|-> [-d N]
    mytools intersect -a <file|-> -b <file>
    mytools subtract  -a <file|-> -b <file>
    mytools --version

Flag names and meanings match bedtools exactly. `-` means stdin; at most one input per
invocation may be `-`. `-b` must be a real file.

## 3. Input format

BED3 through BED6, tab-separated:

    chrom  start  end  [name  score  strand]

- Column count may vary between lines; missing trailing columns are absent, not empty.
- Lines beginning with `#`, `track`, or `browser` are skipped silently.
- Blank lines are skipped.
- `start` and `end` are non-negative integers. `start > end` is an error (§7).
- `start == end` (zero-length) is **legal**. See §4.

## 4. Interval semantics

BED is **0-based, half-open**. `chr1 100 200` covers bases 100..199.

- Overlap predicate: `a.start < b.end AND b.start < a.end`. Strict `<` on both sides.
- Bookended intervals (`a.end == b.start`) do **not** overlap.
- Bookended intervals **do** merge at `-d 0` (the default), because `merge` joins
  features whose gap is `<= d`, and that gap is 0.

### Zero-length intervals — measured, not derived

Zero-length intervals are legal and bedtools handles them in ways the overlap predicate
above does not predict. These are **measured** against bedtools v2.31.1 on our fixtures;
encode them, do not reason about them.

**They overlap things the strict predicate says they should not.** All three zero-length
features in `data/a.bed` are reported as overlapping by `intersect`, and none of them
appear under `intersect -v`:

| Feature              | `-b` feature it "overlaps" | Strict predicate says |
|----------------------|----------------------------|-----------------------|
| `chr1 500 500` (a07) | `chr1 500 500` (b07)       | no overlap            |
| `chr2 0 0` (a12)     | `chr2 0 10` (b10)          | no overlap            |
| `chr2 300 300` (a16) | `chr2 200 300` (b12)       | no overlap            |

Under `intersect` and `intersect -u` all three print as their full original interval.
Under `subtract`, a12 and a16 survive; a07 is removed by the coincident b07.

**`merge` expands them leftward by one base when they touch a neighbour.** A lone
zero-length interval merges to itself unchanged, but one that meets a feature starting at
the same coordinate produces a merged interval starting at `P-1`:

    chr1 500 500 + chr1 500 600   ->  chr1 499 600      (not 500 600)
    chr1 300 300 + chr1 300 400   ->  chr1 299 400      (not 300 400)
    chr1 500 500 alone            ->  chr1 500 500      (unchanged)
    chr1   0   0 alone            ->  chr1   0   0      (no expansion below 0)

This is why `bedtools merge` on sorted `data/a.bed` prints `chr1 499 600`. Reproduce it.

**In `subtract`, a zero-length `-b` feature trims one base off each side of its
position**, including from `-a` features that merely abut it and do not overlap it by the
predicate at all. `b02` is `chr1 100 100`; it is bookended with both `a01` and `a02`:

    a01 = chr1 0 100,  minus b01 (0,50) alone      ->  chr1 50 100
    a01 = chr1 0 100,  minus b01 and b02 (100,100) ->  chr1 50  99   (right end trimmed)
    a02 = chr1 100 200, minus b02 and b03 (180,220) -> chr1 101 180   (left end trimmed)

So `bedtools subtract` on the fixtures prints `chr1 50 99 a01` and `chr1 101 180 a02`.
Both are one base narrower than the predicate in §4 predicts. Encode it; do not correct
it.

### Minimum overlap

One base. `-f` is not in v1.

## 5. Output

- Tab-separated, LF line endings, trailing newline on the final line.
- Empty result: print nothing, exit 0.
- `sort`: all input columns preserved. Order is chrom lexicographic (`chr17` before
  `chr7`), then `start`, then `end`.
- `merge`: BED3 only (`chrom start end`). Input columns are dropped.
- `intersect`: default prints the intersected region carrying `-a`'s trailing columns.
  `-wa` prints `-a`'s original interval, once per overlapping `-b` feature. `-u` prints
  each `-a` feature at most once. `-v` prints `-a` features with no overlap. `-u`, `-v`
  and `-wa` are mutually exclusive.
- `subtract`: `-a` features with `-b` regions removed. One feature may become two, or
  vanish entirely.
- Input order is preserved for `intersect` and `subtract`. bedtools does not sort `-a`
  for you, and neither do we.

## 6. Memory model

Vectorised throughout: read each input fully into a data frame and operate on columns.
R's per-row loops are slow enough that streaming would be the pessimisation, not the
optimisation.

| Subcommand  | Holds in memory        |
|-------------|------------------------|
| `sort`      | whole input            |
| `merge`     | whole input (see §4 and §8) |
| `intersect` | both `-a` and `-b`     |
| `subtract`  | both `-a` and `-b`     |

- Target: inputs up to ~10^6 intervals. `bedtools bamtobed` on the workshop BAM yields
  ~500,000, which is the intended stress case.
- No mmap, no index files, no parallelism, no third-party runtime dependencies.
- Being slower than bedtools is acceptable — it is C and we are not. Being *quadratic*
  is a bug.

## 7. Errors and exit codes

Errors go to **stderr**. stdout carries data only, so it can be piped.

| Situation                           | stderr message                          | exit |
|-------------------------------------|-----------------------------------------|------|
| Success (including empty output)    | —                                       | 0    |
| `--version`, `--help`               | —                                       | 0    |
| Malformed line / non-integer coords | `a.bed:14: malformed BED line`          | 1    |
| `start > end`                       | `a.bed:14: start > end (500 > 400)`     | 1    |
| Negative coordinate                 | `a.bed:14: negative coordinate (-1)`    | 1    |
| Unknown flag                        | `unknown flag: -q`                      | 2    |
| Missing required argument           | `merge: -i is required`                 | 2    |
| Mutually exclusive flags            | `intersect: -u, -v and -wa are exclusive` | 2  |
| Input file does not exist           | `no such file: a.bed`                   | 2    |
| Two inputs given as `-`             | `at most one input may be stdin`        | 2    |

Data problems are 1. Caller problems are 2. Messages name the file and line number.

## 8. Correctness

Real `bedtools` is the oracle. Every case below must produce byte-identical stdout **and
the same exit code** as its bedtools equivalent:

    mytools sort -i data/a.bed                     == bedtools sort -i data/a.bed
    mytools merge -i <sorted a.bed>                == bedtools merge -i <sorted a.bed>
    mytools merge -d 10 -i <sorted a.bed>          == bedtools merge -d 10 -i <sorted a.bed>
    mytools intersect -a data/a.bed -b data/b.bed  == bedtools intersect -a ... -b ...
    mytools intersect -u  ...                      == bedtools intersect -u ...
    mytools intersect -v  ...                      == bedtools intersect -v ...
    mytools intersect -wa ...                      == bedtools intersect -wa ...
    mytools subtract -a data/a.bed -b data/b.bed   == bedtools subtract -a ... -b ...

Plus two stdin cases (`sort -i -` and `intersect -a -`) and `mytools --version`
printing a version and exiting 0. Ten cases in total.

### Accepted deviations from bedtools

**`merge` sorts unsorted input instead of erroring.** bedtools refuses unsorted input to
`merge` and exits 1; we sort silently and exit 0. Chosen for a more forgiving CLI.

Consequences, both deliberate:

- `merge` cannot stream, and holds its whole input in memory (§6).
- The unsorted case has **no oracle** — bedtools fails where we succeed, so it cannot be
  a golden test. Golden tests for `merge` therefore run on pre-sorted input, where the
  two agree exactly. Cover the deviation with a unit test asserting that
  `merge` on unsorted input equals `merge` on the same input pre-sorted.

No other deviations. If you find one you cannot fix, write it down here with the reason.

## 9. Language and layout

R, per `CLAUDE.md`. Every subcommand and every test.

    mytools           dispatcher, --version, usage
    R/bed.R           shared: parse, write, overlap predicate
    R/sort.R          cmd_sort()
    R/merge.R         cmd_merge()
    R/intersect.R     cmd_intersect()
    R/subtract.R      cmd_subtract()
    tests/            golden + unit tests, driven by tests/run_golden.sh

- `mytools` is executable (`#!/usr/bin/env Rscript`) and on `PATH`, or invoked through a
  wrapper the tests set.
- One file per subcommand so that parallel work does not collide. `R/bed.R` is the one
  shared file — write it first, and agree on it before fanning out.
- Not an R package: no `DESCRIPTION`, no `NAMESPACE`, no install step before tests.
- Base R only. No third-party runtime dependencies.
