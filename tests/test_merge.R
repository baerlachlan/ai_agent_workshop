#!/usr/bin/env Rscript
# Unit tests for R/merge.R. No bedtools required; runs in seconds.
# Usage: Rscript tests/test_merge.R
#
# Golden tests (tests/golden_merge.sh) prove we agree with bedtools. These pin
# the cases that needed thinking about -- the zero-length rules of SPEC.md
# section 4, tie order, negative -d -- so a refactor cannot quietly undo them.
# Every expected value here was measured against bedtools v2.31.1, except the
# unsorted-input case, which by construction has no oracle (SPEC.md section 8).

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
HERE <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
ROOT <- normalizePath(file.path(HERE, ".."))
DATA <- file.path(ROOT, "data")

source(file.path(ROOT, "R", "bed.R"))
source(file.path(ROOT, "R", "merge.R"))

pass <- 0
fail <- 0

ok <- function(name, cond) {
  if (isTRUE(cond)) {
    cat("ok   ", name, "\n", sep = "")
    pass <<- pass + 1
  } else {
    cat("FAIL ", name, "\n", sep = "")
    fail <<- fail + 1
  }
}

# Takes and returns whitespace-separated BED3 rows, so a case reads as the
# intervals going in and the intervals coming out: two rows 100-200 and 150-250
# merge, at the default d, to the single row 100-250.
merged <- function(..., d = 0L) {
  rows <- c(...)
  fields <- strsplit(trimws(rows), "[ \t]+")
  out <- merge_intervals(vapply(fields, `[`, "", 1L),
                         as.integer(vapply(fields, `[`, "", 2L)),
                         as.integer(vapply(fields, `[`, "", 3L)),
                         d)
  if (nrow(out) == 0L) return(character())
  paste(out$chrom, out$start, out$end)
}

# --- output shape ----------------------------------------------------------
# SPEC.md section 5: merge prints BED3 only, input columns are dropped.

ok("output is BED3, input columns dropped",
   identical(names(merge_intervals("chr1", 100L, 200L, 0L)),
             c("chrom", "start", "end")))
ok("empty input gives no rows",
   identical(nrow(merge_intervals(character(), integer(), integer(), 0L)), 0L))
empty_merge <- merge_intervals(character(), integer(), integer(), 0L)
ok("write_bed of an empty merge prints nothing",
   identical(capture.output(write_bed(empty_merge)), character()))

# --- the gap rule ----------------------------------------------------------
# Features merge when the gap between them is <= d. Default d is 0.

ok("overlapping features merge",
   identical(merged("chr1 100 200", "chr1 150 250"), "chr1 100 250"))
ok("disjoint features do not merge",
   identical(merged("chr1 100 200", "chr1 300 400"),
             c("chr1 100 200", "chr1 300 400")))
ok("a gap of 1 does not merge at d = 0",
   identical(merged("chr1 100 200", "chr1 201 300"),
             c("chr1 100 200", "chr1 201 300")))
ok("a gap of 1 merges at d = 1",
   identical(merged("chr1 100 200", "chr1 201 300", d = 1L), "chr1 100 300"))
ok("a gap of 2 does not merge at d = 1",
   identical(merged("chr1 100 200", "chr1 202 300", d = 1L),
             c("chr1 100 200", "chr1 202 300")))

# Bookended: a.end == b.start. They do NOT overlap by the SPEC.md section 4
# predicate, and they DO merge here, because their gap is 0. Both are true.
ok("bookended features merge at d = 0",
   identical(merged("chr1 0 100", "chr1 100 200"), "chr1 0 200"))
ok("bookended features do not merge at d = -1",
   identical(merged("chr1 0 100", "chr1 100 200", d = -1L),
             c("chr1 0 100", "chr1 100 200")))
ok("the predicate and the merge rule disagree on bookended, as designed",
   !overlaps(0L, 100L, 100L, 200L) &&
     identical(merged("chr1 0 100", "chr1 100 200"), "chr1 0 200"))

# Nested and identical features.
ok("a nested feature does not extend its container",
   identical(merged("chr1 300 400", "chr1 320 350"), "chr1 300 400"))
ok("identical features collapse to one",
   identical(merged("chr1 700 800", "chr1 700 800"), "chr1 700 800"))

# Position 0 is an ordinary coordinate, not a sentinel.
ok("a feature at position 0 merges normally",
   identical(merged("chr1 0 100", "chr1 50 150"), "chr1 0 150"))

# Chromosomes are independent, whatever the coordinates say.
ok("features on different chromosomes never merge",
   identical(merged("chr1 100 200", "chr2 100 200"),
             c("chr1 100 200", "chr2 100 200")))
ok("chromosomes come out in lexicographic order, chr17 before chr7",
   identical(merged("chr7 100 200", "chr17 100 200"),
             c("chr17 100 200", "chr7 100 200")))

# --- zero-length features (SPEC.md section 4) ------------------------------
# bedtools holds a zero-length feature [p, p) as [p - 1, p + 1) and clusters on
# that, printing the original coordinates only when nothing merged with it.
# Measured, not derived -- do not "correct" these.

ok("a lone zero-length feature is unchanged",
   identical(merged("chr1 500 500"), "chr1 500 500"))
ok("a lone zero-length feature at 0 is unchanged, no expansion below 0",
   identical(merged("chr1 0 0"), "chr1 0 0"))
ok("a zero-length feature expands leftward when a neighbour starts at it",
   identical(merged("chr1 500 500", "chr1 500 600"), "chr1 499 600"))
ok("the same, at another coordinate",
   identical(merged("chr1 300 300", "chr1 300 400"), "chr1 299 400"))
ok("the expansion is not clamped at 0",
   identical(merged("chr1 0 0", "chr1 0 10"), "chr1 -1 10"))
ok("a zero-length feature also reaches one base right: gap to p + 1 is 0",
   identical(merged("chr1 500 500", "chr1 501 600"), "chr1 499 600"))
ok("its right edge shows up in the output too",
   identical(merged("chr1 400 499", "chr1 500 500"), "chr1 400 501"))
ok("two coincident zero-length features expand both ways",
   identical(merged("chr1 500 500", "chr1 500 500"), "chr1 499 501"))
ok("a zero-length feature inside an interval changes nothing",
   identical(merged("chr1 100 200", "chr1 150 150"), "chr1 100 200"))
ok("a lone zero-length feature stays lone when -d is too small to reach",
   identical(merged("chr1 0 0", "chr1 50 150", d = 10L),
             c("chr1 0 0", "chr1 50 150")))
ok("and expands once -d reaches its neighbour",
   identical(merged("chr1 0 0", "chr1 50 150", d = 100L), "chr1 -1 150"))

# --- tie order is part of the answer ---------------------------------------
# bedtools sort is stable: features sharing a chrom and start keep their input
# order. merge then takes the FIRST feature's start, so the two orders below
# give different output. This is why R/merge.R sorts on (chrom, start) only --
# adding `end` to the key silently reorders these and changes the result.

ok("tie order kept: the wider feature first hides the expansion",
   identical(merged("chr1 54 74", "chr1 54 54"), "chr1 54 74"))
ok("tie order kept: the zero-length feature first exposes it",
   identical(merged("chr1 54 54", "chr1 54 74"), "chr1 53 74"))

# --- negative -d -----------------------------------------------------------
# Negative d demands that many bases of overlap. A nested feature can then open
# a cluster that ends before its container does, so the cluster end has to be
# tracked per cluster and not as a running maximum over the chromosome.

ok("negative d requires overlap: 50 bases is not enough for -d -100",
   identical(merged("chr1 100 200", "chr1 150 250", d = -100L),
             c("chr1 100 200", "chr1 150 250")))
ok("negative d: 100 bases of overlap is enough for -d -100",
   identical(merged("chr1 100 200", "chr1 100 250", d = -100L), "chr1 100 250"))
ok("negative d: a nested feature keeps its own end, not its container's",
   identical(merged("chr1 300 400", "chr1 320 350", d = -100L),
             c("chr1 300 400", "chr1 320 350")))
ok("negative d: the feature after a nested one is compared to the right end",
   identical(merged("chr1 38 43", "chr1 43 63", "chr1 43 44", d = -3L),
             c("chr1 38 43", "chr1 43 63")))

# --- unsorted input, the deliberate deviation (SPEC.md section 8) ----------
# bedtools refuses unsorted input to merge and exits 1; we sort and exit 0, so
# this case has no oracle. data/a.bed is deliberately unsorted. The expected
# rows are bedtools' output for `bedtools sort -i a.bed | bedtools merge`,
# which is the answer our sorted-input golden test already pins.

a <- read_bed(file.path(DATA, "a.bed"))
a_merged <- merge_intervals(a$chrom, a$start, a$end, 0L)
want_sorted <- c("chr1 0 250", "chr1 300 400", "chr1 499 600", "chr1 700 800",
                 "chr1 900 1000", "chr2 0 0", "chr2 50 150", "chr2 200 400",
                 "chr2 1000 2000", "chrX 10 20", "chrX 30 40", "chrX 1000 1100")
ok("merge on unsorted a.bed equals merge on the same input pre-sorted",
   identical(paste(a_merged$chrom, a_merged$start, a_merged$end), want_sorted))

# The sort has to be stable to be equivalent to pre-sorting: a.bed holds
# `chr1 500 500` before `chr1 500 600`, and that order is what produces the
# `chr1 499 600` above rather than `chr1 500 600`.
ok("a.bed's tie at start 500 is not reordered",
   identical(merged("chr1 500 500", "chr1 500 600"), "chr1 499 600"))

# --- CLI: flags and exit codes (SPEC.md section 7) -------------------------
# Run in a subprocess, since a usage error calls quit(). Data problems are 1,
# caller problems are 2, and stdout carries data only.

run_mytools <- function(...) {
  out <- tempfile()
  err <- tempfile()
  st <- system2(file.path(ROOT, "mytools"), c("merge", ...),
                stdout = out, stderr = err)
  list(status = st,
       stdout = readLines(out, warn = FALSE),
       stderr = paste(readLines(err, warn = FALSE), collapse = "\n"))
}

r <- run_mytools("-i", file.path(DATA, "a.bed"))
ok("CLI: merge of a.bed exits 0", identical(r$status, 0L))
ok("CLI: merge of a.bed prints the merged rows",
   identical(r$stdout, sub(" ", "\t", sub(" ", "\t", want_sorted))))
ok("CLI: merge of a.bed says nothing on stderr", identical(r$stderr, ""))

r <- run_mytools("-d", "10", "-i", file.path(DATA, "a.bed"))
ok("CLI: -d 10 exits 0", identical(r$status, 0L))

r <- run_mytools()
ok("CLI: no -i is a usage error", identical(r$status, 2L))
ok("CLI: no -i names the missing flag on stderr",
   grepl("-i is required", r$stderr, fixed = TRUE))
ok("CLI: no -i keeps stdout clean", identical(r$stdout, character()))

r <- run_mytools("-q", "-i", file.path(DATA, "a.bed"))
ok("CLI: unknown flag is a usage error", identical(r$status, 2L))
ok("CLI: unknown flag names itself",
   grepl("unknown flag: -q", r$stderr, fixed = TRUE))

r <- run_mytools("-d", "abc", "-i", file.path(DATA, "a.bed"))
ok("CLI: non-integer -d is a usage error", identical(r$status, 2L))

r <- run_mytools("-i", file.path(DATA, "a.bed"), "-d")
ok("CLI: -d with no value is a usage error", identical(r$status, 2L))

r <- run_mytools("-i", file.path(DATA, "nope.bed"))
ok("CLI: a missing input file is a usage error", identical(r$status, 2L))
ok("CLI: a missing input file is named",
   grepl("no such file", r$stderr, fixed = TRUE))

cat("---\n")
cat(pass, " passed, ", fail, " failed\n", sep = "")
quit(save = "no", status = if (fail > 0) 1 else 0)
