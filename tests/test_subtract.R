#!/usr/bin/env Rscript
# Unit tests for R/subtract.R. No bedtools required; runs in seconds.
# Usage: Rscript tests/test_subtract.R
#
# Every expectation below was measured against bedtools v2.31.1 first. Where an
# expectation contradicts the SPEC.md section 4 overlap predicate, that is the
# point of the test -- see the comments.

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
HERE <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
ROOT <- normalizePath(file.path(HERE, ".."))

source(file.path(ROOT, "R", "bed.R"))
source(file.path(ROOT, "R", "subtract.R"))

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

# Runs cmd_subtract on two inline fixtures and returns its stdout as lines.
sub_out <- function(a, b) {
  fa <- tempfile(); writeLines(a, fa)
  fb <- tempfile(); writeLines(b, fb)
  capture.output(cmd_subtract(c("-a", fa, "-b", fb)))
}

same <- function(name, a, b, want) ok(name, identical(sub_out(a, b), want))

# Runs a whole invocation in a subprocess, since usage errors call quit().
run_subprocess <- function(argv) {
  expr <- sprintf('source("%s"); source("%s"); cmd_subtract(c(%s))',
                  file.path(ROOT, "R", "bed.R"),
                  file.path(ROOT, "R", "subtract.R"),
                  paste(sprintf('"%s"', argv), collapse = ", "))
  err <- tempfile()
  st <- system2("Rscript", c("-e", shQuote(expr)), stdout = NULL, stderr = err)
  list(status = st, stderr = paste(readLines(err, warn = FALSE), collapse = "\n"))
}

# --- the ordinary cases ----------------------------------------------------

same("trim from the left",
     "chr1\t100\t200\ta", "chr1\t50\t150\tb",
     "chr1\t150\t200\ta")

same("trim from the right",
     "chr1\t100\t200\ta", "chr1\t150\t250\tb",
     "chr1\t100\t150\ta")

# A -b feature landing in the middle splits one feature into two, and the
# trailing columns are carried onto both fragments (SPEC.md section 5).
same("split in two carries the trailing columns",
     "chr1\t100\t200\ta\t50\t+", "chr1\t140\t160\tb",
     c("chr1\t100\t140\ta\t50\t+", "chr1\t160\t200\ta\t50\t+"))

same("fully covered vanishes",
     "chr1\t100\t200\ta", "chr1\t100\t200\tb",
     character(0))

same("nested -b splits, nested -a survives whole",
     c("chr1\t100\t200\ta1", "chr1\t120\t130\ta2"), "chr1\t300\t400\tb",
     c("chr1\t100\t200\ta1", "chr1\t120\t130\ta2"))

# Overlapping -b features are one region, not two subtractions in sequence.
same("overlapping -b features merge into one removed region",
     "chr1\t300\t400\ta", c("chr1\t320\t350\tb1", "chr1\t340\t360\tb2"),
     c("chr1\t300\t320\ta", "chr1\t360\t400\ta"))

same("empty -b leaves -a untouched",
     "chr1\t100\t200\ta", character(0),
     "chr1\t100\t200\ta")

same("-b on another chromosome is ignored",
     "chr1\t100\t200\ta", "chr2\t100\t200\tb",
     "chr1\t100\t200\ta")

# SPEC.md section 5: bedtools does not sort -a for you, and neither do we.
same("input order of -a is preserved",
     c("chr2\t50\t60\ta1", "chr1\t10\t20\ta2", "chr2\t10\t20\ta3"), character(0),
     c("chr2\t50\t60\ta1", "chr1\t10\t20\ta2", "chr2\t10\t20\ta3"))

# --- position 0 ------------------------------------------------------------

same("feature at position 0 trims correctly",
     "chr1\t0\t100\ta", "chr1\t0\t50\tb",
     "chr1\t50\t100\ta")

# No test for a zero-length -b at position 0: it widens to start -1, where
# bedtools aborts ("illegal bin number -1") and exits 1. There is no oracle to
# encode, so nothing is asserted about it.

# --- bookended: no overlap, and no trim ------------------------------------
# SPEC.md section 4: a.end == b.start is not an overlap. A non-zero-length -b
# that merely abuts a feature leaves it alone.

same("bookended -b on the right does not trim",
     "chr1\t100\t200\ta", "chr1\t200\t300\tb",
     "chr1\t100\t200\ta")

same("bookended -b on the left does not trim",
     "chr1\t100\t200\ta", "chr1\t0\t100\tb",
     "chr1\t100\t200\ta")

# --- zero-length -b: the oracle quirk --------------------------------------
# MEASURED, not derived. bedtools widens a zero-length record to
# [P - 1, P + 1) before searching, so it removes bases P-1 and P -- including
# from features it only abuts, which the section 4 predicate says it misses.
# These three are the exact cases in SPEC.md section 4 and issue #8.

same("zero-length -b bookended on the right trims one base (a01)",
     "chr1\t0\t100\ta01\t10\t+",
     c("chr1\t0\t50\tb01", "chr1\t100\t100\tb02"),
     "chr1\t50\t99\ta01\t10\t+")

same("zero-length -b bookended on the left trims one base (a02)",
     "chr1\t100\t200\ta02\t20\t-",
     c("chr1\t100\t100\tb02", "chr1\t180\t220\tb03"),
     "chr1\t101\t180\ta02\t20\t-")

same("zero-length -b alone, no widening, leaves the predicate answer",
     "chr1\t0\t100\ta01\t10\t+", "chr1\t0\t50\tb01",
     "chr1\t50\t100\ta01\t10\t+")

same("zero-length -b mid-feature removes two bases, not zero",
     "chr1\t100\t200\ta", "chr1\t150\t150\tb",
     c("chr1\t100\t149\ta", "chr1\t151\t200\ta"))

same("zero-length -b two bases clear of the feature does nothing",
     "chr1\t100\t200\ta", "chr1\t201\t201\tb",
     "chr1\t100\t200\ta")

same("zero-length -b just past the end still trims (widened reach)",
     "chr1\t100\t200\ta", "chr1\t200\t200\tb",
     "chr1\t100\t199\ta")

# --- zero-length -a --------------------------------------------------------
# A zero-length -a feature widens too, and if anything survives it prints once,
# back at its original zero-length coordinates.

same("zero-length -a is removed by a coincident zero-length -b (a07)",
     "chr1\t500\t500\ta07\t0\t+", "chr1\t500\t500\tb07\t0\t-",
     character(0))

same("zero-length -a survives a -b that starts at the same base (a12)",
     "chr2\t0\t0\ta12\t0\t+", "chr2\t0\t10\tb10\t18\t-",
     "chr2\t0\t0\ta12\t0\t+")

same("zero-length -a survives a -b bookended on its left (a16)",
     "chr2\t300\t300\ta16\t0\t+", "chr2\t200\t300\tb12\t20\t+",
     "chr2\t300\t300\ta16\t0\t+")

same("zero-length -a is removed by a -b that covers its widened span",
     "chr1\t100\t100\ta", "chr1\t99\t200\tb",
     character(0))

same("zero-length -a at position 0 survives",
     "chr1\t0\t0\ta", "chr1\t0\t1\tb",
     "chr1\t0\t0\ta")

# --- usage errors ----------------------------------------------------------
# SPEC.md section 7: caller problems exit 2, data problems exit 1.

fa <- tempfile(); writeLines("chr1\t100\t200\ta", fa)
fb <- tempfile(); writeLines("chr1\t100\t200\tb", fb)

r <- run_subprocess(c("-a", fa))
ok("error: missing -b exits 2", r$status == 2)
ok("error: missing -b is reported", grepl("-b is required", r$stderr))

r <- run_subprocess(c("-b", fb))
ok("error: missing -a exits 2", r$status == 2)
ok("error: missing -a is reported", grepl("-a is required", r$stderr))

r <- run_subprocess(c("-a", fa, "-b", fb, "-q"))
ok("error: unknown flag exits 2", r$status == 2)
ok("error: unknown flag is reported", grepl("unknown flag: -q", r$stderr))

r <- run_subprocess(c("-a", "-", "-b", "-"))
ok("error: two stdin inputs exits 2", r$status == 2)
ok("error: two stdin inputs is reported", grepl("at most one input may be stdin", r$stderr))

r <- run_subprocess(c("-a", fa, "-b", file.path(tempdir(), "definitely-not-here.bed")))
ok("error: missing -b file exits 2", r$status == 2)
ok("error: missing -b file is reported", grepl("no such file", r$stderr))

bad <- tempfile(); writeLines("chr1\t500\t400\tbad", bad)
r <- run_subprocess(c("-a", bad, "-b", fb))
ok("error: start > end in -a exits 1 (data, not caller)", r$status == 1)

# --- report ----------------------------------------------------------------

cat("---\n")
cat(pass, " passed, ", fail, " failed\n", sep = "")
quit(save = "no", status = if (fail > 0) 1 else 0)
