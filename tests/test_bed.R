#!/usr/bin/env Rscript
# Unit tests for R/bed.R. No bedtools required; runs in seconds.
# Usage: Rscript tests/test_bed.R

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
HERE <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
ROOT <- normalizePath(file.path(HERE, ".."))
DATA <- file.path(ROOT, "data")

source(file.path(ROOT, "R", "bed.R"))

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

# Runs read_bed in a subprocess, since a data error calls quit().
# Returns list(status, stderr).
read_bed_subprocess <- function(path) {
  expr <- sprintf('source("%s"); read_bed("%s")',
                  file.path(ROOT, "R", "bed.R"), path)
  err <- tempfile()
  st <- system2("Rscript", c("-e", shQuote(expr)),
                stdout = NULL, stderr = err)
  list(status = st, stderr = paste(readLines(err, warn = FALSE), collapse = "\n"))
}

# --- the overlap predicate -------------------------------------------------
# SPEC.md section 4: a.start < b.end AND b.start < a.end, strict < both sides.

ok("overlaps: plain overlap", overlaps(100, 200, 150, 250))
ok("overlaps: no overlap, disjoint", !overlaps(0, 100, 200, 300))
ok("overlaps: nested is an overlap", overlaps(300, 400, 320, 350))
ok("overlaps: identical intervals overlap", overlaps(700, 800, 700, 800))
ok("overlaps: is vectorised",
   identical(overlaps(c(0, 0), c(100, 10), c(50, 200), c(150, 300)),
             c(TRUE, FALSE)))

# Bookended: a.end == b.start. Strict < means these do NOT overlap.
# They DO merge at -d 0, which is merge's business, not the predicate's.
ok("overlaps: bookended do not overlap", !overlaps(0, 100, 100, 200))
ok("overlaps: bookended the other way", !overlaps(100, 200, 0, 100))

# Zero-length, by the predicate alone. Note bedtools disagrees in practice --
# see SPEC.md section 4. The subcommands encode the oracle; this function
# stays honest about the predicate.
ok("overlaps: zero-length against itself", !overlaps(500, 500, 500, 500))
# Not uniformly false: a zero-length point strictly inside an interval does
# satisfy the predicate (150 < 200 and 100 < 150). It only goes false at the
# boundaries, where one of the strict comparisons collapses.
ok("overlaps: zero-length strictly inside an interval", overlaps(150, 150, 100, 200))
ok("overlaps: zero-length at an interval's start boundary", !overlaps(100, 100, 100, 200))
ok("overlaps: zero-length at an interval's end boundary", !overlaps(200, 200, 100, 200))

# Position 0 is an ordinary coordinate, not a sentinel.
ok("overlaps: at position 0", overlaps(0, 100, 0, 50))
ok("overlaps: zero-length at position 0", !overlaps(0, 0, 0, 10))

# --- read_bed / write_bed round-trip ---------------------------------------

roundtrip_identical <- function(path) {
  df <- read_bed(path)
  tmp <- tempfile()
  con <- file(tmp, "w"); sink(con); write_bed(df); sink(); close(con)
  identical(readBin(tmp, "raw", file.size(tmp)),
            readBin(path, "raw", file.size(path)))
}

ok("round-trip: a.bed is byte-identical", roundtrip_identical(file.path(DATA, "a.bed")))
ok("round-trip: b.bed is byte-identical", roundtrip_identical(file.path(DATA, "b.bed")))

a <- read_bed(file.path(DATA, "a.bed"))
ok("read_bed: a.bed has 22 rows", nrow(a) == 22)
ok("read_bed: a.bed is BED6", attr(a, "ncol_bed") == 6)
ok("read_bed: preserves input order (a05 first, unsorted)", a$name[1] == "a05")
ok("read_bed: zero-length interval survives parsing",
   any(a$start == a$end))
ok("read_bed: position-0 interval survives parsing", any(a$start == 0))

# Comments, track/browser lines and blanks are skipped silently (SPEC.md s3).
skips <- tempfile()
writeLines(c("# a comment", "track name=x", "browser position chr1",
             "", "chr1\t10\t20\tn1\t0\t+", "chr1\t30\t40\tn2\t0\t-"), skips)
s <- read_bed(skips)
ok("read_bed: skips comment/track/browser/blank lines", nrow(s) == 2)

# Varying column counts: BED3 lines alongside BED6.
mixed <- tempfile()
writeLines(c("chr1\t10\t20", "chr1\t30\t40\tn2\t0\t-"), mixed)
m <- read_bed(mixed)
ok("read_bed: ragged columns widen to the maximum", attr(m, "ncol_bed") == 6)
ok("read_bed: missing trailing fields become '.'", m$name[1] == ".")

bed3 <- tempfile()
writeLines(c("chr1\t10\t20", "chr1\t30\t40"), bed3)
b3 <- read_bed(bed3)
ok("read_bed: BED3 stays BED3", attr(b3, "ncol_bed") == 3)

empty <- tempfile(); invisible(file.create(empty))
ok("read_bed: empty input gives zero rows", nrow(read_bed(empty)) == 0)

# --- errors and exit codes (SPEC.md section 7) -----------------------------

e <- tempfile(); writeLines("chr1\t500\t400\tbad\t0\t+", e)
r <- read_bed_subprocess(e)
ok("error: start > end exits 1", r$status == 1)
ok("error: start > end names file and line", grepl(":1: start > end \\(500 > 400\\)", r$stderr))

e <- tempfile(); writeLines("chr1\tnotanumber\t400", e)
r <- read_bed_subprocess(e)
ok("error: non-integer coordinate exits 1", r$status == 1)
ok("error: non-integer coordinate is reported", grepl("malformed BED line", r$stderr))

e <- tempfile(); writeLines("chr1\t100", e)
r <- read_bed_subprocess(e)
ok("error: too few columns exits 1", r$status == 1)

e <- tempfile(); writeLines("chr1\t-5\t100", e)
r <- read_bed_subprocess(e)
ok("error: negative coordinate exits 1", r$status == 1)
ok("error: negative coordinate is reported", grepl("negative coordinate", r$stderr))

r <- read_bed_subprocess(file.path(tempdir(), "definitely-not-here.bed"))
ok("error: missing file exits 2 (caller error, not data)", r$status == 2)
ok("error: missing file is reported", grepl("no such file", r$stderr))

# --- report ----------------------------------------------------------------

cat("---\n")
cat(pass, " passed, ", fail, " failed\n", sep = "")
quit(save = "no", status = if (fail > 0) 1 else 0)
