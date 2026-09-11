#!/usr/bin/env Rscript
# Unit tests for R/sort.R. No bedtools required; runs in seconds.
# Usage: Rscript tests/test_sort.R

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
HERE <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
ROOT <- normalizePath(file.path(HERE, ".."))
DATA <- file.path(ROOT, "data")
MYTOOLS <- file.path(ROOT, "mytools")

source(file.path(ROOT, "R", "bed.R"))
source(file.path(ROOT, "R", "sort.R"))

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

# Builds a data frame the way read_bed would, without touching the disk.
bed <- function(chrom, start, end, name = NULL) {
  df <- data.frame(chrom = chrom, start = as.integer(start),
                   end = as.integer(end), stringsAsFactors = FALSE)
  if (!is.null(name)) df$name <- name
  df
}

# Runs the CLI in a subprocess. Returns list(status, stdout, stderr).
run_mytools <- function(..., stdin = NULL) {
  out <- tempfile(); err <- tempfile()
  st <- system2(MYTOOLS, c(...), stdout = out, stderr = err,
                stdin = if (is.null(stdin)) "" else stdin)
  list(status = st,
       stdout = readLines(out, warn = FALSE),
       stderr = paste(readLines(err, warn = FALSE), collapse = "\n"))
}

# --- sort order (SPEC.md section 5) ----------------------------------------

s <- sort_bed(bed(c("chr1", "chr2", "chr1"), c(50, 10, 10), c(60, 20, 20)))
ok("sort: chrom first, then start", identical(s$chrom, c("chr1", "chr1", "chr2")))
ok("sort: start within a chrom", identical(s$start, c(10L, 50L, 10L)))

s <- sort_bed(bed(rep("chr1", 3), c(100, 100, 100), c(300, 150, 200)))
ok("sort: end breaks a start tie", identical(s$end, c(150L, 200L, 300L)))

# --- lexicographic, and locale-independently so -----------------------------
# The classic bite: R's default collation follows LC_COLLATE, which puts
# "chr7" before "chr17" under a UTF-8 locale. bedtools sorts bytes.

s <- sort_bed(bed(c("chr7", "chr17", "chr2", "chr10", "chr1", "chrX"),
                  rep(0, 6), rep(10, 6)))
ok("sort: chr17 sorts before chr7 (lexicographic, not numeric)",
   identical(s$chrom, c("chr1", "chr10", "chr17", "chr2", "chr7", "chrX")))

# Case is where the two collations disagree outright: C puts "chrB" before
# "chra" on the byte values, en_US folds case and puts "chra" first. Setting
# the locale is a no-op if it is not installed, in which case this still
# passes -- it just stops proving anything.
old_collate <- Sys.getlocale("LC_COLLATE")
invisible(suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8")))
s <- sort_bed(bed(c("chra", "chrB"), c(0, 0), c(10, 10)))
ok("sort: byte order survives a UTF-8 LC_COLLATE",
   identical(s$chrom, c("chrB", "chra")))
invisible(suppressWarnings(Sys.setlocale("LC_COLLATE", old_collate)))

# --- ties ------------------------------------------------------------------
# bedtools v2.31.1 orders by (chrom, start) alone and does not promise an
# order within a tie -- it is stable on small inputs but reorders 2000
# equal-start features. We pin (chrom, start, end) and a stable sort, which
# agrees with bedtools on every fixture and is reproducible where bedtools
# is not. a03/a04 and a09/a10 are the fixture's tie cases.

s <- sort_bed(bed(rep("chr1", 4), c(150, 150, 700, 700), c(250, 250, 800, 800),
                  c("a03", "a04", "a09", "a10")))
ok("sort: identical coordinates keep input order (stable)",
   identical(s$name, c("a03", "a04", "a09", "a10")))

a <- sort_bed(read_bed(file.path(DATA, "a.bed")))
ok("sort: a03 before a04 in the fixture",
   which(a$name == "a03") < which(a$name == "a04"))
ok("sort: a09 before a10 in the fixture",
   which(a$name == "a09") < which(a$name == "a10"))

# --- edge cases ------------------------------------------------------------
# Position 0 is an ordinary coordinate and sorts first, not last.
s <- sort_bed(bed(rep("chr2", 3), c(50, 0, 200), c(150, 0, 300)))
ok("sort: position 0 sorts first", identical(s$start, c(0L, 50L, 200L)))

# Zero-length intervals sort by end like anything else: chr1 500 500 (a07)
# comes before chr1 500 600 (a08), which is what bedtools prints.
s <- sort_bed(bed(rep("chr1", 2), c(500, 500), c(600, 500), c("a08", "a07")))
ok("sort: zero-length sorts before a longer interval at the same start",
   identical(s$name, c("a07", "a08")))
ok("sort: a07 before a08 in the fixture",
   which(a$name == "a07") < which(a$name == "a08"))

# Nested intervals are just a start tie that isn't: a06 (320) after a05 (300).
ok("sort: nested interval follows its container",
   which(a$name == "a05") < which(a$name == "a06"))

blank <- tempfile(); invisible(file.create(blank))
ok("sort: empty input gives zero rows", nrow(sort_bed(read_bed(blank))) == 0)

# --- columns are preserved (SPEC.md section 5) -----------------------------

ok("sort: BED6 in, BED6 out", identical(names(a), BED_COLS))
ok("sort: row count is unchanged", nrow(a) == 22)
ok("sort: no feature is lost or duplicated",
   identical(sort(a$name, method = "radix"),
             sort(read_bed(file.path(DATA, "a.bed"))$name, method = "radix")))

bed3 <- tempfile()
writeLines(c("chr1\t30\t40", "chr1\t10\t20"), bed3)
s <- sort_bed(read_bed(bed3))
ok("sort: BED3 stays BED3", identical(names(s), c("chrom", "start", "end")))

# --- the CLI (SPEC.md section 7) -------------------------------------------
# These are the deliberate deviations from bedtools, so they have no oracle
# and cannot be golden tests: real bedtools exits 0 with no -i, and 1 on a
# missing file. SPEC.md section 7 makes both caller errors, exit 2.

r <- run_mytools("sort")
ok("cli: no -i exits 2", r$status == 2)
ok("cli: no -i explains itself on stderr", grepl("-i is required", r$stderr))
ok("cli: no -i prints nothing on stdout", length(r$stdout) == 0)

ok("cli: -i with no value exits 2", run_mytools("sort", "-i")$status == 2)

r <- run_mytools("sort", "-i", file.path(tempdir(), "definitely-not-here.bed"))
ok("cli: missing file exits 2", r$status == 2)
ok("cli: missing file is reported", grepl("no such file", r$stderr))

r <- run_mytools("sort", "-q", "-i", file.path(DATA, "a.bed"))
ok("cli: unknown flag exits 2", r$status == 2)
ok("cli: unknown flag is named", grepl("unknown flag: -q", r$stderr))

r <- run_mytools("sort", file.path(DATA, "a.bed"))
ok("cli: bare positional argument exits 2", r$status == 2)

r <- run_mytools("sort", "-i", file.path(DATA, "a.bed"))
ok("cli: success exits 0", r$status == 0)
ok("cli: success prints every row", length(r$stdout) == 22)

r <- run_mytools("sort", "-i", "-", stdin = file.path(DATA, "a.bed"))
ok("cli: '-' reads stdin", r$status == 0 && length(r$stdout) == 22)
ok("cli: stdin gives the same answer as the file",
   identical(r$stdout, run_mytools("sort", "-i", file.path(DATA, "a.bed"))$stdout))

empty <- tempfile(); invisible(file.create(empty))
r <- run_mytools("sort", "-i", empty)
ok("cli: empty input prints nothing and exits 0",
   r$status == 0 && length(r$stdout) == 0)

# --- report ----------------------------------------------------------------

cat("---\n")
cat(pass, " passed, ", fail, " failed\n", sep = "")
quit(save = "no", status = if (fail > 0) 1 else 0)
