#!/usr/bin/env Rscript
# Unit tests for R/intersect.R. No bedtools required; runs in seconds.
# Usage: Rscript tests/test_intersect.R
#
# Expected values here were measured against bedtools v2.31.1 (see SPEC.md s4 and the
# comments in R/intersect.R). Where one looks wrong, it is the oracle's answer, not a
# guess -- the golden tests in tests/golden_intersect.sh are what keep it honest.

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
HERE <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
ROOT <- normalizePath(file.path(HERE, ".."))
DATA <- file.path(ROOT, "data")

source(file.path(ROOT, "R", "bed.R"))
source(file.path(ROOT, "R", "intersect.R"))

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

# Output lines of one intersect run, as a character vector.
run <- function(...) capture.output(cmd_intersect(c(...)))

# Runs mytools in a subprocess, since a usage or data error calls quit().
# Returns list(status, stdout, stderr).
run_subprocess <- function(...) {
  out <- tempfile()
  err <- tempfile()
  st <- system2(file.path(ROOT, "mytools"), c("intersect", ...),
                stdout = out, stderr = err)
  list(status = st,
       stdout = readLines(out, warn = FALSE),
       stderr = paste(readLines(err, warn = FALSE), collapse = "\n"))
}

bed <- function(...) {
  path <- tempfile(fileext = ".bed")
  lines <- c(...)
  writeLines(if (is.null(lines)) character(0) else lines, path)
  path
}

A <- file.path(DATA, "a.bed")
B <- file.path(DATA, "b.bed")

# --- the four modes are four different questions ----------------------------
# a05 (chr1 300 400) overlaps two -b features, b05 (320 350) and b06 (340 360).

ok("default: prints the intersected region, not the -a interval",
   identical(grep("a05", run("-a", A, "-b", B), value = TRUE),
             c("chr1\t320\t350\ta05\t40\t+", "chr1\t340\t360\ta05\t40\t+")))
ok("-wa: prints the -a interval once per overlapping -b feature",
   identical(grep("a05", run("-wa", "-a", A, "-b", B), value = TRUE),
             rep("chr1\t300\t400\ta05\t40\t+", 2)))
ok("-u: prints the -a interval at most once",
   identical(grep("a05", run("-u", "-a", A, "-b", B), value = TRUE),
             "chr1\t300\t400\ta05\t40\t+"))
ok("-v: an overlapping feature is absent",
   length(grep("a05", run("-v", "-a", A, "-b", B))) == 0)

# -a's trailing columns are carried through, and its input order is kept: a.bed is
# deliberately unsorted and starts at a05.
ok("default: -a's input order is preserved (a05 before a01)",
   identical(sub("\t.*", "", sub("^\\S+\t\\d+\t\\d+\t", "",
                                 run("-a", A, "-b", B)))[1:3],
             c("a05", "a05", "a01")))

# --- bookended features do not overlap (SPEC.md s4, strict <) ---------------
# a01 chr1 0 100 and a02 chr1 100 200 are bookended at 100. Tested on their own so
# that b02 (chr1 100 100, zero-length) cannot muddy the result -- its widening is a
# separate case, below.

bookend_a <- bed("chr1\t0\t100\ta01\t10\t+")
bookend_b <- bed("chr1\t100\t200\tbook\t0\t+")
ok("bookended: no overlap reported",
   length(run("-a", bookend_a, "-b", bookend_b)) == 0)
ok("bookended: the -a feature survives -v",
   identical(run("-v", "-a", bookend_a, "-b", bookend_b),
             "chr1\t0\t100\ta01\t10\t+"))
ok("bookended the other way round: no overlap",
   length(run("-a", bookend_b, "-b", bookend_a)) == 0)

# --- zero-length -a features: oracle behaviour ------------------------------
# a07 (chr1 500 500), a12 (chr2 0 0) and a16 (chr2 300 300) are each reported as
# overlapping, and each prints as its full original interval, even though the s4
# predicate says none of them overlaps anything. Measured; do not "correct" it.

zero_lines <- c("chr1\t500\t500\ta07\t0\t+",
                "chr2\t0\t0\ta12\t0\t+",
                "chr2\t300\t300\ta16\t0\t+")
ok("zero-length -a: all three are reported under default, as themselves",
   all(zero_lines %in% run("-a", A, "-b", B)))
ok("zero-length -a: all three are reported under -u, as themselves",
   all(zero_lines %in% run("-u", "-a", A, "-b", B)))
ok("zero-length -a: none of the three appears under -v",
   !any(zero_lines %in% run("-v", "-a", A, "-b", B)))
ok("zero-length -a: each is reported exactly once under -wa",
   identical(as.integer(table(factor(
     grep("a07|a12|a16", run("-wa", "-a", A, "-b", B), value = TRUE),
     levels = zero_lines))), c(1L, 1L, 1L)))

# A zero-length -a feature at position 0 is legal and is not a special case: a12 is
# chr2 0 0 and hits b10 (chr2 0 10), which the strict predicate also denies.
ok("zero-length -a at position 0 is reported",
   identical(run("-a", bed("chr2\t0\t0\ta12\t0\t+"),
                 "-b", bed("chr2\t0\t10\tb10\t18\t-")),
             "chr2\t0\t0\ta12\t0\t+"))

# --- zero-length -b features widen by a base on each side -------------------
# b02 is chr1 100 100. It is bookended with both a01 and a02 and overlaps neither by
# the predicate, yet bedtools reports both -- and the widening shows up in the
# coordinates it prints, one base outside the zero-length feature on each side.

ok("zero-length -b: widens leftward into a bookended -a (a01 -> 99 100)",
   "chr1\t99\t100\ta01\t10\t+" %in% run("-a", A, "-b", B))
ok("zero-length -b: widens rightward into a bookended -a (a02 -> 100 101)",
   "chr1\t100\t101\ta02\t20\t-" %in% run("-a", A, "-b", B))
ok("zero-length -b: a zero-length -a at the same position still prints as itself",
   identical(run("-a", bed("chr1\t500\t500\ta07\t0\t+"),
                 "-b", bed("chr1\t500\t500\tb07\t0\t-")),
             "chr1\t500\t500\ta07\t0\t+"))
ok("zero-length -b: reported region is the overlap of -a with the widened -b",
   identical(run("-a", bed("chr1\t0\t10\ta\t1\t+"),
                 "-b", bed("chr1\t5\t5\tz\t0\t+")),
             "chr1\t4\t6\ta\t1\t+"))

# Widening a zero-length -b feature at position 0 gives it a start of -1, which
# bedtools' bin index refuses: it rejects the whole file and exits 1 with empty
# stdout, in every mode and whatever chromosome the feature is on. We match that.
r <- run_subprocess("-a", A, "-b", bed("chr9\t0\t0\tz\t0\t+"))
ok("zero-length -b at position 0: exits 1 like bedtools", r$status == 1)
ok("zero-length -b at position 0: stdout stays empty", length(r$stdout) == 0)
ok("zero-length -b at position 0: says why on stderr",
   grepl("position 0", r$stderr))

# --- nesting is an overlap, not a containment test --------------------------
# a06 (chr1 320 350) is inside a05; b08 (chr1 750 760) is inside a09 (700 800).

ok("nested -b inside -a: reported as the -b span",
   "chr1\t750\t760\ta09\t60\t+" %in% run("-a", A, "-b", B))
ok("nested -a inside -b: reported as the -a span",
   identical(run("-a", bed("chr1\t320\t350\ta06\t50\t-"),
                 "-b", bed("chr1\t300\t400\touter\t0\t+")),
             "chr1\t320\t350\ta06\t50\t-"))

# --- identical -a features are independent ----------------------------------
# a09 and a10 are both chr1 700 800 and both hit b08. Each is answered on its own;
# nothing is deduplicated across features.

ok("identical -a features are each reported (default)",
   identical(grep("a09|a10", run("-a", A, "-b", B), value = TRUE),
             c("chr1\t750\t760\ta09\t60\t+", "chr1\t750\t760\ta10\t60\t+")))
ok("identical -a features are each reported (-u)",
   length(grep("a09|a10", run("-u", "-a", A, "-b", B))) == 2)

# --- chromosomes present in one input only ----------------------------------
# chr3 exists in -b (b19) and not in -a. chrX 1000 1100 (a22) is the reverse case.

out_all <- c(run("-a", A, "-b", B), run("-u", "-a", A, "-b", B),
             run("-v", "-a", A, "-b", B), run("-wa", "-a", A, "-b", B))
ok("a chromosome only in -b contributes nothing to the output",
   !any(grepl("^chr3\t|b19", out_all)))
ok("an -a feature on a chromosome with no -b features appears under -v",
   "chrX\t1000\t1100\ta22\t95\t+" %in% run("-v", "-a", A, "-b", B))

# --- hit order follows bedtools' bin index, not -b file order ---------------
# Finest bin level first, then bin number, then file order. Two -b features over one
# -a feature: the narrow one is reported first even though the wide one is listed
# first in the file, because the wide one sits in a coarser bin.

wide_first <- bed("chr1\t0\t20000\twide\t0\t+", "chr1\t100\t155\tnarrow\t0\t+")
ok("hit order: narrow bin before wide bin, whatever the -b file order",
   identical(run("-a", bed("chr1\t150\t160\ta1\t1\t+"), "-b", wide_first),
             c("chr1\t150\t155\ta1\t1\t+", "chr1\t150\t160\ta1\t1\t+")))
ok("hit order: within one bin it is -b file order",
   identical(run("-a", bed("chr1\t150\t160\ta1\t1\t+"),
                 "-b", bed("chr1\t140\t158\tsecond\t0\t+",
                           "chr1\t100\t155\tfirst\t0\t+")),
             c("chr1\t150\t158\ta1\t1\t+", "chr1\t150\t155\ta1\t1\t+")))

# --- empty inputs -----------------------------------------------------------

none <- bed()
ok("empty -b: nothing overlaps", length(run("-a", A, "-b", none)) == 0)
ok("empty -b: every -a feature appears under -v",
   length(run("-v", "-a", A, "-b", none)) == 22)
ok("empty -a: no output, and none under -v",
   length(run("-a", none, "-b", B)) == 0 &&
     length(run("-v", "-a", none, "-b", B)) == 0)

# --- usage errors (SPEC.md s7: exit 2, stderr, stdout clean) ----------------

r <- run_subprocess("-u", "-v", "-a", A, "-b", B)
ok("-u -v together exits 2", r$status == 2)
ok("-u -v together says the flags are exclusive",
   grepl("-u, -v and -wa are exclusive", r$stderr))
ok("-u -wa together exits 2", run_subprocess("-u", "-wa", "-a", A, "-b", B)$status == 2)
ok("-v -wa together exits 2", run_subprocess("-v", "-wa", "-a", A, "-b", B)$status == 2)
ok("all three together exits 2",
   run_subprocess("-u", "-v", "-wa", "-a", A, "-b", B)$status == 2)

r <- run_subprocess("-q", "-a", A, "-b", B)
ok("unknown flag exits 2", r$status == 2)
ok("unknown flag is named", grepl("unknown flag: -q", r$stderr))

ok("missing -a exits 2", run_subprocess("-b", B)$status == 2)
ok("missing -b exits 2", run_subprocess("-a", A)$status == 2)
ok("-a without a value exits 2", run_subprocess("-b", B, "-a")$status == 2)

r <- run_subprocess("-a", "-", "-b", "-")
ok("two stdin inputs exits 2", r$status == 2)
ok("two stdin inputs is reported",
   grepl("at most one input may be stdin", r$stderr))

# -b must be a real file (SPEC.md s2).
r <- run_subprocess("-a", A, "-b", "-")
ok("-b as stdin exits 2", r$status == 2)

r <- run_subprocess("-a", A, "-b", file.path(tempdir(), "definitely-not-here.bed"))
ok("missing -b file exits 2 (caller error, not data)", r$status == 2)

# A data error in either input is exit 1, not 2.
r <- run_subprocess("-a", A, "-b", bed("chr1\t500\t400\tbad\t0\t+"))
ok("start > end in -b exits 1", r$status == 1)
r <- run_subprocess("-a", bed("chr1\tx\t400\tbad\t0\t+"), "-b", B)
ok("malformed -a line exits 1", r$status == 1)

# --- report ----------------------------------------------------------------

cat("---\n")
cat(pass, " passed, ", fail, " failed\n", sep = "")
quit(save = "no", status = if (fail > 0) 1 else 0)
