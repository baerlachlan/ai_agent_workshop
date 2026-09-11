# mytools intersect -- report overlaps between -a and -b.
#
# SPEC.md s1, s2, s4, s5. Four output modes:
#   (default)  the intersected region, carrying -a's trailing columns
#   -wa        -a's original interval, once per overlapping -b feature
#   -u         each -a feature at most once, if it overlaps anything
#   -v         -a features with no overlap
# -a's input order is preserved; neither input is sorted.

# --- zero-length intervals: oracle behaviour, measured not derived -----------
#
# The predicate in SPEC.md s4 (a.start < b.end AND b.start < a.end) does not
# describe what bedtools does with zero-length features. Measured against
# bedtools v2.31.1, a zero-length -b feature behaves as though it were one base
# wider on each side, and that widening is visible in the coordinates printed:
#
#   a01 chr1 0 100   x  b02 chr1 100 100  ->  chr1 99  100
#   a02 chr1 100 200 x  b02 chr1 100 100  ->  chr1 100 101
#
# Both -a features are merely bookended with b02 and do not overlap it by the
# predicate at all. SPEC.md s4 records the same widening in subtract.
b_effective <- function(b) {
  zero <- b$start == b$end
  list(start = ifelse(zero, b$start - 1L, b$start),
       end = ifelse(zero, b$end + 1L, b$end))
}

# A zero-length -a feature widens the same way, but for the overlap test only:
# all three zero-length features in data/a.bed are reported as overlapping and
# none appear under -v, yet each prints as its full original interval.
#
#   a07 chr1 500 500 x  b07 chr1 500 500  ->  chr1 500 500
#   a12 chr2 0   0   x  b10 chr2 0   10   ->  chr2 0   0
#   a16 chr2 300 300 x  b12 chr2 200 300  ->  chr2 300 300
#
# So the region printed is max/min over -a's *original* coordinates and -b's
# effective ones. Encode it; do not "correct" the predicate.
a_test_coords <- function(a) {
  zero <- a$start == a$end
  list(start = ifelse(zero, a$start - 1L, a$start),
       end = ifelse(zero, a$end + 1L, a$end))
}

# --- hit order --------------------------------------------------------------
#
# bedtools emits the hits for one -a feature in the order its bin index yields
# them, which is not -b file order: finest bin level first, then bin number,
# then position in the file. Every coordinate in data/a.bed and data/b.bed
# lands in the same finest-level bin, where that reduces to file order, but the
# two diverge as soon as a -b feature is wide enough to sit in a coarser bin:
#
#   -b = [chr1 0 20000, chr1 100 155], -a = chr1 150 160
#   bedtools prints  chr1 150 155  (the narrow feature)  before  chr1 150 160
#
# UCSC binning as bedtools builds it: 8 levels, finest 2^14 wide, 8x per level.
# A feature sits at the finest level that contains it whole.
BIN_FIRST_SHIFT <- 14L
BIN_NEXT_SHIFT <- 3L
BIN_LEVELS <- 8L

# Returns the bin level and bin number of each interval. Arithmetic is in
# doubles, not bitwShiftR: a zero-length feature at position 0 widens to a
# start of -1, and bitwShiftR would read that as a large unsigned value.
bin_keys <- function(start, end) {
  s <- as.numeric(start)
  e <- as.numeric(end) - 1 # bedtools bins on the last base, not the end
  level <- rep(BIN_LEVELS - 1L, length(s))
  bin <- numeric(length(s))
  todo <- rep(TRUE, length(s))
  for (k in seq_len(BIN_LEVELS) - 1L) {
    width <- 2^(BIN_FIRST_SHIFT + BIN_NEXT_SHIFT * k)
    sb <- floor(s / width)
    at_level <- todo & sb == floor(e / width)
    level[at_level] <- k
    bin[at_level] <- sb[at_level]
    todo <- todo & !at_level
    if (!any(todo)) break
  }
  # Coordinates too large for the coarsest bin: bedtools refuses the file
  # outright. Keep a deterministic key rather than inventing an error here.
  if (any(todo)) {
    bin[todo] <- floor(s[todo] / 2^(BIN_FIRST_SHIFT +
                                      BIN_NEXT_SHIFT * (BIN_LEVELS - 1L)))
  }
  list(level = level, bin = bin)
}

# --- overlap search ---------------------------------------------------------
#
# Every overlapping pair has either b's start inside a, or a's start inside b.
# Both are contiguous-range queries against a start-sorted vector, so each side
# is one binary search per feature: O((n + m) log) plus the number of pairs
# found. Neither input is scanned against the other (SPEC.md s6: quadratic is a
# bug).

# Expands per-feature index ranges (lo, hi], half-open over positions in a
# sorted vector, into pairs.
range_pairs <- function(outer_idx, lo, hi, inner_sorted) {
  n <- hi - lo
  keep <- n > 0L
  if (!any(keep)) return(list(outer = integer(0), inner = integer(0)))
  n <- n[keep]
  list(outer = rep(outer_idx[keep], n),
       inner = inner_sorted[sequence(n, from = lo[keep] + 1L)])
}

# All overlapping (a, b) index pairs, in the order bedtools emits them.
intersect_pairs <- function(a, b, at, be) {
  pa <- integer(0)
  pb <- integer(0)
  for (chrom in unique(a$chrom)) {
    ai <- which(a$chrom == chrom)
    bi <- which(b$chrom == chrom)
    if (length(bi) == 0L) next

    # b features whose start lies in [a.start, a.end)
    b_ord <- bi[order(be$start[bi], method = "radix")]
    b_starts <- be$start[b_ord]
    p <- range_pairs(ai,
                     findInterval(at$start[ai] - 1L, b_starts),
                     findInterval(at$end[ai] - 1L, b_starts),
                     b_ord)
    pa <- c(pa, p$outer)
    pb <- c(pb, p$inner)

    # the rest: a features whose start lies in (b.start, b.end)
    a_ord <- ai[order(at$start[ai], method = "radix")]
    a_starts <- at$start[a_ord]
    p <- range_pairs(bi,
                     findInterval(be$start[bi], a_starts),
                     findInterval(be$end[bi] - 1L, a_starts),
                     a_ord)
    pa <- c(pa, p$inner)
    pb <- c(pb, p$outer)
  }

  keys <- bin_keys(be$start[pb], be$end[pb])
  ord <- order(pa, keys$level, keys$bin, pb, method = "radix")
  list(a = pa[ord], b = pb[ord])
}

# TRUE for each -a feature that overlaps at least one -b feature. -u and -v
# only need the answer, not the pairs, so this stays linear in the inputs even
# where every feature overlaps every other.
any_overlap <- function(a, b, at, be) {
  hit <- logical(nrow(a))
  for (chrom in unique(a$chrom)) {
    ai <- which(a$chrom == chrom)
    bi <- which(b$chrom == chrom)
    if (length(bi) == 0L) next
    b_ord <- bi[order(be$start[bi], method = "radix")]
    # b features starting before a ends; a overlaps one of them iff the widest
    # end among them reaches past a's start.
    n <- findInterval(at$end[ai] - 1L, be$start[b_ord])
    furthest <- cummax(be$end[b_ord])
    hit[ai] <- n > 0L & furthest[pmax(n, 1L)] > at$start[ai]
  }
  hit
}

# --- arguments --------------------------------------------------------------

parse_intersect_args <- function(argv) {
  opt <- list(a = NULL, b = NULL, u = FALSE, v = FALSE, wa = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[i]
    if (arg %in% c("-a", "-b")) {
      if (i == length(argv)) die_usage("intersect: ", arg, " requires a value")
      opt[[substring(arg, 2L)]] <- argv[i + 1L]
      i <- i + 2L
      next
    }
    switch(arg,
           "-u" = opt$u <- TRUE,
           "-v" = opt$v <- TRUE,
           "-wa" = opt$wa <- TRUE,
           die_usage("unknown flag: ", arg))
    i <- i + 1L
  }

  if (sum(opt$u, opt$v, opt$wa) > 1L) {
    die_usage("intersect: -u, -v and -wa are exclusive")
  }
  if (is.null(opt$a)) die_usage("intersect: -a is required")
  if (is.null(opt$b)) die_usage("intersect: -b is required")
  if (identical(opt$a, "-") && identical(opt$b, "-")) {
    die_usage("at most one input may be stdin")
  }
  # SPEC.md s2: -b must be a real file.
  if (identical(opt$b, "-")) die_usage("intersect: -b must be a file")
  opt
}

# --- entry point ------------------------------------------------------------

cmd_intersect <- function(argv) {
  opt <- parse_intersect_args(argv)
  a <- read_bed(opt$a)
  b <- read_bed(opt$b)
  be <- b_effective(b)

  # bedtools indexes the whole -b file by bin before it looks at -a, and
  # widening a zero-length -b feature at position 0 (above) gives it a start of
  # -1, which its bin index refuses: it rejects the file, prints nothing and
  # exits 1 -- in every mode, whatever chromosome the feature is on, and even
  # when -a is empty. Measured on v2.31.1; match the oracle rather than
  # inventing a kinder answer it would not agree with. Zero-length features are
  # otherwise legal in -b, and legal at position 0 in -a (a12 in data/a.bed).
  if (any(be$start < 0L)) {
    die_data(basename(opt$b),
             ": zero-length feature at position 0 cannot be indexed")
  }

  if (nrow(a) == 0L) return(invisible(NULL))
  if (nrow(b) == 0L) {
    if (opt$v) write_bed(a)
    return(invisible(NULL))
  }

  at <- a_test_coords(a)

  if (opt$u || opt$v) {
    hit <- any_overlap(a, b, at, be)
    write_bed(a[if (opt$v) !hit else hit, , drop = FALSE])
    return(invisible(NULL))
  }

  pairs <- intersect_pairs(a, b, at, be)
  out <- a[pairs$a, , drop = FALSE]
  if (!opt$wa) {
    out$start <- pmax(a$start[pairs$a], be$start[pairs$b])
    out$end <- pmin(a$end[pairs$a], be$end[pairs$b])
  }
  write_bed(out)
  invisible(NULL)
}
