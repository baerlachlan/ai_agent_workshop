# mytools merge -- see SPEC.md sections 1, 2, 4, 5 and 8.
#
# Merges features whose gap is <= d (default 0), so bookended features merge
# even though they do not overlap by the SPEC.md section 4 predicate. Output is
# BED3 only; input columns are dropped.
#
# Deviation from bedtools (SPEC.md section 8): bedtools refuses unsorted input
# and exits 1, we sort silently and exit 0. Hence the whole input is held in
# memory, and the unsorted case has no oracle -- tests/test_merge.R covers it.

# Zero-length intervals, measured against bedtools v2.31.1 (SPEC.md section 4).
# bedtools stores a zero-length feature [p, p) as [p - 1, p + 1) and clusters on
# those coordinates; the original coordinates are printed again only when the
# feature ends up in a cluster of its own. Every measured case falls out of that
# one rule:
#
#   chr1 500 500 + chr1 500 600  ->  chr1 499 600   ([499,501) u [500,600))
#   chr1 500 500 + chr1 501 600  ->  chr1 499 600   (gap 501 - 501 = 0)
#   chr1 400 499 + chr1 500 500  ->  chr1 400 501   (the +1 end shows up too)
#   chr1 500 500 + chr1 500 500  ->  chr1 499 501   (cluster of two, expanded)
#   chr1 500 500 alone           ->  chr1 500 500   (lone, so printed as read)
#   chr1   0   0 + chr1   0  10  ->  chr1  -1  10   (yes, a negative start)
#   chr1   0   0 alone           ->  chr1   0   0
#
# The expansion is not clamped at 0: `bedtools merge -d 100` on sorted a.bed
# prints `chr2 -1 400`. Encode it; do not "fix" it.
zero_length_span <- function(start, end) {
  zero <- start == end
  list(start = ifelse(zero, start - 1L, start),
       end = ifelse(zero, end + 1L, end),
       zero = zero)
}

# For each feature: the running end of the cluster it lands in, and whether it
# opens that cluster. A feature joins the cluster in progress when its gap to
# that cluster's end is <= d.
cluster_scan <- function(s, e, d, first_of_chrom) {
  n <- length(s)

  if (d >= 0L) {
    # With d >= 0 the running cluster end is just the per-chromosome prefix
    # maximum: a feature only opens a new cluster when s[i] exceeds that
    # maximum by more than d, so its own end exceeds every end before it and
    # the prefix maximum restarts from the new cluster anyway. One loop
    # iteration per chromosome, not per feature (SPEC.md section 6).
    chrom_first <- which(first_of_chrom)
    chrom_last <- c(chrom_first[-1L] - 1L, n)
    run_end <- integer(n)
    for (k in seq_along(chrom_first)) {
      idx <- chrom_first[k]:chrom_last[k]
      run_end[idx] <- cummax(e[idx])
    }
    prev_end <- c(0L, run_end[-n])
    return(list(run_end = run_end,
                new_cluster = first_of_chrom | s - prev_end > d))
  }

  # Negative d demands a given number of overlapping bases, and then the prefix
  # maximum is no longer the cluster's end: a nested feature can open a cluster
  # that ends before the feature containing it (`-d -100` on sorted a.bed keeps
  # `chr1 320 350` inside `chr1 300 400`). That running maximum has to be reset
  # per cluster, which is sequential -- ~0.1s per 500k features, measured.
  run_end <- integer(n)
  new_cluster <- logical(n)
  cur <- e[1L]
  new_cluster[1L] <- TRUE
  run_end[1L] <- cur
  for (i in seq_len(n)[-1L]) {
    if (first_of_chrom[i] || s[i] - cur > d) {
      cur <- e[i]
      new_cluster[i] <- TRUE
    } else if (e[i] > cur) {
      cur <- e[i]
    }
    run_end[i] <- cur
  }
  list(run_end = run_end, new_cluster = new_cluster)
}

# Clusters intervals and returns a BED3 data frame. Sorts first, because our
# input may be unsorted (SPEC.md section 8) and the clustering below walks one
# chromosome at a time in ascending start order.
merge_intervals <- function(chrom, start, end, d = 0L) {
  n <- length(chrom)
  if (n == 0L) {
    return(data.frame(chrom = character(), start = integer(), end = integer(),
                      stringsAsFactors = FALSE))
  }

  # Key is (chrom, start) and the sort must be stable, because `bedtools sort`
  # is: it leaves features sharing a chrom and start in input order, and merge
  # keeps the *first* feature's start, so tie order changes the output.
  # `chr1 54 74` then `chr1 54 54` merges to `chr1 54 74`; the same two the
  # other way round give `chr1 53 74`. order(method = "radix") is stable.
  o <- order(chrom, start, method = "radix")
  chrom <- chrom[o]
  start <- start[o]
  end <- end[o]

  span <- zero_length_span(start, end)
  first_of_chrom <- c(TRUE, chrom[-1L] != chrom[-n])
  scan <- cluster_scan(span$start, span$end, d, first_of_chrom)

  cluster_first <- which(scan$new_cluster)
  cluster_last <- c(cluster_first[-1L] - 1L, n)

  # The cluster's start is its first feature's start, not the smallest start in
  # it: a zero-length feature shifts to p - 1, and when it sorts behind another
  # feature starting at p that lower value is simply not reported (measured --
  # see the stable-sort note above).
  out_start <- span$start[cluster_first]
  out_end <- scan$run_end[cluster_last]

  # A lone zero-length feature prints its original coordinates: the [p-1, p+1)
  # span is bedtools' internal form, not its output.
  lone_zero <- cluster_first == cluster_last & span$zero[cluster_first]
  out_start[lone_zero] <- start[cluster_first][lone_zero]
  out_end[lone_zero] <- end[cluster_first][lone_zero]

  data.frame(chrom = chrom[cluster_first], start = out_start, end = out_end,
             stringsAsFactors = FALSE)
}

cmd_merge <- function(args) {
  input <- NULL
  d <- 0L

  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a == "-i") {
      if (i == length(args)) die_usage("merge: -i requires a file")
      input <- args[i + 1L]
      i <- i + 2L
    } else if (a == "-d") {
      if (i == length(args)) die_usage("merge: -d requires an integer")
      # suppressWarnings: a non-integer -d becomes NA and is reported as a
      # usage error below, which beats R's coercion warning.
      d <- suppressWarnings(as.integer(args[i + 1L]))
      if (is.na(d)) die_usage("merge: -d requires an integer")
      i <- i + 2L
    } else if (a != "-" && startsWith(a, "-")) {
      die_usage("unknown flag: ", a)
    } else {
      die_usage("merge: unexpected argument: ", a)
    }
  }

  if (is.null(input)) die_usage("merge: -i is required")

  bed <- read_bed(input)
  write_bed(merge_intervals(bed$chrom, bed$start, bed$end, d))
}
