# mytools subtract -a <file|-> -b <file>
#
# Removes every -b region from each -a feature. A feature may be trimmed at
# either end, split into two, or vanish entirely. Trailing columns ride along on
# every fragment. -a input order is preserved; neither input is sorted (SPEC.md
# section 5).
#
# Zero-length records are the whole difficulty, and the rule below is MEASURED
# against bedtools v2.31.1, not derived from the section 4 overlap predicate.
# bedtools widens every zero-length record to [start - 1, end + 1) before it
# looks for overlaps, and prints a surviving zero-length -a record back at its
# original coordinates. That one rule reproduces every oracle result:
#
#   a01 chr1 0 100   minus b02 chr1 100 100  ->  chr1 0 99   (b02 acts as 99..101,
#                                                 though it merely abuts a01)
#   a02 chr1 100 200 minus b02 chr1 100 100  ->  chr1 101 200
#   a08 chr1 500 600 minus b07 chr1 500 500  ->  chr1 501 600
#   a07 chr1 500 500 minus b07 chr1 500 500  ->  gone (499..501 covers 499..501)
#   a12 chr2 0 0     minus b10 chr2 0 10     ->  chr2 0 0    (survives: -1..0 is
#                                                 uncovered; printed unwidened)
#   a16 chr2 300 300 minus b12 chr2 200 300  ->  chr2 300 300
#
# This is oracle behaviour. Do not "correct" it -- see SPEC.md section 4.
#
# The one place the widening has no valid answer is a zero-length -b record at
# position 0, which widens to a start of -1. bedtools refuses the whole -b file
# there; cmd_subtract() below matches that refusal.

parse_subtract_args <- function(args) {
  a_path <- NULL
  b_path <- NULL
  i <- 1L
  while (i <= length(args)) {
    flag <- args[i]
    if (flag != "-a" && flag != "-b") die_usage("unknown flag: ", flag)
    if (i == length(args)) die_usage("subtract: ", flag, " requires an argument")
    if (flag == "-a") a_path <- args[i + 1L] else b_path <- args[i + 1L]
    i <- i + 2L
  }
  if (is.null(a_path)) die_usage("subtract: -a is required")
  if (is.null(b_path)) die_usage("subtract: -b is required")
  if (a_path == "-" && b_path == "-") die_usage("at most one input may be stdin")
  if (b_path == "-") die_usage("subtract: -b must be a real file")
  list(a = a_path, b = b_path)
}

# Union of one chromosome's intervals as disjoint ranges sorted by start.
# Subtracting the union is the same as subtracting each range, and it lets the
# per-feature gap walk below stay vectorised.
merge_ranges <- function(s, e) {
  if (length(s) == 0L) return(list(s = integer(0), e = integer(0)))
  o <- order(s, method = "radix")
  s <- s[o]
  e <- e[o]
  reach <- cummax(e)
  opens <- c(TRUE, s[-1L] > reach[-length(reach)])
  closes <- c(opens[-1L], TRUE)
  list(s = s[opens], e = reach[closes])
}

cmd_subtract <- function(args) {
  opt <- parse_subtract_args(args)
  # -b first: it is always a real file, so a bad path is reported before we
  # block reading -a from stdin.
  b <- read_bed(opt$b)
  a <- read_bed(opt$a)

  # Widen zero-length records by one base each side (see header comment).
  b_zero <- b$start == b$end
  b_start <- b$start - b_zero
  b_end <- b$end + b_zero

  # bedtools indexes the whole -b file by bin before it looks at -a, and
  # widening a zero-length -b feature at position 0 gives it a start of -1,
  # which its bin index refuses: it rejects the file, prints nothing and exits
  # 1 -- whatever chromosome the feature sits on, even when it shares the file
  # with perfectly good features, and even when -a is empty. Measured on
  # v2.31.1; match the oracle rather than inventing a kinder answer it would
  # not agree with. Zero-length features are otherwise legal in -b, and legal
  # at position 0 in -a, where a12 (chr2 0 0) survives. R/intersect.R refuses
  # the same case in the same words.
  if (any(b_start < 0L)) {
    die_data(basename(opt$b),
             ": zero-length feature at position 0 cannot be indexed")
  }

  if (nrow(a) == 0L) return(invisible(NULL))

  a_zero <- a$start == a$end
  a_start <- a$start - a_zero
  a_end <- a$end + a_zero

  b_by_chrom <- split(seq_len(nrow(b)), b$chrom)

  # Fragments accumulate per chromosome, tagged with their -a row, then get
  # reordered into -a input order at the end.
  rows_out <- list()
  starts_out <- list()
  ends_out <- list()
  # nolint start: assignment_linter.
  # <<- is deliberate and contained: emit() accumulates into the three lists in
  # this function's own frame, never a global. The alternative is rebuilding the
  # lists on every call, which is the quadratic growth SPEC.md section 6 rules out.
  emit <- function(rows, s, e) {
    k <- length(rows_out) + 1L
    rows_out[[k]] <<- rows
    starts_out[[k]] <<- s
    ends_out[[k]] <<- e
  }
  # nolint end

  for (ch in unique(a$chrom)) {
    rows <- which(a$chrom == ch)
    as_ <- a_start[rows]
    ae_ <- a_end[rows]

    slot <- match(ch, names(b_by_chrom))
    if (is.na(slot)) {
      emit(rows, as_, ae_)
      next
    }
    m <- merge_ranges(b_start[b_by_chrom[[slot]]], b_end[b_by_chrom[[slot]]])
    ms <- m$s
    me <- m$e

    # Ranges j0..j1 are the ones this feature actually meets: the first ending
    # strictly after it starts, through the last starting strictly before it
    # ends. Strict on both sides, so bookended ranges are not in the window.
    j0 <- findInterval(as_, me) + 1L
    j1 <- findInterval(ae_ - 1L, ms)
    cnt <- pmax(0L, j1 - j0 + 1L)

    untouched <- cnt == 0L
    if (any(untouched)) emit(rows[untouched], as_[untouched], ae_[untouched])

    total <- sum(cnt)
    if (total > 0L) {
      feat <- rep.int(seq_along(rows), cnt)
      nth <- seq_len(total) - rep.int(cumsum(cnt) - cnt, cnt)
      j <- rep.int(j0, cnt) + nth - 1L

      # Gap in front of each met range: from the feature start for the first
      # one, from the previous range's end for the rest.
      lo <- me[pmax(j - 1L, 1L)]
      lo[nth == 1L] <- as_[feat[nth == 1L]]
      hi <- ms[j]
      keep <- lo < hi
      emit(rows[feat[keep]], lo[keep], hi[keep])

      # Gap after the last range it meets.
      tail_of <- which(cnt > 0L)
      lo <- me[j1[tail_of]]
      hi <- ae_[tail_of]
      keep <- lo < hi
      emit(rows[tail_of[keep]], lo[keep], hi[keep])
    }
  }

  idx <- unlist(rows_out, use.names = FALSE)
  frag_start <- unlist(starts_out, use.names = FALSE)
  frag_end <- unlist(ends_out, use.names = FALSE)
  if (length(idx) == 0L) return(invisible(NULL))

  ord <- order(idx, frag_start, method = "radix")
  idx <- idx[ord]
  frag_start <- frag_start[ord]
  frag_end <- frag_end[ord]

  # A zero-length -a feature that survives prints once, at its original
  # coordinates rather than the widened ones.
  zl <- a_zero[idx]
  drop <- zl & duplicated(idx)
  if (any(drop)) {
    idx <- idx[!drop]
    frag_start <- frag_start[!drop]
    frag_end <- frag_end[!drop]
    zl <- zl[!drop]
  }
  frag_start[zl] <- a$start[idx[zl]]
  frag_end[zl] <- a$end[idx[zl]]

  out <- a[idx, , drop = FALSE]
  out$start <- frag_start
  out$end <- frag_end
  write_bed(out)
  invisible(NULL)
}
