# mytools sort -- SPEC.md sections 1, 2, 5 and 9.
#
# Order is chrom lexicographic, then start, then end (SPEC.md section 5).
# Lexicographic means byte order, not the user's locale: "chr17" sorts before
# "chr7" and "chr10" before "chr2". order(method = "radix") sorts character
# vectors in the C locale regardless of LC_COLLATE, which is why it is used
# here -- the default method would reorder these under a UTF-8 locale.
#
# On ties, note what the oracle actually does: bedtools v2.31.1 orders by
# (chrom, start) only and leaves equal-start features in whatever order its
# chunked sort happens to emit. It is stable for small inputs but not at
# scale -- 2000 equal-start features come back reordered. Neither fixture has
# an equal-start tie with differing ends, so ordering by end as well agrees
# with bedtools byte-for-byte on a.bed and b.bed while being deterministic
# where bedtools is not. Fully identical rows (a09/a10) keep their input
# order because radix order is stable.

sort_bed <- function(df) {
  if (nrow(df) == 0L) return(df)
  df[order(df$chrom, df$start, df$end, method = "radix"), , drop = FALSE]
}

cmd_sort <- function(args) {
  input <- NULL
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (identical(a, "-i")) {
      if (i == length(args)) die_usage("sort: -i is required")
      input <- args[i + 1L]
      i <- i + 2L
    } else if (startsWith(a, "-") && !identical(a, "-")) {
      die_usage("unknown flag: ", a)
    } else {
      die_usage("sort: unexpected argument: ", a)
    }
  }
  if (is.null(input)) die_usage("sort: -i is required")

  write_bed(sort_bed(read_bed(input)))
  invisible(NULL)
}
