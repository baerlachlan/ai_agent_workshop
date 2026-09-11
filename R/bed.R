BED_COLS <- c("chrom", "start", "end", "name", "score", "strand")

die_data <- function(...) {
  cat(paste0(..., "\n"), file = stderr())
  quit(save = "no", status = 1)
}

die_usage <- function(...) {
  cat(paste0(..., "\n"), file = stderr())
  quit(save = "no", status = 2)
}

read_bed <- function(path) {
  if (identical(path, "-")) {
    lines <- readLines(file("stdin"), warn = FALSE)
    label <- "(stdin)"
  } else {
    if (!file.exists(path)) die_usage("no such file: ", path)
    lines <- readLines(path, warn = FALSE)
    label <- basename(path)
  }

  lineno <- seq_along(lines)
  keep <- !grepl("^\\s*$|^#|^track|^browser", lines)
  lines <- lines[keep]
  lineno <- lineno[keep]

  if (length(lines) == 0) {
    empty <- data.frame(chrom = character(), start = integer(), end = integer(),
                        stringsAsFactors = FALSE)
    attr(empty, "ncol_bed") <- 3L
    return(empty)
  }

  fields <- strsplit(lines, "\t", fixed = TRUE)
  widths <- lengths(fields)

  bad <- which(widths < 3L)
  if (length(bad)) die_data(label, ":", lineno[bad[1]], ": malformed BED line")

  ncol_bed <- min(max(widths), 6L)
  padded <- vapply(fields, function(f) {
    length(f) <- ncol_bed
    f
  }, character(ncol_bed))
  m <- matrix(padded, nrow = ncol_bed)

  df <- data.frame(chrom = m[1, ], stringsAsFactors = FALSE)

  # suppressWarnings: non-numeric coordinates become NA and are reported below
  # with a line number, which is more useful than R's coercion warning.
  start <- suppressWarnings(as.integer(m[2, ]))
  end <- suppressWarnings(as.integer(m[3, ]))

  bad <- which(is.na(start) | is.na(end))
  if (length(bad)) die_data(label, ":", lineno[bad[1]], ": malformed BED line")

  bad <- which(start < 0L | end < 0L)
  if (length(bad)) {
    neg <- if (start[bad[1]] < 0L) start[bad[1]] else end[bad[1]]
    die_data(label, ":", lineno[bad[1]], ": negative coordinate (", neg, ")")
  }

  bad <- which(start > end)
  if (length(bad)) {
    die_data(label, ":", lineno[bad[1]], ": start > end (",
             start[bad[1]], " > ", end[bad[1]], ")")
  }

  df$start <- start
  df$end <- end
  if (ncol_bed >= 4L) df$name <- ifelse(is.na(m[4, ]), ".", m[4, ])
  if (ncol_bed >= 5L) df$score <- ifelse(is.na(m[5, ]), ".", m[5, ])
  if (ncol_bed >= 6L) df$strand <- ifelse(is.na(m[6, ]), ".", m[6, ])

  attr(df, "ncol_bed") <- ncol_bed
  df
}

write_bed <- function(df) {
  if (nrow(df) == 0L) return(invisible(NULL))
  cols <- intersect(BED_COLS, names(df))
  out <- do.call(paste, c(lapply(cols, function(cl) as.character(df[[cl]])),
                          list(sep = "\t")))
  writeLines(out)
  invisible(NULL)
}

overlaps <- function(a_start, a_end, b_start, b_end) {
  a_start < b_end & b_start < a_end
}

chrom_order <- function(chrom) order(chrom, method = "radix")
