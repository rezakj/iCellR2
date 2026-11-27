#' Merge peaks from several samples into one shared set
#'
#' Each sample's peak caller writes slightly different coordinates. This
#' reduces overlapping (and, by default, adjacent) intervals on the same
#' chromosome into one consensus peak, so every sample can use the same
#' feature ids. Peak names use \code{chr_start_end}.
#' \code{chr:start-end} (10x) and \code{chr.start.end} are accepted too.
#'
#' The returned table is the shared peak set. \code{attr(x, "map")} links
#' every input peak to its consensus id, which is what a later count merge
#' uses. This function does not read fragment files and does not build a
#' matrix.
#'
#' @param peaks A character vector of peak ids or BED paths; a named list of
#'   those (names are sample ids); or a data.frame with chromosome, start,
#'   and end columns and an optional sample column.
#' @param gap Extra bases, beyond a shared edge, that still count as one
#'   peak. Default \code{0}: overlapping or book-ended peaks merge
#'   (\code{chr1:1-10} with \code{chr1:11-20}). Use \code{-1} to merge only
#'   real overlaps. A positive value also bridges small gaps.
#' @param min.width Drop consensus peaks shorter than this. Default \code{1}.
#' @return A data.frame with columns \code{chromosome}, \code{start},
#'   \code{end}, \code{peak_id}, \code{n_samples}, \code{n_peaks}. The
#'   \code{map} attribute has \code{sample}, \code{peak_id}, \code{consensus_id}.
#'
#' @examples
#' \dontrun{
#' shared <- merge.peaks(list(
#'   s1 = c("chr1:100-200", "chr2:1-50"),
#'   s2 = c("chr1:180-250", "chr2:80-90")
#' ))
#' shared
#' attr(shared, "map")
#' }
#'
#' @export
merge.peaks <- function(peaks, gap = 0L, min.width = 1L) {
  gap <- as.integer(gap)[1L]
  min.width <- as.integer(min.width)[1L]
  if (is.na(gap)) stop("'gap' must be an integer.")
  if (is.na(min.width) || min.width < 1L) stop("'min.width' must be a positive integer.")

  raw <- .ic2_peaks_collect(peaks)
  if (nrow(raw) < 1L) stop("No peaks to merge.")

  parts <- split(raw, raw$chromosome, drop = TRUE)
  cons_rows <- vector("list", length(parts))
  map_rows <- vector("list", length(parts))
  for (i in seq_along(parts)) {
    one <- .ic2_peaks_reduce_chrom(parts[[i]], gap)
    keep <- (one$consensus$end - one$consensus$start + 1L) >= min.width
    if (!any(keep)) next
    kept_id <- one$consensus$peak_id[keep]
    cons_rows[[i]] <- one$consensus[keep, , drop = FALSE]
    map_i <- one$map[one$map$consensus_id %in% kept_id, , drop = FALSE]
    map_rows[[i]] <- map_i
  }
  cons <- do.call(rbind, cons_rows)
  map <- do.call(rbind, map_rows)
  if (is.null(cons) || nrow(cons) < 1L) stop("No peaks left after min.width.")
  rownames(cons) <- NULL
  rownames(map) <- NULL
  cons <- cons[order(cons$chromosome, cons$start, cons$end), , drop = FALSE]
  map <- map[order(map$sample, map$chromosome, map$start), , drop = FALSE]
  attr(cons, "map") <- map
  attr(cons, "gap") <- gap
  class(cons) <- c("iCellR2_merged_peaks", "data.frame")
  cons
}

#' Parse one peak id into chromosome, start, end.
#' Names use underscores (\code{chr_start_end}) or 10x \code{chr:start-end}.
#' @keywords internal
.ic2_parse_peak_id <- function(id) {
  id <- as.character(id)
  n <- length(id)
  chrom <- rep(NA_character_, n)
  start <- rep(NA_integer_, n)
  end <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    s <- id[[i]]
    if (is.na(s) || !nzchar(s)) next
    hit <- .ic2_parse_one_peak(s)
    if (is.null(hit)) next
    chrom[[i]] <- hit[[1L]]
    start[[i]] <- hit[[2L]]
    end[[i]] <- hit[[3L]]
  }
  data.frame(
    chromosome = chrom, start = start, end = end,
    stringsAsFactors = FALSE
  )
}

.ic2_parse_one_peak <- function(s) {
  s <- trimws(s)
  if (grepl(":", s, fixed = TRUE)) {
    bits <- strsplit(s, ":", fixed = TRUE)[[1L]]
    if (length(bits) < 2L) return(NULL)
    chrom <- paste(bits[-length(bits)], collapse = ":")
    se <- bits[[length(bits)]]
    se <- sub("^\\s+", "", se)
    mm <- regexec("^([0-9]+)\\s*[-]\\s*([0-9]+)$", se)
    rr <- regmatches(se, mm)[[1L]]
    if (length(rr) != 3L) return(NULL)
    return(.ic2_peak_span(chrom, rr[[2L]], rr[[3L]]))
  }
  norm <- gsub(".", "_", s, fixed = TRUE)
  bits <- strsplit(norm, "_", fixed = TRUE)[[1L]]
  if (length(bits) >= 3L &&
      grepl("^[0-9]+$", bits[[length(bits)]]) &&
      grepl("^[0-9]+$", bits[[length(bits) - 1L]])) {
    chrom <- paste(bits[seq_len(length(bits) - 2L)], collapse = "_")
    return(.ic2_peak_span(chrom, bits[[length(bits) - 1L]], bits[[length(bits)]]))
  }
  bits <- strsplit(s, "-", fixed = TRUE)[[1L]]
  if (length(bits) >= 3L &&
      grepl("^[0-9]+$", bits[[length(bits)]]) &&
      grepl("^[0-9]+$", bits[[length(bits) - 1L]])) {
    chrom <- paste(bits[seq_len(length(bits) - 2L)], collapse = "-")
    return(.ic2_peak_span(chrom, bits[[length(bits) - 1L]], bits[[length(bits)]]))
  }
  NULL
}

.ic2_peak_span <- function(chrom, start, end) {
  if (!nzchar(chrom)) return(NULL)
  a <- suppressWarnings(as.integer(start))
  b <- suppressWarnings(as.integer(end))
  if (is.na(a) || is.na(b)) return(NULL)
  if (b < a) {
    tmp <- a
    a <- b
    b <- tmp
  }
  list(chrom, a, b)
}

.ic2_peaks_collect <- function(peaks) {
  if (is.null(peaks)) stop("'peaks' is required.")
  rows <- list()
  add <- function(sample, tab, ids) {
    if (nrow(tab) < 1L) return()
    ok <- !is.na(tab$chromosome) & !is.na(tab$start) & !is.na(tab$end)
    if (!any(ok)) return()
    tab <- tab[ok, , drop = FALSE]
    ids <- ids[ok]
    rows[[length(rows) + 1L]] <<- data.frame(
      sample = sample,
      chromosome = tab$chromosome,
      start = as.integer(tab$start),
      end = as.integer(tab$end),
      peak_id = as.character(ids),
      stringsAsFactors = FALSE
    )
  }

  if (is.data.frame(peaks)) {
    tab <- .ic2_peaks_from_table(peaks)
    sample <- if ("sample" %in% names(tab)) tab$sample else rep("peaks", nrow(tab))
    add_df <- split(tab, sample, drop = TRUE)
    for (nm in names(add_df)) {
      d <- add_df[[nm]]
      add(nm, d, .ic2_peak_ids(d$chromosome, d$start, d$end))
    }
  } else if (is.list(peaks) && !is.character(peaks)) {
    nms <- names(peaks)
    if (is.null(nms) || any(!nzchar(nms))) {
      nms <- paste0("sample", seq_along(peaks))
    }
    for (i in seq_along(peaks)) {
      one <- .ic2_peaks_collect(peaks[[i]])
      if (nrow(one) < 1L) next
      one$sample <- nms[[i]]
      rows[[length(rows) + 1L]] <- one
    }
  } else {
    ids <- as.character(peaks)
    is_file <- !is.na(ids) & nzchar(ids) & file.exists(ids)
    if (any(is_file)) {
      for (path in ids[is_file]) {
        tab <- .ic2_peaks_read_bed(path)
        sample <- .ic2_peak_sample_from_path(path)
        add(
          sample, tab,
          if (!is.null(tab$peak_id)) tab$peak_id else .ic2_peak_ids(tab$chromosome, tab$start, tab$end)
        )
      }
      ids <- ids[!is_file]
    }
    if (length(ids)) {
      parsed <- .ic2_parse_peak_id(ids)
      add("peaks", parsed, ids)
    }
  }
  if (length(rows) < 1L) {
    return(data.frame(
      sample = character(), chromosome = character(),
      start = integer(), end = integer(), peak_id = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.ic2_peak_ids <- function(chrom, start, end) {
  paste(chrom, start, end, sep = "_")
}

.ic2_peak_sample_from_path <- function(path) {
  b <- basename(path)
  b <- sub("\\.(bed|narrowPeak|broadPeak|tsv|txt|csv)(\\.(gz))?$", "", b, ignore.case = TRUE)
  if (!nzchar(b)) "peaks" else b
}

.ic2_peaks_from_table <- function(df) {
  nm <- tolower(names(df))
  names(df) <- nm
  chrom_col <- intersect(c("chromosome", "chrom", "chr", "seqnames", "seqname"), nm)[1L]
  start_col <- intersect(c("start", "chromstart", "start_pos"), nm)[1L]
  end_col <- intersect(c("end", "chromend", "stop", "end_pos"), nm)[1L]
  if (is.na(chrom_col) || is.na(start_col) || is.na(end_col)) {
    stop("A peak table needs chromosome, start, and end columns.")
  }
  sample_col <- intersect(c("sample", "sample_id", "sample.id"), nm)[1L]
  out <- data.frame(
    chromosome = as.character(df[[chrom_col]]),
    start = as.integer(df[[start_col]]),
    end = as.integer(df[[end_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.na(sample_col)) out$sample <- as.character(df[[sample_col]])
  bad <- !is.na(out$start) & !is.na(out$end) & out$end < out$start
  if (any(bad)) {
    tmp <- out$start[bad]
    out$start[bad] <- out$end[bad]
    out$end[bad] <- tmp
  }
  out
}

.ic2_peaks_read_bed <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else path
  lines <- readLines(con, warn = FALSE)
  if (inherits(con, "connection") && isOpen(con)) close(con)
  lines <- lines[!grepl("^\\s*(#|track|browser)", lines)]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 1L) stop("BED file is empty: ", path)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(parts) < 3L)) {
    parts <- strsplit(lines, "\\s+")
  }
  if (any(lengths(parts) < 3L)) stop("BED file needs chromosome, start, and end: ", path)
  chrom <- vapply(parts, `[[`, "", 1L)
  start <- suppressWarnings(as.integer(vapply(parts, `[[`, "", 2L)))
  end <- suppressWarnings(as.integer(vapply(parts, `[[`, "", 3L)))
  ## Header row
  if (is.na(start[[1L]]) || is.na(end[[1L]])) {
    chrom <- chrom[-1L]
    start <- start[-1L]
    end <- end[-1L]
  }
  ## BED is half-open and 0-based. 10x peak names are closed 1-based.
  ## A file of intervals is stored as closed 1-based: start + 1, end unchanged.
  start <- start + 1L
  bad <- !is.na(start) & !is.na(end) & end < start
  if (any(bad)) {
    tmp <- start[bad]
    start[bad] <- end[bad]
    end[bad] <- tmp
  }
  data.frame(
    chromosome = chrom, start = start, end = end,
    peak_id = .ic2_peak_ids(chrom, start, end),
    stringsAsFactors = FALSE
  )
}

#' Reduce one chromosome. gap 0 merges book-ended closed intervals.
#' @keywords internal
.ic2_peaks_reduce_chrom <- function(d, gap) {
  o <- order(d$start, d$end, d$peak_id)
  d <- d[o, , drop = FALSE]
  n <- nrow(d)
  which_k <- integer(n)
  cs <- d$start[[1L]]
  ce <- d$end[[1L]]
  k <- 1L
  starts <- ce_vec <- integer(n)
  starts[[1L]] <- cs
  ce_vec[[1L]] <- ce
  slack <- 1L + gap
  for (i in seq_len(n)) {
    if (d$start[[i]] <= ce + slack) {
      if (d$end[[i]] > ce) ce <- d$end[[i]]
      ce_vec[[k]] <- ce
    } else {
      k <- k + 1L
      cs <- d$start[[i]]
      ce <- d$end[[i]]
      starts[[k]] <- cs
      ce_vec[[k]] <- ce
    }
    which_k[[i]] <- k
  }
  chrom <- d$chromosome[[1L]]
  consensus <- data.frame(
    chromosome = chrom,
    start = starts[seq_len(k)],
    end = ce_vec[seq_len(k)],
    stringsAsFactors = FALSE
  )
  consensus$peak_id <- .ic2_peak_ids(consensus$chromosome, consensus$start, consensus$end)
  consensus$n_peaks <- as.integer(tabulate(which_k, nbins = k))
  samples <- split(d$sample, which_k)
  consensus$n_samples <- vapply(samples, function(s) length(unique(s)), integer(1))
  map <- data.frame(
    sample = d$sample,
    peak_id = d$peak_id,
    chromosome = d$chromosome,
    start = d$start,
    end = d$end,
    consensus_id = consensus$peak_id[which_k],
    stringsAsFactors = FALSE
  )
  list(consensus = consensus, map = map)
}
