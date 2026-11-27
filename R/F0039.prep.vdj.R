#' Prepare V(D)J contigs for an RNA object
#'
#' Accepts a 10x \code{filtered_contig_annotations.csv} (or TSV) or a
#' data.frame, including AIRR names (\code{cell_id}, \code{locus},
#' \code{junction_aa}, \code{clone_id}). Barcodes are not rewritten.
#' Clonotype frequency is the number of cells, not contig rows divided by 2.
#'
#' @param vdj.data Data.frame, or a path to CSV/TSV.
#' @param cond.name Optional sample id used only to match prefixed RNA
#'   barcodes (\code{sample_barcode}). It is not pasted onto the barcode.
#' @param productive.only Keep productive contigs when that column exists.
#' @return The filtered contig table, with a \code{cells} attribute
#'   (one row per barcode).
#' @examples
#' \dontrun{
#' vdj <- prep.vdj("filtered_contig_annotations.csv", cond.name = "pbmc")
#' }
#' @export
prep.vdj <- function(vdj.data = NULL, cond.name = NULL, productive.only = TRUE) {
  df <- .ic2_vdj_read(vdj.data)
  df <- .ic2_vdj_alias_columns(df)
  if (!"barcode" %in% colnames(df)) stop("VDJ table needs a barcode (or cell_id) column.")
  idcol <- if ("raw_clonotype_id" %in% colnames(df)) "raw_clonotype_id" else "clonotype_id"
  if (!idcol %in% colnames(df)) stop("VDJ table needs raw_clonotype_id or clonotype_id (or clone_id).")

  df$barcode <- as.character(df$barcode)
  df <- df[!is.na(df$barcode) & nzchar(df$barcode), , drop = FALSE]
  if ("is_cell" %in% colnames(df)) {
    ok <- tolower(as.character(df$is_cell)) %in% c("true", "t", "yes", "1")
    if (any(ok)) df <- df[ok, , drop = FALSE]
  }
  if (isTRUE(productive.only) && "productive" %in% colnames(df)) {
    ok <- tolower(as.character(df$productive)) %in% c("true", "t", "yes", "1")
    if (any(ok)) df <- df[ok, , drop = FALSE]
  }
  clon <- as.character(df[[idcol]])
  bad <- is.na(clon) | !nzchar(clon) | tolower(clon) %in% c("none", "nan", "na")
  df <- df[!bad, , drop = FALSE]
  if (nrow(df) < 1L) stop("No productive contigs with a clonotype id.")
  if (!identical(idcol, "raw_clonotype_id")) df$raw_clonotype_id <- as.character(df[[idcol]])
  df$raw_clonotype_id <- as.character(df$raw_clonotype_id)

  if (!"umis" %in% colnames(df)) df$umis <- 1L
  df$umis <- suppressWarnings(as.integer(df$umis))
  df$umis[is.na(df$umis)] <- 1L
  if (!"chain" %in% colnames(df)) df$chain <- NA_character_
  if (!"cdr3" %in% colnames(df)) df$cdr3 <- NA_character_
  df$chain <- as.character(df$chain)
  df$cdr3 <- as.character(df$cdr3)

  pair <- unique(df[, c("barcode", "raw_clonotype_id"), drop = FALSE])
  freq <- table(pair$raw_clonotype_id)
  n_cells <- length(unique(pair$barcode))
  df$clonotype_freq <- as.integer(freq[df$raw_clonotype_id])
  df$proportion <- df$clonotype_freq / n_cells

  ord <- order(df$barcode, -df$umis, method = "radix")
  df <- df[ord, , drop = FALSE]
  first <- !duplicated(df$barcode)
  cells <- df[first, c("barcode", "raw_clonotype_id", "clonotype_freq", "proportion"), drop = FALSE]
  colnames(cells)[2] <- "clonotype_id"

  sub <- df[!duplicated(paste(df$barcode, df$chain, sep = "\r")), , drop = FALSE]
  sub <- sub[nzchar(sub$chain) & !is.na(sub$chain), , drop = FALSE]
  if (nrow(sub) > 0L) {
    lab <- paste0(sub$chain, ":", ifelse(is.na(sub$cdr3), "", sub$cdr3))
    cdr <- tapply(lab, sub$barcode, paste, collapse = ";")
    ch <- tapply(sub$chain, sub$barcode, function(z) paste(unique(z), collapse = ";"))
    cells$vdj_chain <- unname(ch[cells$barcode])
    cells$vdj_cdr3 <- unname(cdr[cells$barcode])
  } else {
    cells$vdj_chain <- NA_character_
    cells$vdj_cdr3 <- NA_character_
  }
  cells$vdj_paired <- .ic2_vdj_paired(cells$vdj_chain)
  cells$vdj_expand <- .ic2_vdj_expand(cells$clonotype_freq)
  if (!is.null(cond.name) && length(cond.name) == 1L && nzchar(as.character(cond.name)) &&
      !identical(as.character(cond.name), "NULL")) {
    cells$sample_id <- as.character(cond.name)
    attr(df, "sample_id") <- as.character(cond.name)
  }
  rownames(cells) <- NULL
  attr(df, "cells") <- cells
  attr(df, "n_cells") <- n_cells
  df
}

#' @keywords internal
.ic2_vdj_read <- function(vdj.data) {
  if (is.null(vdj.data) || identical(vdj.data, "data.frame")) {
    stop("Provide vdj.data: a data.frame or a CSV/TSV path.")
  }
  if (is.character(vdj.data) && length(vdj.data) == 1L) {
    if (!file.exists(vdj.data)) stop("VDJ file not found: ", vdj.data)
    line <- readLines(vdj.data, n = 1L, warn = FALSE)
    sep <- if (grepl("\t", line, fixed = TRUE)) "\t" else ","
    hdr <- strsplit(gsub("\"", "", line, fixed = TRUE), sep, fixed = TRUE)[[1]]
    want <- c(
      "barcode", "cell_id", "chain", "locus", "productive",
      "raw_clonotype_id", "clonotype_id", "clone_id",
      "cdr3", "junction_aa", "cdr3_nt", "junction",
      "v_gene", "v_call", "j_gene", "j_call", "c_gene", "c_call",
      "umis", "duplicate_count", "is_cell", "high_confidence"
    )
    classes <- rep("NULL", length(hdr))
    classes[hdr %in% want] <- "character"
    if ("umis" %in% hdr) classes[hdr == "umis"] <- "integer"
    if ("duplicate_count" %in% hdr) classes[hdr == "duplicate_count"] <- "integer"
    return(utils::read.table(
      vdj.data, header = TRUE, sep = sep, quote = "\"", comment.char = "",
      colClasses = classes, stringsAsFactors = FALSE, check.names = FALSE,
      fill = TRUE, numerals = "no.loss"
    ))
  }
  if (is.data.frame(vdj.data)) return(vdj.data)
  stop("vdj.data must be a data.frame or a CSV/TSV path.")
}

#' @keywords internal
.ic2_vdj_alias_columns <- function(df) {
  rename <- c(
    cell_id = "barcode", locus = "chain", junction_aa = "cdr3",
    junction = "cdr3_nt", clone_id = "raw_clonotype_id",
    v_call = "v_gene", j_call = "j_gene", c_call = "c_gene",
    duplicate_count = "umis"
  )
  for (src in names(rename)) {
    dest <- unname(rename[[src]])
    if (src %in% colnames(df) && !dest %in% colnames(df)) {
      colnames(df)[colnames(df) == src] <- dest
    }
  }
  df
}

#' @keywords internal
.ic2_vdj_paired <- function(chain) {
  chain <- as.character(chain)
  parts <- strsplit(chain, ";", fixed = TRUE)
  vapply(parts, function(z) {
    z <- unique(z)
    (all(c("TRA", "TRB") %in% z)) ||
      ("IGH" %in% z && any(c("IGK", "IGL") %in% z))
  }, logical(1))
}

#' @keywords internal
.ic2_vdj_expand <- function(freq) {
  freq <- as.integer(freq)
  out <- rep("Mid", length(freq))
  out[is.na(freq)] <- NA_character_
  out[freq <= 1L] <- "singleton"
  out[freq == 2L] <- "doubleton"
  big <- !is.na(freq) & freq > 2L
  if (any(big)) {
    u <- unique(freq[big])
    q1 <- as.numeric(stats::quantile(u, 0.25, names = FALSE))
    q3 <- as.numeric(stats::quantile(u, 0.75, names = FALSE))
    out[big & freq < q1] <- "Q1"
    out[big & freq >= q1 & freq <= q3] <- "Mid"
    out[big & freq > q3] <- "Q3"
  }
  out
}

#' Write a data.frame as TSV without inventing a cell_id column from rownames.
#' @keywords internal
.ic2_vdj_write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  attributes(out) <- attributes(out)[c("names", "row.names", "class")]
  utils::write.table(out, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  invisible(path)
}

#' Match VDJ barcodes onto RNA cell ids. Exact, then sample_barcode, then a unique suffix.
#' @keywords internal
.ic2_vdj_match_cells <- function(barcodes, cell_ids, sample.id = NULL) {
  m <- match(barcodes, cell_ids)
  if (!is.null(sample.id) && length(sample.id) == 1L && nzchar(sample.id)) {
    m2 <- match(paste0(sample.id, "_", barcodes), cell_ids)
    m[is.na(m)] <- m2[is.na(m)]
  }
  if (anyNA(m)) {
    suf <- sub("^[^_]+_", "", cell_ids)
    tab <- table(suf, useNA = "no")
    uniq <- names(tab)[tab == 1L]
    suf[!suf %in% uniq] <- NA_character_
    m3 <- match(barcodes, suf)
    m[is.na(m)] <- m3[is.na(m)]
  }
  m
}
