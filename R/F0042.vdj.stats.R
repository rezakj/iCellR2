#' Clonotype expansion summary for TCR and BCR
#'
#' Counts cells (not contig rows) in singleton, doubleton, and quartile bins,
#' separately for each chain and for paired TRA+TRB or IGH+light. Writes
#' \code{modalities/vdj/stats.tsv} when given an iCellR2 object.
#'
#' @param x An iCellR2 object after \code{\link{add.vdj}}, or a contig data.frame.
#' @param my.vdj Alias of a contig data.frame when \code{x} is missing.
#' @param plot If \code{TRUE}, draw a bar of expansion bins (not one point per CDR3).
#' @param out.name PNG path when \code{plot = TRUE}.
#' @return The iCellR2 object when \code{x} is one (stats path on
#'   \code{modalities$vdj$repertoire}). A data.frame when \code{x} is a table.
#' @examples
#' \dontrun{
#' obj <- vdj.stats(obj, plot = TRUE, out.name = "vdj_stats.png")
#' }
#' @export
vdj.stats <- function(x = NULL, my.vdj = NULL, plot = FALSE, out.name = NULL) {
  src <- NULL
  root <- NULL
  if (inherits(x, "iCellR2")) {
    root <- x@storage$root_dir
    path <- x@modalities$vdj$cells %||% file.path(root, "modalities", "vdj", "cells.tsv")
    if (!is.character(path) || !file.exists(path)) {
      stop("No VDJ cell table. Run add.vdj() first.")
    }
    cells <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    contigs <- NULL
    cpath <- x@modalities$vdj$contigs
    if (is.character(cpath) && length(cpath) == 1L && file.exists(cpath)) {
      contigs <- utils::read.delim(cpath, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    }
    src <- "object"
  } else {
    tab <- if (!is.null(my.vdj)) my.vdj else x
    contigs <- prep.vdj(tab)
    cells <- attr(contigs, "cells")
  }
  if (!"vdj_chain" %in% colnames(cells)) stop("VDJ cell table is missing vdj_chain.")
  stats <- .ic2_vdj_chain_stats(cells, contigs)
  if (!is.null(root)) {
    stats_path <- .ic2_vdj_write_tsv(stats, file.path(root, "modalities", "vdj", "stats.tsv"))
    if (inherits(x, "iCellR2")) {
      x@modalities$vdj$repertoire <- stats_path
      try(.ic2_log_command(x, "vdj.stats", list(n_chains = nrow(stats)), list(stats = stats_path)), silent = TRUE)
      try(.ic2_soda_write_manifest(x), silent = TRUE)
    }
  }
  if (isTRUE(plot)) {
    bins <- c("singleton", "doubleton", "Q1", "Mid", "Q3")
    long <- do.call(rbind, lapply(bins, function(b) {
      data.frame(chain = stats$chain, bin = b, n_cells = stats[[b]], stringsAsFactors = FALSE)
    }))
    long$bin <- factor(long$bin, levels = bins)
    p <- ggplot2::ggplot(long, ggplot2::aes(x = chain, y = n_cells, fill = bin)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::labs(x = NULL, y = "cells", title = "Clonotype expansion") +
      ggplot2::theme_bw()
    if (!is.null(out.name)) ggplot2::ggsave(out.name, p, width = 8, height = 4, dpi = 150)
  }
  if (inherits(x, "iCellR2")) return(x)
  stats
}

#' @keywords internal
.ic2_vdj_chain_stats <- function(cells, contigs = NULL) {
  bins <- c("singleton", "doubleton", "Q1", "Mid", "Q3")
  chains <- c("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")
  rows <- list()
  add_row <- function(chain, clonotype_id, expand) {
    if (length(clonotype_id) < 1L) return()
    tab <- table(factor(expand, levels = bins))
    rows[[length(rows) + 1L]] <<- data.frame(
      chain = chain,
      n_cells = length(clonotype_id),
      n_clonotypes = length(unique(clonotype_id)),
      singleton = as.integer(tab[["singleton"]]),
      doubleton = as.integer(tab[["doubleton"]]),
      Q1 = as.integer(tab[["Q1"]]),
      Mid = as.integer(tab[["Mid"]]),
      Q3 = as.integer(tab[["Q3"]]),
      stringsAsFactors = FALSE
    )
  }
  present <- chains[vapply(chains, function(ch) any(grepl(ch, cells$vdj_chain, fixed = TRUE)), logical(1))]
  if (!is.null(contigs) && "chain" %in% colnames(contigs)) {
    for (ch in intersect(chains, unique(as.character(contigs$chain)))) {
      bc <- unique(contigs$barcode[contigs$chain == ch])
      sub <- cells[cells$barcode %in% bc, , drop = FALSE]
      add_row(ch, sub$clonotype_id, sub$vdj_expand)
    }
  } else {
    for (ch in present) {
      sub <- cells[grepl(ch, cells$vdj_chain, fixed = TRUE), , drop = FALSE]
      add_row(ch, sub$clonotype_id, sub$vdj_expand)
    }
  }
  tcr <- grepl("TRA", cells$vdj_chain, fixed = TRUE) | grepl("TRB", cells$vdj_chain, fixed = TRUE)
  bcr <- grepl("IGH", cells$vdj_chain, fixed = TRUE)
  if (any(tcr, na.rm = TRUE)) add_row("TCR", cells$clonotype_id[tcr], cells$vdj_expand[tcr])
  if (any(bcr, na.rm = TRUE)) add_row("BCR", cells$clonotype_id[bcr], cells$vdj_expand[bcr])
  paired_flag <- tolower(as.character(cells$vdj_paired)) %in% c("true", "t", "1")
  paired <- cells[paired_flag, , drop = FALSE]
  if (nrow(paired) > 0L) add_row("paired", paired$clonotype_id, paired$vdj_expand)
  if (length(rows) < 1L) {
    return(data.frame(
      chain = character(), n_cells = integer(), n_clonotypes = integer(),
      singleton = integer(), doubleton = integer(), Q1 = integer(),
      Mid = integer(), Q3 = integer()
    ))
  }
  do.call(rbind, rows)
}
