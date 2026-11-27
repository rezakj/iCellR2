#' Add V(D)J contigs onto an iCellR2 RNA object
#'
#' Writes \code{modalities/vdj/contigs.tsv} and \code{cells.tsv}, and adds
#' one clonotype per cell to \code{cell_data}. RNA assays and reductions
#' are not modified. Contig tables are columnar (not a genes-by-cells matrix),
#' so a million cells stays a few million short rows.
#'
#' @param x An iCellR2 object.
#' @param input CSV/TSV path or data.frame. \code{vdj.data} is the same argument.
#' @param vdj.data Alias of \code{input}.
#' @param cond.name Optional sample id for matching \code{sample_barcode} cell ids.
#' @param productive.only Drop non-productive contigs when that column exists.
#' @return The iCellR2 object.
#' @examples
#' \dontrun{
#' obj <- add.vdj(obj, input = "filtered_contig_annotations.csv", cond.name = "pbmc")
#' }
#' @export
add.vdj <- function(
    x = NULL,
    input = NULL,
    vdj.data = NULL,
    cond.name = NULL,
    productive.only = TRUE
) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  src <- if (!is.null(input)) input else vdj.data
  contigs <- prep.vdj(src, cond.name = cond.name, productive.only = productive.only)
  cells <- attr(contigs, "cells")
  sample.id <- attr(contigs, "sample_id") %||% cond.name
  cell_ids <- .ic2_adt_cell_ids(x)
  hit <- .ic2_vdj_match_cells(cells$barcode, cell_ids, sample.id)
  n_hit <- sum(!is.na(hit))
  if (n_hit < 1L) stop("No VDJ barcodes matched cell_data.")
  if (anyNA(hit)) {
    message("add.vdj: ", sum(is.na(hit)), " VDJ barcodes did not match an RNA cell.")
  }

  root <- x@storage$root_dir
  vdir <- file.path(root, "modalities", "vdj")
  dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
  contig_path <- .ic2_vdj_write_tsv(contigs, file.path(vdir, "contigs.tsv"))
  cell_path <- .ic2_vdj_write_tsv(cells, file.path(vdir, "cells.tsv"))

  cd <- x@cell_data
  cd$clonotype_id <- NA_character_
  cd$vdj_chain <- NA_character_
  cd$vdj_cdr3 <- NA_character_
  cd$vdj_freq <- NA_integer_
  cd$vdj_prop <- NA_real_
  cd$vdj_expand <- NA_character_
  cd$vdj_paired <- NA
  ok <- !is.na(hit)
  cd$clonotype_id[hit[ok]] <- cells$clonotype_id[ok]
  cd$vdj_chain[hit[ok]] <- cells$vdj_chain[ok]
  cd$vdj_cdr3[hit[ok]] <- cells$vdj_cdr3[ok]
  cd$vdj_freq[hit[ok]] <- cells$clonotype_freq[ok]
  cd$vdj_prop[hit[ok]] <- cells$proportion[ok]
  cd$vdj_expand[hit[ok]] <- cells$vdj_expand[ok]
  cd$vdj_paired[hit[ok]] <- cells$vdj_paired[ok]
  x@cell_data <- cd
  x <- .ic2_persist_cell_data(x)

  if (is.null(x@modalities$vdj) || !is.list(x@modalities$vdj)) {
    x@modalities$vdj <- list()
  }
  x@modalities$vdj$contigs <- contig_path
  x@modalities$vdj$cells <- cell_path
  x@modalities$vdj$clonotypes <- cell_path
  x@storage$modality_paths$vdj <- vdir
  x@params$vdj <- list(
    n_contigs = nrow(contigs),
    n_vdj_cells = n_hit,
    n_clonotypes = length(unique(cells$clonotype_id[ok]))
  )
  x <- .ic2_log_command(
    x, "add.vdj",
    list(n_contigs = nrow(contigs), n_vdj_cells = n_hit),
    list(contigs = contig_path, cells = cell_path)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}
