#' CLR-normalize ADT counts
#'
#' For each cell, divide each antibody by the geometric mean of \code{log1p}
#' counts across antibodies, then \code{log1p}. Streams cell chunks from the
#' on-disk counts assay and writes \code{modalities/adt/assays/data}.
#' RNA is not touched.
#'
#' @param x An iCellR2 object.
#' @param assay.name Source ADT assay. Default \code{"counts"}.
#' @param new.assay.name Destination assay. Default \code{"data"}.
#' @param chunk.cells Cells per read.
#' @param overwrite Replace an existing normalized assay.
#' @return The iCellR2 object.
#' @export
norm.adt <- function(
    x = NULL,
    assay.name = "counts",
    new.assay.name = "data",
    chunk.cells = 20000L,
    overwrite = TRUE
) {
  if (!inherits(x, "iCellR2")) {
    stop("x should be an object of class iCellR2")
  }
  if (is.null(x@modalities$adt$assays[[assay.name]])) {
    stop("No ADT assay '", assay.name, "'. Run add.adt() first.")
  }
  handles <- .ic2_get_assay_handles(x, "adt", assay.name)
  pieces <- list()
  chunk.cells <- as.integer(chunk.cells)
  for (h in handles) {
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    n <- h$n_cells
    starts <- seq.int(1L, n, by = chunk.cells)
    for (s in starts) {
      e <- min(s + chunk.cells - 1L, n)
      block <- .ic2_read_cols(h, s:e)
      if (is.null(rownames(block))) rownames(block) <- h$feature_ids
      if (is.null(colnames(block))) colnames(block) <- h$cell_ids[s:e]
      pieces[[length(pieces) + 1L]] <- .ic2_adt_clr(block)
    }
  }
  mat <- if (length(pieces) == 1L) pieces[[1L]] else do.call(cbind, pieces)
  cell_ids <- .ic2_adt_cell_ids(x)
  if (!identical(colnames(mat), cell_ids)) {
    if (!all(cell_ids %in% colnames(mat))) {
      stop("Normalized ADT columns do not cover cell_data.")
    }
    mat <- mat[, cell_ids, drop = FALSE]
  }
  dest <- .ic2_soda_assay_path(x, "adt", new.assay.name, ext = ".h5")
  handle <- write_sparse_hdf5(mat, filepath = dest, group = "counts", overwrite = overwrite)
  handle$modality <- "adt"
  x <- .ic2_assign_assay(
    x, modality = "adt", assay.name = new.assay.name,
    handles = list(handle), virtual = FALSE, set.active = TRUE
  )
  x@params$adt$normalization <- list(
    method = "clr",
    assay.name = assay.name,
    new.assay.name = new.assay.name,
    n_features = nrow(mat),
    n_cells = ncol(mat)
  )
  x <- .ic2_log_command(
    x, "norm.adt",
    list(assay.name = assay.name, new.assay.name = new.assay.name, method = "clr"),
    list(data = dest)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}
