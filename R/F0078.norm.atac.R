#' TF-IDF of an ATAC peak assay
#'
#' Signac method 1, streamed by cell chunk. For each peak
#' \code{idf = log(1 + n_cells / total_counts)} and each cell is divided by
#' its own total, then multiplied by \code{idf}. Zeros stay implicit. The
#' result is \code{modalities/atac/assays/tfidf}. Do not log this assay again
#' in \code{\link{run.pca}}: \code{run.pca(..., modality = "atac", assay.name = "tfidf")}
#' already turns scaling off and skips the RNA log transform.
#'
#' @param x An iCellR2 object.
#' @param assay.name Source assay. Default \code{"counts"}.
#' @param new.assay.name Output assay. Default \code{"tfidf"}.
#' @param scale.factor Multiplier after TF-IDF. Default \code{1}.
#' @param chunk.size Cells per read.
#' @param overwrite Replace an existing TF-IDF file.
#' @param verbose Progress style.
#' @return The iCellR2 object.
#'
#' @examples
#' \dontrun{
#' obj <- norm.atac(obj)
#' obj <- run.pca(obj, modality = "atac", assay.name = "tfidf", npcs = 50)
#' obj <- run.umap(obj, modality = "atac", dims = 2:30)
#' }
#'
#' @export
norm.atac <- function(
    x = NULL,
    assay.name = "counts",
    new.assay.name = "tfidf",
    scale.factor = 1,
    chunk.size = 20000L,
    overwrite = TRUE,
    verbose = TRUE
) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  if (is.null(x@modalities$atac$assays[[assay.name]])) {
    stop("No ATAC assay '", assay.name, "'. Run add.atac() first.")
  }
  scale.factor <- as.numeric(scale.factor)[1L]
  if (!is.finite(scale.factor) || scale.factor <= 0) stop("'scale.factor' must be > 0.")
  chunk.size <- as.integer(chunk.size)
  rep <- .ic2_new_reporter("norm.atac", verbose)
  handles <- .ic2_get_assay_handles(x, "atac", assay.name)
  handles <- lapply(handles, .ic2_enrich_handle, load_dimnames = TRUE)
  feat <- as.character(handles[[1L]]$feature_ids)
  p <- length(feat)
  rs <- numeric(p)
  n_cells <- 0L
  .ic2_step(rep, "Pass 1: peak totals")
  for (h in handles) {
    starts <- seq.int(1L, h$n_cells, by = chunk.size)
    for (s in starts) {
      e <- min(s + chunk.size - 1L, h$n_cells)
      block <- .ic2_to_dgCMatrix(.ic2_read_cols(h, s:e))
      rs <- rs + Matrix::rowSums(block)
      n_cells <- n_cells + ncol(block)
    }
  }
  if (n_cells < 1L) stop("ATAC assay has no cells.")
  idf <- numeric(p)
  pos <- rs > 0
  idf[pos] <- log1p(n_cells / rs[pos])

  cell_ids <- .ic2_adt_cell_ids(x)
  pieces <- list()
  .ic2_step(rep, "Pass 2: TF-IDF")
  for (h in handles) {
    starts <- seq.int(1L, h$n_cells, by = chunk.size)
    for (s in starts) {
      e <- min(s + chunk.size - 1L, h$n_cells)
      block <- .ic2_atac_tfidf_block(.ic2_to_dgCMatrix(.ic2_read_cols(h, s:e)), idf, scale.factor)
      colnames(block) <- h$cell_ids[s:e]
      pieces[[length(pieces) + 1L]] <- block
    }
  }
  mat <- do.call(cbind, pieces)
  if (!is.null(colnames(mat)) && all(cell_ids %in% colnames(mat))) {
    mat <- mat[, cell_ids, drop = FALSE]
  }
  rownames(mat) <- feat
  dest <- .ic2_soda_assay_path(x, "atac", new.assay.name, ext = ".h5")
  handle <- write_sparse_hdf5(mat, filepath = dest, group = "counts", overwrite = overwrite)
  handle$modality <- "atac"
  x <- .ic2_assign_assay(
    x, modality = "atac", assay.name = new.assay.name,
    handles = list(handle), virtual = FALSE, set.active = TRUE
  )
  x@params$atac$tfidf <- list(
    source = assay.name, scale.factor = scale.factor, n_cells = n_cells
  )
  x <- .ic2_log_command(
    x, "norm.atac",
    list(assay.name = assay.name, new.assay.name = new.assay.name, scale.factor = scale.factor),
    list(tfidf = dest)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  .ic2_finish(rep, extra_lines = paste0("Wrote ", new.assay.name, "   |   peaks: ", nrow(mat)))
  x
}

#' TF * log(1 + N / peak total) on one cell chunk. Zeros stay implicit.
#' @keywords internal
.ic2_atac_tfidf_block <- function(mat, idf, scale.factor = 1) {
  mat <- .ic2_to_dgCMatrix(mat)
  cs <- Matrix::colSums(mat)
  p <- mat@p
  for (j in seq_len(ncol(mat))) {
    a <- p[[j]] + 1L
    b <- p[[j + 1L]]
    if (a > b || !is.finite(cs[[j]]) || cs[[j]] <= 0) next
    rows <- mat@i[a:b] + 1L
    mat@x[a:b] <- (mat@x[a:b] / cs[[j]]) * idf[rows] * scale.factor
  }
  mat
}
