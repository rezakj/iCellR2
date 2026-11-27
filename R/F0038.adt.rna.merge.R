#' Join RNA PCA and ADT PCA under multimodal/rna_adt
#'
#' Does not paste antibodies into the RNA assay and does not replace
#' \code{reductions/pca}. RNA principal components (already computed) are
#' concatenated with a small PCA of the CLR ADT matrix. Each block is
#' scaled so the two modalities contribute equal total variance.
#' The ADT PCA is a covariance eigen-decomposition on the antibody axis
#' (two passes over ADT chunks), not a second RNA PCA.
#'
#' @param x An iCellR2 object.
#' @param adt.data Ignored. Kept so classic calls do not error. The join
#'   always uses the normalized ADT assay (\code{data}).
#' @param rna.modality RNA modality. Default \code{"rna"}.
#' @param rna.dims Number of RNA PCs to include.
#' @param adt.dims Number of ADT PCs to include.
#' @param adt.assay Normalized ADT assay. Default \code{"data"}.
#' @param chunk.cells Cells per ADT read.
#' @return The iCellR2 object, with \code{@@multimodal$rna_adt} set.
#' @export
adt.rna.merge <- function(
    x = NULL,
    adt.data = "data",
    rna.modality = "rna",
    rna.dims = 30L,
    adt.dims = 10L,
    adt.assay = "data",
    chunk.cells = 20000L
) {
  if (!inherits(x, "iCellR2")) {
    stop("x should be an object of class iCellR2")
  }
  if (is.null(x@modalities$adt$assays[[adt.assay]])) {
    stop("No ADT assay '", adt.assay, "'. Run norm.adt() first.")
  }
  rna_path <- x@modalities[[rna.modality]]$reductions$pca$cell_embeddings
  if (is.null(rna_path) || !file.exists(rna_path)) {
    alt <- file.path(x@storage$root_dir, "reductions", "pca", "cell_embeddings.h5")
    if (file.exists(alt)) rna_path <- alt else {
      stop("RNA PCA is missing. Run run.pca() before adt.rna.merge().")
    }
  }
  rna <- .ic2_h5_read_dense(rna_path)
  cell_ids <- .ic2_adt_cell_ids(x)
  if (is.null(rownames(rna)) || !all(cell_ids %in% rownames(rna))) {
    stop("RNA PCA embeddings are not aligned to cell_data.")
  }
  rna <- rna[cell_ids, , drop = FALSE]
  rna.dims <- min(as.integer(rna.dims), ncol(rna))
  rna <- rna[, seq_len(rna.dims), drop = FALSE]

  adt_pcs <- .ic2_adt_pca_stream(
    x, assay.name = adt.assay, n_pcs = as.integer(adt.dims),
    chunk.cells = as.integer(chunk.cells), cell_ids = cell_ids
  )
  joint <- cbind(
    .ic2_adt_block_scale(rna),
    .ic2_adt_block_scale(adt_pcs$scores)
  )
  rownames(joint) <- cell_ids
  colnames(joint) <- c(
    paste0("RNA_", seq_len(ncol(rna))),
    paste0("ADT_", seq_len(ncol(adt_pcs$scores)))
  )

  root <- x@storage$root_dir
  red_dir <- file.path(root, "multimodal", "rna_adt", "joint_pca")
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  .ic2_h5_write_dense(joint, emb_path)
  info <- list(
    name = "rna_adt",
    method = "scaled_pca_concat",
    rna_modality = rna.modality,
    rna_reduction = "pca",
    rna_dims = rna.dims,
    adt_assay = adt.assay,
    adt_dims = ncol(adt_pcs$scores),
    embeddings = emb_path,
    created = as.character(Sys.time())
  )
  meta_path <- file.path(root, "multimodal", "rna_adt", "multimodal.json")
  .ic2_soda_write_json(info, meta_path)
  if (is.null(x@multimodal)) x@multimodal <- .ic2_empty_multimodal()
  x@multimodal$rna_adt <- info
  x <- .ic2_log_command(
    x, "adt.rna.merge",
    list(rna.dims = rna.dims, adt.dims = ncol(adt_pcs$scores), adt.assay = adt.assay),
    list(multimodal = meta_path, joint_pca = emb_path)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}

#' Two-pass PCA of a short ADT assay (features x cells covariance).
#' @keywords internal
.ic2_adt_pca_stream <- function(x, assay.name, n_pcs, chunk.cells, cell_ids) {
  handles <- .ic2_get_assay_handles(x, "adt", assay.name)
  handles <- lapply(handles, .ic2_enrich_handle, load_dimnames = TRUE)
  feat <- handles[[1L]]$feature_ids
  p <- length(feat)
  n <- length(cell_ids)
  sums <- numeric(p)
  seen <- 0L
  for (h in handles) {
    starts <- seq.int(1L, h$n_cells, by = chunk.cells)
    for (s in starts) {
      e <- min(s + chunk.cells - 1L, h$n_cells)
      block <- as.matrix(.ic2_read_cols(h, s:e))
      sums <- sums + rowSums(block)
      seen <- seen + ncol(block)
    }
  }
  if (seen < 2L) stop("ADT assay has fewer than 2 cells.")
  means <- sums / seen
  covm <- matrix(0, p, p)
  for (h in handles) {
    starts <- seq.int(1L, h$n_cells, by = chunk.cells)
    for (s in starts) {
      e <- min(s + chunk.cells - 1L, h$n_cells)
      block <- as.matrix(.ic2_read_cols(h, s:e))
      block <- block - means
      covm <- covm + tcrossprod(block)
    }
  }
  covm <- covm / (seen - 1L)
  eig <- eigen(covm, symmetric = TRUE)
  k <- min(as.integer(n_pcs), p, seen - 1L)
  k <- max(k, 1L)
  loadings <- eig$vectors[, seq_len(k), drop = FALSE]
  scores <- matrix(NA_real_, n, k)
  rownames(scores) <- cell_ids
  colnames(scores) <- paste0("ADT_PC", seq_len(k))
  for (h in handles) {
    starts <- seq.int(1L, h$n_cells, by = chunk.cells)
    for (s in starts) {
      e <- min(s + chunk.cells - 1L, h$n_cells)
      ids <- h$cell_ids[s:e]
      block <- as.matrix(.ic2_read_cols(h, s:e))
      block <- block - means
      sc <- t(crossprod(loadings, block))
      scores[ids, ] <- sc
    }
  }
  if (anyNA(scores)) {
    stop("ADT PCA could not be aligned to every cell_data id.")
  }
  list(scores = scores, loadings = loadings)
}

#' Scale a cells-by-dims block to unit column sd, then equalize block energy.
#' @keywords internal
.ic2_adt_block_scale <- function(m) {
  m <- scale(m)
  m[is.nan(m)] <- 0
  m / sqrt(ncol(m))
}
