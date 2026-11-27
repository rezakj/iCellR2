#' Join RNA PCA and ATAC PCA under multimodal/rna_atac
#'
#' Does not paste peaks into the RNA assay and does not replace
#' \code{reductions/pca}. The two embeddings are scaled so each block has
#' equal total variance, then concatenated. ATAC components start at 2
#' because the first component usually tracks depth.
#'
#' Run \code{\link{run.pca}} on RNA and
#' \code{run.pca(x, modality = "atac", assay.name = "tfidf")} first.
#' The ATAC embedding is read from \code{modalities/atac/reductions/pca}.
#'
#' @param x An iCellR2 object.
#' @param rna.dims RNA components to keep. Default \code{1:30}.
#' @param atac.dims ATAC components to keep. Default \code{2:30}.
#' @param rna.modality RNA modality. Default \code{"rna"}.
#' @return The iCellR2 object, with \code{@@multimodal$rna_atac} set.
#'
#' @examples
#' \dontrun{
#' obj <- atac.rna.merge(obj, rna.dims = 1:30, atac.dims = 2:30)
#' }
#'
#' @export
atac.rna.merge <- function(
    x = NULL,
    rna.dims = 1:30,
    atac.dims = 2:30,
    rna.modality = "rna"
) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  rna_path <- x@modalities[[rna.modality]]$reductions$pca$cell_embeddings
  if (is.null(rna_path) || !file.exists(rna_path)) {
    alt <- file.path(x@storage$root_dir, "reductions", "pca", "cell_embeddings.h5")
    if (file.exists(alt)) rna_path <- alt else {
      stop("RNA PCA is missing. Run run.pca() before atac.rna.merge().")
    }
  }
  atac_path <- x@modalities$atac$reductions$pca$cell_embeddings
  if (is.null(atac_path) || !file.exists(atac_path)) {
    alt <- file.path(x@storage$root_dir, "modalities", "atac", "reductions", "pca", "cell_embeddings.h5")
    if (file.exists(alt)) atac_path <- alt else {
      stop(
        "ATAC PCA is missing. Run run.pca(x, modality = \"atac\", assay.name = \"tfidf\") first."
      )
    }
  }
  if (identical(normalizePath(rna_path, winslash = "/", mustWork = FALSE),
                normalizePath(atac_path, winslash = "/", mustWork = FALSE))) {
    stop("ATAC PCA points at the RNA embedding. Re-run ATAC PCA so it is stored under modalities/atac.")
  }

  cell_ids <- .ic2_adt_cell_ids(x)
  rna <- .ic2_h5_read_dense(rna_path)
  atac <- .ic2_h5_read_dense(atac_path)
  if (is.null(rownames(rna)) || !all(cell_ids %in% rownames(rna))) {
    stop("RNA PCA embeddings are not aligned to cell_data.")
  }
  if (is.null(rownames(atac)) || !all(cell_ids %in% rownames(atac))) {
    stop("ATAC PCA embeddings are not aligned to cell_data.")
  }
  rna <- rna[cell_ids, , drop = FALSE]
  atac <- atac[cell_ids, , drop = FALSE]
  rna.dims <- .ic2_atac_dim_index(rna.dims, ncol(rna))
  atac.dims <- .ic2_atac_dim_index(atac.dims, ncol(atac))
  rna <- rna[, rna.dims, drop = FALSE]
  atac <- atac[, atac.dims, drop = FALSE]
  joint <- cbind(.ic2_adt_block_scale(rna), .ic2_adt_block_scale(atac))
  rownames(joint) <- cell_ids
  colnames(joint) <- c(
    paste0("RNA_", seq_len(ncol(rna))),
    paste0("ATAC_", seq_len(ncol(atac)))
  )

  root <- x@storage$root_dir
  red_dir <- file.path(root, "multimodal", "rna_atac", "joint_pca")
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  .ic2_h5_write_dense(joint, emb_path)
  info <- list(
    name = "rna_atac",
    method = "scaled_pca_concat",
    rna_modality = rna.modality,
    rna_reduction = "pca",
    rna_dims = as.integer(rna.dims),
    atac_dims = as.integer(atac.dims),
    embeddings = emb_path,
    created = as.character(Sys.time())
  )
  meta_path <- file.path(root, "multimodal", "rna_atac", "multimodal.json")
  .ic2_soda_write_json(info, meta_path)
  if (is.null(x@multimodal)) x@multimodal <- .ic2_empty_multimodal()
  x@multimodal$rna_atac <- info
  x <- .ic2_log_command(
    x, "atac.rna.merge",
    list(rna.dims = as.integer(rna.dims), atac.dims = as.integer(atac.dims)),
    list(multimodal = meta_path, joint_pca = emb_path)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}

#' @keywords internal
.ic2_atac_dim_index <- function(dims, n) {
  dims <- as.integer(dims)
  dims <- dims[is.finite(dims) & dims >= 1L & dims <= n]
  if (!length(dims)) stop("No PCA dimensions in range 1..", n, ".")
  unique(dims)
}
