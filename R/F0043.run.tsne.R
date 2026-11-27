#' tSNE on an existing PCA
#'
#' Barnes-Hut tSNE (\pkg{Rtsne}) of the PCA already stored on the object
#' (\code{aligned_pca} when \code{use = "auto"} and it exists). The gene
#' matrix is not read. Exact tSNE (\code{theta = 0}) is refused above
#' 20,000 cells. A second 3D map is written only when \code{add.3d = TRUE}.
#'
#' @param x An iCellR2 object after \code{\link{run.pca}}.
#' @param use \code{"auto"}, \code{"pca"}, or \code{"aligned_pca"}.
#' @param dims PCA columns. Default \code{1:30}.
#' @param perplexity Rtsne perplexity. Default \code{30}.
#' @param theta Barnes-Hut tradeoff. Default \code{0.5}.
#' @param max_iter Iterations. Default \code{750}.
#' @param add.3d Also write a 3-column map (\code{tsne_3d}). Default \code{FALSE}.
#' @param check_duplicates Passed to Rtsne. Default \code{FALSE} (required for large data).
#' @param n.threads BLAS threads inside Rtsne when that argument exists.
#' @param my.seed Seed.
#' @param modality Modality. Default \code{"rna"}.
#' @param pca,initial_dims,verbose,clust.method,top.rank,gene.list Classic
#'   arguments. Gene ranking is ignored; tSNE always starts from PCA.
#' @return The object, with \code{reductions/tsne} or \code{aligned_tsne}.
#' @examples
#' \dontrun{
#' obj <- run.tsne(obj, perplexity = 30, theta = 0.5)
#' }
#' @export
run.tsne <- function(
    x = NULL,
    use = "auto",
    dims = 1:30,
    perplexity = 30,
    theta = 0.5,
    max_iter = 750L,
    add.3d = FALSE,
    check_duplicates = FALSE,
    n.threads = 1L,
    my.seed = 1,
    modality = "rna",
    pca = FALSE,
    initial_dims = NULL,
    verbose = TRUE,
    clust.method = NULL,
    top.rank = NULL,
    gene.list = NULL
) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package 'Rtsne' is required for run.tsne().")
  }
  if (theta <= 0 && .ic2_n_cells(x) > 20000L) {
    stop("theta = 0 is exact tSNE. Use theta = 0.5 above 20,000 cells.")
  }
  modality <- tolower(modality)
  src <- .ic2_resolve_use(x, modality, "pca", use)
  emb <- .ic2_get_plot_embedding(x, modality, src$name)
  dims <- as.integer(dims)
  dims <- dims[is.finite(dims) & dims >= 1L & dims <= ncol(emb)]
  if (length(dims) < 2L) stop("Need at least 2 PCA columns.")
  emb <- emb[, dims, drop = FALSE]
  n <- nrow(emb)
  perplexity <- min(as.numeric(perplexity), floor((n - 1) / 3))
  if (perplexity < 5) stop("Too few cells for tSNE.")
  out_name <- .ic2_aligned_out_name("tsne", src$aligned)
  set.seed(as.integer(my.seed))
  lay <- .ic2_rtsne_run(emb, 2L, perplexity, theta, as.integer(max_iter), check_duplicates, as.integer(n.threads))
  rownames(lay) <- rownames(emb)
  colnames(lay) <- c("tSNE_1", "tSNE_2")
  x <- .ic2_tsne_store(x, modality, lay, out_name, src, dims, perplexity, theta)
  if (isTRUE(add.3d)) {
    lay3 <- .ic2_rtsne_run(emb, 3L, perplexity, theta, as.integer(max_iter), check_duplicates, as.integer(n.threads))
    rownames(lay3) <- rownames(emb)
    colnames(lay3) <- c("tSNE_1", "tSNE_2", "tSNE_3")
    x <- .ic2_tsne_store(
      x, modality, lay3, .ic2_aligned_out_name("tsne_3d", src$aligned),
      src, dims, perplexity, theta
    )
  }
  x
}

#' @keywords internal
.ic2_n_cells <- function(x) {
  if (!is.null(x@cell_data)) nrow(x@cell_data) else NA_integer_
}

#' @keywords internal
.ic2_rtsne_run <- function(emb, n_comp, perplexity, theta, max_iter, check_duplicates, n.threads) {
  args <- list(
    X = emb, dims = n_comp, pca = FALSE, perplexity = perplexity,
    theta = theta, max_iter = max_iter, check_duplicates = isTRUE(check_duplicates),
    verbose = FALSE
  )
  if ("num_threads" %in% names(formals(Rtsne::Rtsne))) args$num_threads <- n.threads
  out <- do.call(Rtsne::Rtsne, args)
  as.matrix(out$Y)
}

#' @keywords internal
.ic2_tsne_store <- function(x, modality, lay, name, src, dims, perplexity, theta) {
  root <- x@storage$root_dir
  red_dir <- file.path(root, "reductions", name)
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  meta_path <- file.path(red_dir, "reduction.json")
  .ic2_h5_write_dense(lay, emb_path)
  meta <- list(
    name = name, mark = src$mark, aligned = isTRUE(src$aligned),
    source = src$name, based_on = src$name, modality = modality,
    dims = as.integer(dims),
    perplexity = perplexity, theta = theta,
    paths = list(embeddings = emb_path)
  )
  .ic2_soda_write_json(meta, meta_path)
  x@modalities[[modality]]$reductions[[name]] <- list(
    cell_embeddings = emb_path, meta_path = meta_path,
    aligned = isTRUE(src$aligned), mark = src$mark,
    source = src$name, based_on = src$name
  )
  x@modalities[[modality]]$active_reduction <- name
  x@storage$reduction_paths[[name]] <- red_dir
  x <- .ic2_log_command(
    x, "run.tsne",
    list(source = src$name, perplexity = perplexity, theta = theta),
    list(reduction = red_dir)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}
