#' Reorder / relabel clusters
#'
#' Relabels \code{cluster_id} either randomly or by hierarchical clustering
#' of per-cluster average expression (pseudotime-ish order along the
#' dendrogram). Requires \code{\link{clust.avg.exp}} for
#' \code{how.to.order = "distance"} (runs it automatically if missing).
#'
#' @param x An \code{iCellR2} object with \code{cluster_id}.
#' @param top.rank Genes (by mean across clusters) used for the dendrogram.
#'   Default \code{500}.
#' @param clust.method \code{hclust} method. Default \code{"complete"}.
#' @param dist.method Distance for cluster-average profiles.
#' @param how.to.order \code{"distance"} (default) or \code{"random"}.
#' @param modality Modality for averages / metadata. Default \code{"rna"}.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} with remapped \code{cluster_id}.
#'
#' @export
clust.ord <- function(
    x = NULL,
    top.rank = 500L,
    dist.method = "euclidean",
    clust.method = "complete",
    how.to.order = c("distance", "random"),
    modality = "rna",
    verbose = TRUE
) {
  .ic2_require_object(x)
  how.to.order <- match.arg(how.to.order)
  modality <- tolower(modality)
  top.rank <- as.integer(top.rank)

  if (!"cluster_id" %in% colnames(x@cell_data)) {
    stop("No 'cluster_id' in cell_data. Run run.clustering() first.")
  }

  rep <- .ic2_new_reporter("clust.ord", verbose, soda_root = x@storage$root_dir)
  old <- as.integer(x@cell_data$cluster_id)
  if (anyNA(old)) stop("'cluster_id' contains NA.")
  if (!is.numeric(old)) stop("Clusters must be numeric.")

  uniq <- sort(unique(old))
  rep <- .ic2_banner_start(rep, c(
    paste0("Clusters: ", length(uniq), "   |   order: ", how.to.order),
    paste0("Cells: ", .ic2_fmt_num(length(old)))
  ))

  .ic2_step(rep, "[1/2] Computing new label order")
  if (identical(how.to.order, "random")) {
    set.seed(1)
    BestOrd <- sample(uniq)
  } else {
    avg <- x@misc$clust.avg
    if (is.null(avg) && !is.null(x@storage$metadata_paths$clust_avg) &&
        file.exists(x@storage$metadata_paths$clust_avg)) {
      avg <- .ic2_soda_read_table(x@storage$metadata_paths$clust_avg)
      x@misc$clust.avg <- avg
    }
    if (is.null(avg)) {
      .ic2_step(rep, "  No clust.avg yet  -  running clust.avg.exp()")
      x <- clust.avg.exp(x, modality = modality, verbose = verbose)
      avg <- x@misc$clust.avg
    }
    mat <- as.matrix(avg[, setdiff(colnames(avg), "gene"), drop = FALSE])
    rownames(mat) <- avg$gene
    ## column names cluster_N
    cn <- colnames(mat)
    cl_nums <- as.integer(sub("^cluster_", "", cn))
    mat <- mat[, order(cl_nums), drop = FALSE]
    cl_nums <- sort(cl_nums)

    mat <- mat[order(rowMeans(mat), decreasing = TRUE), , drop = FALSE]
    mat <- utils::head(mat, min(top.rank, nrow(mat)))
    if (ncol(mat) < 2L) {
      BestOrd <- cl_nums
    } else {
      d <- stats::dist(scale(t(mat)), method = dist.method)
      hc <- stats::hclust(d, method = clust.method)
      BestOrd <- cl_nums[hc$order]
    }
  }
  rep <- .ic2_progress_bar(rep, 1L, 2L, prefix = "  order")

  .ic2_step(rep, "[2/2] Remapping cluster_id")
  ## map old label -> new 1..K in BestOrd sequence
  new_map <- setNames(seq_along(BestOrd), as.character(BestOrd))
  ## also map any leftover labels not in BestOrd
  leftover <- setdiff(uniq, BestOrd)
  if (length(leftover)) {
    new_map <- c(new_map, setNames(
      length(BestOrd) + seq_along(leftover),
      as.character(leftover)
    ))
  }
  new_id <- as.integer(new_map[as.character(old)])
  x@cell_data$cluster_id <- new_id
  .ic2_persist_cell_data(x)

  best <- data.frame(
    clusters = new_id,
    row.names = rownames(x@cell_data),
    stringsAsFactors = FALSE
  )
  x@misc$best.clust <- best

  cl_dir <- file.path(x@storage$root_dir, "metadata", "clustering")
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
  cl_stem <- file.path(cl_dir, paste0(modality, "_cluster_id"))
  cl_df <- data.frame(
    cell_id = rownames(x@cell_data),
    cluster_id = as.character(new_id),
    stringsAsFactors = FALSE
  )
  cl_path <- .ic2_soda_write_table(cl_df, cl_stem)
  if (is.null(x@storage$metadata_paths)) x@storage$metadata_paths <- list()
  x@storage$metadata_paths$clustering <- cl_path

  ## Invalidate stale averages (labels changed)
  x@misc$clust.avg <- NULL
  x@modalities[[modality]]$clust.avg <- NULL

  x <- .ic2_log_command(
    x, "clust.ord",
    arguments = list(
      how.to.order = how.to.order, top.rank = top.rank,
      dist.method = dist.method, clust.method = clust.method
    ),
    outputs = list(
      order = BestOrd,
      n_clusters = length(unique(new_id)),
      path = cl_path
    )
  )

  rep <- .ic2_progress_bar(rep, 2L, 2L, prefix = "  remap")
  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  .ic2_finish(rep, extra_lines = paste0(
    "New order: ", paste(BestOrd, collapse = "  ->  ")
  ))
  x
}
