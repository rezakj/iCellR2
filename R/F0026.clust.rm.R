#' Remove cells in selected clusters
#'
#' Drops cells whose \code{cluster_id} is in \code{clust.to.rm} from
#' \code{x@cell_data} only. Raw / normalized assays in SODA are left
#' untouched (same idea as keeping full counts after \code{cell.filter}).
#' Plotters already intersect embeddings with \code{cell_data}, so this
#' stays memory-light even at tens of millions of cells.
#'
#' @param x An \code{iCellR2} object with \code{cluster_id}.
#' @param clust.to.rm Cluster label(s) to remove.
#' @param modality Modality for clustering sidecar. Default \code{"rna"}.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} with fewer rows in \code{cell_data}.
#'
#' @examples
#' \dontrun{
#' obj <- clust.rm(obj, clust.to.rm = c(8, 12))
#' }
#'
#' @export
clust.rm <- function(
    x,
    clust.to.rm = NULL,
    modality = "rna",
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  rep <- .ic2_new_reporter("clust.rm", verbose, soda_root = x@storage$root_dir)

  if (!"cluster_id" %in% colnames(x@cell_data)) {
    stop("No 'cluster_id' in cell_data. Run iclust() / run.clustering() first.")
  }
  if (is.null(clust.to.rm) || length(clust.to.rm) < 1L ||
      identical(clust.to.rm, "numeric")) {
    stop("Provide cluster label(s) in 'clust.to.rm'.")
  }

  n0 <- nrow(x@cell_data)
  drop_lab <- as.character(clust.to.rm)
  keep <- !as.character(x@cell_data$cluster_id) %in% drop_lab
  n_drop <- sum(!keep)
  if (n_drop < 1L) {
    warning("No cells matched clust.to.rm = ", paste(drop_lab, collapse = ", "))
    .ic2_finish(rep, extra_lines = "No cells removed")
    return(x)
  }
  if (sum(keep) < 2L) {
    stop("Removing these clusters would leave fewer than 2 cells.")
  }

  .ic2_step(rep, paste0(
    "Dropping ", .ic2_fmt_num(n_drop), " cells in cluster(s) ",
    paste(drop_lab, collapse = ", ")
  ))

  ## Snapshot original labels for the kept cells only (aligned)
  if (!is.null(x@misc$cluster_id_original)) {
    orig <- x@misc$cluster_id_original
    nm <- names(orig)
    if (!is.null(nm)) {
      x@misc$cluster_id_original <- orig[rownames(x@cell_data)[keep]]
    } else if (length(orig) == n0) {
      x@misc$cluster_id_original <- orig[keep]
      names(x@misc$cluster_id_original) <- rownames(x@cell_data)[keep]
    }
  }

  x@cell_data <- x@cell_data[keep, , drop = FALSE]
  x <- .ic2_sync_cluster_id(x, modality = modality)
  x@misc$clust.avg <- NULL
  if (!is.null(x@modalities[[modality]])) {
    x@modalities[[modality]]$clust.avg <- NULL
  }

  x <- .ic2_log_command(
    x, "clust.rm",
    arguments = list(clust.to.rm = drop_lab, modality = modality),
    outputs = list(
      n_removed = n_drop,
      n_cells = nrow(x@cell_data)
    )
  )
  .ic2_finish(rep, extra_lines = paste0(
    "Removed: ", .ic2_fmt_num(n_drop),
    "   |   remaining: ", .ic2_fmt_num(nrow(x@cell_data)),
    "   |   assays untouched (cell_data only)"
  ))
  x
}
