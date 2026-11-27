#' Change / rename cluster labels
#'
#' Edits \code{x@cell_data$cluster_id} only (no assay I/O)  -  safe at
#' tens of millions of cells. Optionally reset to the snapshot taken
#' before the first rename in this session (\code{misc$cluster_id_original}).
#'
#' @param x An \code{iCellR2} object with \code{cluster_id}.
#' @param change.clust Cluster label(s) to replace.
#' @param to.clust New label (length 1, or same length as
#'   \code{change.clust} for parallel remap).
#' @param clust.reset If \code{TRUE}, restore \code{cluster_id_original}
#'   when available.
#' @param modality Used only for SODA clustering sidecar path.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} object.
#'
#' @examples
#' \dontrun{
#' obj <- change.clust(obj, change.clust = 3, to.clust = 1)  # merge 3 -> 1
#' obj <- change.clust(obj, clust.reset = TRUE)
#' }
#'
#' @export
change.clust <- function(
    x,
    change.clust = 0,
    to.clust = 0,
    clust.reset = FALSE,
    modality = "rna",
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  rep <- .ic2_new_reporter("change.clust", verbose, soda_root = x@storage$root_dir)

  if (!"cluster_id" %in% colnames(x@cell_data)) {
    stop("No 'cluster_id' in cell_data. Run iclust() / run.clustering() first.")
  }

  x <- .ic2_ensure_cluster_original(x)
  n <- nrow(x@cell_data)

  if (isTRUE(clust.reset)) {
    orig <- x@misc$cluster_id_original
    if (is.null(orig) || length(orig) != n) {
      stop(
        "No cluster_id_original snapshot found. ",
        "Reset only works after at least one prior change.clust / gate.to.clust."
      )
    }
    .ic2_step(rep, paste0("Resetting cluster_id for ", .ic2_fmt_num(n), " cells"))
    x@cell_data$cluster_id <- as.integer(orig[rownames(x@cell_data)])
  } else {
    from <- as.character(change.clust)
    to <- as.character(to.clust)
    if (length(to) == 1L && length(from) > 1L) to <- rep(to, length(from))
    if (length(from) != length(to)) {
      stop("'change.clust' and 'to.clust' must have the same length (or to.clust length 1).")
    }
    if (any(!nzchar(from)) || any(!nzchar(to))) {
      stop("Cluster labels must be non-empty.")
    }
    cur <- as.character(x@cell_data$cluster_id)
    n_hit <- 0L
    for (i in seq_along(from)) {
      hit <- which(cur == from[i])
      if (length(hit)) {
        cur[hit] <- to[i]
        n_hit <- n_hit + length(hit)
      }
    }
    ## Prefer integer labels when all look numeric
    num <- suppressWarnings(as.integer(cur))
    x@cell_data$cluster_id <- if (!anyNA(num) && all(as.character(num) == cur)) {
      num
    } else {
      cur
    }
    .ic2_step(rep, paste0(
      "Remapped ", .ic2_fmt_num(n_hit), " / ", .ic2_fmt_num(n), " cells   |   ",
      paste(paste0(from, " -> ", to), collapse = ", ")
    ))
  }

  x <- .ic2_sync_cluster_id(x, modality = modality)
  x@misc$clust.avg <- NULL
  if (!is.null(x@modalities[[modality]])) {
    x@modalities[[modality]]$clust.avg <- NULL
  }

  x <- .ic2_log_command(
    x, "change.clust",
    arguments = list(
      change.clust = change.clust, to.clust = to.clust,
      clust.reset = isTRUE(clust.reset), modality = modality
    ),
    outputs = list(n_cells = n, n_clusters = length(unique(x@cell_data$cluster_id)))
  )
  .ic2_finish(rep, extra_lines = paste0(
    "Clusters: ", length(unique(x@cell_data$cluster_id)),
    "   |   cells: ", .ic2_fmt_num(n)
  ))
  x
}

#' Snapshot original cluster_id once (for change.clust reset)
#' @keywords internal
.ic2_ensure_cluster_original <- function(x) {
  if (is.null(x@misc$cluster_id_original) &&
      "cluster_id" %in% colnames(x@cell_data)) {
    v <- x@cell_data$cluster_id
    names(v) <- rownames(x@cell_data)
    x@misc$cluster_id_original <- v
  }
  x
}

#' Persist cluster_id to cell_data + best.clust + SODA sidecar
#' @keywords internal
.ic2_sync_cluster_id <- function(x, modality = "rna") {
  modality <- tolower(modality)
  .ic2_persist_cell_data(x)
  cd_ids <- rownames(x@cell_data)
  best <- data.frame(
    clusters = x@cell_data$cluster_id,
    row.names = cd_ids,
    stringsAsFactors = FALSE
  )
  x@misc$best.clust <- best
  if (!is.null(x@storage$root_dir) && nzchar(x@storage$root_dir)) {
    cl_dir <- file.path(x@storage$root_dir, "metadata", "clustering")
    dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
    cl_stem <- file.path(cl_dir, paste0(modality, "_cluster_id"))
    cl_df <- data.frame(
      cell_id = cd_ids,
      cluster_id = as.character(x@cell_data$cluster_id),
      stringsAsFactors = FALSE
    )
    cl_path <- .ic2_soda_write_table(cl_df, cl_stem)
    if (is.null(x@storage$metadata_paths)) x@storage$metadata_paths <- list()
    x@storage$metadata_paths$clustering <- cl_path
  }
  x
}
