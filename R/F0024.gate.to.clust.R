#' Assign a cluster id to gated cells
#'
#' Sets \code{cluster_id} for selected cell barcodes. Metadata-only  - 
#' safe at tens of millions of cells.
#'
#' @param x An \code{iCellR2} object with \code{cluster_id}.
#' @param my.gate Character vector of cell ids, or a path to a text file
#'   with one cell id per line.
#' @param to.clust Cluster label to assign.
#' @param modality SODA clustering sidecar modality. Default \code{"rna"}.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} object.
#'
#' @examples
#' \dontrun{
#' obj <- gate.to.clust(obj, my.gate = ids, to.clust = 5)
#' obj <- gate.to.clust(obj, my.gate = "selected_cells.txt", to.clust = 5)
#' }
#'
#' @export
gate.to.clust <- function(
    x,
    my.gate = NULL,
    to.clust = 0,
    modality = "rna",
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  rep <- .ic2_new_reporter("gate.to.clust", verbose, soda_root = x@storage$root_dir)

  if (!"cluster_id" %in% colnames(x@cell_data)) {
    stop("No 'cluster_id' in cell_data. Run iclust() / run.clustering() first.")
  }
  if (is.null(my.gate) || length(my.gate) < 1L) {
    stop("'my.gate' must be a character vector of cell ids or a file path.")
  }

  ## File path (single existing path) or character vector of ids
  if (is.character(my.gate) && length(my.gate) == 1L &&
      file.exists(my.gate) && !my.gate %in% rownames(x@cell_data)) {
    gates <- unique(trimws(readLines(my.gate, warn = FALSE)))
    gates <- gates[nzchar(gates)]
  } else {
    gates <- unique(as.character(my.gate))
  }
  if (length(gates) < 1L) stop("No cell ids found in 'my.gate'.")

  x <- .ic2_ensure_cluster_original(x)
  cd_ids <- rownames(x@cell_data)
  hit <- which(cd_ids %in% gates)
  if (length(hit) < 1L) {
    stop("None of the gated cell ids match rownames(x@cell_data).")
  }
  miss <- setdiff(gates, cd_ids)
  if (length(miss) > 0L) {
    warning(
      length(miss), " gated ids not in cell_data (first: ",
      paste(utils::head(miss, 3L), collapse = ", "), ")."
    )
  }

  to.clust <- to.clust[1L]
  cur <- as.character(x@cell_data$cluster_id)
  cur[hit] <- as.character(to.clust)
  num <- suppressWarnings(as.integer(cur))
  x@cell_data$cluster_id <- if (!anyNA(num) && all(as.character(num) == cur)) {
    num
  } else {
    cur
  }

  .ic2_step(rep, paste0(
    "Assigned cluster ", to.clust, " to ", .ic2_fmt_num(length(hit)),
    " / ", .ic2_fmt_num(length(cd_ids)), " cells"
  ))

  x <- .ic2_sync_cluster_id(x, modality = modality)
  x@misc$clust.avg <- NULL
  if (!is.null(x@modalities[[modality]])) {
    x@modalities[[modality]]$clust.avg <- NULL
  }

  x <- .ic2_log_command(
    x, "gate.to.clust",
    arguments = list(to.clust = to.clust, n_gate = length(gates), modality = modality),
    outputs = list(n_assigned = length(hit))
  )
  .ic2_finish(rep, extra_lines = paste0(
    "Assigned: ", .ic2_fmt_num(length(hit)), " cells  ->  cluster ", to.clust
  ))
  x
}
