#' Average expression per cluster (chunked, on-disk)
#'
#' Streams the chosen assay from SODA and computes per-gene means for each
#' \code{cluster_id} without loading the full matrix into memory.
#'
#' @param x An \code{iCellR2} object with \code{cluster_id} in
#'   \code{cell_data} (from \code{\link{run.clustering}}).
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Assay to average; default active assay.
#' @param conds.to.avg Optional character vector of \code{condition_id}
#'   values to include. Default \code{NULL} (all cells).
#' @param low.cell.filt Drop clusters with fewer than this many cells.
#'   Default \code{5}.
#' @param rounding.digits Digits when \code{round.num = TRUE}.
#' @param round.num Round means. Default \code{FALSE}.
#' @param chunk.size Cells per HDF5 read. Default \code{200000}.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} with cluster averages in
#'   \code{x@misc$clust.avg} and SODA metadata.
#'
#' @examples
#' \dontrun{
#' obj <- clust.avg.exp(obj)
#' head(obj@misc$clust.avg)
#' }
#'
#' @export
clust.avg.exp <- function(
    x = NULL,
    modality = "rna",
    assay.name = NULL,
    conds.to.avg = NULL,
    low.cell.filt = 5L,
    rounding.digits = 4L,
    round.num = FALSE,
    chunk.size = 200000L,
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  chunk.size <- as.integer(chunk.size)
  low.cell.filt <- as.integer(low.cell.filt)

  if (!"cluster_id" %in% colnames(x@cell_data)) {
    stop("No 'cluster_id' in cell_data. Run run.clustering() first.")
  }

  rep <- .ic2_new_reporter("clust.avg.exp", verbose, soda_root = x@storage$root_dir)

  meta <- x@cell_data
  keep <- !is.na(meta$cluster_id)
  if (!is.null(conds.to.avg)) {
    if (!"condition_id" %in% colnames(meta)) {
      stop("'conds.to.avg' set but no 'condition_id' in cell_data.")
    }
    keep <- keep & meta$condition_id %in% conds.to.avg
  }
  meta <- meta[keep, , drop = FALSE]
  if (nrow(meta) < 1L) stop("No cells left after filtering.")

  ## Drop tiny clusters
  tab <- table(meta$cluster_id)
  drop_cl <- as.integer(names(tab)[as.integer(tab) < low.cell.filt])
  if (length(drop_cl)) {
    .ic2_step(rep, paste0(
      "  Excluding ", length(drop_cl), " cluster(s) with <",
      low.cell.filt, " cells: ", paste(drop_cl, collapse = ",")
    ))
    meta <- meta[!meta$cluster_id %in% drop_cl, , drop = FALSE]
  }
  clusters <- sort(unique(as.integer(meta$cluster_id)))
  if (length(clusters) < 1L) stop("No clusters remaining after low.cell.filt.")

  cell_to_cl <- setNames(as.integer(meta$cluster_id), rownames(meta))
  n_cl <- length(clusters)

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  ## Feature ids from first handle
  h0 <- .ic2_enrich_handle(handles[[1]], load_dimnames = TRUE)
  genes <- h0$feature_ids
  n_feat <- length(genes)
  if (n_feat < 1L) stop("Assay has no features.")

  sum_mat <- matrix(0, nrow = n_feat, ncol = n_cl)
  count_vec <- integer(n_cl)
  colnames(sum_mat) <- paste0("cluster_", clusters)

  rep <- .ic2_banner_start(rep, c(
    paste0("Assay: ", assay.name, "   |   modality: ", modality),
    paste0(
      "Clusters: ", n_cl,
      "   |   cells: ", .ic2_fmt_num(nrow(meta)),
      if (!is.null(conds.to.avg)) {
        paste0("   |   conditions: ", paste(conds.to.avg, collapse = ","))
      } else ""
    )
  ))

  .ic2_step(rep, "[1/2] Streaming assay and accumulating cluster sums")
  total_cells_plan <- nrow(meta)
  done_cells <- 0L

  for (hi in seq_along(handles)) {
    handle <- .ic2_enrich_handle(handles[[hi]], load_dimnames = TRUE)
    if (!identical(handle$feature_ids, genes)) {
      ## align rows if needed
      if (!all(genes %in% handle$feature_ids)) {
        stop("Feature ids differ across assay handles.")
      }
      row_map <- match(genes, handle$feature_ids)
    } else {
      row_map <- NULL
    }

    idx <- which(handle$cell_ids %in% names(cell_to_cl))
    if (length(idx) < 1L) next
    runs <- .ic2_contig_runs(sort(idx))

    for (run in runs) {
      starts <- seq.int(1L, length(run), by = chunk.size)
      for (s in starts) {
        e <- min(s + chunk.size - 1L, length(run))
        cols <- run[s:e]
        block <- .ic2_read_cols(handle, cols)
        if (!is.null(row_map)) block <- block[row_map, , drop = FALSE]
        cids <- handle$cell_ids[cols]
        cl <- cell_to_cl[cids]
        ## accumulate per cluster in this chunk
        for (cj in seq_along(clusters)) {
          hit <- which(cl == clusters[cj])
          if (length(hit) == 0L) next
          sub <- block[, hit, drop = FALSE]
          sum_mat[, cj] <- sum_mat[, cj] + as.numeric(Matrix::rowSums(sub))
          count_vec[cj] <- count_vec[cj] + length(hit)
        }
        done_cells <- done_cells + length(cols)
        rep <- .ic2_progress_bar(rep, done_cells, total_cells_plan, prefix = "  avg")
      }
    }
  }

  .ic2_step(rep, "[2/2] Writing cluster averages")
  if (any(count_vec == 0L)) {
    warning("Some clusters had zero accumulated cells; check assay cell ids.")
  }
  mean_mat <- sweep(sum_mat, 2L, pmax(count_vec, 1L), `/`)
  if (isTRUE(round.num)) {
    mean_mat <- round(mean_mat, digits = as.integer(rounding.digits))
  }
  out <- data.frame(
    gene = genes,
    mean_mat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  avg_dir <- file.path(x@storage$root_dir, "metadata", "clustering")
  dir.create(avg_dir, recursive = TRUE, showWarnings = FALSE)
  avg_stem <- file.path(avg_dir, paste0(modality, "_clust_avg"))
  avg_path <- .ic2_soda_write_table(out, avg_stem)

  x@misc$clust.avg <- out
  x@modalities[[modality]]$clust.avg <- list(
    path = avg_path,
    assay.name = assay.name,
    n_clusters = n_cl,
    n_genes = n_feat,
    cluster_sizes = setNames(as.integer(count_vec), paste0("cluster_", clusters))
  )
  if (is.null(x@storage$metadata_paths)) x@storage$metadata_paths <- list()
  x@storage$metadata_paths$clust_avg <- avg_path

  x <- .ic2_log_command(
    x, "clust.avg.exp",
    arguments = list(
      modality = modality, assay.name = assay.name,
      low.cell.filt = low.cell.filt,
      conds.to.avg = conds.to.avg
    ),
    outputs = list(
      n_clusters = n_cl, n_genes = n_feat, path = avg_path
    )
  )

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  rep <- .ic2_add_output(rep, "Cluster avg", avg_path)
  .ic2_finish(rep, extra_lines = paste0(
    "Means: ", .ic2_fmt_num(n_feat), " genes x ", n_cl, " clusters"
  ))
  x
}
