#' Compute QC statistics for an iCellR2 object
#'
#' Computes per-cell QC from on-disk assays by streaming column chunks.
#' Works with iCellR2 CSC HDF5, RDS handles, and native \code{tenx_h5}
#' (from \code{load10x(..., output = "virtual")}). Virtual multi-sample
#' aggregations are streamed sample-by-sample  -  no merged matrix.
#'
#' @param x An \code{iCellR2} object.
#' @param modality Modality name. Default \code{"rna"}.
#' @param assay.name Assay name; default active assay.
#' @param mito.genes Optional mitochondrial feature names.
#' @param s.phase.genes Optional S-phase genes.
#' @param g2m.phase.genes Optional G2/M genes.
#' @param chunk.size Cells per chunk. Default \code{200000} (fast default for
#'   workstations with tens of GB RAM).
#' @param use.rcpp Use Rcpp sparse kernel when available. Default \code{TRUE}.
#' @param update.object Update \code{x@cell_data}. Default \code{TRUE}.
#' @param verbose Progress style: \code{TRUE}/\code{"full"}, \code{"compact"},
#'   or \code{FALSE}/\code{"silent"}.
#'
#' @return Updated \code{iCellR2} or a QC data frame.
#'
#' @importFrom Matrix colSums
#' @export
qc.stats <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    mito.genes = NULL,
    s.phase.genes = NULL,
    g2m.phase.genes = NULL,
    chunk.size = 200000L,
    use.rcpp = TRUE,
    update.object = TRUE,
    verbose = TRUE
) {
  .ic2_require_object(x)
  if (!is.character(modality) || length(modality) != 1L || !nzchar(modality)) {
    stop("'modality' must be a non-empty character string.")
  }
  modality <- tolower(modality)
  if (is.null(x@modalities[[modality]])) {
    stop("Modality '", modality, "' was not found in 'x@modalities'.")
  }
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  if (!is.numeric(chunk.size) || length(chunk.size) != 1L || is.na(chunk.size) || chunk.size < 1) {
    stop("'chunk.size' must be a positive number.")
  }
  chunk.size <- as.integer(chunk.size)
  rep <- .ic2_new_reporter("qc.stats", verbose, soda_root = x@storage$root_dir)

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  if (!length(handles)) {
    stop("Could not resolve on-disk assay handles for QC.")
  }
  rep <- .ic2_register_assay_inputs(rep, x, modality, assay.name, label = "Assay")

  n_cells_guess <- sum(vapply(handles, function(h) {
    h2 <- .ic2_enrich_handle(h, load_dimnames = FALSE)
    as.integer(h2$n_cells %||% 0L)
  }, integer(1)))
  n_feat_guess <- {
    h0 <- .ic2_enrich_handle(handles[[1L]], load_dimnames = FALSE)
    as.integer(h0$n_features %||% NA_integer_)
  }
  backends_used <- unique(vapply(handles, function(h) {
    as.character(h$backend %||% "unknown")
  }, character(1)))

  rep <- .ic2_banner_start(rep, c(
    paste0("Modality: ", modality, " | Assay: ", assay.name),
    paste0(
      "Cells: ", .ic2_fmt_num(n_cells_guess),
      if (is.finite(n_feat_guess)) paste0(" | Features: ", .ic2_fmt_num(n_feat_guess)) else "",
      " | Samples: ", length(handles),
      if (.ic2_assay_is_virtual(x, modality, assay.name)) " (virtual)" else ""
    ),
    paste0(
      "chunk.size: ", .ic2_fmt_num(chunk.size),
      " | Rcpp: ", if (isTRUE(use.rcpp)) "on" else "off",
      " | backend: ", paste(backends_used, collapse = ",")
    )
  ))

  if (is.null(s.phase.genes)) {
    s.phase.file <- system.file("extdata", "s.phase.txt", package = "iCellR2")
    s.phase.genes <- if (nzchar(s.phase.file) && file.exists(s.phase.file)) {
      scan(s.phase.file, what = character(), quiet = TRUE)
    } else character(0)
  }
  if (is.null(g2m.phase.genes)) {
    g2m.phase.file <- system.file("extdata", "g2m.phase.txt", package = "iCellR2")
    g2m.phase.genes <- if (nzchar(g2m.phase.file) && file.exists(g2m.phase.file)) {
      scan(g2m.phase.file, what = character(), quiet = TRUE)
    } else character(0)
  }

  match_idx <- function(gene.set, feature_ids) {
    idx <- match(unique(toupper(as.character(gene.set))), toupper(feature_ids))
    idx[!is.na(idx)]
  }

  ## Preallocate (avoids O(n^2) growth on multi-sample atlas runs)
  n_alloc <- max(0L, as.integer(n_cells_guess))
  all_cell_ids <- character(n_alloc)
  UMIs_all <- numeric(n_alloc)
  nGenes_all <- numeric(n_alloc)
  mito_all <- numeric(n_alloc)
  s_all <- numeric(n_alloc)
  g2m_all <- numeric(n_alloc)
  fill_pos <- 0L
  n_features_total <- NA_integer_
  mito.idx.n <- 0L
  s.idx.n <- 0L
  g2m.idx.n <- 0L

  for (hi in seq_along(handles)) {
    handle <- .ic2_enrich_handle(handles[[hi]], load_dimnames = TRUE)
    feature_ids <- handle$feature_ids
    cell_ids <- handle$cell_ids
    n_features_total <- length(feature_ids)

    if (is.null(mito.genes)) {
      mito.idx <- grep("^mt[-.]", feature_ids, ignore.case = TRUE)
    } else {
      mito.idx <- match(unique(as.character(mito.genes)), feature_ids)
      mito.idx <- mito.idx[!is.na(mito.idx)]
    }
    s.idx <- match_idx(s.phase.genes, feature_ids)
    g2m.idx <- match_idx(g2m.phase.genes, feature_ids)
    mito.idx.n <- length(mito.idx)
    s.idx.n <- length(s.idx)
    g2m.idx.n <- length(g2m.idx)

    n_cells <- length(cell_ids)
    UMIs <- numeric(n_cells)
    nGenes <- numeric(n_cells)
    mito_pct <- numeric(n_cells)
    s_score <- numeric(n_cells)
    g2m_score <- numeric(n_cells)

    starts <- seq.int(1L, n_cells, by = chunk.size)
    .ic2_step(rep, paste0(
      "[", hi, "/", length(handles), "] Streaming",
      if (!is.null(handle$sample_id)) paste0(" sample ", handle$sample_id) else "",
      " (", .ic2_fmt_num(n_cells), " cells, ", handle$backend %||% "?", ")"
    ))

    for (ci in seq_along(starts)) {
      start <- starts[ci]
      end <- min(start + chunk.size - 1L, n_cells)
      cols <- start:end
      block <- .ic2_read_cols(handle, cols)

      if (inherits(block, "dgCMatrix") && isTRUE(use.rcpp) &&
          exists("cpp_sparse_col_stats", mode = "function")) {
        base_stats <- cpp_sparse_col_stats(block)
        u <- base_stats$col_sums
        ng <- base_stats$nnz
        mito_counts <- if (length(mito.idx) > 0L) {
          cpp_sparse_col_stats(block[mito.idx, , drop = FALSE])$col_sums
        } else rep(0, length(cols))
        s_counts <- if (length(s.idx) > 0L) {
          cpp_sparse_col_stats(block[s.idx, , drop = FALSE])$col_sums
        } else rep(0, length(cols))
        g2m_counts <- if (length(g2m.idx) > 0L) {
          cpp_sparse_col_stats(block[g2m.idx, , drop = FALSE])$col_sums
        } else rep(0, length(cols))
      } else {
        u <- Matrix::colSums(block)
        ng <- if (inherits(block, "dgCMatrix")) diff(block@p) else Matrix::colSums(block != 0)
        mito_counts <- if (length(mito.idx) > 0L) Matrix::colSums(block[mito.idx, , drop = FALSE]) else rep(0, length(cols))
        s_counts <- if (length(s.idx) > 0L) Matrix::colSums(block[s.idx, , drop = FALSE]) else rep(0, length(cols))
        g2m_counts <- if (length(g2m.idx) > 0L) Matrix::colSums(block[g2m.idx, , drop = FALSE]) else rep(0, length(cols))
      }

      UMIs[cols] <- as.numeric(u)
      nGenes[cols] <- as.numeric(ng)
      mito_pct[cols] <- ifelse(u > 0, mito_counts / u, 0)
      s_score[cols] <- ifelse(u > 0, s_counts / u, 0)
      g2m_score[cols] <- ifelse(u > 0, g2m_counts / u, 0)

      rm(block)
      rep <- .ic2_progress_bar(
        rep, end, n_cells,
        prefix = paste0("  sample ", hi, "/", length(handles))
      )
    }

    need <- fill_pos + n_cells
    if (need > length(all_cell_ids)) {
      extra <- max(need - length(all_cell_ids), n_cells)
      all_cell_ids <- c(all_cell_ids, character(extra))
      UMIs_all <- c(UMIs_all, numeric(extra))
      nGenes_all <- c(nGenes_all, numeric(extra))
      mito_all <- c(mito_all, numeric(extra))
      s_all <- c(s_all, numeric(extra))
      g2m_all <- c(g2m_all, numeric(extra))
    }
    idx <- (fill_pos + 1L):(fill_pos + n_cells)
    all_cell_ids[idx] <- cell_ids
    UMIs_all[idx] <- UMIs
    nGenes_all[idx] <- nGenes
    mito_all[idx] <- mito_pct
    s_all[idx] <- s_score
    g2m_all[idx] <- g2m_score
    fill_pos <- fill_pos + n_cells

    rm(handle, UMIs, nGenes, mito_pct, s_score, g2m_score, cell_ids)
  }

  if (fill_pos < length(all_cell_ids)) {
    keep <- seq_len(fill_pos)
    all_cell_ids <- all_cell_ids[keep]
    UMIs_all <- UMIs_all[keep]
    nGenes_all <- nGenes_all[keep]
    mito_all <- mito_all[keep]
    s_all <- s_all[keep]
    g2m_all <- g2m_all[keep]
  }

  if (!all(all_cell_ids %in% rownames(x@cell_data))) {
    stop("Assay cell ids must be present in 'x@cell_data'.")
  }

  stats.df <- data.frame(
    cell_id = all_cell_ids,
    nGenes = nGenes_all,
    UMIs = UMIs_all,
    mito.percent = mito_all,
    S.phase.probability = s_all,
    g2m.phase.probability = g2m_all,
    stringsAsFactors = FALSE,
    row.names = all_cell_ids
  )

  if (!update.object) {
    .ic2_finish(rep)
    return(stats.df)
  }

  x@cell_data[all_cell_ids, "nGenes"] <- stats.df[all_cell_ids, "nGenes"]
  x@cell_data[all_cell_ids, "UMIs"] <- stats.df[all_cell_ids, "UMIs"]
  x@cell_data[all_cell_ids, "mito.percent"] <- stats.df[all_cell_ids, "mito.percent"]
  x@cell_data[all_cell_ids, "S.phase.probability"] <- stats.df[all_cell_ids, "S.phase.probability"]
  x@cell_data[all_cell_ids, "g2m.phase.probability"] <- stats.df[all_cell_ids, "g2m.phase.probability"]

  ## persist updated cell_data if SODA path known
  x <- .ic2_persist_cell_data(x)

  x@params$qc[[modality]] <- list(
    assay.name = assay.name,
    chunk.size = chunk.size,
    use.rcpp = use.rcpp,
    mito.genes.provided = !is.null(mito.genes),
    n_mito_features = mito.idx.n,
    n_s_phase_genes = s.idx.n,
    n_g2m_phase_genes = g2m.idx.n,
    backend = "on_disk_chunked"
  )

  qc_command_args <- list(
    modality = modality, assay.name = assay.name,
    chunk.size = chunk.size, use.rcpp = use.rcpp,
    update.object = update.object
  )
  qc_command_outs <- list(
    n_cells = length(all_cell_ids),
    n_features = n_features_total
  )
  x <- .ic2_log_command(x, "qc.stats", qc_command_args, qc_command_outs)

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
    rep <- .ic2_add_output(rep, "Log", file.path(x@storage$root_dir, "iCellR2.log"))
  }
  rep <- .ic2_add_output(rep, "QC cols", "nGenes, UMIs, mito.percent, S/G2M scores  ->  cell_data")
  if (!is.null(x@storage$metadata_paths$cell_data)) {
    rep <- .ic2_add_output(rep, "Saved", x@storage$metadata_paths$cell_data)
  }
  rep <- .ic2_add_scale(rep, c(
    "Mode:     chunked stream over on-disk assays (no full matrix in RAM)",
    paste0(
      "Backend:  ", paste(backends_used, collapse = "+"),
      if (isTRUE(use.rcpp)) " + Rcpp column kernels" else ""
    ),
    "vs RAM:   QC does not require a resident multi-GB sparse object",
    paste0("Chunks:   ", .ic2_fmt_num(chunk.size), " cells/chunk across ",
           length(handles), " handle(s)")
  ))
  wall <- as.numeric(difftime(Sys.time(), rep$wall0, units = "secs"))
  if (is.finite(wall) && wall > 0) {
    rep$stats$cells_per_sec <- length(all_cell_ids) / wall
  }
  .ic2_finish(rep)
  x
}
