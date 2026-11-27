#' Find marker genes for each cluster (chunked, on-disk)
#'
#' One-vs-rest markers without loading the full assay. Streams the chosen
#' assay once, accumulates per-cluster sums / sum-of-squares, then computes
#' fold-change and Welch \code{t.test} p-values from those aggregates.
#' \code{wilcox.test} is supported via a second pass on FC-passing genes
#' with optional cell downsampling.
#'
#' @param x An \code{iCellR2} object with \code{cluster_id} in
#'   \code{cell_data}.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Assay to use; default active / \code{"data"}.
#' @param transform \code{"auto"} (default): \code{log2p1} unless assay is
#'   \code{scale.data}; else \code{"none"}, \code{"log2p1"}, \code{"log1p"}.
#' @param fold.change Minimum fold-change (cluster / other). Default \code{2}.
#' @param pval.test \code{"t.test"} (default, aggregate Welch) or
#'   \code{"wilcox.test"} (FC-filtered genes + downsampled cells).
#' @param p.adjust.method Passed to \code{p.adjust}. Default \code{"hochberg"}.
#' @param padjval Maximum adjusted p-value to keep. Default \code{0.1}.
#' @param low.cell.filt Drop clusters with fewer than this many cells.
#'   Default \code{5}.
#' @param Inf.FCs Keep infinite fold-changes. Default \code{FALSE}.
#' @param uniq Keep only genes that mark a single cluster. Default \code{FALSE}.
#' @param positive Keep only up-regulated markers. Default \code{TRUE}.
#' @param max.cells Max cells for \code{wilcox.test} second pass.
#'   Default \code{5000000}. Ignored for \code{t.test}.
#' @param chunk.size Cells per HDF5 read. Default \code{200000}.
#' @param verbose Progress style.
#'
#' @return A markers \code{data.frame} (same spirit as classic iCellR).
#'
#' @examples
#' \dontrun{
#' markers <- findMarkers(obj)
#' top.markers(markers, topde = 10)
#' }
#'
#' @export
findMarkers <- function(
    x = NULL,
    modality = "rna",
    assay.name = NULL,
    transform = c("auto", "none", "log2p1", "log1p"),
    fold.change = 2,
    pval.test = c("t.test", "wilcox.test"),
    p.adjust.method = "hochberg",
    padjval = 0.1,
    low.cell.filt = 5L,
    Inf.FCs = FALSE,
    uniq = FALSE,
    positive = TRUE,
    max.cells = 5000000L,
    chunk.size = 200000L,
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  transform <- match.arg(transform)
  pval.test <- match.arg(pval.test)
  chunk.size <- as.integer(chunk.size)
  low.cell.filt <- as.integer(low.cell.filt)
  max.cells <- as.integer(max.cells)
  fold.change <- as.numeric(fold.change)[1L]
  padjval <- as.numeric(padjval)[1L]
  if (!is.finite(fold.change) || fold.change <= 0) {
    stop("'fold.change' must be > 0.")
  }

  cl_col <- "cluster_id"
  alt_col <- paste0(modality, "_cluster_id")
  if (!identical(modality, "rna") && alt_col %in% colnames(x@cell_data)) {
    cl_col <- alt_col
  }
  if (!cl_col %in% colnames(x@cell_data)) {
    stop("No '", cl_col, "' in cell_data. Run run.clustering() first.")
  }
  if (!identical(cl_col, "cluster_id")) {
    x@cell_data$cluster_id <- x@cell_data[[cl_col]]
  }

  assay.names <- names(x@modalities[[modality]]$assays %||% list())
  if (is.null(assay.name)) {
    assay.name <- if ("data" %in% assay.names) "data" else NULL
  }
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  if (identical(transform, "auto")) {
    transform <- if (assay.name %in% c("scale.data", "tfidf")) "none" else "log2p1"
  }

  rep <- .ic2_new_reporter("findMarkers", verbose, soda_root = x@storage$root_dir)

  meta <- x@cell_data
  keep <- !is.na(meta$cluster_id)
  meta <- meta[keep, , drop = FALSE]
  if (nrow(meta) < 2L) stop("Need at least 2 cells with cluster_id.")

  tab <- table(meta$cluster_id)
  drop_cl <- names(tab)[as.integer(tab) < low.cell.filt]
  if (length(drop_cl)) {
    .ic2_step(rep, paste0(
      "  Excluding ", length(drop_cl), " cluster(s) with <",
      low.cell.filt, " cells: ", paste(drop_cl, collapse = ",")
    ))
    meta <- meta[!as.character(meta$cluster_id) %in% as.character(drop_cl), , drop = FALSE]
  }

  clusters <- sort(unique(as.integer(meta$cluster_id)))
  if (length(clusters) < 2L) {
    stop("Need at least 2 clusters after low.cell.filt.")
  }
  if (anyNA(clusters)) stop("'cluster_id' must be numeric (no NA after filter).")

  cell_to_cl <- setNames(as.integer(meta$cluster_id), rownames(meta))
  n_cl <- length(clusters)

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  size_factors <- .ic2_assay_size_factors(x, modality, assay.name)
  h0 <- .ic2_enrich_handle(handles[[1]], load_dimnames = TRUE)
  genes <- as.character(h0$feature_ids %||% character(0))
  ## Lazy-norm assay may store the universe even when source handles are thin
  if (length(genes) < 1L) {
    aref <- x@modalities[[modality]]$assays[[assay.name]]
    if (is.list(aref) && length(aref$feature_universe)) {
      genes <- as.character(aref$feature_universe)
    }
  }
  if (length(genes) < 1L) {
    fd <- x@modalities[[modality]]$feature_data
    if (is.data.frame(fd) && nrow(fd) > 0L) {
      genes <- as.character(fd[[1L]])
    }
  }
  n_feat <- length(genes)
  if (n_feat < 1L) stop("Assay has no features.")

  sum_mat <- matrix(0, nrow = n_feat, ncol = n_cl)
  sumsq_mat <- matrix(0, nrow = n_feat, ncol = n_cl)
  count_vec <- integer(n_cl)
  colnames(sum_mat) <- colnames(sumsq_mat) <- paste0("cluster_", clusters)

  rep <- .ic2_banner_start(rep, c(
    paste0(
      "Assay: ", assay.name, "   |   transform: ", transform,
      if (!is.null(size_factors)) "   |   lazy_norm" else ""
    ),
    paste0(
      "Clusters: ", n_cl,
      "   |   cells: ", .ic2_fmt_num(nrow(meta)),
      "   |   test: ", pval.test
    ),
    paste0(
      "FC >= ", fold.change,
      "   |   padj < ", padjval,
      "   |   positive=", positive
    )
  ))

  .ic2_step(rep, "[1/3] Streaming assay (per-cluster sum / sumsq)")
  total_cells_plan <- nrow(meta)
  done_cells <- 0L

  for (hi in seq_along(handles)) {
    handle <- .ic2_enrich_handle(handles[[hi]], load_dimnames = TRUE)
    if (!identical(handle$feature_ids, genes)) {
      if (!all(genes %in% handle$feature_ids)) {
        ## Allow genes from feature_universe / feature_data that map into handle
        if (is.null(handle$feature_ids) || length(handle$feature_ids) < 1L) {
          stop("Feature ids missing on assay handle (tenx_meta?). Re-open SODA or re-run make.obj.")
        }
        row_map <- match(genes, handle$feature_ids)
        if (all(is.na(row_map))) {
          stop("Feature ids differ across assay handles.")
        }
      } else {
        row_map <- match(genes, handle$feature_ids)
      }
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
        if (!is.null(row_map)) {
          ## NA rows  ->  zeros for genes absent from this handle
          out <- matrix(0, nrow = length(genes), ncol = ncol(block),
                        dimnames = list(genes, colnames(block)))
          hit <- which(!is.na(row_map))
          if (length(hit)) {
            sub <- block[row_map[hit], , drop = FALSE]
            if (inherits(sub, "Matrix")) sub <- as.matrix(.ic2_to_dgCMatrix(sub))
            else sub <- as.matrix(sub)
            storage.mode(sub) <- "double"
            out[hit, ] <- sub
          }
          block <- out
        }
        if (!is.null(size_factors)) {
          block <- .ic2_scale_cols_by_sf(block, size_factors)
        }
        block <- .ic2_assay_transform_block(block, transform = transform)
        cids <- handle$cell_ids[cols]
        cl <- cell_to_cl[cids]
        for (cj in seq_along(clusters)) {
          hit <- which(cl == clusters[cj])
          if (length(hit) == 0L) next
          sub <- block[, hit, drop = FALSE]
          if (inherits(sub, "Matrix")) {
            rs <- as.numeric(Matrix::rowSums(sub))
            rssq <- as.numeric(Matrix::rowSums(sub * sub))
          } else {
            sub <- as.matrix(sub)
            storage.mode(sub) <- "double"
            rs <- rowSums(sub)
            rssq <- rowSums(sub * sub)
          }
          sum_mat[, cj] <- sum_mat[, cj] + rs
          sumsq_mat[, cj] <- sumsq_mat[, cj] + rssq
          count_vec[cj] <- count_vec[cj] + length(hit)
        }
        done_cells <- done_cells + length(cols)
        rep <- .ic2_progress_bar(rep, done_cells, total_cells_plan, prefix = "  markers")
      }
    }
  }

  if (any(count_vec < 2L)) {
    bad <- clusters[count_vec < 2L]
    stop(
      "Clusters with <2 accumulated cells (check assay cell ids): ",
      paste(bad, collapse = ",")
    )
  }

  .ic2_step(rep, "[2/3] Computing FC + p-values per cluster")
  mean_mat <- sweep(sum_mat, 2L, count_vec, `/`)
  ## unbiased sample variance per gene x cluster: (sumsq - sum^2/n)/(n-1)
  var_mat <- (sumsq_mat - sweep(sum_mat^2, 2L, count_vec, `/`)) /
    matrix(pmax(count_vec - 1L, 1L), nrow = n_feat, ncol = n_cl, byrow = TRUE)
  var_mat[!is.finite(var_mat) | var_mat < 0] <- 0

  total_sum <- rowSums(sum_mat)
  total_sumsq <- rowSums(sumsq_mat)
  total_n <- sum(count_vec)
  base_mean <- total_sum / total_n
  base_var <- (total_sumsq - (total_sum^2) / total_n) / max(total_n - 1L, 1L)
  base_var[!is.finite(base_var) | base_var < 0] <- 0
  base_sd <- sqrt(base_var)

  avg_df <- data.frame(
    gene = genes,
    mean_mat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  log2_fc_cut <- log2(fold.change)
  out_list <- vector("list", n_cl)

  ## Optional wilcox cell sample (shared across clusters)
  wilcox_cells <- NULL
  if (identical(pval.test, "wilcox.test")) {
    all_ids <- names(cell_to_cl)
    if (length(all_ids) > max.cells) {
      ## Stratified sample by cluster
      set.seed(1L)
      parts <- lapply(clusters, function(cl) {
        ids <- names(cell_to_cl)[cell_to_cl == cl]
        n_take <- max(
          2L,
          as.integer(round(
            as.numeric(max.cells) * as.numeric(length(ids)) /
              as.numeric(length(all_ids))
          ))
        )
        n_take <- min(n_take, length(ids))
        if (length(ids) <= n_take) ids else sample(ids, n_take)
      })
      wilcox_cells <- unique(unlist(parts, use.names = FALSE))
      if (length(wilcox_cells) > max.cells) {
        wilcox_cells <- sample(wilcox_cells, max.cells)
      }
      .ic2_step(rep, paste0(
        "  wilcox downsample: ", .ic2_fmt_num(length(wilcox_cells)),
        " / ", .ic2_fmt_num(length(all_ids)), " cells"
      ))
    } else {
      wilcox_cells <- all_ids
    }
  }

  for (cj in seq_along(clusters)) {
    cl <- clusters[cj]
    n1 <- as.numeric(count_vec[cj])
    n2 <- as.numeric(total_n - n1)
    m1 <- mean_mat[, cj]
    m2 <- (total_sum - sum_mat[, cj]) / n2
    v1 <- var_mat[, cj]
    ## rest variance from pooled aggregates excluding cluster cj
    sum2 <- total_sum - sum_mat[, cj]
    sumsq2 <- total_sumsq - sumsq_mat[, cj]
    v2 <- (sumsq2 - (sum2^2) / n2) / max(n2 - 1, 1)
    v2[!is.finite(v2) | v2 < 0] <- 0

    fc <- m1 / m2
    fc[!is.finite(fc) & m1 > 0 & m2 == 0] <- Inf
    fc[!is.finite(fc) & m1 == 0 & m2 == 0] <- NaN
    lfc <- log2(fc)

    if (isTRUE(positive)) {
      keep_fc <- is.finite(lfc) & lfc > log2_fc_cut
      if (isTRUE(Inf.FCs)) keep_fc <- keep_fc | (lfc == Inf)
    } else {
      keep_fc <- is.finite(lfc) & abs(lfc) > log2_fc_cut
      if (isTRUE(Inf.FCs)) keep_fc <- keep_fc | !is.finite(lfc)
    }
    ## When Inf.FCs=FALSE, drop Inf/-Inf (classic)
    if (!isTRUE(Inf.FCs)) {
      keep_fc <- keep_fc & is.finite(lfc)
    }

    idx_fc <- which(keep_fc)
    if (length(idx_fc) < 1L) {
      out_list[[cj]] <- NULL
      next
    }

    if (identical(pval.test, "t.test")) {
      se2 <- v1[idx_fc] / n1 + v2[idx_fc] / n2
      se2[se2 <= 0 | !is.finite(se2)] <- .Machine$double.eps
      tstat <- (m1[idx_fc] - m2[idx_fc]) / sqrt(se2)
      ## Welch-Satterthwaite df
      a <- v1[idx_fc] / n1
      b <- v2[idx_fc] / n2
      df <- (a + b)^2 / (a^2 / max(n1 - 1, 1) + b^2 / max(n2 - 1, 1))
      df[!is.finite(df) | df < 1] <- 1
      pval <- 2 * stats::pt(-abs(tstat), df = df)
    } else {
      ## Load FC-passing genes for downsampled cells, then wilcox
      cand <- genes[idx_fc]
      expr <- .ic2_load_features_for_cells(
        handles = handles,
        features = cand,
        cell_ids = wilcox_cells,
        chunk.size = chunk.size,
        transform = transform,
        size_factors = size_factors
      )
      cl_s <- cell_to_cl[wilcox_cells]
      in_cl <- cl_s == cl
      pval <- vapply(cand, function(g) {
        if (!g %in% colnames(expr)) return(NA_real_)
        x1 <- expr[wilcox_cells[in_cl], g]
        x2 <- expr[wilcox_cells[!in_cl], g]
        if (length(x1) < 1L || length(x2) < 1L) return(NA_real_)
        tryCatch(
          stats::wilcox.test(x = x1, y = x2)$p.value,
          error = function(e) NA_real_
        )
      }, numeric(1))
    }

    padj <- stats::p.adjust(pval, method = p.adjust.method)
    keep_p <- which(is.finite(padj) & padj < padjval)
    if (length(keep_p) < 1L) {
      out_list[[cj]] <- NULL
      next
    }
    ii <- idx_fc[keep_p]

    stats_df <- data.frame(
      gene = genes[ii],
      baseMean = base_mean[ii],
      baseSD = base_sd[ii],
      AvExpInCluster = m1[ii],
      AvExpInOtherClusters = m2[ii],
      foldChange = fc[ii],
      log2FoldChange = lfc[ii],
      pval = pval[keep_p],
      padj = padj[keep_p],
      clusters = cl,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    ## Attach full cluster-average columns (classic behavior)
    stats_df <- merge(stats_df, avg_df, by = "gene", all.x = TRUE, sort = FALSE)
    out_list[[cj]] <- stats_df
    .ic2_step(rep, paste0(
      "  cluster ", cl, ": ", .ic2_fmt_num(nrow(stats_df)), " markers"
    ))
  }

  .ic2_step(rep, "[3/3] Binding + filtering results")
  out_list <- Filter(Negate(is.null), out_list)
  if (length(out_list) < 1L) {
    .ic2_finish(rep, extra_lines = "No markers passed filters.")
    return(data.frame())
  }
  df <- do.call(rbind, out_list)
  rownames(df) <- NULL

  if (isTRUE(uniq)) {
    freq <- table(df$gene)
    keep_g <- names(freq)[freq == 1L]
    df <- df[df$gene %in% keep_g, , drop = FALSE]
  }
  if (!isTRUE(Inf.FCs)) {
    df <- df[is.finite(df$log2FoldChange), , drop = FALSE]
  }

  ## Sort: by cluster asc, then log2FC desc (classic)
  ord <- order(df$clusters, -df$log2FoldChange, na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  rownames(df) <- make.names(df$gene, unique = TRUE)

  ## Persist for reproducibility
  out_path <- NULL
  if (!is.null(x@storage$root_dir) && nzchar(x@storage$root_dir)) {
    out_dir <- file.path(x@storage$root_dir, "metadata", "markers")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_path <- .ic2_soda_write_table(df, file.path(out_dir, paste0(modality, "_findMarkers")))
    invisible(.ic2_log_command(
      x, "findMarkers",
      arguments = list(
        modality = modality, assay.name = assay.name, transform = transform,
        fold.change = fold.change, pval.test = pval.test,
        p.adjust.method = p.adjust.method, padjval = padjval,
        low.cell.filt = low.cell.filt, Inf.FCs = Inf.FCs,
        uniq = uniq, positive = positive, max.cells = max.cells
      ),
      outputs = list(n_markers = nrow(df), path = out_path)
    ))
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
    rep <- .ic2_add_output(rep, "Markers", out_path)
  }

  .ic2_finish(rep, extra_lines = paste0(
    "Markers: ", .ic2_fmt_num(nrow(df)), " rows   |   clusters: ", n_cl
  ))
  df
}
