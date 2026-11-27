#' Cluster cells on a dimensionality reduction
#'
#' Graph clustering (Louvain / PhenoGraph / k-means) on PCA/UMAP embeddings
#' stored in the SODA archive. \code{iclust} is the same function (classic
#' iCellR name)  -  both names share one implementation.
#'
#' By default clusters all cells with weighted kNN Louvain when
#' \code{n <= max.cells.direct} (5M). Set \code{use.sketch = TRUE} for
#' sketch -> assign on larger objects.
#'
#' Labels are written to \code{x@cell_data$cluster_id} (and mirrored under
#' SODA). Use \code{\link{cluster.plot}} to visualize.
#'
#' @param x An \code{iCellR2} object (after \code{run.pca} / \code{run.umap}).
#' @param modality Modality. Default \code{"rna"}.
#' @param reduction \code{"pca"} (default) or \code{"umap"}.
#'   Alias: \code{data.type} (classic \code{iclust} name).
#' @param use Which stored reduction to read when several exist.
#'   \code{"auto"} uses an aligned embedding when present, otherwise
#'   \code{reduction}.
#' @param data.type Alias for \code{reduction}.
#' @param dims Dimensions to use. Default: \code{opt.pcs} for PCA, all
#'   UMAP components otherwise.
#' @param method \code{"louvain"} (default), \code{"phenograph"}, or
#'   \code{"kmeans"}.
#' @param k Neighbors for Louvain / PhenoGraph. Default \code{30}.
#'   Alias: \code{sensitivity} (classic \code{iclust} name).
#' @param sensitivity Alias for \code{k}.
#' @param n.clust Number of clusters when \code{method = "kmeans"}.
#'   Default \code{NULL} (uses \code{max.clust} search on a sketch sample).
#' @param max.clust Upper bound for optional k-means k search.
#' @param use.sketch If \code{TRUE}, prefer sketch -> assign when n is large.
#'   Default \code{FALSE} (cluster all cells directly when
#'   \code{n <= max.cells.direct}).
#' @param sketch.name Sketch name from \code{\link{run.sketch}}.
#' @param max.cells.direct Cluster all cells directly below this size;
#'   above it, sketch -> assign when \code{use.sketch = TRUE}.
#'   Default \code{5000000}.
#' @param assign.k Neighbors used when assigning non-sketch cells.
#' @param return.graph If \code{TRUE}, return the igraph object instead of
#'   updating \code{x} (direct Louvain path only  -  not sketch -> assign).
#' @param dist.method Distance for the Louvain kNN graph when
#'   \code{return.graph = TRUE} or \code{method = "louvain"}.
#' @param seed RNG seed.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} object with \code{cluster_id} in
#'   \code{cell_data}, or an igraph object if \code{return.graph = TRUE}.
#'
#' @examples
#' \dontrun{
#' obj <- iclust(obj)                              # Louvain on PCA
#' obj <- run.clustering(obj, method = "louvain", k = 40)
#' obj <- iclust(obj, data.type = "umap", dims = 1:2, sensitivity = 30)
#' cluster.plot(obj)
#' }
#'
#' @export
run.clustering <- function(
    x,
    modality = "rna",
    reduction = c("pca", "umap"),
    use = "auto",
    data.type = NULL,
    dims = NULL,
    method = c("louvain", "phenograph", "kmeans"),
    k = 30L,
    sensitivity = NULL,
    n.clust = NULL,
    max.clust = 25L,
    use.sketch = FALSE,
    sketch.name = "sketch1",
    max.cells.direct = 5000000L,
    assign.k = 5L,
    return.graph = FALSE,
    dist.method = "euclidean",
    seed = 1,
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  if (!is.null(data.type)) {
    dt <- tolower(sub("^reductions/", "", as.character(data.type)[1L]))
    if (dt %in% c("knetl", "aligned_knetl")) {
      reduction <- dt
    } else if (dt %in% c("aligned_pca", "aligned_umap")) {
      use <- dt
      reduction <- if (grepl("umap", dt)) "umap" else "pca"
    } else {
      reduction <- match.arg(dt, c("pca", "umap"))
    }
  } else {
    reduction <- match.arg(reduction)
  }
  if (reduction %in% c("knetl", "aligned_knetl")) {
    if (is.null(x@modalities[[modality]]$reductions[[reduction]])) {
      stop("No '", reduction, "' reduction. Run run.knetl() first.")
    }
    src <- list(
      name = reduction,
      aligned = startsWith(reduction, "aligned_"),
      mark = if (startsWith(reduction, "aligned_")) "Aligned" else NULL,
      family = "knetl"
    )
  } else {
    fam <- if (identical(reduction, "umap")) "umap" else "pca"
    src <- .ic2_resolve_use(x, modality, fam, use)
    reduction <- src$name
  }
  method <- match.arg(method)
  if (!is.null(sensitivity)) {
    k <- as.integer(sensitivity)
  } else {
    k <- as.integer(k)
  }
  max.cells.direct <- as.integer(max.cells.direct)
  assign.k <- as.integer(assign.k)

  if (isTRUE(return.graph)) {
    emb <- .ic2_get_plot_embedding(x, modality, reduction)
    if (is.null(dims)) {
      dims <- seq_len(min(10L, ncol(emb)))
    }
    dims <- as.integer(dims)
    dims <- dims[dims <= ncol(emb)]
    emb <- emb[, dims, drop = FALSE]
    lv <- .ic2_louvain_knn(emb, k = k, dist.method = dist.method)
    return(lv$graph)
  }

  if (!requireNamespace("igraph", quietly = TRUE) &&
      method %in% c("phenograph", "louvain")) {
    stop("Package 'igraph' is required for method = '", method, "'.")
  }

  rep <- .ic2_new_reporter("run.clustering", verbose, soda_root = x@storage$root_dir)

  .ic2_step(rep, paste0("[1/4] Loading ", reduction, " embeddings"))
  emb <- .ic2_get_plot_embedding(x, modality, reduction)
  ## Log concrete embedding path when available
  red_ref <- x@modalities[[modality]]$reductions[[reduction]]
  if (!is.null(red_ref$cell_embeddings) && nzchar(as.character(red_ref$cell_embeddings)[1L])) {
    rep <- .ic2_add_input(rep, paste0(toupper(reduction), " emb"), red_ref$cell_embeddings)
  }
  cell_ids <- rownames(emb)
  n_cells <- nrow(emb)
  if (n_cells < 3L) stop("Need at least 3 cells to cluster.")

  ## Resolve dims
  if (is.null(dims)) {
    if (reduction %in% c("pca", "aligned_pca")) {
      pca <- x@modalities[[modality]]$reductions[[reduction]]
      opt <- x@misc$opt.pcs %||% pca$opt.pcs %||% min(10L, ncol(emb))
      dims <- seq_len(min(as.integer(opt), ncol(emb)))
    } else {
      dims <- seq_len(ncol(emb))
    }
  }
  dims <- as.integer(dims)
  if (max(dims) > ncol(emb)) {
    stop("Requested dims exceed available ", reduction, " columns (", ncol(emb), ").")
  }
  emb <- emb[, dims, drop = FALSE]

  ## Sketch decision
  .ic2_step(rep, "[2/4] Choosing cells to cluster")
  fit_ids <- cell_ids
  sketch_mode <- FALSE
  if (isTRUE(use.sketch) && n_cells > max.cells.direct) {
    sk_ids <- NULL
    sk <- x@modalities[[modality]]$sketches$sketch_models[[sketch.name]]
    if (!is.null(sk$path) && file.exists(sk$path)) {
      meta <- readRDS(sk$path)
      sk_ids <- intersect(meta$cell_ids, cell_ids)
    } else if (!is.null(x@modalities[[modality]]$sketches$sketch_ids)) {
      sk_ids <- intersect(x@modalities[[modality]]$sketches$sketch_ids, cell_ids)
    }
    if (length(sk_ids) >= max(k + 2L, 50L) && length(sk_ids) < n_cells) {
      fit_ids <- sk_ids
      sketch_mode <- TRUE
    } else {
      set.seed(seed)
      nsk <- min(max.cells.direct, n_cells)
      fit_ids <- cell_ids[sort(sample.int(n_cells, nsk))]
      sketch_mode <- TRUE
      .ic2_step(rep, paste0(
        "  Auto-sketching ", .ic2_fmt_num(nsk),
        " / ", .ic2_fmt_num(n_cells), " cells for clustering"
      ))
    }
  } else if (!isTRUE(use.sketch) && n_cells > max.cells.direct) {
    warning(
      "n_cells (", n_cells, ") > max.cells.direct (", max.cells.direct,
      "). Consider use.sketch = TRUE."
    )
  }

  emb_fit <- emb[fit_ids, , drop = FALSE]
  rep <- .ic2_banner_start(rep, c(
    paste0("Reduction: ", reduction, "   |   dims: ", paste(range(dims), collapse = ":")),
    paste0("Method: ", method, "   |   k=", k),
    paste0(
      "Fit on: ", .ic2_fmt_num(length(fit_ids)),
      if (sketch_mode) " (sketch  ->  assign)" else " (direct)",
      "   |   total: ", .ic2_fmt_num(n_cells)
    )
  ))

  .ic2_step(rep, paste0("[3/4] Clustering (", method, ")"))
  set.seed(seed)
  if (identical(method, "phenograph")) {
    pg <- Rphenograph(emb_fit, k = k, verbose = FALSE)
    memb <- as.integer(igraph::membership(pg$community))
    names(memb) <- fit_ids
    rep <- .ic2_progress_bar(rep, 1L, 2L, prefix = "  cluster")
  } else if (identical(method, "louvain")) {
    lv <- .ic2_louvain_knn(emb_fit, k = k, dist.method = dist.method)
    memb <- as.integer(igraph::membership(lv$community))
    names(memb) <- fit_ids
    rep <- .ic2_progress_bar(rep, 1L, 2L, prefix = "  cluster")
  } else {
    ## kmeans
    kk <- n.clust
    if (is.null(kk)) {
      ## silhouette-ish search on a subsample of the fit set
      max.clust <- as.integer(min(as.integer(max.clust), nrow(emb_fit) - 1L))
      max.clust <- max(2L, max.clust)
      search_n <- min(nrow(emb_fit), 5000L)
      set.seed(seed)
      sub <- if (nrow(emb_fit) > search_n) {
        emb_fit[sort(sample.int(nrow(emb_fit), search_n)), , drop = FALSE]
      } else emb_fit
      best_k <- 2L
      best_ss <- -Inf
      for (cand in 2:max.clust) {
        km <- stats::kmeans(sub, centers = cand, nstart = 5L, iter.max = 50L)
        ## simple score: between/total
        score <- km$betweenss / max(km$totss, .Machine$double.eps)
        if (score > best_ss) {
          best_ss <- score
          best_k <- cand
        }
        rep <- .ic2_progress_bar(rep, cand - 1L, max.clust - 1L, prefix = "  k-search")
      }
      kk <- best_k
      .ic2_step(rep, paste0("  Selected kmeans k=", kk))
    }
    km <- stats::kmeans(emb_fit, centers = as.integer(kk), nstart = 20L, iter.max = 100L)
    memb <- as.integer(km$cluster)
    names(memb) <- fit_ids
  }

  n_clust_fit <- length(unique(memb))
  .ic2_step(rep, paste0(
    "  Found ", n_clust_fit, " clusters on ",
    .ic2_fmt_num(length(fit_ids)), " cells"
  ))

  .ic2_step(rep, "[4/4] Writing cluster_id")
  if (sketch_mode) {
    labels <- .ic2_assign_clusters_knn(
      emb, fit_ids = fit_ids, fit_labels = memb, k = assign.k
    )
    rep <- .ic2_progress_bar(rep, 2L, 2L, prefix = "  assign")
  } else {
    labels <- as.character(memb)
    names(labels) <- fit_ids
  }

  ## Align to cell_data order
  cd_ids <- rownames(x@cell_data)
  common <- intersect(cd_ids, names(labels))
  if (length(common) < 2L) stop("Cluster labels do not match cell_data.")
  clust_num <- as.integer(as.character(labels[cd_ids]))
  ## Relabel to 1..N contiguous by frequency
  tab <- sort(table(clust_num[!is.na(clust_num)]), decreasing = TRUE)
  map <- setNames(seq_along(tab), names(tab))
  cluster_id <- as.integer(map[as.character(clust_num)])
  names(cluster_id) <- cd_ids

  ## RNA clusters stay in cluster_id. A later ATAC clustering must not replace them.
  cl_col <- if (identical(modality, "rna")) "cluster_id" else paste0(modality, "_cluster_id")
  write_main <- identical(cl_col, "cluster_id") || !"cluster_id" %in% colnames(x@cell_data)
  x@cell_data[[cl_col]] <- cluster_id[cd_ids]
  if (isTRUE(write_main)) {
    x@cell_data$cluster_id <- cluster_id[cd_ids]
    x@misc$cluster_id_original <- cluster_id[cd_ids]
    x@cell_data$cluster_space <- if (isTRUE(src$aligned)) "Aligned" else "PCA"
  }
  .ic2_persist_cell_data(x)

  ## Compat + SODA metadata
  best <- data.frame(
    clusters = as.integer(x@cell_data[[cl_col]]),
    row.names = cd_ids,
    stringsAsFactors = FALSE
  )
  if (isTRUE(write_main)) x@misc$best.clust <- best

  cl_dir <- file.path(x@storage$root_dir, "metadata", "clustering")
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
  cl_stem <- file.path(cl_dir, paste0(modality, "_cluster_id"))
  cl_path <- .ic2_soda_write_table(
    data.frame(
      cell_id = cd_ids,
      clusters = as.integer(x@cell_data[[cl_col]]),
      stringsAsFactors = FALSE
    ),
    cl_stem
  )
  unlink(c(
    paste0(cl_stem, ".rds"),
    file.path(cl_dir, paste0(modality, "_clustering_meta.rds"))
  ), force = TRUE)

  meta <- list(
    modality = modality,
    reduction = reduction,
    based_on = reduction,
    aligned = isTRUE(src$aligned),
    mark = src$mark,
    dims = dims,
    method = method,
    k = k,
    sketch_mode = sketch_mode,
    n_fit = length(fit_ids),
    n_cells = n_cells,
    n_clusters = length(unique(stats::na.omit(cluster_id))),
    seed = seed,
    path = cl_path,
    created = as.character(Sys.time())
  )
  meta_path <- file.path(cl_dir, "clustering.json")
  .ic2_soda_write_json(meta, meta_path)
  x@modalities[[modality]]$clustering <- meta
  if (is.null(x@storage$metadata_paths)) x@storage$metadata_paths <- list()
  x@storage$metadata_paths$clustering <- cl_path
  x@storage$metadata_paths$clustering_meta <- meta_path

  x <- .ic2_log_command(
    x, "run.clustering",
    arguments = list(
      modality = modality, reduction = reduction, dims = dims,
      method = method, k = k, use.sketch = use.sketch,
      sketch_mode = sketch_mode
    ),
    outputs = list(
      n_clusters = meta$n_clusters,
      n_cells = n_cells,
      path = cl_path
    )
  )

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  rep <- .ic2_add_output(rep, "Clusters", cl_path)
  .ic2_finish(rep, extra_lines = paste0(
    "Clusters: ", meta$n_clusters,
    "   |   cells: ", .ic2_fmt_num(n_cells),
    if (sketch_mode) "   |   sketch -> assign" else ""
  ))
  x
}

#' @rdname run.clustering
#' @export
iclust <- run.clustering
