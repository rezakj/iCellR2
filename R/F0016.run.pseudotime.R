#' Pseudotime via graph diffusion (DPT) with soft branches
#'
#' Builds a kNN affinity graph on PCA (or UMAP) embeddings, computes
#' diffusion / graph pseudotime from a user-chosen root, and optionally
#' soft branch probabilities toward terminal tips (absorbing random walk).
#'
#' For large objects (\code{n > max.cells.direct}), the graph is fit on a
#' sketch (or random subsample) and remaining cells inherit pseudotime /
#' branch probs from nearest sketch neighbors.
#'
#' @param x An \code{iCellR2} object after \code{\link{run.pca}}
#'   (and ideally \code{\link{run.clustering}} if using cluster roots/tips).
#' @param modality Modality. Default \code{"rna"}.
#' @param reduction Embedding to use: \code{"pca"} (default) or \code{"umap"}.
#' @param use Which stored reduction to read. \code{"auto"} uses an aligned
#'   embedding when present, otherwise \code{reduction}.
#' @param dims Dimensions of \code{reduction} to use. Default \code{NULL}:
#'   for PCA uses \code{1:opt.pcs} (or \code{1:10}); for UMAP uses all.
#' @param root.cells Character cell ids that define the trajectory start.
#'   Pseudotime 0 is the mean embedding of these cells' nearest sketch /
#'   fit root. Mutually exclusive with \code{root.cluster}.
#' @param root.cluster Cluster id(s) in \code{x@cell_data$cluster_id} used
#'   as the root (centroid of those cells). Mutually exclusive with
#'   \code{root.cells}.
#' @param tip.cells Optional character cell ids for terminal fates
#'   (soft branches). If \code{NULL} and \code{soft.branches = TRUE},
#'   tips are auto-picked as clusters farthest from the root in
#'   pseudotime (requires \code{cluster_id}).
#' @param tip.clusters Optional cluster id(s) used as tips.
#' @param n.neighbors k for the kNN graph. Default \code{30}.
#' @param n.dcs Number of diffusion components (when eigen-decomposition
#'   is available). Default \code{10}.
#' @param soft.branches If \code{TRUE} (default), compute absorption
#'   probabilities toward tips and write \code{branch_id} +
#'   \code{branch_prob_*}.
#' @param n.tips Number of auto tips when none are supplied.
#'   Default \code{2}.
#' @param use.sketch Prefer an existing sketch for the fit graph when
#'   present. Default \code{TRUE}.
#' @param sketch.name Sketch name from \code{\link{run.sketch}}.
#' @param max.cells.direct Fit on all cells up to this size; above it,
#'   sketch / subsample. Default \code{100000}.
#' @param assign.k Neighbors used to project non-fit cells.
#'   Default \code{15}.
#' @param seed RNG seed for subsampling. Default \code{1}.
#' @param verbose Progress style.
#'
#' @return Updated object. Writes:
#'   \itemize{
#'     \item \code{x@cell_data$Pseudotime} (0-1 scaled)
#'     \item \code{x@cell_data$branch_id} (if soft branches)
#'     \item \code{x@cell_data$branch_prob_<tip>} columns
#'     \item SODA \code{reductions/pseudotime/} embeddings + JSON
#'   }
#'
#' @examples
#' \dontrun{
#' obj <- run.pca(obj)
#' obj <- run.clustering(obj)
#' obj <- run.pseudotime(obj, root.cluster = 1)
#' obj <- run.pseudotime(obj, root.cells = "AAACCTGAGAAACCAT-1",
#'                       tip.clusters = c(5, 8))
#' }
#'
#' @export
run.pseudotime <- function(
    x,
    modality = "rna",
    reduction = c("pca", "umap"),
    use = "auto",
    dims = NULL,
    root.cells = NULL,
    root.cluster = NULL,
    tip.cells = NULL,
    tip.clusters = NULL,
    n.neighbors = 30L,
    n.dcs = 10L,
    soft.branches = TRUE,
    n.tips = 2L,
    use.sketch = TRUE,
    sketch.name = "sketch1",
    max.cells.direct = 100000L,
    assign.k = 15L,
    seed = 1L,
    verbose = TRUE
) {
  .ic2_require_object(x)
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required for run.pseudotime().")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for run.pseudotime().")
  }
  modality <- tolower(modality)
  reduction <- match.arg(reduction)
  src <- .ic2_resolve_use(x, modality, reduction, use)
  reduction <- src$name
  pt_name <- .ic2_aligned_out_name("pseudotime", src$aligned)
  n.neighbors <- as.integer(n.neighbors)
  n.dcs <- as.integer(n.dcs)
  n.tips <- as.integer(n.tips)
  max.cells.direct <- as.integer(max.cells.direct)
  assign.k <- as.integer(assign.k)

  if (is.null(root.cells) && is.null(root.cluster)) {
    stop("Provide root.cells or root.cluster to define trajectory start.")
  }
  if (!is.null(root.cells) && !is.null(root.cluster)) {
    stop("Provide only one of root.cells or root.cluster.")
  }

  rep <- .ic2_new_reporter("run.pseudotime", verbose, soda_root = x@storage$root_dir)

  ## ---- embeddings ----
  emb <- .ic2_pseudotime_load_embedding(x, modality, reduction, dims)
  cell_ids <- rownames(emb)
  cd_ids <- rownames(x@cell_data)
  common <- intersect(cd_ids, cell_ids)
  if (length(common) < 10L) {
    stop("Fewer than 10 cells share embedding and cell_data ids.")
  }
  emb <- emb[common, , drop = FALSE]
  cell_ids <- common
  n_cells <- length(cell_ids)
  .ic2_step(rep, paste0(
    "Embedding: ", reduction, " (", ncol(emb), " dims)  |  ",
    .ic2_fmt_num(n_cells), " cells"
  ))

  ## ---- resolve root / tips on full id space ----
  root_ids <- .ic2_pseudotime_resolve_ids(
    x, root.cells = root.cells, clusters = root.cluster, cell_ids = cell_ids,
    label = "root"
  )
  tip_ids <- character(0)
  if (isTRUE(soft.branches)) {
    if (!is.null(tip.cells) || !is.null(tip.clusters)) {
      tip_ids <- .ic2_pseudotime_resolve_ids(
        x, root.cells = tip.cells, clusters = tip.clusters, cell_ids = cell_ids,
        label = "tip"
      )
    }
  }

  ## ---- choose fit set (sketch / subsample / all) ----
  fit_ids <- .ic2_pseudotime_fit_ids(
    x, modality = modality, cell_ids = cell_ids,
    use.sketch = use.sketch, sketch.name = sketch.name,
    max.cells.direct = max.cells.direct, seed = seed,
    must_include = unique(c(root_ids, tip_ids))
  )
  sketch_mode <- !setequal(fit_ids, cell_ids)
  .ic2_step(rep, paste0(
    "Fit graph on ", .ic2_fmt_num(length(fit_ids)),
    if (sketch_mode) " cells (sketch/subsample)" else " cells (full)"
  ))

  emb_fit <- emb[fit_ids, , drop = FALSE]
  root_fit <- intersect(root_ids, fit_ids)
  if (!length(root_fit)) {
    ## map each root to nearest fit cell
    root_fit <- .ic2_pseudotime_nearest_fit(emb, root_ids, fit_ids)
  }
  tip_fit <- intersect(tip_ids, fit_ids)
  if (length(tip_ids) && !length(tip_fit)) {
    tip_fit <- .ic2_pseudotime_nearest_fit(emb, tip_ids, fit_ids)
  }

  ## ---- kNN affinity + transition ----
  .ic2_step(rep, paste0("Building kNN affinity graph (k=", n.neighbors, ")"))
  Aff <- .ic2_knn_affinity(emb_fit, k = n.neighbors)
  Trans <- .ic2_row_stochastic(Aff)

  root_idx <- match(root_fit[1L], fit_ids)
  if (is.na(root_idx)) stop("Internal error: root index missing.")

  ## ---- diffusion / geodesic pseudotime on fit ----
  .ic2_step(rep, "Computing diffusion / graph pseudotime from root")
  pt_fit <- .ic2_graph_dpt(Aff, Trans, root_idx = root_idx, n_dcs = n.dcs)
  names(pt_fit) <- fit_ids

  ## ---- auto tips if needed ----
  if (isTRUE(soft.branches) && !length(tip_fit)) {
    tip_fit <- .ic2_pseudotime_auto_tips(
      x, fit_ids = fit_ids, pt_fit = pt_fit, root_fit = root_fit,
      n.tips = n.tips
    )
    .ic2_step(rep, paste0(
      "Auto tips: ", paste(tip_fit, collapse = ", ")
    ))
  }

  ## ---- soft branch absorption probs ----
  branch_fit <- NULL
  tip_names <- character(0)
  if (isTRUE(soft.branches) && length(tip_fit) >= 1L) {
    .ic2_step(rep, paste0(
      "Soft branches: absorption probs  ->  ", length(tip_fit), " tip(s)"
    ))
    tip_idx <- match(tip_fit, fit_ids)
    tip_idx <- tip_idx[!is.na(tip_idx)]
    if (length(tip_idx) >= 1L) {
      branch_fit <- .ic2_absorption_probs(Trans, tip_idx)
      rownames(branch_fit) <- fit_ids
      tip_names <- fit_ids[tip_idx]
      colnames(branch_fit) <- paste0("branch_prob_", seq_along(tip_names))
    }
  }

  ## ---- project to all cells ----
  .ic2_step(rep, "Projecting pseudotime to all cells")
  if (sketch_mode) {
    pt_all <- .ic2_pseudotime_project_knn(
      emb, fit_ids = fit_ids, values = pt_fit, k = assign.k
    )
    branch_all <- NULL
    if (!is.null(branch_fit)) {
      branch_all <- matrix(NA_real_, nrow = n_cells, ncol = ncol(branch_fit))
      rownames(branch_all) <- cell_ids
      colnames(branch_all) <- colnames(branch_fit)
      for (j in seq_len(ncol(branch_fit))) {
        branch_all[, j] <- .ic2_pseudotime_project_knn(
          emb, fit_ids = fit_ids, values = branch_fit[, j], k = assign.k
        )
      }
      ## renormalize rows
      rs <- rowSums(branch_all)
      rs[rs <= 0 | !is.finite(rs)] <- 1
      branch_all <- sweep(branch_all, 1L, rs, "/")
    }
  } else {
    pt_all <- pt_fit[cell_ids]
    branch_all <- if (is.null(branch_fit)) NULL else branch_fit[cell_ids, , drop = FALSE]
  }

  ## scale PT to [0, 1]
  pt_all <- as.numeric(pt_all)
  names(pt_all) <- cell_ids
  ok <- is.finite(pt_all)
  if (any(ok)) {
    rng <- range(pt_all[ok])
    if (diff(rng) > 0) {
      pt_all[ok] <- (pt_all[ok] - rng[1L]) / diff(rng)
    } else {
      pt_all[ok] <- 0
    }
  }
  pt_all[!ok] <- NA_real_

  ## branch_id = argmax
  branch_id <- rep(NA_integer_, n_cells)
  names(branch_id) <- cell_ids
  if (!is.null(branch_all)) {
    branch_id[] <- max.col(branch_all, ties.method = "first")
  }

  ## ---- write cell_data ----
  .ic2_step(rep, "Writing Pseudotime / branch columns to cell_data")
  x@cell_data$Pseudotime <- NA_real_
  x@cell_data[cell_ids, "Pseudotime"] <- pt_all[cell_ids]
  if (!is.null(branch_all)) {
    ## drop old branch_prob_* columns
    old <- grep("^branch_prob_", colnames(x@cell_data), value = TRUE)
    if (length(old)) {
      x@cell_data <- x@cell_data[, setdiff(colnames(x@cell_data), old), drop = FALSE]
    }
    x@cell_data$branch_id <- NA_integer_
    x@cell_data[cell_ids, "branch_id"] <- as.integer(branch_id[cell_ids])
    for (j in seq_len(ncol(branch_all))) {
      cn <- colnames(branch_all)[j]
      x@cell_data[[cn]] <- NA_real_
      x@cell_data[cell_ids, cn] <- as.numeric(branch_all[cell_ids, j])
    }
  }
  x <- .ic2_persist_cell_data(x)

  ## ---- SODA reductions/pseudotime ----
  red_dir <- file.path(x@storage$root_dir, "reductions", pt_name)
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_mat <- cbind(Pseudotime = pt_all[cell_ids])
  if (!is.null(branch_all)) {
    emb_mat <- cbind(emb_mat, branch_all[cell_ids, , drop = FALSE])
  }
  rownames(emb_mat) <- cell_ids
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  .ic2_h5_write_dense(emb_mat, emb_path)
  unlink(file.path(red_dir, "cell_embeddings.rds"), force = TRUE)

  tip_table <- NULL
  if (length(tip_fit)) {
    tip_table <- data.frame(
      tip_index = seq_along(tip_fit),
      tip_cell = tip_fit,
      stringsAsFactors = FALSE
    )
    utils::write.table(
      tip_table,
      file = file.path(red_dir, "tips.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }

  meta <- list(
    name = pt_name,
    aligned = isTRUE(src$aligned),
    mark = src$mark,
    modality = modality,
    reduction = reduction,
    based_on = reduction,
    dims = as.integer(attr(emb, "dims_used") %||% seq_len(ncol(emb))),
    n_neighbors = n.neighbors,
    n_dcs = n.dcs,
    soft_branches = isTRUE(soft.branches) && !is.null(branch_all),
    n_tips = length(tip_fit),
    root_cells = root_fit,
    tip_cells = tip_fit,
    sketch_mode = sketch_mode,
    n_fit = length(fit_ids),
    n_cells = n_cells,
    method = "graph_dpt",
    created = as.character(Sys.time()),
    paths = list(
      embeddings = emb_path,
      tips = if (length(tip_fit)) file.path(red_dir, "tips.tsv") else NULL
    )
  )
  meta_path <- file.path(red_dir, "reduction.json")
  .ic2_soda_write_json(meta, meta_path)

  x@modalities[[modality]]$reductions[[pt_name]] <- list(
    cell_embeddings = emb_path,
    meta_path = meta_path,
    root_cells = root_fit,
    tip_cells = tip_fit,
    sketch_mode = sketch_mode
  )
  x@storage$reduction_paths[[pt_name]] <- red_dir
  try(.ic2_soda_write_manifest(x), silent = TRUE)

  if (is.null(x@params$trajectory)) x@params$trajectory <- list()
  x@params$trajectory[[modality]] <- list(
    method = "graph_dpt",
    reduction = reduction,
    soft_branches = isTRUE(soft.branches) && !is.null(branch_all),
    n_tips = length(tip_fit),
    sketch_mode = sketch_mode,
    n_fit = length(fit_ids)
  )

  x <- .ic2_log_command(
    x, "run.pseudotime",
    arguments = list(
      modality = modality, reduction = reduction,
      root.cells = root.cells, root.cluster = root.cluster,
      tip.clusters = tip.clusters, n.neighbors = n.neighbors,
      soft.branches = soft.branches, sketch_mode = sketch_mode
    ),
    outputs = list(
      embeddings = emb_path,
      n_cells = n_cells,
      n_tips = length(tip_fit)
    )
  )

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  rep <- .ic2_add_output(rep, "Pseudotime", emb_path)
  .ic2_finish(rep, extra_lines = c(
    paste0("Cells: ", .ic2_fmt_num(n_cells),
           "  |  fit: ", .ic2_fmt_num(length(fit_ids))),
    paste0("Root: ", paste(utils::head(root_fit, 3L), collapse = ", "),
           if (length(root_fit) > 3L) "..." else ""),
    if (length(tip_fit)) {
      paste0("Tips: ", length(tip_fit), "  |  soft branches written")
    } else {
      "Soft branches: off / none"
    }
  ))
  x
}

# ---- internals -------------------------------------------------------------

#' @keywords internal
.ic2_pseudotime_load_embedding <- function(x, modality, reduction, dims) {
  if (grepl("pca", reduction)) {
    pca <- x@modalities[[modality]]$reductions[[reduction]]
    if (is.null(pca)) stop("No '", reduction, "' reduction. Run run.pca() or iba() first.")
    emb <- .ic2_get_plot_embedding(x, modality, reduction)
    npcs <- ncol(emb)
    if (is.null(dims)) {
      opt <- x@misc$opt.pcs %||% pca$opt.pcs %||% min(10L, npcs)
      dims <- seq_len(min(as.integer(opt), npcs))
    } else {
      dims <- as.integer(dims)
      dims <- dims[dims >= 1L & dims <= npcs]
    }
    if (length(dims) < 2L) stop("'dims' must select at least 2 PCs.")
    out <- emb[, dims, drop = FALSE]
    attr(out, "dims_used") <- dims
    return(out)
  }
  ## umap
  uref <- x@modalities[[modality]]$reductions[[reduction]]
  if (is.null(uref)) stop("No '", reduction, "' found. Run run.umap() or use reduction = \"pca\".")
  path <- uref$cell_embeddings %||% uref$paths$embeddings
  if (is.null(path) || !file.exists(path)) {
    stop("UMAP embeddings not found on disk.")
  }
  emb <- if (grepl("\\.(h5|hdf5)$", path, ignore.case = TRUE)) {
    .ic2_h5_read_dense(path)
  } else {
    readRDS(path)
  }
  if (!is.matrix(emb)) emb <- as.matrix(emb)
  if (is.null(dims)) {
    dims <- seq_len(ncol(emb))
  } else {
    dims <- as.integer(dims)
    dims <- dims[dims >= 1L & dims <= ncol(emb)]
  }
  out <- emb[, dims, drop = FALSE]
  attr(out, "dims_used") <- dims
  out
}

#' @keywords internal
.ic2_pseudotime_resolve_ids <- function(x, root.cells = NULL, clusters = NULL,
                                        cell_ids, label = "root") {
  if (!is.null(root.cells)) {
    ids <- intersect(as.character(root.cells), cell_ids)
    if (!length(ids)) {
      stop("No ", label, ".cells found among embedding cell ids.")
    }
    return(unique(ids))
  }
  if (is.null(clusters)) {
    stop("Internal: clusters required for ", label)
  }
  if (is.null(x@cell_data$cluster_id)) {
    stop(label, ".cluster requires x@cell_data$cluster_id. Run run.clustering() first.")
  }
  cl <- as.character(clusters)
  hit <- rownames(x@cell_data)[as.character(x@cell_data$cluster_id) %in% cl]
  ids <- intersect(hit, cell_ids)
  if (!length(ids)) {
    stop("No cells for ", label, ".cluster = ", paste(cl, collapse = ","))
  }
  unique(ids)
}

#' @keywords internal
.ic2_pseudotime_fit_ids <- function(x, modality, cell_ids, use.sketch,
                                    sketch.name, max.cells.direct, seed,
                                    must_include = character(0)) {
  must_include <- intersect(as.character(must_include), cell_ids)
  fit_ids <- NULL
  if (isTRUE(use.sketch)) {
    sk <- x@modalities[[modality]]$sketches$sketch_models[[sketch.name]]
    if (!is.null(sk$path) && file.exists(sk$path)) {
      meta <- readRDS(sk$path)
      fit_ids <- intersect(meta$cell_ids, cell_ids)
    } else if (!is.null(x@modalities[[modality]]$sketches$sketch_ids)) {
      fit_ids <- intersect(x@modalities[[modality]]$sketches$sketch_ids, cell_ids)
    }
  }
  if (is.null(fit_ids) || length(fit_ids) < 50L) {
    if (length(cell_ids) <= max.cells.direct) {
      fit_ids <- cell_ids
    } else {
      set.seed(seed)
      fit_ids <- cell_ids[sort(sample.int(length(cell_ids), max.cells.direct))]
    }
  }
  unique(c(must_include, fit_ids))
}

#' @keywords internal
.ic2_pseudotime_nearest_fit <- function(emb, query_ids, fit_ids) {
  q <- intersect(query_ids, rownames(emb))
  if (!length(q)) return(fit_ids[1L])
  nn <- FNN::get.knnx(emb[fit_ids, , drop = FALSE], emb[q, , drop = FALSE], k = 1L)
  unique(fit_ids[as.integer(nn$nn.index[, 1L])])
}

#' Adaptive Gaussian kNN affinity (symmetric)
#' @keywords internal
.ic2_knn_affinity <- function(X, k = 30L) {
  k <- as.integer(min(k, nrow(X) - 1L))
  if (k < 1L) stop("Need at least 2 cells for kNN graph.")
  kn <- FNN::get.knn(X, k = k)
  n <- nrow(X)
  ii <- rep(seq_len(n), each = k)
  jj <- as.vector(t(kn$nn.index))
  dd <- as.vector(t(kn$nn.dist))
  sigma <- apply(kn$nn.dist, 1L, stats::median)
  sigma[!is.finite(sigma) | sigma <= 0] <- 1e-6
  ww <- exp(-(dd * dd) / (sigma[ii] * sigma[ii]))
  A <- Matrix::sparseMatrix(i = ii, j = jj, x = ww, dims = c(n, n))
  A <- (A + Matrix::t(A)) / 2
  ## keep sparse: zero diagonal without densifying
  if (any(Matrix::diag(A) != 0)) {
    A <- A - Matrix::Diagonal(x = Matrix::diag(A))
  }
  A
}

#' @keywords internal
.ic2_row_stochastic <- function(A) {
  deg <- Matrix::rowSums(A)
  deg[!is.finite(deg) | deg <= 0] <- 1
  Matrix::Diagonal(x = 1 / as.numeric(deg)) %*% A
}

#' Graph / diffusion pseudotime from root index (1-based)
#' @keywords internal
.ic2_graph_dpt <- function(Aff, Trans, root_idx, n_dcs = 10L) {
  n <- nrow(Aff)
  root_idx <- as.integer(root_idx)[1L]
  ## Prefer diffusion distance via leading eigenvectors of symmetric operator
  pt <- NULL
  if (requireNamespace("irlba", quietly = TRUE) && n >= 20L) {
    deg <- Matrix::rowSums(Aff)
    deg[!is.finite(deg) | deg <= 0] <- 1
    dinv_sqrt <- 1 / sqrt(as.numeric(deg))
    ## Symmetric normalized affinity: D^{-1/2} A D^{-1/2}
    S <- Matrix::Diagonal(x = dinv_sqrt) %*% Aff %*% Matrix::Diagonal(x = dinv_sqrt)
    n_dcs <- as.integer(min(n_dcs, n - 2L, 50L))
    if (n_dcs >= 2L) {
      ## irlba SVD ~ eigen for symmetric PSD S
      ev <- tryCatch(
        irlba::irlba(S, nv = n_dcs + 1L, nu = n_dcs + 1L),
        error = function(e) NULL
      )
      if (!is.null(ev) && !is.null(ev$u) && ncol(ev$u) >= 2L) {
        ## drop leading component; diffusion distance to root
        vecs <- ev$u[, -1L, drop = FALSE]
        vals <- ev$d[-1L]
        vals <- pmax(as.numeric(vals), 0)
        scaled <- sweep(vecs, 2L, vals, `*`)
        root_v <- scaled[root_idx, , drop = FALSE]
        diff <- sweep(scaled, 2L, as.numeric(root_v), `-`)
        pt <- sqrt(rowSums(diff * diff))
      }
    }
  }
  ## Fallback: weighted shortest-path distance on affinity graph
  if (is.null(pt) || !all(is.finite(pt))) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package 'igraph' is required for geodesic pseudotime fallback.")
    }
    ## edge weight = 1/affinity for distance
    A2 <- Aff
    A2@x <- 1 / pmax(A2@x, 1e-12)
    g <- igraph::graph_from_adjacency_matrix(A2, mode = "undirected", weighted = TRUE, diag = FALSE)
    pt <- as.numeric(igraph::distances(g, v = root_idx, to = igraph::V(g), weights = igraph::E(g)$weight))
  }
  pt[!is.finite(pt)] <- max(pt[is.finite(pt)], na.rm = TRUE)
  pt
}

#' Absorbing Markov chain probabilities toward tip indices
#' @keywords internal
.ic2_absorption_probs <- function(Trans, tip_idx, max_iter = 400L, tol = 1e-7) {
  n <- nrow(Trans)
  tip_idx <- unique(as.integer(tip_idx))
  tip_idx <- tip_idx[tip_idx >= 1L & tip_idx <= n]
  n_tips <- length(tip_idx)
  if (n_tips < 1L) return(NULL)

  ## Build absorbing transition: tips are absorbing states
  Tm <- Trans
  Tm <- as(Tm, "dgCMatrix")
  for (t in tip_idx) {
    ## zero row then put 1 on diagonal
    Tm[t, ] <- 0
    Tm[t, t] <- 1
  }
  ## Fixed point: B = T_abs %*% B with tip rows = I
  B <- matrix(0, nrow = n, ncol = n_tips)
  for (j in seq_len(n_tips)) B[tip_idx[j], j] <- 1
  for (it in seq_len(max_iter)) {
    B_new <- as.matrix(Tm %*% B)
    B_new[tip_idx, ] <- 0
    for (j in seq_len(n_tips)) B_new[tip_idx[j], j] <- 1
    if (max(abs(B_new - B)) < tol) {
      B <- B_new
      break
    }
    B <- B_new
  }
  ## Normalize remaining rows
  rs <- rowSums(B)
  miss <- which(rs <= 1e-12 | !is.finite(rs))
  if (length(miss)) {
    B[miss, ] <- 1 / n_tips
    rs[miss] <- 1
  }
  B <- sweep(B, 1L, rs, "/")
  B
}

#' Project scalar values from fit cells to all cells via kNN
#' @keywords internal
.ic2_pseudotime_project_knn <- function(emb, fit_ids, values, k = 15L) {
  all_ids <- rownames(emb)
  out <- rep(NA_real_, length(all_ids))
  names(out) <- all_ids
  values <- values[fit_ids]
  out[fit_ids] <- as.numeric(values)
  query <- setdiff(all_ids, fit_ids)
  if (!length(query)) return(out)
  k <- as.integer(min(k, length(fit_ids)))
  nn <- FNN::get.knnx(
    emb[fit_ids, , drop = FALSE],
    emb[query, , drop = FALSE],
    k = k
  )
  ## inverse-distance weights
  d <- nn$nn.dist
  d[!is.finite(d) | d < 1e-9] <- 1e-9
  w <- 1 / d
  w <- w / rowSums(w)
  idx <- nn$nn.index
  pred <- rowSums(w * matrix(as.numeric(values)[idx], nrow = length(query), ncol = k))
  out[query] <- pred
  out
}

#' Auto-select tip cells: farthest clusters from root by mean PT
#' @keywords internal
.ic2_pseudotime_auto_tips <- function(x, fit_ids, pt_fit, root_fit, n.tips = 2L) {
  n.tips <- as.integer(max(1L, n.tips))
  if (!is.null(x@cell_data$cluster_id)) {
    cl <- as.character(x@cell_data[fit_ids, "cluster_id"])
    root_cl <- unique(as.character(x@cell_data[root_fit, "cluster_id"]))
    tab <- tapply(pt_fit[fit_ids], cl, mean, na.rm = TRUE)
    tab <- tab[is.finite(tab)]
    tab <- tab[setdiff(names(tab), root_cl)]
    if (length(tab)) {
      ord <- names(sort(tab, decreasing = TRUE))
      pick_cl <- utils::head(ord, n.tips)
      tips <- character(0)
      for (pc in pick_cl) {
        ids <- fit_ids[cl == pc]
        tips <- c(tips, ids[which.max(pt_fit[ids])])
      }
      return(unique(tips))
    }
  }
  ## no clusters: farthest cells by PT, with simple diversity
  ord <- fit_ids[order(pt_fit[fit_ids], decreasing = TRUE)]
  utils::head(ord, n.tips)
}
