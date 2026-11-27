#' KNetL map (kNN graph embedding)
#'
#' Classic iCellR KNetL math: build a kNN graph on a source embedding
#' (default PCA / \code{aligned_pca}), lay it out with igraph
#' (\code{layout_nicely} by default), then optionally unpack that 2D layout
#' with UMAP or PCA (\code{do.redux = TRUE}).
#'
#' \code{zoom} is honored as-is (e.g. \code{zoom = 400} -> \eqn{k \approx 399}
#' neighbors). There is no auto-cap, no \code{force.zoom}, and no HNSW / ANN
#' substitute for classic exact neighbors.
#'
#' Neighbors:
#' \itemize{
#'   \item \strong{Small \eqn{n}} (\code{n <= max.cells.full.dist}, default
#'     5000): classic full \eqn{n \times n} \code{stats::dist} with
#'     \code{order(...)[2:zoom]} and linear-index weights \code{D[MyOrd]}.
#'   \item \strong{Larger \eqn{n}} + Euclidean: chunked \pkg{FNN}
#'     \code{get.knnx} (exact Euclidean kNN; no \eqn{n \times n} matrix).
#'     Neighbor index/distances may be streamed to SODA
#'     \code{graphs/knetl/}. Edge weights are true Euclidean distances
#'     (not the classic \code{D[MyOrd]} quirk).
#'   \item Non-Euclidean distances above the full-dist threshold are refused.
#' }
#'
#' Same-math speedups (do not change the algorithm): parallel chunked FNN
#' via \code{n_threads}, preallocated edge buffers, integer-vertex igraph
#' construction, SODA knn reuse when n/k match. At atlas scale
#' (\eqn{n \sim 10^6}, \code{zoom = 400}) \code{layout_nicely} on
#' \eqn{\sim n\cdot(\mathrm{zoom}-1)} edges remains day-scale -- that cost is
#' inherent to classic math, not an ANN tradeoff.
#'
#' Embeddings are written to \code{reductions/knetl}, or
#' \code{reductions/aligned_knetl} when the source reduction is aligned.
#'
#' @param x An \code{iCellR2} object after \code{\link{run.pca}} (or \code{\link{iba}}).
#' @param dist.method Distance method. Default \code{"euclidean"}. Full
#'   \code{stats::dist} is only used when \code{n <= max.cells.full.dist}.
#' @param zoom Neighborhood size, classic default \code{300}
#'   (non-joint edges are \code{2:zoom}). Honored as-is
#'   (\code{zoom = 400} uses \eqn{k \approx 399}).
#' @param data.type Source family: \code{"pca"} (default) or \code{"umap"}.
#'   Ignored when \code{use} names a reduction.
#' @param use \code{"auto"} (default) prefers \code{aligned_pca} /
#'   \code{aligned_umap}. Or \code{"pca"}, \code{"aligned_pca"},
#'   \code{"reductions/aligned_pca"}.
#' @param dims Columns of the source embedding. Default \code{1:20}.
#' @param joint If \code{TRUE}, take \code{zoom} neighbors inside every
#'   condition (\code{condition.by}), including self, as in classic CPCA.
#'   Default \code{FALSE}.
#' @param my.seed RNG seed for the igraph layout (classic). The nested UMAP
#'   unpack uses seed \code{0} to match classic \code{run.umap} defaults.
#' @param n_components UMAP components when \code{do.redux = TRUE}.
#'   Default \code{2}. \code{add.3d} stores a separate 3D igraph layout.
#' @param add.3d Store a 3-column KNetL from \code{layout.3d} on the same
#'   graph (classic), independent of the 2D / UMAP unpack.
#' @param do.redux If \code{TRUE} (default), UMAP (or PCA) on the scaled
#'   igraph layout. If \code{FALSE}, keep the scaled layout as KNetL.
#' @param min_dist,spread Passed to \pkg{uwot} when unpacking with UMAP
#'   (classic defaults \code{0.01} / \code{1}).
#' @param n_threads Threads for chunked FNN (parallel query chunks) and
#'   UMAP unpack. Default \code{4}.
#' @param max.graph.edges Cap on graph edges. Default \code{Inf} (classic
#'   has no refuse). A warning is issued when edges exceed \code{2e6}.
#' @param max.cells.full.dist Max cells for classic full \code{dist}.
#'   Default \code{5000}. Above this, Euclidean uses chunked exact FNN.
#' @param chunk.size Query cells per FNN chunk when streaming. Default
#'   \code{25000}.
#' @param return.graph If \code{TRUE}, return the igraph instead of writing
#'   a reduction.
#' @param col.by Vertex color when \code{return.graph = TRUE}:
#'   \code{"clusters"} or \code{"conditions"}.
#' @param run.iclust If \code{TRUE}, Louvain-cluster the KNetL map
#'   (\code{k = 30}). Default \code{FALSE}.
#' @param condition.by Condition column for \code{joint = TRUE}.
#' @param modality Modality. Default \code{"rna"}.
#' @param verbose Progress style.
#' @param layout.2d,layout.3d,dim.redux Classic layout / unpack choices.
#'   \code{dim.redux = "umap"} (default) or \code{"pca"}.
#'
#' @return Updated \code{iCellR2} with a \code{knetl} (or \code{aligned_knetl})
#'   reduction, or an igraph if \code{return.graph = TRUE}.
#'
#' @examples
#' \dontrun{
#' obj <- run.knetl(obj)                         # auto PCA, zoom 300
#' obj <- run.knetl(obj, zoom = 400, n_threads = 4L)  # classic math, k~=399
#' obj <- run.knetl(obj, use = "pca")
#' cluster.plot(obj, plot.type = "knetl")
#' }
#'
#' @export
run.knetl <- function(
    x = NULL,
    dist.method = "euclidean",
    zoom = 300,
    data.type = "pca",
    use = "auto",
    dims = 1:20,
    joint = FALSE,
    col.by = c("clusters", "conditions"),
    my.seed = 1,
    layout.2d = "layout_nicely",
    layout.3d = "layout_with_fr",
    add.3d = FALSE,
    dim.redux = c("umap", "pca"),
    do.redux = TRUE,
    run.iclust = FALSE,
    return.graph = FALSE,
    min_dist = 0.01,
    spread = 1,
    n_components = 2L,
    n_threads = 4L,
    max.graph.edges = Inf,
    max.cells.full.dist = 5000L,
    chunk.size = 25000L,
    condition.by = "condition_id",
    modality = "rna",
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  col.by <- match.arg(col.by)
  dim.redux <- match.arg(dim.redux)
  dist.method <- match.arg(dist.method, c(
    "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  ))
  zoom <- as.integer(zoom)
  if (!is.finite(zoom) || zoom < 3L) {
    stop("'zoom' must be >= 3 (classic neighbors are 2:zoom).")
  }
  n_components <- as.integer(n_components)
  if (!n_components %in% c(2L, 3L)) stop("'n_components' must be 2 or 3.")
  max.cells.full.dist <- as.integer(max.cells.full.dist)
  if (!is.finite(max.cells.full.dist) || max.cells.full.dist < 3L) {
    max.cells.full.dist <- 5000L
  }
  chunk.size <- as.integer(chunk.size)
  if (!is.finite(chunk.size) || chunk.size < 1L) chunk.size <- 25000L
  n_threads <- as.integer(n_threads)
  if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 4L
  if (is.null(max.graph.edges) || length(max.graph.edges) < 1L ||
      is.na(max.graph.edges[1L])) {
    max.graph.edges <- Inf
  } else {
    max.graph.edges <- as.numeric(max.graph.edges)[1L]
    if (!is.finite(max.graph.edges) && !is.infinite(max.graph.edges)) {
      max.graph.edges <- Inf
    }
  }

  family <- if (identical(tolower(data.type), "umap")) "umap" else "pca"
  src <- .ic2_resolve_use(x, modality, family, use)
  out_name <- .ic2_aligned_out_name("knetl", src$aligned)
  out_3d <- .ic2_aligned_out_name("knetl_3d", src$aligned)

  rep <- .ic2_new_reporter("run.knetl", verbose, soda_root = x@storage$root_dir)

  .ic2_step(rep, paste0("[1/4] Loading ", src$name))
  emb <- .ic2_get_plot_embedding(x, modality, src$name)
  dims <- as.integer(dims)
  dims <- dims[is.finite(dims) & dims >= 1L & dims <= ncol(emb)]
  if (length(dims) < 2L) stop("Need at least 2 source dimensions.")
  emb <- emb[, dims, drop = FALSE]
  cell_ids <- rownames(emb)
  n_cells <- nrow(emb)
  if (zoom >= n_cells) {
    stop("'zoom' (", zoom, ") must be < number of cells (", n_cells, ").")
  }
  X <- as.matrix(emb)
  if (any(!is.finite(X))) stop("Source embedding has non-finite values.")

  est_edges <- as.numeric(n_cells) * as.numeric(zoom - 1L)
  if (isTRUE(verbose) && est_edges > 2e6) {
    message(
      "run.knetl: classic path will build ~", .ic2_fmt_num(est_edges),
      " edges (n=", .ic2_fmt_num(n_cells), ", zoom=", zoom,
      "). Exact kNN + igraph layout_nicely at this scale is typically ",
      "multi-hour to day-scale; ANN/HNSW is intentionally not used."
    )
  }

  pipeline <- "classic_layout_then_redux"
  use_full_dist <- n_cells <= max.cells.full.dist
  nn_label <- if (use_full_dist) {
    paste0("classic full ", dist.method, " dist")
  } else {
    paste0("chunked exact FNN (", dist.method, ")")
  }

  rep <- .ic2_banner_start(rep, c(
    paste0("Source: ", src$name, if (isTRUE(src$aligned)) "   |   Aligned" else ""),
    paste0(
      "zoom = ", zoom, " (k \u2248 ", zoom - 1L, ")   |   joint = ", isTRUE(joint),
      "   |   pipeline = ", pipeline
    )
  ))

  .ic2_step(rep, paste0(
    "[2/4] kNN via ", nn_label, "   |   ", .ic2_fmt_num(n_cells), " cells"
  ))
  conditions <- if (isTRUE(joint)) {
    .ic2_iba_conditions(x, cell_ids, condition.by)
  } else {
    NULL
  }
  graph_dir <- file.path(x@storage$root_dir, "graphs", "knetl")
  nn <- .ic2_knetl_neighbors(
    X,
    zoom = zoom,
    joint = isTRUE(joint),
    dist.method = dist.method,
    conditions = conditions,
    max.cells.full.dist = max.cells.full.dist,
    chunk.size = chunk.size,
    soda_graph_dir = graph_dir,
    cell_ids = cell_ids,
    n_threads = n_threads,
    verbose = isTRUE(verbose)
  )

  .ic2_step(rep, paste0(
    "[3/4] Building kNN graph (", .ic2_fmt_num(length(nn$from)), " edges)"
  ))
  g <- .ic2_knetl_graph(nn$from, nn$to, nn$weight, cell_ids, max.graph.edges)
  if (isTRUE(return.graph)) {
    g <- .ic2_knetl_color_graph(g, x, cell_ids, col.by)
    .ic2_finish(rep, extra_lines = "Returned igraph (not written to SODA).")
    return(g)
  }

  .ic2_step(rep, paste0("[4/4] ", layout.2d, if (isTRUE(do.redux)) {
    paste0(" \u2192 ", dim.redux, " unpack")
  } else {
    " (no redux)"
  }))
  set.seed(as.integer(my.seed))
  lay2 <- .ic2_knetl_igraph_layout(g, layout.2d, dim = 2L)
  rownames(lay2) <- cell_ids
  colnames(lay2) <- paste0("V", seq_len(ncol(lay2)))
  ## Classic scales the igraph layout before (optional) UMAP/PCA unpack.
  lay2 <- as.matrix(scale(lay2))
  rownames(lay2) <- cell_ids
  colnames(lay2) <- paste0("V", seq_len(ncol(lay2)))

  if (isTRUE(do.redux)) {
    lay <- .ic2_knetl_unpack(
      lay2,
      dim.redux = dim.redux,
      n_components = n_components,
      min_dist = min_dist,
      spread = spread,
      n_threads = n_threads
    )
    rownames(lay) <- cell_ids
    colnames(lay) <- paste0("V", seq_len(ncol(lay)))
  } else {
    lay <- lay2
  }

  x <- .ic2_knetl_store(
    x, modality, lay, out_name, src,
    zoom = zoom, k = zoom - 1L, dims = dims, joint = isTRUE(joint),
    neighbor_method = nn$method, pipeline = pipeline
  )

  if (isTRUE(add.3d)) {
    ## Classic: separate 3D igraph layout on the same graph (not UMAP unpack).
    set.seed(as.integer(my.seed))
    lay3 <- .ic2_knetl_igraph_layout(g, layout.3d, dim = 3L)
    lay3 <- as.matrix(scale(lay3))
    rownames(lay3) <- cell_ids
    colnames(lay3) <- paste0("V", seq_len(ncol(lay3)))
    x <- .ic2_knetl_store(
      x, modality, lay3, out_3d, src,
      zoom = zoom, k = zoom - 1L, dims = dims, joint = isTRUE(joint),
      neighbor_method = nn$method, pipeline = pipeline
    )
  }

  if (isTRUE(run.iclust)) {
    .ic2_step(rep, "Clustering on KNetL")
    x <- run.clustering(
      x, modality = modality, data.type = out_name, dims = 1:2,
      k = 30L, seed = my.seed, verbose = verbose
    )
  }

  x <- .ic2_log_command(
    x, "run.knetl",
    arguments = list(
      zoom = as.integer(zoom),
      source = src$name, dims = as.integer(dims),
      joint = isTRUE(joint), do.redux = isTRUE(do.redux),
      dim.redux = dim.redux, layout.2d = layout.2d,
      neighbor_method = nn$method,
      pipeline = pipeline,
      max.cells.full.dist = as.integer(max.cells.full.dist),
      n_threads = as.integer(n_threads)
    ),
    outputs = list(reduction = file.path(x@storage$root_dir, "reductions", out_name))
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  .ic2_finish(rep, extra_lines = paste0(
    out_name, "   |   ", .ic2_fmt_num(n_cells), " cells   |   zoom = ", zoom,
    "   |   neighbors = ", nn$method,
    "   |   ", pipeline,
    if (isTRUE(src$aligned)) "   |   Aligned" else ""
  ))
  x
}


#' KNetL neighbor edges: full dist (small n) or chunked exact FNN (large n).
#'
#' @keywords internal
.ic2_knetl_neighbors <- function(
    X,
    zoom,
    joint = FALSE,
    dist.method = "euclidean",
    conditions = NULL,
    max.cells.full.dist = 5000L,
    chunk.size = 25000L,
    soda_graph_dir = NULL,
    cell_ids = NULL,
    n_threads = 4L,
    verbose = FALSE
) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  zoom <- as.integer(zoom)
  if (zoom >= n) stop("'zoom' must be < number of cells.")
  max.cells.full.dist <- as.integer(max.cells.full.dist)
  if (!is.finite(max.cells.full.dist) || max.cells.full.dist < 3L) {
    max.cells.full.dist <- 5000L
  }

  ## Try reuse of prior SODA knn sidecar (same n / k).
  reused <- .ic2_knetl_try_reuse_knn(
    soda_graph_dir, n = n, k = as.integer(zoom - 1L), joint = isTRUE(joint),
    need_edges = TRUE
  )
  if (!is.null(reused)) {
    if (isTRUE(verbose)) message("run.knetl: reusing graphs/knetl knn sidecar.")
    reused$method <- "reuse_soda"
    return(reused)
  }

  if (n <= max.cells.full.dist) {
    out <- .ic2_knetl_neighbors_full_dist(
      X, zoom = zoom, joint = isTRUE(joint),
      dist.method = dist.method, conditions = conditions
    )
    out$method <- "full_dist"
    return(out)
  }

  if (!identical(dist.method, "euclidean")) {
    stop(
      "dist.method '", dist.method, "' needs a full n\u00d7n distance matrix. ",
      "Use 'euclidean' for large objects, or set max.cells.full.dist >= n (",
      n, ")."
    )
  }

  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required for large-n run.knetl() neighbor search.")
  }

  out <- .ic2_knetl_neighbors_fnn(
    X,
    zoom = zoom,
    joint = isTRUE(joint),
    conditions = conditions,
    chunk.size = as.integer(chunk.size),
    soda_graph_dir = soda_graph_dir,
    cell_ids = cell_ids,
    n_threads = as.integer(n_threads),
    verbose = isTRUE(verbose)
  )
  out$method <- "fnn_chunked"
  out
}


#' Classic full-dist KNetL edges (small n only).
#' @keywords internal
.ic2_knetl_neighbors_full_dist <- function(
    X, zoom, joint = FALSE, dist.method = "euclidean", conditions = NULL
) {
  n <- nrow(X)
  D <- as.matrix(stats::dist(X, method = dist.method))
  if (is.null(rownames(D))) rownames(D) <- as.character(seq_len(n))
  if (is.null(colnames(D))) colnames(D) <- rownames(D)
  ids <- rownames(D)

  if (!isTRUE(joint)) {
    ## Classic: MyOrd <- order(My.distances[, findKNN])[2:cell.num]
    ##          MyDist <- My.distances[MyOrd]   # linear index quirk
    from <- integer(n * (zoom - 1L))
    to <- integer(n * (zoom - 1L))
    weight <- numeric(n * (zoom - 1L))
    pos <- 0L
    for (i in seq_len(n)) {
      MyOrd <- order(D[, i])[2:zoom]
      MyDist <- D[MyOrd]
      idx <- (pos + 1L):(pos + zoom - 1L)
      from[idx] <- i
      to[idx] <- as.integer(MyOrd)
      weight[idx] <- as.numeric(MyDist)
      pos <- pos + zoom - 1L
    }
    return(list(from = from, to = to, weight = weight))
  }

  if (is.null(conditions) || length(conditions) != n) {
    stop("'joint = TRUE' requires a condition label per cell.")
  }
  conditions <- as.character(conditions)
  lev <- unique(conditions)
  ## Pre-size upper bound: n * zoom * n_levels; trim after.
  n_lev <- length(lev)
  cap <- as.integer(n) * as.integer(zoom) * as.integer(n_lev)
  from <- integer(cap)
  to <- integer(cap)
  weight <- numeric(cap)
  pos <- 0L
  for (i in seq_len(n)) {
    CellOrd <- ids[order(D[, i])]
    for (j in seq_along(lev)) {
      hit <- CellOrd[conditions[match(CellOrd, ids)] == lev[[j]]]
      if (length(hit) < 1L) next
      take <- hit[seq_len(min(zoom, length(hit)))]
      MyOrd <- match(take, ids)
      m <- length(MyOrd)
      idx <- (pos + 1L):(pos + m)
      from[idx] <- i
      to[idx] <- as.integer(MyOrd)
      weight[idx] <- as.numeric(D[MyOrd])
      pos <- pos + m
    }
  }
  if (pos < 1L) {
    return(list(from = integer(0), to = integer(0), weight = numeric(0)))
  }
  list(
    from = from[seq_len(pos)],
    to = to[seq_len(pos)],
    weight = weight[seq_len(pos)]
  )
}


#' Drop self-neighbor column from knn matrices (vectorized).
#' @keywords internal
.ic2_knetl_drop_self <- function(idx, distm, k) {
  n <- nrow(idx)
  k <- as.integer(k)
  ## Prefer drop first column when it is self (sorted FNN typical).
  if (ncol(idx) >= k + 1L) {
    self_first <- sum(idx[, 1L] == seq_len(n))
    if (self_first >= max(1L, as.integer(0.9 * n))) {
      return(list(
        idx = idx[, 2:(k + 1L), drop = FALSE],
        dist = distm[, 2:(k + 1L), drop = FALSE]
      ))
    }
  }
  ## Fallback: per-row drop (rare ties / missing self).
  out_i <- matrix(0L, n, k)
  out_d <- matrix(0, n, k)
  kk <- ncol(idx)
  for (r in seq_len(n)) {
    hit <- idx[r, ]
    dhit <- distm[r, ]
    self <- which(hit == r | dhit <= 0)[1L]
    if (is.na(self)) {
      keep <- seq_len(k)
    } else {
      keep <- setdiff(seq_len(kk), self)
      if (length(keep) < k) keep <- c(keep, setdiff(seq_len(kk), keep))
      keep <- keep[seq_len(k)]
    }
    out_i[r, ] <- as.integer(hit[keep])
    out_d[r, ] <- as.numeric(dhit[keep])
  }
  list(idx = out_i, dist = out_d)
}


#' Chunked exact FNN KNetL edges; optional SODA knn sidecar; parallel chunks.
#' @keywords internal
.ic2_knetl_neighbors_fnn <- function(
    X,
    zoom,
    joint = FALSE,
    conditions = NULL,
    chunk.size = 25000L,
    soda_graph_dir = NULL,
    cell_ids = NULL,
    n_threads = 4L,
    verbose = FALSE
) {
  n <- nrow(X)
  k <- as.integer(zoom - 1L)
  chunk.size <- as.integer(max(1L, chunk.size))
  n_threads <- as.integer(max(1L, n_threads))
  if (is.null(cell_ids) || length(cell_ids) != n) {
    cell_ids <- rownames(X)
    if (is.null(cell_ids)) cell_ids <- as.character(seq_len(n))
  }

  if (!isTRUE(joint)) {
    if (n <= chunk.size) {
      kn <- FNN::get.knn(data = X, k = k + 1L)
      dropped <- .ic2_knetl_drop_self(kn$nn.index, kn$nn.dist, k)
      idx <- dropped$idx
      distm <- dropped$dist
    } else {
      starts <- seq.int(1L, n, by = chunk.size)
      ends <- pmin(starts + chunk.size - 1L, n)
      if (isTRUE(verbose)) {
        message(
          "run.knetl: exact FNN chunks = ", length(starts),
          " n_threads = ", n_threads, " k = ", k
        )
      }
      chunks <- .ic2_knetl_fnn_query_chunks(
        X, starts = starts, ends = ends, kk = k + 1L, n_threads = n_threads
      )
      idx_raw <- matrix(0L, n, k + 1L)
      dist_raw <- matrix(0, n, k + 1L)
      for (i in seq_along(starts)) {
        s <- starts[[i]]
        e <- ends[[i]]
        idx_raw[s:e, ] <- chunks[[i]]$nn.index
        dist_raw[s:e, ] <- chunks[[i]]$nn.dist
      }
      dropped <- .ic2_knetl_drop_self(idx_raw, dist_raw, k)
      idx <- dropped$idx
      distm <- dropped$dist
    }
    .ic2_knetl_write_knn_soda(
      soda_graph_dir, idx, distm, cell_ids, zoom, joint = FALSE,
      method = "fnn_chunked"
    )
    from <- rep(seq_len(n), each = k)
    to <- as.integer(as.vector(t(idx)))
    weight <- as.numeric(as.vector(t(distm)))
    return(list(from = from, to = to, weight = weight, idx = idx, dist = distm))
  }

  if (is.null(conditions) || length(conditions) != n) {
    stop("'joint = TRUE' requires a condition label per cell.")
  }
  conditions <- as.character(conditions)
  lev <- unique(conditions)
  ## Preallocate: each cell contributes zoom neighbors per condition.
  n_lev <- length(lev)
  cap <- as.integer(n) * as.integer(zoom) * as.integer(n_lev)
  frm <- integer(cap)
  too <- integer(cap)
  ww <- numeric(cap)
  pos <- 0L
  for (j in seq_along(lev)) {
    sel <- which(conditions == lev[[j]])
    kk <- min(as.integer(zoom), length(sel))
    if (kk < 1L) next
    data <- X[sel, , drop = FALSE]
    starts <- seq.int(1L, n, by = chunk.size)
    ends <- pmin(starts + chunk.size - 1L, n)
    for (ii in seq_along(starts)) {
      s <- starts[[ii]]
      e <- ends[[ii]]
      kn <- FNN::get.knnx(data = data, query = X[s:e, , drop = FALSE], k = kk)
      m <- (e - s + 1L) * kk
      idx <- (pos + 1L):(pos + m)
      frm[idx] <- rep(s:e, each = kk)
      too[idx] <- as.integer(sel[as.vector(t(kn$nn.index))])
      ww[idx] <- as.numeric(as.vector(t(kn$nn.dist)))
      pos <- pos + m
    }
  }
  if (pos < 1L) {
    return(list(from = integer(0), to = integer(0), weight = numeric(0)))
  }
  list(
    from = frm[seq_len(pos)],
    to = too[seq_len(pos)],
    weight = ww[seq_len(pos)]
  )
}


#' Parallel (or serial) FNN get.knnx over query row ranges.
#' Exact same calls as serial chunked FNN; only scheduling differs.
#' @keywords internal
.ic2_knetl_fnn_query_chunks <- function(X, starts, ends, kk, n_threads = 1L) {
  n_chunks <- length(starts)
  n_threads <- as.integer(max(1L, n_threads))
  run_one <- function(i) {
    s <- starts[[i]]
    e <- ends[[i]]
    FNN::get.knnx(data = X, query = X[s:e, , drop = FALSE], k = as.integer(kk))
  }
  can_fork <- .Platform$OS.type != "windows" && n_threads > 1L && n_chunks > 1L
  if (can_fork && requireNamespace("parallel", quietly = TRUE)) {
    n_workers <- min(n_threads, n_chunks)
    return(parallel::mclapply(
      seq_len(n_chunks), run_one, mc.cores = n_workers
    ))
  }
  lapply(seq_len(n_chunks), run_one)
}


#' Reuse SODA knn sidecar when n/k/joint match.
#' @keywords internal
.ic2_knetl_try_reuse_knn <- function(
    soda_graph_dir, n, k, joint = FALSE, need_edges = TRUE
) {
  if (is.null(soda_graph_dir) || !nzchar(as.character(soda_graph_dir)[1L])) {
    return(NULL)
  }
  meta_path <- file.path(soda_graph_dir, "knetl.json")
  knn_path <- file.path(soda_graph_dir, "knn_index.h5")
  dist_path <- file.path(soda_graph_dir, "knn_dist.h5")
  if (!file.exists(meta_path) || !file.exists(knn_path) || !file.exists(dist_path)) {
    return(NULL)
  }
  meta <- tryCatch(
    jsonlite::fromJSON(meta_path, simplifyVector = TRUE, simplifyDataFrame = FALSE),
    error = function(e) NULL
  )
  if (is.null(meta)) return(NULL)
  if (isTRUE(as.logical(meta$joint)) != isTRUE(joint)) return(NULL)
  if (as.integer(meta$n_cells)[1L] != as.integer(n)) return(NULL)
  meta_k <- as.integer(meta$k)[1L]
  if (!is.finite(meta_k) || meta_k < as.integer(k)) return(NULL)
  ## Only reuse exact FNN / full_dist sidecars -- never ANN/HNSW.
  meth <- as.character(meta$method)[1L]
  if (!is.na(meth) && grepl("hnsw|ann|pynn", meth, ignore.case = TRUE)) {
    return(NULL)
  }
  idx <- tryCatch(.ic2_h5_read_dense(knn_path), error = function(e) NULL)
  distm <- tryCatch(.ic2_h5_read_dense(dist_path), error = function(e) NULL)
  if (is.null(idx) || is.null(distm)) return(NULL)
  if (nrow(idx) != n || ncol(idx) < k) return(NULL)
  idx <- idx[, seq_len(k), drop = FALSE]
  distm <- distm[, seq_len(k), drop = FALSE]
  storage.mode(idx) <- "integer"
  if (isTRUE(need_edges)) {
    from <- rep(seq_len(n), each = k)
    to <- as.integer(as.vector(t(idx)))
    weight <- as.numeric(as.vector(t(distm)))
  } else {
    from <- integer(0)
    to <- integer(0)
    weight <- numeric(0)
  }
  list(from = from, to = to, weight = weight, idx = idx, dist = distm)
}


#' Persist knn idx/dist under SODA graphs/knetl/ (audit / lower peak reuse).
#' @keywords internal
.ic2_knetl_write_knn_soda <- function(
    soda_graph_dir, idx, distm, cell_ids, zoom, joint = FALSE,
    method = "fnn_chunked"
) {
  if (is.null(soda_graph_dir) || !nzchar(as.character(soda_graph_dir)[1L])) {
    return(invisible(NULL))
  }
  dir.create(soda_graph_dir, recursive = TRUE, showWarnings = FALSE)
  knn_path <- file.path(soda_graph_dir, "knn_index.h5")
  dist_path <- file.path(soda_graph_dir, "knn_dist.h5")
  tryCatch(
    {
      .ic2_h5_write_dense(idx, knn_path, group = "data")
      .ic2_h5_write_dense(distm, dist_path, group = "data")
      utils::write.table(
        data.frame(cell_id = cell_ids, stringsAsFactors = FALSE),
        file = file.path(soda_graph_dir, "cell_ids.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
      )
      .ic2_soda_write_json(
        list(
          method = as.character(method)[1L],
          zoom = as.integer(zoom),
          k = as.integer(ncol(idx)),
          n_cells = as.integer(nrow(idx)),
          joint = isTRUE(joint),
          knn_path = basename(knn_path),
          dist_path = basename(dist_path)
        ),
        file.path(soda_graph_dir, "knetl.json")
      )
    },
    error = function(e) invisible(NULL)
  )
  invisible(NULL)
}


#' Integer-vertex igraph (same topology as classic named graph_from_data_frame).
#' @keywords internal
.ic2_knetl_graph <- function(from, to, weight, cell_ids, max.graph.edges) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for run.knetl().")
  }
  n_edges <- length(from)
  n_vert <- length(cell_ids)
  if (is.finite(max.graph.edges) && n_edges > max.graph.edges) {
    stop(
      "KNetL graph has ", n_edges, " edges (cap ", max.graph.edges,
      "). Lower zoom or raise max.graph.edges (default Inf matches classic)."
    )
  }
  if (n_edges > 2e6) {
    warning(
      "KNetL graph has ", format(n_edges, scientific = FALSE),
      " edges; classic layout_nicely at this size is typically day-scale.",
      call. = FALSE
    )
  }
  ok <- is.finite(to) & to >= 1L & to <= n_vert &
    is.finite(from) & from >= 1L & from <= n_vert
  from <- as.integer(from[ok])
  to <- as.integer(to[ok])
  weight <- as.numeric(weight[ok])
  ## Integer edgelist avoids materializing 2xn_edges character cell IDs.
  el <- cbind(from, to)
  g <- igraph::graph_from_edgelist(el, directed = FALSE)
  if (igraph::vcount(g) < n_vert) {
    g <- igraph::add_vertices(g, n_vert - igraph::vcount(g))
  }
  igraph::V(g)$name <- cell_ids
  igraph::E(g)$weight <- weight
  g
}


#' @keywords internal
.ic2_knetl_color_graph <- function(g, x, cell_ids, col.by) {
  if (!requireNamespace("igraph", quietly = TRUE)) return(g)
  if (identical(col.by, "clusters") && "cluster_id" %in% colnames(x@cell_data)) {
    labs <- x@cell_data[cell_ids, "cluster_id"]
  } else if ("condition_id" %in% colnames(x@cell_data)) {
    labs <- x@cell_data[cell_ids, "condition_id"]
  } else {
    return(g)
  }
  labs <- as.integer(factor(labs))
  n <- max(labs, na.rm = TRUE)
  hues <- seq(15, 375, length.out = n + 1L)
  cols <- grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
  igraph::V(g)$color <- cols[labs]
  g
}


#' @keywords internal
.ic2_knetl_igraph_layout <- function(g, layout_name, dim = 2L) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for run.knetl().")
  }
  fn <- tryCatch(
    get(layout_name, mode = "function", envir = asNamespace("igraph")),
    error = function(e) NULL
  )
  if (is.null(fn)) {
    fn <- if (as.integer(dim) >= 3L) igraph::layout_with_fr else igraph::layout_nicely
  }
  ## layout_nicely ignores dim; FR / others accept dim when requested.
  lay <- tryCatch(
    fn(g, dim = as.integer(dim)),
    error = function(e) fn(g)
  )
  as.matrix(lay)
}


#' Classic unpack: UMAP/PCA on the scaled 2D igraph layout.
#' UMAP uses classic run.umap defaults (n_neighbors = 15, seed 0).
#' @keywords internal
.ic2_knetl_unpack <- function(
    lay2, dim.redux, n_components, min_dist, spread, n_threads
) {
  if (identical(dim.redux, "pca")) {
    return(stats::prcomp(lay2, center = TRUE, scale. = FALSE)$x[, seq_len(min(as.integer(n_components), ncol(lay2))), drop = FALSE])
  }
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for KNetL UMAP unpack (do.redux = TRUE).")
  }
  ## Classic: y <- run.umap(y, dims = 1:2) with run.umap my.seed = 0 default.
  set.seed(0L)
  n_neighbors <- min(15L, nrow(lay2) - 1L)
  if (n_neighbors < 2L) stop("Need at least 3 cells for UMAP unpack.")
  uwot::umap(
    X = lay2,
    n_neighbors = n_neighbors,
    n_components = as.integer(n_components),
    metric = "euclidean",
    min_dist = min_dist,
    spread = spread,
    scale = FALSE,
    n_threads = n_threads,
    n_sgd_threads = 0L,
    verbose = isTRUE(getOption("iCellR2.uwot.verbose", FALSE))
  )
}


#' @keywords internal
.ic2_knetl_store <- function(
    x, modality, lay, name, src, zoom, k, dims, joint,
    neighbor_method = "full_dist",
    pipeline = "classic_layout_then_redux",
    zoom_requested = NULL
) {
  root <- x@storage$root_dir
  red_dir <- file.path(root, "reductions", name)
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  meta_path <- file.path(red_dir, "reduction.json")
  .ic2_h5_write_dense(lay, emb_path)
  meta <- list(
    name = name,
    mark = src$mark,
    aligned = isTRUE(src$aligned),
    source = src$name,
    based_on = src$name,
    zoom = as.integer(zoom),
    zoom_requested = if (is.null(zoom_requested)) as.integer(zoom) else as.integer(zoom_requested),
    k = as.integer(k),
    dims = as.integer(dims),
    joint = isTRUE(joint),
    n_cells = nrow(lay),
    n_components = ncol(lay),
    pipeline = as.character(pipeline)[1L],
    neighbor_method = as.character(neighbor_method)[1L],
    created = as.character(Sys.time()),
    paths = list(embeddings = emb_path)
  )
  .ic2_soda_write_json(meta, meta_path)
  x@modalities[[modality]]$reductions[[name]] <- list(
    cell_embeddings = emb_path,
    meta_path = meta_path,
    aligned = isTRUE(src$aligned),
    mark = src$mark,
    source = src$name,
    based_on = src$name,
    n_components = ncol(lay)
  )
  x@storage$reduction_paths[[name]] <- red_dir
  x@modalities[[modality]]$active_reduction <- name
  x
}
