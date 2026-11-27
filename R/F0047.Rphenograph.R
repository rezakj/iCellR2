#' PhenoGraph clustering (kNN + Jaccard + Louvain)
#'
#' R implementation of PhenoGraph for embedding matrices (cells  x  dims).
#' Uses \pkg{FNN} for nearest neighbors when available, the package
#' \code{jaccard_coeff} C++ kernel, and \pkg{igraph} Louvain.
#'
#' Prefer \code{\link{run.clustering}} for iCellR2 objects (progress, SODA,
#' sketch -> assign). This function is the low-level engine.
#'
#' @param data Numeric matrix; rows = cells, columns = dimensions.
#' @param k Integer; number of nearest neighbours (default 30).
#' @param verbose If \code{TRUE}, print stage messages.
#'
#' @return A list with \code{graph} (igraph) and \code{community}
#'   (Louvain communities object).
#'
#' @source \url{https://github.com/JinmiaoChenLab/Rphenograph}
#' @references Levine et al. Cell 2015.
#' @importFrom igraph graph_from_data_frame cluster_louvain modularity membership
#' @export
Rphenograph <- function(data, k = 30, verbose = FALSE) {
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data)) stop("'data' must be a matrix (cells x dims).")
  storage.mode(data) <- "double"
  k <- as.integer(k)
  if (k < 1L) stop("'k' must be a positive integer.")
  if (k > nrow(data) - 2L) {
    stop("'k' must be smaller than nrow(data) - 1.")
  }

  if (isTRUE(verbose)) {
    message(
      "Rphenograph: ", nrow(data), " cells x ", ncol(data),
      " dims   |   k=", k
    )
  }

  if (isTRUE(verbose)) message("  Finding nearest neighbors...")
  neighborMatrix <- .ic2_find_neighbors(data, k = k + 1L)[, -1L, drop = FALSE]

  if (isTRUE(verbose)) message("  Jaccard coefficients...")
  links <- jaccard_coeff(neighborMatrix)
  links <- links[links[, 1] > 0, , drop = FALSE]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  if (isTRUE(verbose)) message("  Building graph + Louvain...")
  g <- igraph::graph_from_data_frame(relations, directed = FALSE)
  community <- igraph::cluster_louvain(g)

  if (isTRUE(verbose)) {
    message(
      "  Clusters: ", length(unique(igraph::membership(community))),
      "   |   modularity: ", round(igraph::modularity(community), 4)
    )
  }
  list(graph = g, community = community)
}

#' kNN index matrix (self in column 1), matching classic PhenoGraph layout
#' @param data cells x dims matrix
#' @param k number of columns including self
#' @keywords internal
.ic2_find_neighbors <- function(data, k) {
  data <- as.matrix(data)
  n <- nrow(data)
  k <- as.integer(k)
  if (k < 1L || k > n) {
    stop("Invalid k for neighbor search: ", k, " (n=", n, ").")
  }
  out <- matrix(0L, nrow = n, ncol = k)
  out[, 1L] <- seq_len(n)
  if (k == 1L) return(out)

  knn_k <- k - 1L
  if (requireNamespace("FNN", quietly = TRUE)) {
    kn <- FNN::get.knn(data, k = knn_k)
    out[, -1L] <- kn$nn.index
    return(out)
  }

  ## Fallback exact kNN (prefer FNN for large n)
  if (n > 8000L) {
    warning(
      "Package 'FNN' not installed; using slower exact kNN. ",
      "Install FNN for scalable clustering."
    )
  }
  for (i in seq_len(n)) {
    di <- sqrt(colSums((t(data) - data[i, ])^2))
    di[i] <- Inf
    out[i, -1L] <- order(di)[seq_len(knn_k)]
  }
  out
}

#' Build a distance-weighted kNN graph and Louvain-cluster (classic iclust)
#' @keywords internal
.ic2_louvain_knn <- function(data, k = 30L, dist.method = "euclidean") {
  data <- as.matrix(data)
  n <- nrow(data)
  k <- as.integer(min(as.integer(k), n - 1L))
  if (k < 2L) stop("Need k >= 2 for Louvain kNN graph.")

  if (requireNamespace("FNN", quietly = TRUE) && identical(dist.method, "euclidean")) {
    kn <- FNN::get.knn(data, k = k)
    from <- rep(seq_len(n), each = k)
    to <- as.integer(t(kn$nn.index))
    w <- as.numeric(t(kn$nn.dist))
    w <- 1 / (w + 1e-6)
  } else {
    D <- as.matrix(stats::dist(data, method = dist.method))
    from <- integer(n * k)
    to <- integer(n * k)
    w <- numeric(n * k)
    pos <- 1L
    for (i in seq_len(n)) {
      ord <- order(D[i, ])[2:(k + 1L)]
      from[pos:(pos + k - 1L)] <- i
      to[pos:(pos + k - 1L)] <- ord
      w[pos:(pos + k - 1L)] <- 1 / (D[i, ord] + 1e-6)
      pos <- pos + k
    }
  }
  relations <- data.frame(from = from, to = to, weight = w)
  g <- igraph::graph_from_data_frame(relations, directed = FALSE)
  community <- igraph::cluster_louvain(g)
  list(graph = g, community = community)
}

#' Assign labels from sketch cells to all cells via kNN in embedding space
#' @keywords internal
.ic2_assign_clusters_knn <- function(emb, fit_ids, fit_labels, k = 5L) {
  emb <- as.matrix(emb)
  all_ids <- rownames(emb)
  if (is.null(all_ids)) stop("Embedding must have rownames (cell ids).")
  fit_labels <- as.character(fit_labels)
  names(fit_labels) <- fit_ids
  out <- character(length(all_ids))
  names(out) <- all_ids
  out[fit_ids] <- fit_labels[fit_ids]

  query_ids <- setdiff(all_ids, fit_ids)
  if (length(query_ids) == 0L) return(out)

  ref <- emb[fit_ids, , drop = FALSE]
  qry <- emb[query_ids, , drop = FALSE]
  k <- as.integer(min(as.integer(k), nrow(ref)))

  if (requireNamespace("FNN", quietly = TRUE)) {
    kn <- FNN::get.knnx(ref, qry, k = k)
    idx <- kn$nn.index
  } else {
    ## Exact assign: nearest sketch cell
    idx <- matrix(0L, nrow = length(query_ids), ncol = k)
    for (i in seq_along(query_ids)) {
      di <- sqrt(colSums((t(ref) - qry[i, ])^2))
      idx[i, ] <- order(di)[seq_len(k)]
    }
  }

  for (i in seq_along(query_ids)) {
    labs <- fit_labels[fit_ids[idx[i, ]]]
    tab <- sort(table(labs), decreasing = TRUE)
    out[query_ids[i]] <- names(tab)[1L]
  }
  out
}
