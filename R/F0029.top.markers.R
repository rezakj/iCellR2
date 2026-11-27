#' Choose top marker genes from a markers table
#'
#' Table-only helper (no assay I/O). Input is typically the markers /
#' DE data frame from \code{\link{findMarkers}} (columns include
#' \code{gene}, \code{clusters}, \code{baseMean}).
#'
#' @param x A markers / DE \code{data.frame} (not an iCellR2 object).
#' @param topde Number of top genes per cluster. Default \code{10}.
#' @param min.base.mean Minimum \code{baseMean}. Default \code{0.2}.
#' @param filt.ambig If \code{TRUE}, drop genes listed for more than one
#'   cluster. Default \code{TRUE}.
#' @param cluster Cluster id to restrict to, or \code{0} (default) for all.
#'
#' @return Character vector of gene ids.
#'
#' @export
top.markers <- function(
    x,
    topde = 10,
    min.base.mean = 0.2,
    filt.ambig = TRUE,
    cluster = 0
) {
  if (is.null(x) || !is.data.frame(x)) {
    stop("'x' must be a markers/DE data.frame (with columns gene, clusters, ...).")
  }
  need <- c("gene", "clusters")
  if (!all(need %in% colnames(x))) {
    stop("Markers table must contain columns: ", paste(need, collapse = ", "))
  }
  topde <- as.integer(topde)[1L]
  if (!is.finite(topde) || topde < 1L) stop("'topde' must be >= 1.")

  df <- x
  df$gene <- as.character(df$gene)
  df$clusters <- as.character(df$clusters)

  if (isTRUE(filt.ambig)) {
    tab <- table(df$gene)
    keep_genes <- names(tab)[tab == 1L]
    df <- df[df$gene %in% keep_genes, , drop = FALSE]
  }

  if (!identical(as.character(cluster)[1L], "0") &&
      !(is.numeric(cluster) && length(cluster) == 1L && isTRUE(cluster == 0))) {
    cl <- as.character(cluster)
    df <- df[df$clusters %in% cl, , drop = FALSE]
    if ("baseMean" %in% colnames(df)) {
      df <- df[as.numeric(df$baseMean) >= as.numeric(min.base.mean), , drop = FALSE]
    }
    return(utils::head(unique(df$gene), topde))
  }

  ## All clusters: take topde per cluster (table assumed pre-sorted by significance)
  clusts <- unique(df$clusters)
  ## Numeric-aware order when possible
  num <- suppressWarnings(as.numeric(clusts))
  if (all(is.finite(num))) clusts <- clusts[order(num)] else clusts <- sort(clusts)

  out <- character(0)
  for (cl in clusts) {
    sub <- df[df$clusters == cl, , drop = FALSE]
    if ("baseMean" %in% colnames(sub)) {
      sub <- sub[as.numeric(sub$baseMean) >= as.numeric(min.base.mean), , drop = FALSE]
    }
    out <- c(out, utils::head(unique(as.character(sub$gene)), topde))
  }
  unique(out)
}
