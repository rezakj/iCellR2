#' Cluster-guided UMAP (cgUMAP)
#'
#' PCA dims + \code{cluster_id}  ->  UMAP layout. Stores as \code{cgumap}.
#' Plot with \code{cluster.plot(x, plot.type = "cgumap")}.
#'
#' @inheritParams .ic2_run_cgmap
#' @param dims PCA dimensions. Default \code{1:3}.
#' @param ... Passed to \code{\link{.ic2_run_cgmap}} (uwot knobs, sketch, etc.).
#'
#' @return Updated \code{iCellR2} with a \code{cgumap} reduction.
#'
#' @examples
#' \dontrun{
#' obj <- run.cgumap(obj)
#' cluster.plot(obj, plot.type = "cgumap")
#' }
#'
#' @export
run.cgumap <- function(
    x,
    modality = "rna",
    dims = 1:3,
    scale = TRUE,
    cluster.weight = 1,
    cluster.col = "cluster_id",
    center = TRUE,
    verbose = TRUE,
    ...
) {
  .ic2_run_cgmap(
    x,
    modality = modality,
    input = "pca",
    dims = dims,
    method = "umap",
    scale = scale,
    cluster.weight = cluster.weight,
    cluster.col = cluster.col,
    center = center,
    reduction.name = "cgumap",
    verbose = verbose,
    ...
  )
}
