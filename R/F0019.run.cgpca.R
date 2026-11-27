#' Cluster-guided PCA (cgPCA)
#'
#' PCA dims + \code{cluster_id}  ->  PCA. Stores as \code{cgpca}.
#' Plot with \code{cluster.plot(x, plot.type = "cgpca")}.
#'
#' @inheritParams .ic2_run_cgmap
#' @param dims PCA dimensions. Default \code{1:3}.
#' @param ... Passed to \code{\link{.ic2_run_cgmap}}.
#'
#' @return Updated \code{iCellR2} with a \code{cgpca} reduction.
#'
#' @examples
#' \dontrun{
#' obj <- run.cgpca(obj)
#' cluster.plot(obj, plot.type = "cgpca")
#' }
#'
#' @export
run.cgpca <- function(
    x,
    modality = "rna",
    dims = 1:3,
    npcs = 3L,
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
    method = "pca",
    npcs = npcs,
    scale = scale,
    cluster.weight = cluster.weight,
    cluster.col = cluster.col,
    center = center,
    reduction.name = "cgpca",
    verbose = verbose,
    ...
  )
}
