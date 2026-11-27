#' Pseudotime-guided UMAP (ptUMAP)
#'
#' PCA dims (default 10) + \code{Pseudotime} + \code{branch_id}  ->  write
#' \code{guided_input.tsv} (12 columns by default)  ->  UMAP. Stores as \code{ptumap}.
#' Plot with \code{cluster.plot(x, plot.type = "ptumap")}.
#'
#' @param x An \code{iCellR2} object after \code{\link{run.pseudotime}}.
#' @param modality Modality. Default \code{"rna"}.
#' @param dims PCA dimensions. Default \code{1:10}.
#' @param scale Scale the pseudotime and branch columns before UMAP.
#'   Default \code{TRUE}.
#' @param pt.weight Weight of the pseudotime column. Default \code{1}.
#' @param branch.weight Weight of the branch column. Default \code{1}.
#' @param pseudotime.col Pseudotime column in \code{cell_data}.
#'   Default \code{"Pseudotime"}.
#' @param branch.col Branch column to append. Default \code{"branch_id"}.
#' @param center Center columns before UMAP. Default \code{TRUE}.
#' @param verbose Progress style.
#' @param ... Passed to \code{\link{.ic2_run_ptmap}} (uwot knobs, sketch, etc.).
#'
#' @return Updated \code{iCellR2} with a \code{ptumap} reduction.
#'
#' @examples
#' \dontrun{
#' obj <- run.pca(obj)
#' obj <- run.pseudotime(obj, root.cluster = 1)
#' obj <- run.ptumap(obj)
#' cluster.plot(obj, plot.type = "ptumap", col.by = "branch_id")
#' cluster.plot(obj, plot.type = "ptumap", col.by = "Pseudotime")
#' }
#'
#' @export
run.ptumap <- function(
    x,
    modality = "rna",
    dims = 1:10,
    scale = TRUE,
    pt.weight = 1,
    branch.weight = 1,
    pseudotime.col = "Pseudotime",
    branch.col = "branch_id",
    center = TRUE,
    verbose = TRUE,
    ...
) {
  .ic2_run_ptmap(
    x,
    modality = modality,
    input = "pca",
    dims = dims,
    method = "umap",
    scale = scale,
    pt.weight = pt.weight,
    branch.weight = branch.weight,
    pseudotime.col = pseudotime.col,
    branch.col = branch.col,
    center = center,
    reduction.name = "ptumap",
    verbose = verbose,
    ...
  )
}
