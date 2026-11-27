#' Cell-cycle phase, or extra gene-set scores
#'
#' With no \code{scoring.List}, this is \code{\link{cc}} (S, G2M, and Phase).
#' A named list is scored with \code{\link{i.score}}. If that list contains
#' \code{S} and \code{G2M}, \code{Phase} is set from those two scores.
#' Only the genes in the sets are read.
#'
#' @param object An iCellR2 object. Alias \code{x}.
#' @param x Alias for \code{object}.
#' @param scoring.List Named list of gene sets. \code{NULL} runs S and G2M.
#' @param return.stats If \code{TRUE}, return the score columns. Default
#'   \code{FALSE} returns the object.
#' @param scoring.method Passed to \code{i.score}. Default \code{"tirosh"}.
#' @param assay.name Assay. Default \code{"data"}.
#' @param modality Modality. Default \code{"rna"}.
#' @return The object, or a data.frame when \code{return.stats = TRUE}.
#' @export
cell.cycle <- function(
    object = NULL,
    x = NULL,
    scoring.List = NULL,
    return.stats = FALSE,
    scoring.method = "tirosh",
    assay.name = "data",
    modality = "rna"
) {
  if (is.null(object)) object <- x
  if (!inherits(object, "iCellR2")) stop("object should be an object of class iCellR2")
  if (is.null(scoring.List)) {
    object <- cc(object, assay.name = assay.name, modality = modality)
  } else {
    object <- i.score(
      object,
      scoring.List = scoring.List,
      scoring.method = scoring.method,
      assay.name = assay.name,
      modality = modality
    )
    cd <- object@cell_data
    if (all(c("score_S", "score_G2M") %in% colnames(cd))) {
      object@cell_data$Phase <- .ic2_cycle_phase(cd$score_S, cd$score_G2M)
      object <- .ic2_persist_cell_data(object)
    }
  }
  if (isTRUE(return.stats)) {
    cd <- object@cell_data
    cols <- grep("^(score_|Phase$)", colnames(cd), value = TRUE)
    return(cd[, cols, drop = FALSE])
  }
  object
}
