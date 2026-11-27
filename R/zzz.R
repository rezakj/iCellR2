#' Package attach / load hooks
#' @keywords internal
#' @noRd
NULL

utils::globalVariables(c(
  ".data", "LogSDs", "PCs", "SDs", "assignment.annotation", "bin", "chain",
  "cluster_id", "d1", "d2", "low.cov", "n_cells", "percent.match",
  "spatial_x", "spatial_y", "y"
))

.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion(pkgname)
  width <- getOption("width", 80L)
  if (!is.finite(width) || width < 40L) width <- 80L
  .center <- function(txt) {
    txt <- as.character(txt)[1L]
    pad <- max(0L, as.integer((width - nchar(txt, type = "width")) / 2))
    paste0(strrep(" ", pad), txt)
  }
  packageStartupMessage(
    "\n",
    .center(paste0("iCellR2 v", ver)), "\n\n",
    .center("large-scale single-cell analysis with a familiar R interface."), "\n\n",
    .center("Assays stay on disk in SODA  -  built for millions of cells."), "\n\n",
    .center("Successor to iCellR (since 2018) for large-scale workflows."), "\n\n",
    .center(" | "),
    "\n"
  )
}
