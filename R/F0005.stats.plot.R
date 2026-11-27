#' QC statistics dashboard (HTML / PDF / PNG / SVG)
#'
#' Replaces the older single-panel \code{stats.plot} helpers with the unified
#' QC dashboard: an interactive HTML filter explorer (default) and optional
#' static multi-panel exports matching the publication figure.
#'
#' Run \code{\link{qc.stats}} first. This is a thin wrapper around
#' \code{\link{qc.dashboard}} so existing pipelines can keep calling
#' \code{stats.plot()}.
#'
#' @inheritParams qc.dashboard
#' @param interactive Logical. If \code{TRUE} (default) and \code{format} is
#'   left at its default, write HTML. If \code{FALSE} and \code{format} is
#'   still the default, write PDF instead. Ignored when \code{format} is set
#'   explicitly to something other than \code{"html"}.
#' @param plot.type,cell.color,cell.size,jitter.width,cell.transparency,box.color,box.line.col,back.col,scatter.mode,downsample,save.plot
#'   Ignored (kept for backward compatibility with older \code{stats.plot} calls).
#'
#' @return Invisibly, a named character vector of files written (from
#'   \code{\link{qc.dashboard}}).
#'
#' @examples
#' \dontrun{
#' iCellR2.obj <- qc.stats(iCellR2.obj)
#' stats.plot(iCellR2.obj, out.name = "plot1_QC")                 # HTML
#' stats.plot(iCellR2.obj, out.name = "plot1_QC", format = "pdf")
#' stats.plot(iCellR2.obj, out.name = "plot1_QC", format = c("html", "png"))
#' }
#'
#' @export
stats.plot <- function(
    x,
    out.name = "stats_plot",
    format = "html",
    min.genes = 200,
    min.umis = 500,
    max.mito = 0.1,
    group.by = "condition_id",
    max.cells = 500000L,
    width = 12,
    height = 9,
    dpi = 300,
    seed = 1L,
    verbose = TRUE,
    interactive = TRUE,
    plot.type = NULL,
    cell.color = NULL,
    cell.size = NULL,
    jitter.width = NULL,
    cell.transparency = NULL,
    box.color = NULL,
    box.line.col = NULL,
    back.col = NULL,
    scatter.mode = NULL,
    downsample = NULL,
    save.plot = NULL
) {
  ignored <- c(
    "plot.type", "cell.color", "cell.size", "jitter.width",
    "cell.transparency", "box.color", "box.line.col", "back.col",
    "scatter.mode", "downsample", "save.plot"
  )
  provided <- vapply(
    ignored,
    function(nm) !is.null(eval(as.name(nm), envir = environment())),
    logical(1L)
  )
  if (any(provided)) {
    warning(
      "stats.plot() now writes the QC dashboard; ignoring legacy argument(s): ",
      paste(ignored[provided], collapse = ", "),
      call. = FALSE
    )
  }

  format <- unique(tolower(as.character(format)))
  # Legacy: interactive=FALSE meant a static device plot  ->  PDF snapshot.
  if (isFALSE(interactive) && identical(format, "html")) {
    format <- "pdf"
  }

  qc.dashboard(
    x = x,
    out.name = out.name,
    format = format,
    min.genes = min.genes,
    min.umis = min.umis,
    max.mito = max.mito,
    group.by = group.by,
    max.cells = max.cells,
    width = width,
    height = height,
    dpi = dpi,
    seed = seed,
    verbose = verbose
  )
}
