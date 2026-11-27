
#' Plot one Xenium gene from transcripts.parquet
#'
#' Streams a single gene with \pkg{arrow}. The whole transcript table is not
#' loaded. \code{xlim} and \code{ylim} are microns. At most \code{max.points}
#' transcripts are drawn.
#'
#' @param x An iCellR2 object after \code{\link{add.xenium}}, or a path to
#'   \code{transcripts.parquet}.
#' @param gene Gene name (\code{feature_name}).
#' @param xlim,ylim Optional crop in microns.
#' @param max.points Maximum transcripts drawn. Default \code{200000}.
#' @param plot Draw the points. \code{FALSE} returns the sampled table.
#' @param out.name Optional PNG path.
#' @return A ggplot, or a data.frame when \code{plot = FALSE}.
#' @examples
#' \dontrun{
#' spatial.transcripts(obj, gene = "CD3E", max.points = 200000, out.name = "cd3e_tx.png")
#' }
#' @export
spatial.transcripts <- function(
    x = NULL,
    gene = NULL,
    xlim = NULL,
    ylim = NULL,
    max.points = 200000L,
    plot = TRUE,
    out.name = NULL
) {
  if (is.null(gene) || !nzchar(as.character(gene)[1L])) stop("Provide gene.")
  gene <- as.character(gene)[1L]
  path <- .ic2_xenium_transcript_path(x)
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required to stream transcripts.parquet.")
  }
  ds <- arrow::open_dataset(path)
  cols <- names(ds)
  gcol <- intersect(c("feature_name", "gene"), cols)[1L]
  xcol <- intersect(c("x_location", "x"), cols)[1L]
  ycol <- intersect(c("y_location", "y"), cols)[1L]
  if (anyNA(c(gcol, xcol, ycol))) stop("Transcript file is missing gene or x/y columns: ", path)
  expr <- arrow::Expression$create("==", arrow::Expression$field_ref(gcol), gene)
  if (!is.null(xlim)) {
    expr <- expr &
      arrow::Expression$create(">=", arrow::Expression$field_ref(xcol), xlim[[1L]]) &
      arrow::Expression$create("<=", arrow::Expression$field_ref(xcol), xlim[[2L]])
  }
  if (!is.null(ylim)) {
    expr <- expr &
      arrow::Expression$create(">=", arrow::Expression$field_ref(ycol), ylim[[1L]]) &
      arrow::Expression$create("<=", arrow::Expression$field_ref(ycol), ylim[[2L]])
  }
  reader <- arrow::Scanner$create(
    ds, filter = expr, columns = c(xcol, ycol), batch_size = 65536L
  )$ToRecordBatchReader()
  xs <- numeric(0)
  ys <- numeric(0)
  cap <- as.integer(max.points)
  repeat {
    batch <- reader$read_next_batch()
    if (is.null(batch)) break
    df <- as.data.frame(batch)
    xs <- c(xs, as.numeric(df[[xcol]]))
    ys <- c(ys, as.numeric(df[[ycol]]))
    if (length(xs) >= cap) break
  }
  if (length(xs) < 1L) stop("No transcripts for gene: ", gene)
  if (length(xs) > cap) {
    keep <- sample.int(length(xs), cap)
    xs <- xs[keep]
    ys <- ys[keep]
  }
  tab <- data.frame(x = xs, y = ys, gene = gene, stringsAsFactors = FALSE)
  if (!isTRUE(plot)) return(tab)
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(size = 0.15, alpha = 0.45, colour = "#c23b22") +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::labs(title = gene)
  if (!is.null(out.name)) ggplot2::ggsave(out.name, p, width = 7, height = 7, dpi = 150)
  p
}

#' @keywords internal
.ic2_xenium_transcript_path <- function(x) {
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    return(normalizePath(x, winslash = "/", mustWork = TRUE))
  }
  if (!inherits(x, "iCellR2")) stop("x should be an iCellR2 object or a transcripts.parquet path.")
  tx <- x@misc$spatial$transcripts
  if (is.character(tx) && length(tx) == 1L && !is.na(tx) && file.exists(tx)) {
    return(normalizePath(tx, winslash = "/", mustWork = TRUE))
  }
  root <- x@storage$root_dir
  hits <- list.files(file.path(root, "images"), pattern = "^xenium\\.json$", recursive = TRUE, full.names = TRUE)
  for (js in hits) {
    meta <- jsonlite::fromJSON(js, simplifyVector = TRUE)
    if (is.character(meta$transcripts) && file.exists(meta$transcripts)) {
      return(normalizePath(meta$transcripts, winslash = "/", mustWork = TRUE))
    }
  }
  stop("No transcripts.parquet. Run add.xenium() on a Xenium output folder.")
}
