#' Cluster and condition frequencies
#'
#' Counts cells in each \code{cluster_id} by \code{condition_id}. This is the
#' iCellR2 form of classic \code{clust.cond.info}.
#'
#' \code{normalize.by = "percentage"} is the share of each condition that falls
#' in a cluster (columns sum to 100 within a condition). \code{"sf"} divides
#' raw counts by the size factor of that condition (condition size / smallest
#' condition). \code{normalize.ncell = FALSE} plots the raw counts.
#'
#' @param x An iCellR2 object with \code{cluster_id} and \code{condition_id}.
#' @param plot.type \code{"pie"} (conditions inside each cluster),
#'   \code{"bar"} (same, as bars), \code{"pie.cond"} (clusters inside each
#'   condition), or \code{"bar.cond"}.
#' @param my.out.put \code{"data"} (default) stores the table and returns the
#'   object. \code{"plot"} returns the ggplot.
#' @param normalize.ncell If \code{TRUE} (default), plot the normalized value
#'   chosen by \code{normalize.by}. If \code{FALSE}, plot raw counts.
#' @param normalize.by \code{"percentage"} (default) or \code{"sf"}.
#'
#' @return The object when \code{my.out.put = "data"}, with the table at
#'   \code{x@my.freq} and \code{x@misc$my.freq}. A ggplot when
#'   \code{my.out.put = "plot"}.
#'
#' @examples
#' \dontrun{
#' clust.cond.info(obj, plot.type = "pie", normalize.ncell = TRUE,
#'                 my.out.put = "plot", normalize.by = "percentage")
#' clust.cond.info(obj, plot.type = "bar", normalize.ncell = TRUE,
#'                 my.out.put = "plot", normalize.by = "percentage")
#' clust.cond.info(obj, plot.type = "pie.cond", normalize.ncell = TRUE,
#'                 my.out.put = "plot", normalize.by = "percentage")
#' clust.cond.info(obj, plot.type = "bar.cond", normalize.ncell = TRUE,
#'                 my.out.put = "plot", normalize.by = "percentage")
#' obj <- clust.cond.info(obj)
#' head(obj@my.freq)
#' }
#'
#' @export
clust.cond.info <- function(
    x = NULL,
    plot.type = c("pie", "bar", "pie.cond", "bar.cond"),
    my.out.put = c("data", "plot"),
    normalize.ncell = TRUE,
    normalize.by = c("percentage", "sf")
) {
  .ic2_require_object(x)
  plot.type <- match.arg(plot.type)
  my.out.put <- match.arg(my.out.put)
  normalize.by <- match.arg(normalize.by)

  meta <- x@cell_data
  if (!"cluster_id" %in% colnames(meta)) {
    stop("No cluster_id in cell_data. Run iclust() first.")
  }
  if (!"condition_id" %in% colnames(meta)) {
    stop("No condition_id in cell_data. Aggregate samples before make.obj().")
  }
  meta <- meta[!is.na(meta$cluster_id) & !is.na(meta$condition_id), , drop = FALSE]
  conds <- unique(as.character(meta$condition_id))
  if (length(conds) < 2L) {
    stop("clust.cond.info needs at least two conditions. Found ", length(conds), ".")
  }

  tab <- as.data.frame(table(
    conditions = as.character(meta$condition_id),
    clusters = as.character(meta$cluster_id)
  ), stringsAsFactors = FALSE)
  totals <- as.data.frame(table(conditions = as.character(meta$condition_id)), stringsAsFactors = FALSE)
  colnames(totals)[2L] <- "TC"
  min_n <- min(totals$TC)
  totals$SF <- round(totals$TC / min_n, 3)
  DATA <- merge(totals, tab, by = "conditions")
  DATA$Norm.Freq <- round(DATA$Freq / DATA$SF, 3)
  DATA$percentage <- round((DATA$Freq / DATA$TC) * 100, 2)
  DATA$clusters <- factor(DATA$clusters, levels = unique(DATA$clusters[order(suppressWarnings(as.numeric(DATA$clusters)), DATA$clusters)]))

  y_col <- "Freq"
  if (isTRUE(normalize.ncell)) {
    y_col <- if (identical(normalize.by, "sf")) "Norm.Freq" else "percentage"
  }
  y_lab <- switch(y_col, Freq = "cells", Norm.Freq = "size-factor normalized cells", percentage = "percent of condition")
  DATA$y <- DATA[[y_col]]

  p <- if (identical(plot.type, "bar")) {
    ggplot2::ggplot(DATA, ggplot2::aes(x = .data$conditions, y = .data$y, fill = .data$conditions)) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~clusters, scales = "free") +
      ggplot2::labs(y = y_lab, x = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else if (identical(plot.type, "bar.cond")) {
    ggplot2::ggplot(DATA, ggplot2::aes(x = .data$conditions, y = .data$y, fill = .data$clusters)) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~conditions, scales = "free") +
      ggplot2::labs(y = y_lab, x = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else if (identical(plot.type, "pie")) {
    ggplot2::ggplot(DATA, ggplot2::aes(x = "", y = .data$y, fill = .data$conditions)) +
      ggplot2::geom_col(position = "fill", width = 1) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::facet_wrap(~clusters) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  } else {
    ggplot2::ggplot(DATA, ggplot2::aes(x = "", y = .data$y, fill = .data$clusters)) +
      ggplot2::geom_col(position = "fill", width = 1) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::facet_wrap(~conditions) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  }

  if (identical(my.out.put, "plot")) return(p)

  x@my.freq <- DATA
  x@misc$my.freq <- DATA
  if (!is.null(x@storage$root_dir)) {
    dest <- file.path(x@storage$root_dir, "metadata", "clustering")
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      DATA,
      file.path(dest, "clust_cond_freq.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  x
}
