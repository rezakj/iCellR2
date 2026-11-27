#' Cluster QC / stats plots (metadata only)
#'
#' Violin + box plots of QC columns in \code{cell_data} by \code{cluster_id}.
#' No assay I/O  -  safe at tens of millions of cells with downsampling.
#'
#' @param x An \code{iCellR2} object with \code{cluster_id} (and QC cols from
#'   \code{\link{qc.stats}} for box plots).
#' @param plot.type One of \code{"box.mito"} (default), \code{"box.umi"},
#'   \code{"box.gene"}, \code{"three.in.one"}, \code{"bar.cc"}, \code{"pie.cc"}.
#' @param conds.to.plot Optional \code{condition_id} values to include.
#'   Default \code{NULL} (all).
#' @param cell.color Jitter point color. Default \code{"slategray3"}.
#' @param cell.size Jitter point size. Default \code{0.4}.
#' @param jitter.width Horizontal jitter width. Default \code{0.15}.
#' @param cell.transparency Alpha. Default \code{0.5}.
#' @param box.color Box fill. Default \code{"red"}.
#' @param box.line.col Box outline. Default \code{"green"}.
#' @param back.col Panel background. Default \code{"white"}.
#' @param notch Notch boxplots. Default \code{FALSE}.
#' @param max.cells Max cells to plot. Default \code{500000}.
#' @param seed RNG seed for downsampling.
#' @param interactive If \code{TRUE}, return plotly (single plots only).
#'   Default \code{FALSE}.
#' @param out.name If interactive and non-\code{NULL}, also write HTML.
#' @param verbose Progress style.
#'
#' @return A \code{ggplot}, plotly widget, or \code{gridExtra} grob.
#'
#' @examples
#' \dontrun{
#' clust.stats.plot(obj)                          # mito by cluster
#' clust.stats.plot(obj, plot.type = "box.umi")
#' clust.stats.plot(obj, plot.type = "three.in.one")
#' clust.stats.plot(obj, conds.to.plot = "KO")
#' }
#'
#' @importFrom rlang .data
#' @export
clust.stats.plot <- function(
    x = NULL,
    plot.type = c(
      "box.mito", "box.umi", "box.gene",
      "three.in.one", "bar.cc", "pie.cc"
    ),
    conds.to.plot = NULL,
    cell.color = "slategray3",
    cell.size = 0.4,
    jitter.width = 0.15,
    cell.transparency = 0.5,
    box.color = "red",
    box.line.col = "green",
    back.col = "white",
    notch = FALSE,
    max.cells = 500000L,
    seed = 1L,
    interactive = FALSE,
    out.name = NULL,
    verbose = TRUE
) {
  .ic2_require_object(x)
  plot.type <- match.arg(plot.type)
  max.cells <- as.integer(max.cells)

  rep <- .ic2_new_reporter("clust.stats.plot", verbose, soda_root = x@storage$root_dir)

  meta <- x@cell_data
  if (!is.data.frame(meta) || nrow(meta) < 1L) {
    stop("'x@cell_data' is missing or empty.")
  }
  if (!"cluster_id" %in% colnames(meta)) {
    stop("No 'cluster_id' in cell_data. Run run.clustering() first.")
  }

  keep <- !is.na(meta$cluster_id)
  if (!is.null(conds.to.plot)) {
    if (!"condition_id" %in% colnames(meta)) {
      stop("'conds.to.plot' set but no 'condition_id' in cell_data.")
    }
    keep <- keep & meta$condition_id %in% conds.to.plot
  }
  meta <- meta[keep, , drop = FALSE]
  if (nrow(meta) < 1L) stop("No cells left after filtering.")

  ## Columns needed for this plot
  need <- switch(
    plot.type,
    box.mito = "mito.percent",
    box.umi = "UMIs",
    box.gene = "nGenes",
    three.in.one = c("mito.percent", "UMIs", "nGenes"),
    bar.cc = character(0),
    pie.cc = character(0)
  )
  if (length(need)) {
    missing <- setdiff(need, colnames(meta))
    if (length(missing)) {
      stop(
        "Missing QC columns in cell_data: ", paste(missing, collapse = ", "),
        ". Run qc.stats() first."
      )
    }
  }

  ## Cell-cycle phase for bar/pie
  if (plot.type %in% c("bar.cc", "pie.cc")) {
    meta$Phase <- .ic2_clust_stats_phase(meta)
  }

  n_all <- nrow(meta)
  if (!is.null(max.cells) && n_all > max.cells) {
    set.seed(as.integer(seed))
    meta <- meta[sample.int(n_all, max.cells), , drop = FALSE]
    .ic2_step(rep, paste0(
      "  Downsampled ", .ic2_fmt_num(n_all), "  ->  ",
      .ic2_fmt_num(nrow(meta)), " cells"
    ))
  }

  meta$cluster <- .ic2_factor_sorted(meta$cluster_id)

  rep <- .ic2_banner_start(rep, c(
    paste0("Plot: ", plot.type, "   |   by cluster_id"),
    paste0(
      "Cells plotted: ", .ic2_fmt_num(nrow(meta)),
      if (!is.null(conds.to.plot)) {
        paste0("   |   conditions: ", paste(conds.to.plot, collapse = ","))
      } else ""
    )
  ))

  .ic2_step(rep, "[1/1] Building plot")

  build_box <- function(df, yvar, ylab.txt) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$cluster, y = .data[[yvar]])) +
      ggplot2::geom_jitter(
        height = 0, width = jitter.width, color = cell.color,
        size = cell.size, alpha = cell.transparency
      ) +
      ggplot2::geom_violin(
        trim = FALSE, col = "black", alpha = cell.transparency
      ) +
      ggplot2::geom_boxplot(
        fill = box.color, col = box.line.col, notch = notch,
        outlier.shape = NA, alpha = cell.transparency
      ) +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
      ggplot2::xlab("cluster") +
      ggplot2::ylab(ylab.txt) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = back.col, colour = "grey30"),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  }

  if (identical(plot.type, "box.mito")) {
    p <- build_box(meta, "mito.percent", "percent of mito genes per cell")
  } else if (identical(plot.type, "box.umi")) {
    p <- build_box(meta, "UMIs", "number of UMIs per cell")
  } else if (identical(plot.type, "box.gene")) {
    p <- build_box(meta, "nGenes", "number of genes per cell")
  } else if (identical(plot.type, "three.in.one")) {
    if (isTRUE(interactive)) {
      warning("three.in.one is static; setting interactive = FALSE.")
      interactive <- FALSE
    }
    p1 <- build_box(meta, "nGenes", "nGenes")
    p2 <- build_box(meta, "UMIs", "UMIs")
    p3 <- build_box(meta, "mito.percent", "mito.percent")
    .ic2_finish(rep, extra_lines = "Returned three.in.one panel")
    return(gridExtra::grid.arrange(p1, p2, p3, ncol = 3))
  } else if (identical(plot.type, "bar.cc")) {
    ## Stacked counts by phase within cluster
    tab <- as.data.frame(table(cluster = meta$cluster, Phase = meta$Phase),
                         stringsAsFactors = FALSE)
    colnames(tab)[3] <- "n"
    p <- ggplot2::ggplot(
      tab,
      ggplot2::aes(x = .data$cluster, y = .data$n, fill = .data$Phase)
    ) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::xlab("cluster") +
      ggplot2::ylab("cell count") +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = back.col, colour = "grey30"),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  } else {
    ## pie.cc  -  phase composition per cluster (facet)
    tab <- as.data.frame(table(cluster = meta$cluster, Phase = meta$Phase),
                         stringsAsFactors = FALSE)
    colnames(tab)[3] <- "n"
    p <- ggplot2::ggplot(
      tab,
      ggplot2::aes(x = "", y = .data$n, fill = .data$Phase)
    ) +
      ggplot2::geom_col(position = "fill", width = 1) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::facet_wrap(~cluster) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }

  if (isTRUE(interactive)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive = TRUE.")
    }
    p <- plotly::ggplotly(p)
    if (!is.null(out.name) && nzchar(out.name)) {
      htmlwidgets::saveWidget(p, paste0(out.name, ".html"))
    }
  }

  .ic2_finish(rep, extra_lines = paste0(
    "Plot: ", plot.type, "   |   cells: ", .ic2_fmt_num(nrow(meta))
  ))
  p
}

#' Infer cell-cycle Phase label from cell_data
#' @keywords internal
.ic2_clust_stats_phase <- function(meta) {
  if ("Phase" %in% colnames(meta)) {
    return(as.character(meta$Phase))
  }
  if (all(c("S.phase.probability", "g2m.phase.probability") %in% colnames(meta))) {
    s <- as.numeric(meta$S.phase.probability)
    g <- as.numeric(meta$g2m.phase.probability)
    phase <- rep("G1", length(s))
    phase[is.finite(s) & is.finite(g) & s >= g & s >= 0.5] <- "S"
    phase[is.finite(s) & is.finite(g) & g > s & g >= 0.5] <- "G2M"
    return(phase)
  }
  stop(
    "Cell-cycle plots need 'Phase' or ",
    "'S.phase.probability' + 'g2m.phase.probability' in cell_data. ",
    "Run qc.stats() / cell.cycle() first."
  )
}
