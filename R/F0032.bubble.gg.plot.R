#' Bubble / dot plot of genes by cluster (or condition)
#'
#' Scaled mean expression (color) and percent of cells expressing (size)
#' for selected genes, streamed from the on-disk assay (no full matrix).
#'
#' @param x An \code{iCellR2} object.
#' @param gene Character vector of gene ids to plot.
#' @param group.by \code{"clusters"} (default; by \code{cluster_id}) or
#'   \code{"conditions"} (by \code{condition_id}).
#' @param modality Modality. Default \code{"rna"}.
#' @param data.type Which matrix to plot. Default \code{"main.data"}
#'   (normalized). Same options as \code{\link{gene.plot}}
#'   (\code{"raw.data"}, \code{"imputed.data"}, \code{"scaled.data"}, ...).
#'   Ignored if \code{assay.name} is set.
#' @param assay.name Assay name (overrides \code{data.type}).
#' @param transform \code{"auto"} (default): \code{log2p1} unless assay is
#'   \code{scale.data}; else \code{"none"}, \code{"log2p1"}, \code{"log1p"}.
#' @param conds.to.plot Optional condition ids to include.
#' @param clusts.to.plot Optional cluster ids to include.
#' @param min.scale,max.scale Clip scaled expression. Defaults \code{-2.5},
#'   \code{2.5}.
#' @param colour Aesthetic for point color: \code{"Expression"} (default)
#'   or \code{"Percent.Expressed"}.
#' @param size Aesthetic for point size: \code{"Percent.Expressed"} (default)
#'   or \code{"Expression"}.
#' @param heat.colors Length-2 color vector (low, high) for the color
#'   gradient. Default \code{c("blue", "red")}.
#' @param chunk.size Cells per HDF5 read. Default \code{200000}.
#' @param interactive If \code{TRUE}, return plotly. Default \code{FALSE}.
#' @param out.name If interactive and non-\code{NULL}, also write HTML.
#' @param verbose Progress style.
#'
#' @return A \code{ggplot} or plotly widget.
#'
#' @examples
#' \dontrun{
#' bubble.gg.plot(obj, gene = top.markers(markers, topde = 5))
#' bubble.gg.plot(obj, gene = c("CD3E", "MS4A1", "LYZ"),
#'                data.type = "imputed.data")
#' }
#'
#' @importFrom rlang .data
#' @export
bubble.gg.plot <- function(
    x = NULL,
    gene = NULL,
    group.by = c("clusters", "conditions"),
    modality = "rna",
    data.type = "main.data",
    assay.name = NULL,
    transform = c("auto", "none", "log2p1", "log1p"),
    conds.to.plot = NULL,
    clusts.to.plot = NULL,
    min.scale = -2.5,
    max.scale = 2.5,
    colour = c("Expression", "Percent.Expressed"),
    size = c("Percent.Expressed", "Expression"),
    heat.colors = c("blue", "red"),
    chunk.size = 200000L,
    interactive = FALSE,
    out.name = NULL,
    verbose = TRUE
) {
  .ic2_require_object(x)
  group.by <- match.arg(group.by)
  modality <- tolower(modality)
  transform <- match.arg(transform)
  colour <- match.arg(colour)
  size <- match.arg(size)
  chunk.size <- as.integer(chunk.size)

  if (is.null(gene) || length(gene) < 1L) {
    stop("Provide at least one gene in 'gene'.")
  }
  gene <- unique(as.character(gene))
  if (length(heat.colors) != 2L) {
    stop("'heat.colors' must be length 2 (low, high).")
  }

  assay.name <- .ic2_resolve_plot_assay(
    x, modality, assay.name = assay.name, data.type = data.type, prefer.data = TRUE
  )
  if (identical(transform, "auto")) {
    transform <- if (identical(assay.name, "scale.data")) "none" else "log2p1"
  }

  rep <- .ic2_new_reporter("bubble.gg.plot", verbose, soda_root = x@storage$root_dir)
  meta <- x@cell_data
  if (!is.data.frame(meta) || nrow(meta) < 1L) {
    stop("'x@cell_data' is missing or empty.")
  }

  keep <- rep(TRUE, nrow(meta))
  if (!is.null(conds.to.plot)) {
    if (!"condition_id" %in% colnames(meta)) {
      stop("'conds.to.plot' set but no 'condition_id' in cell_data.")
    }
    keep <- keep & meta$condition_id %in% conds.to.plot
  }
  if (!is.null(clusts.to.plot)) {
    if (!"cluster_id" %in% colnames(meta)) {
      stop("'clusts.to.plot' set but no 'cluster_id' in cell_data.")
    }
    keep <- keep & meta$cluster_id %in% clusts.to.plot
  }
  meta <- meta[keep, , drop = FALSE]
  if (nrow(meta) < 1L) stop("No cells left after filtering.")

  group_col <- if (identical(group.by, "clusters")) "cluster_id" else "condition_id"
  if (!group_col %in% colnames(meta)) {
    stop("No '", group_col, "' in cell_data.")
  }

  cell_to_group <- setNames(
    as.character(meta[[group_col]]),
    rownames(meta)
  )

  .ic2_step(rep, "[1/2] Streaming mean + % expressed by group")
  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  size_factors <- .ic2_assay_size_factors(x, modality, assay.name)
  stats <- .ic2_feature_group_means(
    handles = handles,
    features = gene,
    cell_to_group = cell_to_group,
    transform = transform,
    chunk.size = chunk.size,
    return.pct = TRUE,
    size_factors = size_factors
  )
  mean_mat <- stats$mean
  pct_mat <- stats$pct
  gene <- rownames(mean_mat)
  if (length(gene) < 1L) stop("None of the requested genes were found.")

  ## Numeric-aware column order
  cn <- colnames(mean_mat)
  num <- suppressWarnings(as.numeric(cn))
  if (all(is.finite(num))) {
    ord <- order(num)
    mean_mat <- mean_mat[, ord, drop = FALSE]
    pct_mat <- pct_mat[, ord, drop = FALSE]
  } else {
    ord <- order(cn)
    mean_mat <- mean_mat[, ord, drop = FALSE]
    pct_mat <- pct_mat[, ord, drop = FALSE]
  }

  ## Scale expression across groups (genes  x  groups), then clip
  .ic2_step(rep, "[2/2] Scaling and plotting")
  scaled <- t(scale(t(mean_mat)))
  scaled[!is.finite(scaled)] <- 0
  scaled[scaled > max.scale] <- max.scale
  scaled[scaled < min.scale] <- min.scale

  ## Melt to long form (base R; no reshape dependency)
  groups <- colnames(mean_mat)
  dat <- data.frame(
    gene = factor(rep(gene, times = length(groups)), levels = rev(gene)),
    cluster = factor(rep(groups, each = length(gene)), levels = groups),
    Expression = as.vector(scaled),
    Percent.Expressed = as.vector(pct_mat),
    stringsAsFactors = FALSE
  )

  colour_sym <- rlang::sym(colour)
  size_sym <- rlang::sym(size)
  size_name <- if (identical(size, "Percent.Expressed")) "% expressed" else size
  colour_name <- if (identical(colour, "Expression")) "Expression (scaled)" else colour
  ## Size range matches classic iCellR bubble (visible 20 / 50 / 80% legend)
  size_range <- c(1.2, 9)
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data$cluster,
      y = .data$gene,
      colour = !!colour_sym,
      size = !!size_sym
    )
  ) +
    ggplot2::geom_point(shape = 16, stroke = 0) +
    ggplot2::scale_size_continuous(
      name = size_name,
      range = size_range,
      breaks = c(20, 50, 80),
      limits = c(0, 100),
      labels = c("20% expr", "50% expr", "80% expr")
    ) +
    ggplot2::scale_colour_gradient(
      low = heat.colors[1L],
      high = heat.colors[2L],
      name = colour_name,
      limits = c(min.scale, max.scale)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        order = 1,
        override.aes = list(colour = "grey50")
      ),
      colour = ggplot2::guide_colorbar(order = 2)
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.6),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9)
    ) +
    ggplot2::xlab(if (identical(group.by, "clusters")) "cluster" else "condition") +
    ggplot2::ylab("gene")

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
    "Bubble: ", length(gene), " genes x ", length(groups), " ", group.by
  ))
  invisible(p)
}
