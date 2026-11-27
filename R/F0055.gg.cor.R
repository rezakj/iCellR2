#' Gene-gene correlation scatter (chunked fetch)
#'
#' Pulls two genes for a downsampled cell set from on-disk assays and
#' plots expression vs expression, colored by \code{cluster_id}.
#' Correlation is Pearson on the plotted (transformed) values.
#'
#' @param x An \code{iCellR2} object.
#' @param gene1,gene2 Gene ids (length 1 each, must differ).
#' @param modality Modality. Default \code{"rna"}.
#' @param data.type Which matrix to plot. Default \code{"main.data"}
#'   (normalized). Same options as \code{\link{gene.plot}}.
#'   Ignored if \code{assay.name} is set.
#' @param assay.name Assay name (overrides \code{data.type}).
#' @param transform \code{"auto"} (default): \code{log2p1} unless assay is
#'   \code{scale.data}; else \code{"none"}, \code{"log2p1"}, \code{"log1p"}.
#' @param conds Optional \code{condition_id} values to include (vector OK).
#' @param clusts Optional \code{cluster_id} values to include (vector OK).
#' @param max.cells Max cells to plot / correlate. Default \code{500000}.
#' @param cell.size Point size. Default \code{0.4}.
#' @param cell.transparency Alpha. Default \code{0.5}.
#' @param back.col Panel background. Default \code{"white"}.
#' @param chunk.size Cells per HDF5 read when fetching genes.
#' @param seed RNG seed for downsampling.
#' @param interactive If \code{TRUE}, return plotly. Default \code{FALSE}.
#' @param out.name If interactive and non-\code{NULL}, also write HTML.
#' @param verbose Progress style.
#'
#' @return A \code{ggplot} or plotly widget. Correlation stats are in the
#'   subtitle / attribute \code{attr(,"stats")}.
#'
#' @examples
#' \dontrun{
#' gg.cor(obj, gene1 = "CD3D", gene2 = "MS4A1")
#' gg.cor(obj, gene1 = "CD3D", gene2 = "CD8A", clusts = c(1, 2, 5))
#' gg.cor(obj, gene1 = "CD3D", gene2 = "MS4A1", data.type = "imputed.data")
#' }
#'
#' @importFrom rlang .data
#' @export
gg.cor <- function(
    x = NULL,
    gene1 = NULL,
    gene2 = NULL,
    modality = "rna",
    data.type = "main.data",
    assay.name = NULL,
    transform = c("auto", "none", "log2p1", "log1p"),
    conds = NULL,
    clusts = NULL,
    max.cells = 500000L,
    cell.size = 0.4,
    cell.transparency = 0.5,
    back.col = "white",
    chunk.size = 200000L,
    seed = 1L,
    interactive = FALSE,
    out.name = NULL,
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  transform <- match.arg(transform)
  max.cells <- as.integer(max.cells)
  chunk.size <- as.integer(chunk.size)

  if (is.null(gene1) || is.null(gene2) || length(gene1) != 1L || length(gene2) != 1L) {
    stop("Provide single 'gene1' and 'gene2'.")
  }
  gene1 <- as.character(gene1)
  gene2 <- as.character(gene2)
  if (identical(gene1, gene2)) stop("'gene1' and 'gene2' must be different.")

  assay.name <- .ic2_resolve_plot_assay(
    x, modality, assay.name = assay.name, data.type = data.type, prefer.data = TRUE
  )
  if (identical(transform, "auto")) {
    transform <- if (identical(assay.name, "scale.data")) "none" else "log2p1"
  }

  rep <- .ic2_new_reporter("gg.cor", verbose, soda_root = x@storage$root_dir)

  meta <- x@cell_data
  if (!is.data.frame(meta) || nrow(meta) < 2L) {
    stop("'x@cell_data' needs at least 2 cells.")
  }
  keep <- rep(TRUE, nrow(meta))
  if (!is.null(conds)) {
    if (!"condition_id" %in% colnames(meta)) {
      stop("'conds' set but no 'condition_id' in cell_data.")
    }
    keep <- keep & meta$condition_id %in% conds
  }
  if (!is.null(clusts)) {
    if (!"cluster_id" %in% colnames(meta)) {
      stop("'clusts' set but no 'cluster_id' in cell_data.")
    }
    keep <- keep & meta$cluster_id %in% clusts
  }
  meta <- meta[keep, , drop = FALSE]
  if (nrow(meta) < 2L) stop("Fewer than 2 cells after filtering.")

  cell_ids <- rownames(meta)
  n_all <- length(cell_ids)
  if (!is.null(max.cells) && n_all > max.cells) {
    set.seed(as.integer(seed))
    take <- sample.int(n_all, max.cells)
    meta <- meta[take, , drop = FALSE]
    cell_ids <- rownames(meta)
    .ic2_step(rep, paste0(
      "  Downsampled ", .ic2_fmt_num(n_all), "  ->  ",
      .ic2_fmt_num(length(cell_ids)), " cells"
    ))
  }

  .ic2_step(rep, paste0(
    "[1/2] Fetching ", gene1, " + ", gene2,
    " from assay '", assay.name, "' (transform=", transform, ")"
  ))
  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  expr <- .ic2_load_features_for_cells(
    handles = handles,
    features = c(gene1, gene2),
    cell_ids = cell_ids,
    chunk.size = chunk.size,
    transform = transform
  )
  found <- intersect(c(gene1, gene2), colnames(expr))
  if (length(found) < 2L) {
    miss <- setdiff(c(gene1, gene2), found)
    stop("Gene(s) not found in assay: ", paste(miss, collapse = ", "))
  }

  xval <- as.numeric(expr[cell_ids, gene1])
  yval <- as.numeric(expr[cell_ids, gene2])
  ct <- suppressWarnings(stats::cor.test(xval, yval, method = "pearson"))
  r <- unname(ct$estimate)
  if (!is.finite(r)) r <- NA_real_
  pval <- ct$p.value
  if (!is.finite(pval) || pval == 0) pval_lab <- "< 2.2e-16" else pval_lab <- format(pval, digits = 3)
  sub <- paste0(
    "P value: ", pval_lab,
    " | correlation value: ", round(r, 3)
  )

  has_cl <- "cluster_id" %in% colnames(meta)
  df <- data.frame(
    gene1 = xval,
    gene2 = yval,
    cell = cell_ids,
    stringsAsFactors = FALSE
  )
  if (has_cl) {
    df$cluster <- .ic2_factor_sorted(meta$cluster_id)
  }

  rep <- .ic2_banner_start(rep, c(
    paste0(gene1, " vs ", gene2),
    paste0(
      "Assay: ", assay.name, "   |   transform: ", transform,
      "   |   cells: ", .ic2_fmt_num(nrow(df))
    ),
    sub
  ))

  .ic2_step(rep, "[2/2] Building plot")
  xlab <- paste0(gene1, " (", transform, ")")
  ylab <- paste0(gene2, " (", transform, ")")

  if (has_cl) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$gene1, y = .data$gene2,
        colour = .data$cluster, text = .data$cell
      )
    ) +
      ggplot2::geom_point(size = cell.size, alpha = cell.transparency) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))
      )
  } else {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$gene1, y = .data$gene2, text = .data$cell)
    ) +
      ggplot2::geom_point(size = cell.size, alpha = cell.transparency, colour = "grey30")
  }

  p <- p +
    ggplot2::labs(
      title = "Gene-gene correlation",
      subtitle = sub,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = back.col, colour = "grey30"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (isTRUE(interactive)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive = TRUE.")
    }
    p <- plotly::ggplotly(p)
    if (!is.null(out.name) && nzchar(out.name)) {
      htmlwidgets::saveWidget(p, paste0(out.name, ".html"))
    }
  }

  attr(p, "stats") <- list(
    gene1 = gene1, gene2 = gene2,
    correlation = r, p.value = ct$p.value,
    n_cells = nrow(df), transform = transform, assay.name = assay.name
  )

  .ic2_finish(rep, extra_lines = paste0(
    "r = ", round(r, 3), "   |   cells: ", .ic2_fmt_num(nrow(df))
  ))
  p
}
