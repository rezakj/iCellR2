#' Highlight clonotypes on an existing RNA embedding
#'
#' Colors the UMAP (or PCA / KNetL) already stored on the object.
#' Background cells are down-sampled for plotting; every cell of the
#' requested clonotype is kept. The full embedding is not sent to an
#' interactive widget.
#'
#' @param x An iCellR2 object after \code{\link{add.vdj}} and a reduction.
#' @param clono Clonotype id or ids. If \code{NULL}, color by expansion bin.
#' @param plot.type \code{"umap"} (default), \code{"pca"}, \code{"knetl"},
#'   or \code{"aligned_umap"}.
#' @param max.cells Maximum points drawn. Highlighted clonotypes are extra.
#' @param cell.size,cell.transparency,cell.colors,back.col Passed to ggplot2.
#' @param interactive If \code{TRUE}, write an HTML widget of the sampled plot.
#' @param out.name PNG (or HTML when \code{interactive = TRUE}) path.
#' @param barcode.column,clonotype.column,clust.dim,conds.to.plot,box.cell.col
#'   Classic arguments. Embedding dimension follows \code{clust.dim} when
#'   \code{plot.type} has that many columns. Condition filtering uses
#'   \code{cell_data$condition_id}.
#' @return A ggplot (or the widget path when interactive).
#' @examples
#' \dontrun{
#' clono.plot(obj, clono = "clonotype1", plot.type = "umap", out.name = "clono.png")
#' }
#' @export
clono.plot <- function(
    x = NULL,
    clono = NULL,
    plot.type = "umap",
    max.cells = 200000L,
    cell.size = 0.4,
    cell.transparency = 0.7,
    cell.colors = c("#c23b22", "grey80"),
    back.col = "white",
    interactive = FALSE,
    out.name = NULL,
    barcode.column = NULL,
    clonotype.column = NULL,
    clust.dim = 2L,
    conds.to.plot = NULL,
    box.cell.col = "black"
) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  if (!"clonotype_id" %in% colnames(x@cell_data)) {
    stop("No clonotype_id on cell_data. Run add.vdj() first.")
  }
  modality <- "rna"
  emb <- .ic2_get_plot_embedding(x, modality, plot.type)
  cell_ids <- .ic2_adt_cell_ids(x)
  if (is.null(rownames(emb)) || !any(rownames(emb) %in% cell_ids)) {
    stop("Embedding rownames do not match cell_data.")
  }
  emb <- emb[intersect(rownames(emb), cell_ids), , drop = FALSE]
  meta <- x@cell_data
  rownames(meta) <- cell_ids
  meta <- meta[rownames(emb), , drop = FALSE]
  if (!is.null(conds.to.plot)) {
    if (!"condition_id" %in% colnames(meta)) stop("No condition_id column to filter.")
    keep_c <- meta$condition_id %in% conds.to.plot
    emb <- emb[keep_c, , drop = FALSE]
    meta <- meta[keep_c, , drop = FALSE]
  }
  dims <- seq_len(min(as.integer(clust.dim), ncol(emb)))
  if (length(dims) < 2L) stop("Embedding needs at least 2 columns.")

  if (is.null(clono)) {
    group <- as.character(meta$vdj_expand)
    group[is.na(group)] <- "no VDJ"
  } else {
    group <- ifelse(meta$clonotype_id %in% as.character(clono), "clonotype", "other")
  }
  max.cells <- as.integer(max.cells)
  idx <- seq_len(nrow(emb))
  if (nrow(emb) > max.cells) {
    if (!is.null(clono)) {
      hit <- which(group == "clonotype")
      rest <- setdiff(idx, hit)
      n_rest <- max(0L, max.cells - length(hit))
      if (length(rest) > n_rest) rest <- sample(rest, n_rest)
      idx <- sort(c(hit, rest))
    } else {
      idx <- sort(sample.int(nrow(emb), max.cells))
    }
  }
  emb <- emb[idx, , drop = FALSE]
  group <- group[idx]
  df <- data.frame(
    d1 = emb[, dims[1L]],
    d2 = emb[, dims[2L]],
    group = factor(group),
    stringsAsFactors = FALSE
  )
  if (length(dims) >= 3L) df$d3 <- emb[, dims[3L]]

  title <- if (is.null(clono)) {
    paste(plot.type, " |  clonotype expansion")
  } else {
    paste(plot.type, " | ", paste(as.character(clono), collapse = ", "))
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = d1, y = d2, colour = group)) +
    ggplot2::geom_point(size = cell.size, alpha = cell.transparency) +
    ggplot2::scale_colour_manual(values = .ic2_vdj_plot_colors(levels(df$group), cell.colors)) +
    ggplot2::labs(x = paste0(plot.type, dims[1L]), y = paste0(plot.type, dims[2L]), title = title, colour = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = back.col),
      plot.background = ggplot2::element_rect(fill = back.col, colour = box.cell.col)
    )
  if (!is.null(out.name) && !isTRUE(interactive)) {
    ggplot2::ggsave(out.name, p, width = 7, height = 6, dpi = 150)
  }
  if (isTRUE(interactive)) {
    if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
      stop("Packages 'plotly' and 'htmlwidgets' are required for interactive = TRUE.")
    }
    html <- if (is.null(out.name)) "clono.plot.html" else {
      if (grepl("\\.html$", out.name, ignore.case = TRUE)) out.name else paste0(out.name, ".html")
    }
    w <- if (length(dims) >= 3L) {
      plotly::plot_ly(df, x = ~d1, y = ~d2, z = ~d3, color = ~group, type = "scatter3d",
                      marker = list(size = cell.size + 2))
    } else {
      plotly::ggplotly(p)
    }
    htmlwidgets::saveWidget(w, html)
    return(invisible(html))
  }
  p
}

#' @keywords internal
.ic2_vdj_plot_colors <- function(levels, cell.colors) {
  base <- setNames(rep("grey75", length(levels)), levels)
  if ("clonotype" %in% levels) base[["clonotype"]] <- cell.colors[[1L]]
  if ("other" %in% levels) base[["other"]] <- cell.colors[[min(2L, length(cell.colors))]]
  named <- c(
    "no VDJ" = "grey90", singleton = "#4C78A8", doubleton = "#72B7B2",
    Q1 = "#F2CF5B", Mid = "#F58518", Q3 = "#E45756"
  )
  for (nm in intersect(names(named), levels)) base[[nm]] <- unname(named[[nm]])
  base
}
