#' Create MA and Volcano plots from a DE table
#'
#' Table-only plotting (no assay I/O). Expects a data frame with
#' \code{log2FoldChange} and \code{padj} or \code{pval} (and
#' \code{baseMean} for MA plots). Gene names should be rownames or a
#' \code{gene} column.
#'
#' @param x DE results \code{data.frame}.
#' @param sig.value \code{"padj"} (default) or \code{"pval"}.
#' @param sig.line Significance threshold on the raw p / padj scale.
#'   Default \code{0.1}.
#' @param plot.type \code{"volcano"} (default) or \code{"ma"}.
#' @param x.limit,y.limit Axis limits when \code{limit.force = TRUE}.
#' @param limit.force Force axis limits.
#' @param scale.ax Log1p-scale the secondary axis (y for volcano, x for MA).
#' @param dot.size Point size.
#' @param dot.transparency Point alpha.
#' @param dot.col Colors for up / down / NS.
#' @param interactive If \code{TRUE}, return / write plotly HTML.
#' @param out.name HTML stem when \code{interactive = TRUE}.
#'
#' @return A ggplot (static) or plotly widget / side-effect HTML write.
#'
#' @export
volcano.ma.plot <- function(
    x = NULL,
    sig.value = c("padj", "pval"),
    sig.line = 0.1,
    plot.type = c("volcano", "ma"),
    x.limit = 2,
    y.limit = 2,
    limit.force = FALSE,
    scale.ax = TRUE,
    dot.size = 1.75,
    dot.transparency = 0.5,
    dot.col = c("#E64B35", "#3182bd", "#636363"),
    interactive = FALSE,
    out.name = "plot"
) {
  sig.value <- match.arg(sig.value)
  plot.type <- match.arg(plot.type)
  if (is.null(x) || !is.data.frame(x)) {
    stop("'x' must be a DE data.frame.")
  }
  DATA <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(DATA)) {
    DATA$gene <- rownames(DATA)
  }
  if (!"log2FoldChange" %in% colnames(DATA)) {
    stop("DE table needs a 'log2FoldChange' column.")
  }
  pcol <- sig.value
  if (!pcol %in% colnames(DATA)) {
    stop("DE table needs a '", pcol, "' column.")
  }

  Pval <- -log10(as.numeric(sig.line))
  lfc <- as.numeric(DATA$log2FoldChange)
  pv <- -log10(as.numeric(DATA[[pcol]]))
  pv[!is.finite(pv)] <- 0
  color <- ifelse(lfc > 0 & pv > Pval, "Treated",
                  ifelse(lfc < 0 & pv > Pval, "Untreated", "none"))

  if (identical(plot.type, "volcano")) {
    data <- data.frame(
      gene = as.character(DATA$gene),
      pvalue = pv,
      lfc = lfc,
      color = color,
      stringsAsFactors = FALSE
    )
    Myylab <- if (identical(sig.value, "padj")) {
      "-log10 (adjusted p-value)"
    } else {
      "-log10 (p-value)"
    }
    myPLOT <- ggplot2::ggplot(
      data, ggplot2::aes(x = .data$lfc, y = .data$pvalue, text = .data$gene)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = factor(.data$color)),
        size = dot.size, alpha = dot.transparency, na.rm = TRUE
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(label = "Volcano Plot", subtitle = "Colored by directionality") +
      ggplot2::xlab("log2 (Fold Change)") +
      ggplot2::ylab(Myylab) +
      ggplot2::geom_vline(xintercept = 0, colour = "black") +
      ggplot2::geom_hline(yintercept = Pval, colour = "black") +
      ggplot2::scale_color_manual(values = c(
        Treated = dot.col[1], Untreated = dot.col[2], none = dot.col[3]
      ))
    if (isTRUE(limit.force)) {
      myPLOT <- myPLOT + ggplot2::scale_x_continuous(limits = c(-x.limit, x.limit))
    }
    if (isTRUE(scale.ax)) {
      myPLOT <- myPLOT + ggplot2::scale_y_continuous(trans = "log1p")
    } else if (isTRUE(limit.force)) {
      myPLOT <- myPLOT + ggplot2::scale_y_continuous(limits = c(0, y.limit))
    }
  } else {
    if (!"baseMean" %in% colnames(DATA)) {
      stop("MA plot needs a 'baseMean' column.")
    }
    data <- data.frame(
      gene = as.character(DATA$gene),
      baseMean = as.numeric(DATA$baseMean),
      pvalue = pv,
      lfc = lfc,
      color = color,
      stringsAsFactors = FALSE
    )
    myPLOT <- ggplot2::ggplot(
      data, ggplot2::aes(x = .data$baseMean, y = .data$lfc, text = .data$gene)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = factor(.data$color)),
        size = dot.size, alpha = dot.transparency, na.rm = TRUE
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(label = "MA Plot", subtitle = "Colored by directionality") +
      ggplot2::xlab("Base Mean") +
      ggplot2::ylab("log2 (Fold Change)") +
      ggplot2::geom_hline(yintercept = 0, colour = "black") +
      ggplot2::scale_color_manual(values = c(
        Treated = dot.col[1], Untreated = dot.col[2], none = dot.col[3]
      ))
    if (isTRUE(limit.force)) {
      myPLOT <- myPLOT + ggplot2::scale_y_continuous(limits = c(-y.limit, y.limit))
    }
    if (isTRUE(scale.ax)) {
      myPLOT <- myPLOT + ggplot2::scale_x_continuous(trans = "log1p")
    } else if (isTRUE(limit.force)) {
      myPLOT <- myPLOT + ggplot2::scale_x_continuous(limits = c(0, x.limit))
    }
  }

  if (isTRUE(interactive)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive = TRUE.")
    }
    pw <- plotly::ggplotly(myPLOT, tooltip = c("text", "x", "y"))
    if (!is.null(out.name) && nzchar(out.name)) {
      htmlwidgets::saveWidget(pw, paste0(out.name, ".html"), selfcontained = TRUE)
    }
    return(pw)
  }
  myPLOT
}
