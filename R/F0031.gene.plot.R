#' Plot gene expression (feature plot / violin)
#'
#' Slim plotting helper in the spirit of \code{\link{cluster.plot}}:
#' \itemize{
#'   \item \code{plot.type = "scatter"}  -  color cells on UMAP/PCA/KNetL/... by expression
#'   \item \code{plot.type = "violin"}  -  violin by \code{cluster_id}
#'   \item \code{plot.type = "boxplot"}  -  soft violin with a narrow white box
#'   \item \code{plot.type = "barplot"}  -  mean bars with SD whiskers by cluster
#' }
#' Genes are streamed from the on-disk assay for the plotted cell sample
#' only (\code{max.cells}); the full matrix is never densified.
#'
#' @param x An \code{iCellR2} object.
#' @param gene Gene id(s). Multiple genes  ->  faceted scatter / faceted violin.
#' @param modality Modality. Default \code{"rna"}.
#' @param data.type Which matrix to plot. Default \code{"main.data"}
#'   (normalized). Choose from \code{"raw.data"}, \code{"main.data"},
#'   \code{"imputed.data"}, \code{"scaled.data"} (also accepts short
#'   aliases \code{"raw"}/\code{"counts"}, \code{"main"}/\code{"data"},
#'   \code{"imputed"}, \code{"scale.data"}, plus ATAC
#'   \code{"atac"} / \code{"atac.imputed"}). Ignored if
#'   \code{assay.name} is set.
#' @param assay.name Assay name to read (overrides \code{data.type}).
#'   Internal SODA names: \code{"counts"}, \code{"data"},
#'   \code{"imputed"}, \code{"scale.data"}.
#' @param transform \code{"auto"} (default): \code{log2p1} when assay is
#'   not \code{"scale.data"}, else \code{"none"}. Or set explicitly.
#' @param scale.rows If \code{TRUE} (default), z-score each gene across
#'   plotted cells after \code{transform}. Used for scatter feature plots
#'   only. Violin, box, and bar keep the transformed expression
#'   (\code{log2(x+1)} by default) so the axis does not go below zero.
#'   Ignored when assay is already \code{"scale.data"}.
#' @param reduction Reduction for scatter. Default \code{"auto"} follows
#'   \code{active_reduction} / prefers \code{aligned_umap} when present
#'   (same rule as \code{\link{cluster.plot}}). Explicit:
#'   \code{"umap"}, \code{"aligned_umap"}, \code{"pca"}, \code{"aligned_pca"},
#'   \code{"knetl"}, \code{"aligned_knetl"}, \code{"cgmap"}, \code{"cgpca"},
#'   \code{"cgumap"}, \code{"ptumap"}, \code{"tsne"}, \code{"aligned_tsne"}.
#'   Ignored for violin / box / bar.
#' @param dims Dimensions for scatter (length 2). Default \code{1:2}.
#' @param plot.type \code{"scatter"} (default), \code{"violin"},
#'   \code{"boxplot"} (jitter + violin + box + mean), or \code{"barplot"}
#'   (mean +/- SD whiskers). Aliases: \code{"box"}, \code{"bar"}.
#'   Heatmap-style plots are separate helpers and are unchanged.
#' @param cell.size Point size (scatter / boxplot jitter). Default \code{0.2}
#'   for scatter; boxplot jitter uses the same argument (try \code{0.4}).
#' @param cell.color Jitter point color for \code{boxplot}. Default
#'   \code{"black"}.
#' @param jitter.width Horizontal jitter width for \code{boxplot}.
#'   Default \code{0.15}.
#' @param cell.transparency Point / layer alpha. Default \code{0.7}.
#' @param box.color Box fill for \code{boxplot}. Default \code{"red"}.
#' @param box.line.col Box outline for \code{boxplot}. Default \code{"green"}.
#' @param box.transparency Box/violin fill alpha. Default \code{0.5}.
#' @param notch Notch boxplots. Default \code{FALSE}.
#' @param back.col Panel background. Default \code{"white"}.
#' @param max.cells Max cells to draw / sample for gene fetch.
#'   Default \code{500000}.
#' @param scale.color Clip expression color scale to quantiles
#'   \code{c(low, high)} in \eqn{[0,1]} (scatter only). Default
#'   \code{c(0.01, 0.99)}. Set \code{NULL} to disable.
#' @param colors Two-color gradient for scatter (low, high).
#'   Default \code{c("grey", "red")}.
#' @param interactive If \code{TRUE}, return plotly (scatter only).
#' @param out.name Optional HTML stem when interactive.
#' @param chunk.size Cells per HDF5 read when fetching genes.
#' @param seed Downsample seed.
#' @param verbose Progress style.
#'
#' @return A ggplot or plotly widget.
#'
#' @examples
#' \dontrun{
#' gene.plot(obj, gene = "CD3D")                          # active / aligned UMAP
#' gene.plot(obj, gene = "CD3D", data.type = "raw.data")
#' gene.plot(obj, gene = "CD3D", data.type = "imputed.data")
#' gene.plot(obj, gene = c("CD3D", "MS4A1"), reduction = "umap")
#' gene.plot(obj, gene = "CD3D", reduction = "knetl")
#' gene.plot(obj, gene = "CD3D", plot.type = "violin")
#' gene.plot(obj, gene = "CD3D", plot.type = "boxplot")
#' gene.plot(obj, gene = "CD3D", plot.type = "barplot")
#' }
#'
#' @importFrom rlang .data
#' @export
gene.plot <- function(
    x,
    gene = NULL,
    modality = "rna",
    data.type = "main.data",
    assay.name = NULL,
    transform = c("auto", "none", "log2p1", "log1p"),
    scale.rows = TRUE,
    reduction = "auto",
    dims = 1:2,
    plot.type = c("scatter", "violin", "boxplot", "barplot"),
    cell.size = 0.2,
    cell.color = "black",
    jitter.width = 0.15,
    cell.transparency = 0.7,
    box.color = "red",
    box.line.col = "green",
    box.transparency = 0.5,
    notch = FALSE,
    back.col = "white",
    max.cells = 500000L,
    scale.color = c(0.01, 0.99),
    colors = c("grey", "red"),
    interactive = FALSE,
    out.name = NULL,
    chunk.size = 200000L,
    seed = 1,
    verbose = TRUE
) {
  .ic2_require_object(x)
  transform <- match.arg(transform)
  modality <- tolower(modality)
  dims <- as.integer(dims)
  chunk.size <- as.integer(chunk.size)

  ## Normalize aliases (box / bar / scatterplot)
  pt <- tolower(as.character(plot.type)[1L])
  pt <- switch(
    pt,
    scatterplot = "scatter",
    box = "boxplot",
    bar = "barplot",
    pt
  )
  plot.type <- match.arg(pt, c("scatter", "violin", "boxplot", "barplot"))

  if (is.null(gene) || !length(gene) || !any(nzchar(as.character(gene)))) {
    stop("Provide at least one gene id in 'gene'.")
  }
  gene <- unique(as.character(gene))
  gene <- gene[nzchar(gene)]

  if (identical(plot.type, "scatter") && length(dims) != 2L) {
    stop("Scatter gene.plot requires dims of length 2.")
  }
  if (isTRUE(interactive) && !identical(plot.type, "scatter")) {
    warning(
      "interactive = TRUE is only supported for plot.type = \"scatter\"; ",
      "using static ggplot."
    )
    interactive <- FALSE
  }

  ## Same first-class reductions as cluster.plot (incl. knetl / aligned_*).
  red_choices <- c(
    "auto", "umap", "aligned_umap", "pca", "aligned_pca",
    "tsne", "aligned_tsne", "cgmap", "cgpca", "cgumap", "ptumap",
    "knetl", "aligned_knetl"
  )
  red_arg <- tolower(as.character(reduction %||% "auto")[1L])
  if (identical(plot.type, "scatter")) {
    if (identical(red_arg, "auto") || !nzchar(red_arg)) {
      reduction <- .ic2_default_scatter_reduction(
        x, modality,
        allowed = setdiff(red_choices, "auto")
      )
    } else {
      reduction <- match.arg(red_arg, red_choices)
      if (identical(reduction, "auto")) {
        reduction <- .ic2_default_scatter_reduction(
          x, modality,
          allowed = setdiff(red_choices, "auto")
        )
      }
    }
  } else if (!identical(red_arg, "auto") && nzchar(red_arg)) {
    reduction <- match.arg(red_arg, red_choices)
  }

  rep <- .ic2_new_reporter("gene.plot", verbose, soda_root = x@storage$root_dir)

  ## Assay (data.type / assay.name)
  assay.name <- .ic2_resolve_plot_assay(
    x, modality, assay.name = assay.name, data.type = data.type, prefer.data = TRUE
  )
  if (identical(transform, "auto")) {
    transform <- if (identical(assay.name, "scale.data")) "none" else "log2p1"
  }

  .ic2_step(rep, paste0("[1/3] Preparing cells + ", plot.type))

  meta <- x@cell_data
  cell_ids <- rownames(meta)
  if (length(cell_ids) < 2L) stop("cell_data has fewer than 2 cells.")

  emb <- NULL
  if (identical(plot.type, "scatter")) {
    emb <- .ic2_get_plot_embedding(x, modality, reduction)
    common <- intersect(cell_ids, rownames(emb))
    if (length(common) < 2L) {
      stop("Embedding rownames do not match cell_data.")
    }
    if (max(dims) > ncol(emb)) {
      stop("Requested dims exceed available ", reduction, " columns.")
    }
    emb <- emb[common, dims, drop = FALSE]
    colnames(emb) <- c("Dim1", "Dim2")
    meta <- meta[common, , drop = FALSE]
    cell_ids <- common
  }

  ## Downsample before gene fetch
  keep <- seq_along(cell_ids)
  n_all <- length(cell_ids)
  if (!is.null(max.cells) && is.finite(max.cells) && n_all > as.integer(max.cells)) {
    set.seed(seed)
    keep <- sort(sample.int(n_all, as.integer(max.cells)))
    .ic2_step(rep, paste0(
      "  Downsampled ", .ic2_fmt_num(length(keep)),
      " / ", .ic2_fmt_num(n_all), " cells"
    ))
  }
  cell_ids <- cell_ids[keep]
  meta <- meta[cell_ids, , drop = FALSE]
  if (!is.null(emb)) emb <- emb[cell_ids, , drop = FALSE]

  .ic2_step(rep, paste0(
    "[2/3] Fetching ", length(gene), " gene(s) from assay '", assay.name, "'",
    " (transform=", transform, ")"
  ))
  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  size_factors <- .ic2_assay_size_factors(x, modality, assay.name)
  expr <- .ic2_load_features_for_cells(
    handles = handles,
    features = gene,
    cell_ids = cell_ids,
    chunk.size = chunk.size,
    transform = transform,
    size_factors = size_factors
  )
  found <- intersect(gene, colnames(expr))
  missing <- setdiff(gene, found)
  if (length(missing)) {
    warning(
      "Gene(s) not found in assay: ", paste(missing, collapse = ", ")
    )
  }
  if (length(found) < 1L) stop("None of the requested genes were found.")
  expr <- expr[, found, drop = FALSE]
  gene <- found

  ## Z-score is for scatter contrast only. Violin / box / bar stay on
  ## log2(x+1) (or the chosen transform), which is >= 0.
  do_z <- isTRUE(scale.rows) && !identical(assay.name, "scale.data") &&
    identical(plot.type, "scatter")
  if (do_z && ncol(expr) > 0L && nrow(expr) > 1L) {
    expr <- scale(expr)
    expr[!is.finite(expr)] <- 0
    expr[expr > 2.5] <- 2.5
    expr[expr < -2.5] <- -2.5
  }

  rep <- .ic2_banner_start(rep, c(
    paste0("Genes: ", paste(gene, collapse = ", ")),
    paste0(
      "Assay: ", assay.name, "   |   transform: ", transform,
      if (do_z) " + z" else "",
      "   |   plot: ", plot.type,
      if (identical(plot.type, "scatter")) paste0("   |   ", reduction) else ""
    ),
    paste0("Cells plotted: ", .ic2_fmt_num(length(cell_ids)))
  ))

  .ic2_step(rep, "[3/3] Building plot")

  if (identical(plot.type, "scatter")) {
    ## Long format for facet when multi-gene (preserve request order)
    dfs <- lapply(gene, function(g) {
      data.frame(
        Dim1 = emb[, 1L],
        Dim2 = emb[, 2L],
        cell = cell_ids,
        gene = g,
        expression = as.numeric(expr[, g]),
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, dfs)
    df$gene <- factor(df$gene, levels = gene)
    color_name <- if (do_z) "z" else "expr"
    lims <- NULL
    if (do_z) {
      ## Per-gene positive-z rescale to [0,1] so sparse markers still saturate
      df$expression <- ave(df$expression, df$gene, FUN = function(v) {
        v <- pmax(as.numeric(v), 0)
        hi <- as.numeric(stats::quantile(v, probs = 0.98, na.rm = TRUE, names = FALSE))
        hi <- max(hi, 0.75, na.rm = TRUE)
        pmin(v / hi, 1)
      })
      lims <- c(0, 1)
      color_name <- "z"
    } else if (!is.null(scale.color) && length(scale.color) == 2L) {
      qs <- as.numeric(stats::quantile(
        df$expression, probs = scale.color, na.rm = TRUE, names = FALSE
      ))
      if (all(is.finite(qs)) && qs[1] < qs[2]) {
        lims <- qs
        df$expression <- pmin(pmax(df$expression, lims[1]), lims[2])
      }
    }
    ink <- "#2C2C2C"
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$Dim1, y = .data$Dim2, colour = .data$expression, text = .data$cell)
    ) +
      ggplot2::geom_point(size = cell.size, alpha = cell.transparency) +
      ggplot2::scale_colour_gradientn(
        colours = colors,
        name = color_name,
        limits = lims
      ) +
      ggplot2::labs(
        title = if (length(gene) == 1L) gene else "Gene expression",
        x = paste0(reduction, dims[1]),
        y = paste0(reduction, dims[2])
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = back.col, colour = NA),
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = ink, linewidth = 0.45, lineend = "square"),
        axis.ticks = ggplot2::element_line(colour = ink, linewidth = 0.35),
        axis.title = ggplot2::element_text(colour = ink, size = 13),
        axis.text = ggplot2::element_text(colour = "#6B6B6B", size = 11),
        plot.title = ggplot2::element_text(
          colour = ink, face = "bold", size = 14, hjust = 0,
          margin = ggplot2::margin(b = 8)
        ),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(colour = ink, face = "bold", size = 12),
        panel.spacing = grid::unit(0.7, "lines"),
        legend.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank()
      )
    if (length(gene) > 1L) {
      p <- p + .ic2_facet_wrap(~gene)
    }
  } else {
    if (!"cluster_id" %in% colnames(meta)) {
      stop(plot.type, " needs 'cluster_id' in cell_data. Run run.clustering() first.")
    }
    dfs <- lapply(gene, function(g) {
      data.frame(
        cluster = .ic2_factor_sorted(meta$cluster_id),
        gene = g,
        expression = as.numeric(expr[, g]),
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, dfs)
    df$gene <- factor(df$gene, levels = gene)
    ylab.txt <- switch(
      transform,
      log2p1 = "log2 expression",
      log1p = "log1p expression",
      paste0("expression (", assay.name, ")")
    )
    title.txt <- if (length(gene) == 1L) gene else "Gene expression"
    ink <- "#2C2C2C"
    theme_dist <- ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = back.col, colour = NA),
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = ink, linewidth = 0.45, lineend = "square"),
        axis.ticks = ggplot2::element_line(colour = ink, linewidth = 0.35),
        axis.title = ggplot2::element_text(colour = ink, size = 13),
        axis.text = ggplot2::element_text(colour = "#6B6B6B", size = 11),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = ggplot2::element_text(
          colour = ink, face = "bold", size = 14, hjust = 0,
          margin = ggplot2::margin(b = 8)
        ),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(colour = ink, face = "bold", size = 12),
        panel.spacing = grid::unit(0.7, "lines"),
        legend.position = "none"
      )

    if (identical(plot.type, "violin")) {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$cluster, y = .data$expression, fill = .data$cluster)
      ) +
        ggplot2::geom_violin(scale = "width", trim = TRUE, colour = "grey40", linewidth = 0.25) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0.16)),
          limits = c(0, NA)
        ) +
        ggplot2::labs(title = title.txt, x = "cluster", y = ylab.txt) +
        theme_dist
    } else if (identical(plot.type, "boxplot")) {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$cluster, y = .data$expression, fill = .data$cluster)
      ) +
        ggplot2::geom_violin(
          scale = "width", trim = TRUE, colour = NA, alpha = 0.28, linewidth = 0
        ) +
        ggplot2::geom_boxplot(
          width = 0.16, outlier.shape = NA, fill = "white",
          colour = "#2C2C2C", linewidth = 0.35, notch = isTRUE(notch)
        ) +
        ggplot2::stat_summary(
          fun = mean, geom = "point", size = 1.3, colour = "#2C2C2C"
        ) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0.12)),
          limits = c(0, NA)
        ) +
        ggplot2::labs(title = title.txt, x = "cluster", y = ylab.txt) +
        theme_dist +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_line(colour = "#F0F0F0", linewidth = 0.3)
        )
    } else {
      ## barplot: mean +/- SD whiskers (classic gene.plot)
      keys <- interaction(df$gene, df$cluster, drop = TRUE, lex.order = TRUE)
      summ <- do.call(rbind, lapply(split(df, keys), function(d) {
        data.frame(
          cluster = d$cluster[1L],
          gene = d$gene[1L],
          mean = mean(d$expression, na.rm = TRUE),
          sd = stats::sd(d$expression, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }))
      summ$gene <- factor(summ$gene, levels = gene)
      summ$cluster <- factor(summ$cluster, levels = levels(df$cluster))
      summ$sd[!is.finite(summ$sd)] <- 0
      p <- ggplot2::ggplot(
        summ,
        ggplot2::aes(x = .data$cluster, y = .data$mean, fill = .data$cluster)
      ) +
        ggplot2::geom_col(
          colour = "black", alpha = cell.transparency, width = 0.7
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = .data$mean,
            ymax = .data$mean + .data$sd
          ),
          width = 0.2
        ) +
        ggplot2::labs(
          title = title.txt, x = "cluster",
          y = paste0("mean +/- SD (", assay.name, ")")
        ) +
        theme_dist
    }
    if (length(gene) > 1L) {
      p <- p + .ic2_facet_wrap(~gene, scales = "free_y")
    }
  }

  if (isTRUE(interactive)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive = TRUE.")
    }
    pw <- plotly::ggplotly(p, tooltip = c("text", "colour", "x", "y"))
    if (!is.null(out.name) && nzchar(out.name)) {
      htmlwidgets::saveWidget(pw, paste0(out.name, ".html"), selfcontained = TRUE)
    }
    .ic2_finish(rep, extra_lines = "Returned interactive plotly widget")
    return(pw)
  }

  .ic2_finish(rep, extra_lines = paste0(
    "ggplot   |   ", plot.type, "   |   genes: ", length(gene),
    "   |   cells: ", .ic2_fmt_num(length(cell_ids))
  ))
  invisible(p)
}

#' Load selected features for selected cells from assay handles
#'
#' Returns a dense cells  x  features matrix. Missing genes in a sample
#' are filled with zeros. Values are optionally transformed.
#'
#' @keywords internal
.ic2_load_features_for_cells <- function(
    handles,
    features,
    cell_ids,
    chunk.size = 200000L,
    transform = c("none", "log2p1", "log1p"),
    size_factors = NULL
) {
  transform <- match.arg(transform)
  features <- as.character(features)
  cell_ids <- as.character(cell_ids)
  n_cells <- length(cell_ids)
  n_feat <- length(features)
  if (n_cells < 1L || n_feat < 1L) {
    stop("Need at least one cell and one feature.")
  }

  out <- matrix(0, nrow = n_cells, ncol = n_feat,
                dimnames = list(cell_ids, features))
  remaining <- cell_ids

  for (h in handles) {
    if (length(remaining) < 1L) break
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    local <- intersect(remaining, h$cell_ids)
    if (length(local) < 1L) next
    row_idx <- match(features, h$feature_ids)
    present <- which(!is.na(row_idx))
    if (length(present) < 1L) {
      remaining <- setdiff(remaining, local)
      next
    }
    col_idx <- match(local, h$cell_ids)
    ord <- order(col_idx)
    local <- local[ord]
    col_idx <- col_idx[ord]
    runs <- .ic2_contig_runs(col_idx)
    for (run in runs) {
      starts <- seq.int(1L, length(run), by = as.integer(chunk.size))
      for (s in starts) {
        e <- min(s + as.integer(chunk.size) - 1L, length(run))
        cols <- run[s:e]
        raw <- .ic2_read_cols(h, cols)
        cn <- colnames(raw)
        if (is.null(cn)) cn <- h$cell_ids[cols]
        block <- matrix(0, nrow = length(features), ncol = length(cols),
                        dimnames = list(features, cn))
        sub <- raw[row_idx[present], , drop = FALSE]
        if (inherits(sub, "Matrix")) {
          sub <- as.matrix(.ic2_to_dgCMatrix(sub))
        } else {
          sub <- as.matrix(sub)
        }
        storage.mode(sub) <- "double"
        block[present, ] <- sub
        block[!is.finite(block)] <- 0
        if (!is.null(size_factors)) {
          block <- .ic2_scale_cols_by_sf(block, size_factors)
        }
        block <- .ic2_assay_transform_block(block, transform = transform)
        ## cells  x  features into out
        out[cn, ] <- t(block)
      }
    }
    remaining <- setdiff(remaining, local)
  }
  if (length(remaining) > 0L) {
    warning(
      length(remaining), " cells had no assay data for gene fetch ",
      "(left as 0)."
    )
  }
  out
}
