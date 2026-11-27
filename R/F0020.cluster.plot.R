#' Plot cells on a dimensionality reduction
#'
#' Visualize PCA / UMAP / KNetL (and related reductions) stored in the SODA archive.
#' For large objects, points are **downsampled** (or drawn as hex bins) so
#' plotting stays interactive and memory-light  -  embeddings themselves are
#' never fully densified beyond the plot sample.
#'
#' @param x An \code{iCellR2} object.
#' @param modality Modality. Default \code{"rna"}.
#' @param plot.type Reduction to plot. Default \code{"auto"} follows
#'   \code{active_reduction} (e.g. \code{aligned_umap} after \code{iba()} +
#'   \code{\link{run.umap}}), else prefers \code{aligned_umap} / \code{umap}.
#'   Explicit choices: \code{"umap"}, \code{"aligned_umap"}, \code{"pca"},
#'   \code{"aligned_pca"}, \code{"cgpca"} (\code{\link{run.cgpca}}),
#'   \code{"cgumap"} (\code{\link{run.cgumap}}), \code{"ptumap"}
#'   (\code{\link{run.ptumap}}), \code{"knetl"}, \code{"aligned_knetl"},
#'   \code{"tsne"}, \code{"aligned_tsne"}, or legacy \code{"cgmap"}.
#' @param dims Which dimensions to plot (length 2 or 3). Default \code{1:2}.
#' @param col.by Color mapping: \code{"clusters"} (\code{cluster_id}),
#'   \code{"conditions"} (\code{condition_id}), \code{"monochrome"},
#'   \code{"cc"} (phase from QC), or any column name in \code{x@cell_data}
#'   (e.g. \code{"mito.percent"}, \code{"nGenes"}, \code{"branch_id"},
#'   \code{"Pseudotime"}). Integer ID columns such as \code{branch_id} use
#'   discrete categorical colors; continuous metrics use a heatmap scale.
#' @param cell.size Point size. Default \code{1}.
#' @param cell.transparency Alpha in \code{[0,1]}. Default \code{1}
#'   (full opacity, like classic iCellR).
#' @param cell.color Color when \code{col.by = "monochrome"}.
#' @param back.col Plot background fill (ggplot panel). Default \code{"white"}.
#' @param max.cells Maximum cells to draw. \code{NULL} = all (can be slow).
#'   Default \code{500000}.
#' @param hex If \code{TRUE}, use hexbin density instead of points
#'   (great for millions of cells). Default \code{FALSE}.
#' @param hex.bins Number of hex bins when \code{hex = TRUE}.
#' @param anno.clust Label cluster centroids when coloring by clusters.
#' @param anno.size Annotation text size.
#' @param cond.facet Facet by \code{condition_id} when present.
#' @param interactive If \code{TRUE}, write a Teal-style interactive HTML explorer
#'   (colour toggles, PNG/SVG download). For \code{pca} / \code{cgpca} the HTML
#'   includes both 2D and 3D when a third dimension exists. Default \code{FALSE}
#'   keeps the classic static ggplot.
#' @param out.name Output path prefix when \code{interactive = TRUE}
#'   (writes \code{paste0(out.name, ".html")}). Default \code{NULL} uses
#'   \code{cluster_plot_<plot.type>}.
#' @param seed RNG seed for downsampling.
#' @param verbose Progress style.
#'
#' @return A \code{ggplot} when \code{interactive = FALSE}. When
#'   \code{interactive = TRUE}, invisibly the HTML file path.
#'
#' @examples
#' \dontrun{
#' cluster.plot(obj)                                 # active / aligned UMAP
#' cluster.plot(obj, plot.type = "umap")             # force unaligned umap
#' cluster.plot(obj, plot.type = "knetl")            # KNetL map (run.knetl)
#' cluster.plot(obj, plot.type = "pca", dims = 1:2)
#' cluster.plot(obj, plot.type = "cgmap")
#' cluster.plot(obj, plot.type = "cgpca", dims = 1:2)
#' cluster.plot(obj, plot.type = "cgumap")
#' cluster.plot(obj, plot.type = "ptumap", col.by = "branch_id")
#' cluster.plot(obj, col.by = "mito.percent", hex = TRUE)
#' cluster.plot(obj, interactive = TRUE, out.name = "umap_qc")
#' cluster.plot(obj, plot.type = "cgpca", interactive = TRUE, out.name = "cgpca_qc")
#' }
#'
#' @importFrom rlang .data
#' @importFrom htmlwidgets saveWidget
#' @export
cluster.plot <- function(
    x,
    modality = "rna",
    plot.type = "auto",
    dims = 1:2,
    col.by = "clusters",
    cell.size = 1,
    cell.transparency = 1,
    cell.color = "grey20",
    back.col = "white",
    max.cells = 500000L,
    hex = FALSE,
    hex.bins = 80L,
    anno.clust = TRUE,
    anno.size = 4,
    cond.facet = FALSE,
    interactive = FALSE,
    out.name = NULL,
    seed = 1,
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  plot_choices <- c(
    "auto", "umap", "aligned_umap", "pca", "aligned_pca",
    "tsne", "aligned_tsne", "cgmap", "cgpca", "cgumap", "ptumap",
    "knetl", "aligned_knetl"
  )
  pt <- tolower(as.character(plot.type %||% "auto")[1L])
  if (identical(pt, "auto") || !nzchar(pt)) {
    plot.type <- .ic2_default_scatter_reduction(
      x, modality,
      allowed = setdiff(plot_choices, "auto")
    )
  } else {
    plot.type <- match.arg(pt, plot_choices)
    if (identical(plot.type, "auto")) {
      plot.type <- .ic2_default_scatter_reduction(
        x, modality,
        allowed = setdiff(plot_choices, "auto")
      )
    }
  }
  dims <- as.integer(dims)
  if (!length(dims) %in% c(2L, 3L)) {
    stop("'dims' must have length 2 or 3.")
  }
  if (isTRUE(hex) && length(dims) != 2L) {
    stop("'hex = TRUE' requires 2D dims.")
  }
  if (isTRUE(hex) && isTRUE(interactive)) {
    warning("hex plots are static ggplot; setting interactive = FALSE.")
    interactive <- FALSE
  }

  rep <- .ic2_new_reporter("cluster.plot", verbose, soda_root = x@storage$root_dir)

  ## ---- Load embedding from SODA ----
  .ic2_step(rep, paste0("[1/3] Loading ", plot.type, " embeddings"))
  emb_full <- .ic2_get_plot_embedding(x, modality, plot.type)
  n_all <- nrow(emb_full)
  if (max(dims) > ncol(emb_full)) {
    stop(
      "Requested dims up to ", max(dims),
      " but ", plot.type, " only has ", ncol(emb_full), " columns."
    )
  }

  ## Align metadata on full embedding
  meta <- x@cell_data
  common <- intersect(rownames(emb_full), rownames(meta))
  if (length(common) < 2L) {
    stop("Embedding rownames do not match cell_data.")
  }
  emb_full <- emb_full[common, , drop = FALSE]
  meta <- meta[common, , drop = FALSE]

  ## Friendly defaults before clustering exists
  if (identical(col.by, "clusters") && !"cluster_id" %in% colnames(meta)) {
    if ("condition_id" %in% colnames(meta)) {
      .ic2_step(rep, "  No cluster_id yet  -  coloring by condition_id")
      col.by <- "conditions"
    } else {
      .ic2_step(rep, "  No cluster_id yet  -  monochrome")
      col.by <- "monochrome"
    }
  }

  ## Downsample cell indices once
  .ic2_step(rep, "[2/3] Preparing plot sample")
  keep <- seq_len(nrow(emb_full))
  if (!is.null(max.cells) && is.finite(max.cells) && nrow(emb_full) > as.integer(max.cells)) {
    set.seed(seed)
    keep <- sort(sample.int(nrow(emb_full), as.integer(max.cells)))
    .ic2_step(rep, paste0(
      "  Downsampled ", .ic2_fmt_num(length(keep)),
      " / ", .ic2_fmt_num(n_all), " cells for plotting"
    ))
  }
  emb_s <- emb_full[keep, , drop = FALSE]
  meta_s <- meta[keep, , drop = FALSE]
  rep <- .ic2_progress_bar(rep, 1L, 1L, prefix = "  sample")

  ## ---- Interactive Teal HTML (does not change static ggplot look) ----
  if (isTRUE(interactive)) {
    dims2 <- dims[seq_len(min(2L, length(dims)))]
    if (length(dims2) < 2L) dims2 <- c(1L, 2L)
    emb2 <- emb_s[, dims2, drop = FALSE]

    emb3 <- NULL
    dims3 <- NULL
    ## Only PCA / CGPCA get a 3D panel (they are the reductions with a 3rd dim).
    if (plot.type %in% c("pca", "cgpca") && ncol(emb_s) >= 3L) {
      dims3 <- unique(c(dims2, 3L))
      if (length(dims3) < 3L) dims3 <- 1:3
      dims3 <- dims3[seq_len(3L)]
      if (max(dims3) <= ncol(emb_s)) {
        emb3 <- emb_s[, dims3, drop = FALSE]
      }
    }

    .ic2_step(rep, "[3/3] Writing interactive HTML")
    fhtml <- .ic2_cluster_plot_html(
      emb2 = emb2,
      emb3 = emb3,
      meta = meta_s,
      plot.type = plot.type,
      dims2 = dims2,
      dims3 = dims3,
      col.by = col.by,
      out.name = out.name,
      anno.clust = anno.clust
    )
    .ic2_finish(rep, extra_lines = paste0("Wrote interactive HTML: ", fhtml))
    return(invisible(fhtml))
  }

  ## ---- Static path (unchanged look) ----
  emb <- emb_s[, dims, drop = FALSE]
  colnames(emb) <- paste0("Dim", seq_along(dims))
  color_info <- .ic2_plot_color_vector(meta_s, col.by = col.by)
  color_p <- color_info$values
  color_discrete <- color_info$discrete
  color_title <- color_info$title
  meta_p <- meta_s

  title <- toupper(plot.type)

  rep <- .ic2_banner_start(rep, c(
    paste0("Reduction: ", plot.type, "   |   dims: ", paste(dims, collapse = ",")),
    paste0("Color: ", color_title, "   |   hex=", hex, "   |   interactive=", interactive)
  ))

  .ic2_step(rep, "[3/3] Building plot")
  df <- as.data.frame(emb, stringsAsFactors = FALSE)
  df$cell <- rownames(emb)
  df$.color <- color_p
  if ("condition_id" %in% colnames(meta_p)) {
    df$condition_id <- as.character(meta_p$condition_id)
  }

  ## ---- 3D plotly (static request with dims length 3 only) ----
  if (length(dims) == 3L) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for 3D cluster.plot().")
    }
    p <- plotly::plot_ly(
      df,
      x = ~Dim1, y = ~Dim2, z = ~Dim3,
      color = ~.color,
      type = "scatter3d",
      mode = "markers",
      text = ~cell,
      marker = list(size = cell.size + 1.5, opacity = cell.transparency)
    )
    p <- plotly::layout(
      p,
      title = title,
      scene = list(
        xaxis = list(title = paste0(plot.type, dims[1])),
        yaxis = list(title = paste0(plot.type, dims[2])),
        zaxis = list(title = paste0(plot.type, dims[3]))
      )
    )
    if (!is.null(out.name) && nzchar(out.name)) {
      htmlwidgets::saveWidget(p, paste0(out.name, ".html"), selfcontained = TRUE)
    }
    .ic2_finish(rep, extra_lines = "Returned 3D plotly widget")
    return(p)
  }

  ## ---- 2D ggplot ----
  xlab <- paste0(plot.type, dims[1])
  ylab <- paste0(plot.type, dims[2])

  p <- NULL
  use_hex_plot <- FALSE
  if (isTRUE(hex)) {
    if (identical(col.by, "monochrome") || !isTRUE(color_discrete)) {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2)) +
        ggplot2::geom_hex(bins = as.integer(hex.bins)) +
        ggplot2::scale_fill_viridis_c(name = "count", option = "C")
      use_hex_plot <- TRUE
    } else {
      ## per-group hex is noisy; fall back to points
      .ic2_step(rep, "  hex + discrete color: using points (hex is for density)")
    }
  }

  if (is.null(p)) {
    ## stroke = 0 with very small sizes vanishes on PNG/PDF devices; keep a
    ## hairline stroke so faceted plots with cell.size ~ 0.1 still show cells.
    pt_stroke <- if (isTRUE(as.numeric(cell.size) < 0.5)) 0.15 else 0
    if (identical(col.by, "monochrome")) {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2, text = .data$cell)) +
        ggplot2::geom_point(
          size = cell.size, alpha = cell.transparency, colour = cell.color,
          stroke = pt_stroke, shape = 16
        )
    } else if (isTRUE(color_discrete)) {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$Dim1, y = .data$Dim2, colour = .data$.color, text = .data$cell)
      ) +
        ggplot2::geom_point(
          size = cell.size, alpha = cell.transparency,
          stroke = pt_stroke, shape = 16
        ) +
        ggplot2::scale_colour_discrete(name = color_title) +
        ## Large legend dots, no boxes (classic iCellR look)
        ggplot2::guides(colour = ggplot2::guide_legend(
          override.aes = list(size = 5.5, alpha = 1, stroke = 0, shape = 16),
          title.position = "top",
          keywidth = 0.9,
          keyheight = 0.9
        ))
    } else {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$Dim1, y = .data$Dim2, colour = .data$.color, text = .data$cell)
      ) +
        ggplot2::geom_point(
          size = cell.size, alpha = cell.transparency,
          stroke = pt_stroke, shape = 16
        ) +
        ggplot2::scale_colour_viridis_c(name = color_title, option = "C") +
        ggplot2::guides(colour = ggplot2::guide_colorbar(
          barwidth = 0.55,
          barheight = 6,
          ticks = FALSE,
          frame.colour = NA,
          title.position = "top"
        ))
    }
  } else if (isTRUE(use_hex_plot)) {
    p <- p + ggplot2::guides(fill = ggplot2::guide_colorbar(
      barwidth = 0.55,
      barheight = 6,
      ticks = FALSE,
      frame.colour = NA,
      title.position = "top"
    ))
  }

  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    .ic2_theme_cluster_plot(back.col = back.col)

  if (isTRUE(cond.facet) && "condition_id" %in% colnames(df)) {
    p <- p + .ic2_facet_wrap(~condition_id)
  }

  if (isTRUE(anno.clust) && identical(col.by, "clusters") &&
      isTRUE(color_discrete) && !isTRUE(use_hex_plot)) {
    ## geom_text, not annotate(): faceting repeats annotate data once per
    ## panel and then the label vector no longer matches the layer rows.
    by_cond <- isTRUE(cond.facet) && "condition_id" %in% colnames(df)
    if (by_cond) {
      cents <- stats::aggregate(
        cbind(Dim1, Dim2) ~ .color + condition_id,
        data = df, FUN = mean, na.rm = TRUE
      )
    } else {
      cents <- stats::aggregate(
        cbind(Dim1, Dim2) ~ .color, data = df, FUN = mean, na.rm = TRUE
      )
    }
    p <- p + ggplot2::geom_text(
      data = cents,
      mapping = ggplot2::aes(
        x = .data$Dim1, y = .data$Dim2, label = .data$.color
      ),
      size = anno.size,
      fontface = "bold",
      colour = "black",
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  .ic2_finish(rep, extra_lines = paste0(
    "ggplot   |   plotted ", .ic2_fmt_num(nrow(emb)), " cells"
  ))
  invisible(p)
}

#' Facet wrap with an x/y axis line on every panel when ggplot2 supports it
#'
#' @keywords internal
.ic2_facet_wrap <- function(...) {
  dots <- list(...)
  has_axes <- utils::packageVersion("ggplot2") >= "3.5.0"
  if (isTRUE(has_axes)) {
    if (is.null(dots$axes)) dots$axes <- "all"
  } else {
    dots$axes <- NULL
  }
  keep <- !vapply(dots, is.null, logical(1))
  do.call(ggplot2::facet_wrap, dots[keep])
}

#' Modern open-axis theme for cluster / embedding plots
#'
#' Left + bottom axis lines only (no panel box), blank legend keys so
#' colour swatches show as large points without boxes  -  closer to classic
#' iCellR, with a cleaner contemporary look.
#'
#' @keywords internal
.ic2_theme_cluster_plot <- function(back.col = "white", base_size = 13) {
  ink <- "#2C2C2C"
  mute <- "#6B6B6B"
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = back.col, colour = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = ink, linewidth = 0.45, lineend = "square"),
      axis.ticks = ggplot2::element_line(colour = ink, linewidth = 0.35),
      axis.ticks.length = grid::unit(3.5, "pt"),
      axis.title = ggplot2::element_text(colour = ink, size = base_size),
      axis.text = ggplot2::element_text(colour = mute, size = base_size - 2),
      plot.title = ggplot2::element_text(
        colour = ink, face = "bold", size = base_size + 1, hjust = 0,
        margin = ggplot2::margin(b = 8)
      ),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = grid::unit(1.1, "lines"),
      legend.title = ggplot2::element_text(colour = ink, size = base_size - 1, face = "plain"),
      legend.text = ggplot2::element_text(colour = mute, size = base_size - 1),
      legend.spacing.y = grid::unit(2, "pt"),
      legend.margin = ggplot2::margin(l = 4),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = ink, face = "bold", size = base_size - 1)
    )
}

#' Load reduction embedding matrix for plotting
#' @keywords internal
.ic2_get_plot_embedding <- function(x, modality, plot.type) {
  red <- x@modalities[[modality]]$reductions[[plot.type]]
  if (is.null(red)) {
    hint <- switch(
      plot.type,
      umap = "run.umap()",
      aligned_umap = "run.umap() after iba()",
      pca = "run.pca()",
      aligned_pca = "iba() / run.pca()",
      knetl = "run.knetl()",
      aligned_knetl = "run.knetl() after iba()",
      cgmap = "run.cgumap() / run.cgpca()",
      cgpca = "run.cgpca()",
      cgumap = "run.cgumap()",
      ptumap = "run.ptumap()",
      paste0("run.", plot.type, "()")
    )
    stop(
      "No '", plot.type, "' reduction for modality '", modality,
      "'. Run ", hint, " first."
    )
  }
  path <- red$cell_embeddings
  if (is.null(path) && !is.null(red$meta_path) && file.exists(red$meta_path)) {
    meta <- .ic2_soda_read_meta(red$meta_path)
    path <- meta$paths$embeddings %||% meta$cell_embeddings
  }
  if (is.null(path) || !file.exists(path)) {
    stop("Embedding file for '", plot.type, "' not found on disk.")
  }
  emb <- if (grepl("\\.(h5|hdf5)$", path, ignore.case = TRUE)) {
    .ic2_h5_read_dense(path)
  } else {
    readRDS(path)
  }
  if (!is.matrix(emb)) emb <- as.matrix(emb)
  emb
}

#' Resolve a color vector from cell metadata
#' @keywords internal
.ic2_plot_color_vector <- function(meta, col.by = "clusters") {
  if (identical(col.by, "monochrome")) {
    return(list(values = rep(NA_character_, nrow(meta)), discrete = TRUE, title = ""))
  }
  if (identical(col.by, "clusters")) {
    if (!"cluster_id" %in% colnames(meta)) {
      warning("No 'cluster_id' in cell_data  -  coloring all cells as 0.")
      v <- factor(rep("0", nrow(meta)))
    } else {
      v <- .ic2_factor_sorted(meta$cluster_id)
    }
    return(list(values = v, discrete = TRUE, title = "cluster"))
  }
  if (identical(col.by, "conditions")) {
    if (!"condition_id" %in% colnames(meta)) {
      stop("No 'condition_id' in cell_data.")
    }
    v <- .ic2_factor_sorted(meta$condition_id)
    return(list(values = v, discrete = TRUE, title = "condition"))
  }
  if (identical(col.by, "cc")) {
    ## Prefer explicit Phase; else crude S/G2M score labels from qc.stats
    if ("Phase" %in% colnames(meta)) {
      labs <- unique(as.character(meta$Phase))
      pref <- c("G1", "S", "G2M")
      lev <- c(intersect(pref, labs), setdiff(labs, pref))
      v <- factor(as.character(meta$Phase), levels = lev)
      return(list(values = v, discrete = TRUE, title = "phase"))
    }
    if (all(c("S.phase.probability", "g2m.phase.probability") %in% colnames(meta))) {
      s <- as.numeric(meta$S.phase.probability)
      g <- as.numeric(meta$g2m.phase.probability)
      lab <- ifelse(s >= g & s > 0.2, "S",
                    ifelse(g > s & g > 0.2, "G2M", "G1"))
      return(list(
        values = factor(lab, levels = c("G1", "S", "G2M")),
        discrete = TRUE,
        title = "phase"
      ))
    }
    stop("Cell-cycle columns not found. Run qc.stats() first.")
  }
  if (!col.by %in% colnames(meta)) {
    stop(
      "'col.by' = '", col.by, "' not found in cell_data. ",
      "Available: ", paste(utils::head(colnames(meta), 20L), collapse = ", "), "..."
    )
  }
  v <- meta[[col.by]]
  if (.ic2_plot_col_is_discrete(v, col.by)) {
    list(values = .ic2_factor_sorted(v), discrete = TRUE, title = col.by)
  } else {
    list(values = as.numeric(v), discrete = FALSE, title = col.by)
  }
}

#' Decide discrete vs continuous coloring for a cell_data column
#'
#' Factors / characters are discrete. Known continuous scores
#' (\code{Pseudotime}, QC metrics) stay continuous. Integer-like IDs with
#' few unique values (e.g. \code{branch_id}) use categorical colors.
#'
#' @keywords internal
.ic2_plot_col_is_discrete <- function(v, col.by) {
  if (is.factor(v) || is.character(v) || is.logical(v)) return(TRUE)
  ## Explicit continuous columns (even if few unique values in a demo)
  cont_known <- c(
    "Pseudotime", "mito.percent", "ribo.percent", "nGenes", "nUMI", "UMIs",
    "S.phase.probability", "g2m.phase.probability"
  )
  if (col.by %in% cont_known || grepl("^branch_prob_", col.by)) {
    return(FALSE)
  }
  ## Explicit categorical ID columns
  if (identical(col.by, "branch_id") ||
      identical(col.by, "cluster_id") ||
      identical(col.by, "condition_id") ||
      identical(col.by, "sample_id") ||
      grepl("_id$", col.by)) {
    return(TRUE)
  }
  if (!(is.numeric(v) || is.integer(v))) return(TRUE)
  vv <- v[is.finite(as.numeric(v))]
  if (length(vv) < 1L) return(TRUE)
  n_uniq <- length(unique(vv))
  ## Small-cardinality integer-like IDs  ->  discrete palette like clusters
  is_intish <- all(abs(vv - round(vv)) < 1e-8)
  if (is_intish && n_uniq <= 50L && n_uniq <= max(20L, as.integer(0.2 * length(vv)))) {
    return(TRUE)
  }
  FALSE
}

#' Factor with numeric-aware level order (1,2,...,10 not 1,10,2)
#' @keywords internal
.ic2_factor_sorted <- function(x) {
  labs <- as.character(x)
  uniq <- unique(labs[!is.na(labs)])
  if (length(uniq) < 1L) {
    return(factor(labs))
  }
  num <- suppressWarnings(as.numeric(uniq))
  if (all(is.finite(num))) {
    lev <- uniq[order(num)]
  } else {
    ## Mixed / non-numeric: numeric-looking first in numeric order, then alpha
    is_num <- is.finite(num)
    lev <- c(uniq[is_num][order(num[is_num])], sort(uniq[!is_num]))
  }
  factor(labs, levels = lev)
}
