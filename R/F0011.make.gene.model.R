#' Build a gene model for clustering / dimensionality reduction
#'
#' Selects highly variable / informative features from
#' \code{gene.stats()} output. Default gate is \strong{classic iCellR}
#' (cell-count  x  dispersion via \code{base.mean.rank})  -  that gene set is
#' what PCA uses. An alternate \strong{mean  x  dispersion} view is always
#' drawn alongside for exploration (HTML + static exports).
#'
#' @param x An \code{iCellR2} object with gene stats available
#'   (run \code{\link{gene.stats}} first).
#' @param modality Modality. Default \code{"rna"}.
#' @param gate.mode Gate used to build the PCA gene model:
#'   \code{"classic"} (default) or \code{"mean"}.
#' @param dispersion.limit Dispersion cutoff. Default \code{1.5}.
#' @param mean.limit Minimum mean expression for the mean gate / mean plot.
#'   Default \code{NULL} uses \eqn{2^{0.5}-1} (\eqn{\approx 0.414}).
#' @param base.mean.rank Rank used to set the classic cell-count cut.
#'   Default \code{1000}.
#' @param gene.num.max Maximum genes in the model. Default \code{2000}.
#' @param no.mito.model Exclude mitochondrial genes. Default \code{TRUE}.
#' @param no.cell.cycle Exclude S / G2M cell-cycle genes. Default \code{TRUE}.
#' @param mark.mito Highlight / label mito genes on plots. Default \code{TRUE}.
#' @param require.all.groups If multiple conditions exist in gene stats,
#'   keep genes selected in all groups. Default \code{TRUE}.
#' @param my.out.put \code{"data"}, \code{"plot"}, or \code{"both"}.
#' @param interactive If \code{TRUE} and \code{format} is unset, write HTML.
#' @param out.name Output path prefix when writing plots. Default \code{"plot"}.
#' @param format \code{NULL} or character vector: \code{"html"}, \code{"pdf"},
#'   \code{"png"}, \code{"svg"}. HTML opens on the PCA gate first; static
#'   exports are a two-panel figure (classic | mean).
#' @param width,height Device size in inches \emph{per panel} for PDF/PNG/SVG.
#'   Dual export uses \code{width * 2}  x  \code{height}. Default \code{6}.
#' @param dpi DPI for PNG. Default \code{300}.
#' @param downsample.plot Downsample points for static plots. Default \code{TRUE}.
#' @param max.points Max points on static plots. Default \code{500000}.
#' @param selected.col,other.col,mito.col,disp.line.col,mean.line.col,rank.line.col
#'   Plot colors.
#' @param cell.size,cell.transparency Point aesthetics.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} (and/or a ggplot), depending on
#'   \code{my.out.put}. Written file paths are attached as attribute
#'   \code{files} when \code{format} is set.
#'
#' @examples
#' \dontrun{
#' obj <- gene.stats(obj, assay.name = "data")
#' obj <- make.gene.model(obj, dispersion.limit = 1.5, base.mean.rank = 1000)
#' obj <- make.gene.model(
#'   obj, out.name = "plot2_GM", format = c("html", "png")
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme ggtitle xlab ylab
#'   geom_vline geom_hline scale_color_manual scale_y_continuous
#' @export
make.gene.model <- function(
    x,
    modality = "rna",
    gate.mode = c("classic", "mean"),
    dispersion.limit = 1.5,
    mean.limit = NULL,
    base.mean.rank = 1000,
    gene.num.max = 2000,
    no.mito.model = TRUE,
    no.cell.cycle = TRUE,
    mark.mito = TRUE,
    require.all.groups = TRUE,
    my.out.put = c("data", "plot", "both"),
    interactive = FALSE,
    out.name = "plot",
    format = NULL,
    width = 6,
    height = 6,
    dpi = 300,
    downsample.plot = TRUE,
    max.points = 500000,
    selected.col = "#0F766E",
    other.col = "#94A3B8",
    mito.col = "#BE123C",
    disp.line.col = "#0F766E",
    mean.line.col = "#0E7490",
    rank.line.col = "#BE123C",
    cell.size = 1.75,
    cell.transparency = 0.55,
    verbose = TRUE
) {
  .ic2_require_object(x)
  my.out.put <- match.arg(my.out.put)
  gate.mode <- match.arg(gate.mode)
  modality <- tolower(modality)

  write_formats <- character(0)
  if (!is.null(format)) {
    write_formats <- unique(tolower(as.character(format)))
    write_formats <- write_formats[nzchar(write_formats)]
  } else if (isTRUE(interactive)) {
    write_formats <- "html"
  }
  ok_fmt <- c("html", "pdf", "png", "svg")
  bad_fmt <- setdiff(write_formats, ok_fmt)
  if (length(bad_fmt)) {
    stop(
      "Unknown format(s): ", paste(bad_fmt, collapse = ", "),
      ". Use: ", paste(ok_fmt, collapse = ", ")
    )
  }
  need_plot <- my.out.put %in% c("plot", "both") || length(write_formats) > 0L

  rep <- .ic2_new_reporter("make.gene.model", verbose, soda_root = x@storage$root_dir)

  DATA <- x@modalities[[modality]]$gene.data
  if (is.null(DATA) || !is.data.frame(DATA) || nrow(DATA) < 1L) {
    sp <- x@storage$metadata_paths$gene_stats
    cand <- c(
      sp,
      file.path(x@storage$root_dir, "metadata", paste0(modality, "_gene_stats.tsv")),
      file.path(x@storage$root_dir, "metadata", paste0(modality, "_gene_stats.rds"))
    )
    cand <- unique(cand[!vapply(cand, function(z) {
      is.null(z) || !nzchar(as.character(z)[1L])
    }, logical(1))])
    for (p in cand) {
      if (!file.exists(p)) next
      DATA <- tryCatch(.ic2_soda_read_table(p), error = function(e) NULL)
      if (is.null(DATA) && grepl("\\.rds$", p, ignore.case = TRUE)) {
        DATA <- tryCatch(readRDS(p), error = function(e) NULL)
      }
      if (is.data.frame(DATA) && nrow(DATA) > 0L) {
        x@modalities[[modality]]$gene.data <- DATA
        break
      }
    }
  }
  if (is.null(DATA) || !is.data.frame(DATA) || nrow(DATA) < 1L) {
    stop("No gene stats found. Run gene.stats() first.")
  }

  needed <- c("genes", "numberOfCells", "meanExp", "SDs", "condition")
  if (!all(needed %in% colnames(DATA))) {
    stop("gene.data is missing required columns: ", paste(needed, collapse = ", "))
  }

  rep <- .ic2_banner_start(rep, c(
    paste0("Modality: ", modality, "   |   gate: ", gate.mode),
    paste0(
      "dispersion.limit: ", dispersion.limit,
      " | base.mean.rank: ", base.mean.rank,
      " | gene.num.max: ", gene.num.max
    )
  ))

  .ic2_load_cc <- function() {
    s.file <- system.file("extdata", "S.txt", package = "iCellR2")
    g.file <- system.file("extdata", "G2M.txt", package = "iCellR2")
    if (!nzchar(s.file) || !file.exists(s.file)) {
      s.file <- system.file("extdata", "s.phase.txt", package = "iCellR2")
    }
    if (!nzchar(g.file) || !file.exists(g.file)) {
      g.file <- system.file("extdata", "g2m.phase.txt", package = "iCellR2")
    }
    s.genes <- if (nzchar(s.file) && file.exists(s.file)) {
      scan(s.file, what = character(), quiet = TRUE)
    } else character(0)
    g.genes <- if (nzchar(g.file) && file.exists(g.file)) {
      scan(g.file, what = character(), quiet = TRUE)
    } else character(0)
    unique(c(s.genes, g.genes))
  }

  .ic2_select_from_table <- function(data, gate.mode, mean.limit) {
    me <- as.numeric(data$meanExp)
    nc <- as.numeric(data$numberOfCells)
    sd <- as.numeric(data$SDs)
    ord <- order(me, decreasing = TRUE, na.last = TRUE)
    k <- min(max(1L, as.integer(base.mean.rank)), length(ord))
    cellCountLimit <- nc[ord[k]]
    if (!is.finite(cellCountLimit)) cellCountLimit <- Inf
    top.rank.line <- log2(pmax(cellCountLimit, 1))

    if (is.null(mean.limit) || !is.finite(as.numeric(mean.limit)[1L])) {
      ## Design default: log2(mean+1) = 0.5   ->   mean = 2^0.5 - 1
      mean.limit <- 2^0.5 - 1
    } else {
      mean.limit <- as.numeric(mean.limit)[1L]
    }

    mito.genes <- grep("^mt[-.]", data$genes, value = TRUE, ignore.case = TRUE)
    cc.genes <- .ic2_load_cc()
    Cell.Cycle <- character(0)
    if (length(cc.genes) > 0L) {
      Cell.Cycle <- data$genes[toupper(data$genes) %in% toupper(cc.genes)]
    }

    data$color <- "other"
    if (identical(gate.mode, "mean")) {
      pass <- is.finite(sd) & is.finite(me) &
        sd >= dispersion.limit & me >= mean.limit
      data$color[pass] <- "selected"
    } else {
      ## Classic iCellR: high detection OR (high dispersion & low detection)
      data$color[nc > cellCountLimit] <- "selected"
      data$color[sd > dispersion.limit & nc < cellCountLimit] <- "selected"
    }

    sel <- as.character(data$genes[data$color == "selected"])
    if (isTRUE(no.mito.model)) sel <- setdiff(sel, mito.genes)
    if (isTRUE(no.cell.cycle)) sel <- setdiff(sel, Cell.Cycle)

    ## Rank remaining by dispersion and cap
    if (length(sel) > 0L) {
      hit <- match(sel, data$genes)
      o <- order(sd[hit], decreasing = TRUE, na.last = TRUE)
      sel <- utils::head(sel[o], as.integer(gene.num.max))
    }

    list(
      data = data,
      selected = unique(sel),
      mito.genes = mito.genes,
      mean.limit = mean.limit,
      cellCountLimit = cellCountLimit,
      top.rank.line = top.rank.line,
      SDlimit = dispersion.limit
    )
  }

  .ic2_step(rep, paste0("[1/2] Selecting genes (", gate.mode, " gate)"))
  conds <- unique(as.character(DATA$condition))
  data_all <- DATA[DATA$condition == "all", , drop = FALSE]
  if (nrow(data_all) < 1L) {
    stop("gene.data has no condition == 'all' rows.")
  }

  sel_all <- .ic2_select_from_table(data_all, gate.mode, mean.limit)
  mean.limit <- sel_all$mean.limit
  my.clust.genes <- sel_all$selected

  other_conds <- setdiff(conds, "all")
  if (length(other_conds) > 0L && isTRUE(require.all.groups)) {
    gene_lists <- lapply(other_conds, function(cn) {
      d <- DATA[DATA$condition == cn, , drop = FALSE]
      .ic2_select_from_table(d, gate.mode, mean.limit)$selected
    })
    inter <- Reduce(intersect, gene_lists)
    if (length(inter) > 0L) {
      my.clust.genes <- inter
      .ic2_step(rep, paste0(
        "  require.all.groups: ", .ic2_fmt_num(length(my.clust.genes)),
        " genes shared across ", length(other_conds), " groups"
      ))
    }
  }

  my.clust.genes <- utils::head(unique(my.clust.genes), as.integer(gene.num.max))
  .ic2_step(rep, paste0(
    "  Gene model size: ", .ic2_fmt_num(length(my.clust.genes)),
    if (identical(gate.mode, "mean")) {
      paste0("   |   mean.limit=", signif(mean.limit, 4))
    } else {
      paste0("   |   cell-count cut=", signif(sel_all$cellCountLimit, 4))
    }
  ))

  ## Mark final selected on plot table
  pdata_full <- sel_all$data
  pdata_full$color <- "other"
  pdata_full$color[pdata_full$genes %in% my.clust.genes] <- "selected"
  if (isTRUE(mark.mito) && length(sel_all$mito.genes)) {
    pdata_full$color[pdata_full$genes %in% sel_all$mito.genes &
                       pdata_full$color != "selected"] <- "mito"
  }

  plot_obj <- NULL
  written <- character(0)
  if (isTRUE(need_plot)) {
    .ic2_step(rep, "[2/2] Building classic + mean-dispersion plots")
    pdata <- pdata_full
    if (isTRUE(downsample.plot) && nrow(pdata) > max.points) {
      set.seed(1)
      keep <- unique(c(
        which(pdata$genes %in% my.clust.genes),
        which(pdata$genes %in% sel_all$mito.genes),
        sample.int(nrow(pdata), as.integer(max.points), replace = FALSE)
      ))
      pdata <- pdata[keep, , drop = FALSE]
    }

    pca_tag <- "  |  used for PCA"
    p_classic <- .ic2_gene_model_panel(
      pdata = pdata,
      view = "classic",
      dispersion.limit = dispersion.limit,
      mean.limit = mean.limit,
      cellCountLimit = sel_all$cellCountLimit,
      top.rank.line = sel_all$top.rank.line,
      base.mean.rank = base.mean.rank,
      gene.num.max = gene.num.max,
      title_suffix = if (identical(gate.mode, "classic")) pca_tag else "",
      selected.col = selected.col,
      other.col = other.col,
      mito.col = mito.col,
      disp.line.col = disp.line.col,
      mean.line.col = mean.line.col,
      rank.line.col = rank.line.col,
      cell.size = cell.size,
      cell.transparency = cell.transparency,
      mark.mito = mark.mito,
      mito.genes = sel_all$mito.genes
    )
    p_mean <- .ic2_gene_model_panel(
      pdata = pdata,
      view = "mean",
      dispersion.limit = dispersion.limit,
      mean.limit = mean.limit,
      cellCountLimit = sel_all$cellCountLimit,
      top.rank.line = sel_all$top.rank.line,
      base.mean.rank = base.mean.rank,
      gene.num.max = gene.num.max,
      title_suffix = if (identical(gate.mode, "mean")) pca_tag else "",
      selected.col = selected.col,
      other.col = other.col,
      mito.col = mito.col,
      disp.line.col = disp.line.col,
      mean.line.col = mean.line.col,
      rank.line.col = rank.line.col,
      cell.size = cell.size,
      cell.transparency = cell.transparency,
      mark.mito = mark.mito,
      mito.genes = sel_all$mito.genes
    )

    ## Primary panel = PCA gate; combined figure for export / device print
    primary <- if (identical(gate.mode, "mean")) p_mean else p_classic
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      ## Classic left (PCA default), mean right
      dual <- gridExtra::arrangeGrob(
        p_classic, p_mean,
        ncol = 2L,
        top = grid::textGrob(
          paste0(
            "Gene model   |   PCA gate: ", gate.mode,
            "   |   ", .ic2_fmt_num(length(my.clust.genes)), " genes"
          ),
          gp = grid::gpar(fontsize = 12, fontface = "bold", col = "#14242A"),
          x = 0.01, hjust = 0
        )
      )
      plot_obj <- primary
      attr(plot_obj, "dual") <- dual
      attr(plot_obj, "classic") <- p_classic
      attr(plot_obj, "mean") <- p_mean
    } else {
      plot_obj <- primary
      attr(plot_obj, "classic") <- p_classic
      attr(plot_obj, "mean") <- p_mean
      dual <- NULL
    }

    if (length(write_formats) > 0L) {
      out_stem <- .ic2_plot_out_stem(out.name, default_stem = "gene_model")
      static_formats <- intersect(write_formats, c("pdf", "png", "svg"))
      for (fmt in static_formats) {
        fpath <- paste0(out_stem, ".", fmt)
        grob_out <- if (!is.null(dual)) dual else primary
        w_out <- if (!is.null(dual)) width * 2 else width
        if (identical(fmt, "pdf")) {
          ggplot2::ggsave(
            fpath, grob_out, width = w_out, height = height,
            device = grDevices::pdf
          )
        } else if (identical(fmt, "svg")) {
          ggplot2::ggsave(
            fpath, grob_out, width = w_out, height = height,
            device = grDevices::svg
          )
        } else {
          ggplot2::ggsave(
            fpath, grob_out, width = w_out, height = height,
            dpi = dpi, device = "png"
          )
        }
        written <- c(written, stats::setNames(fpath, fmt))
        .ic2_step(rep, paste0("  Wrote: ", basename(fpath)))
      }
      if ("html" %in% write_formats) {
        fhtml <- paste0(out_stem, ".html")
        ok <- .ic2_gene_model_write_html(
          data = pdata_full,
          file = fhtml,
          gate.mode = gate.mode,
          dispersion.limit = dispersion.limit,
          mean.limit = mean.limit,
          base.mean.rank = base.mean.rank,
          cell.count.limit = sel_all$cellCountLimit,
          gene.num.max = gene.num.max,
          no.mito.model = no.mito.model,
          mark.mito = mark.mito,
          selected = my.clust.genes,
          n_all = nrow(data_all)
        )
        if (isTRUE(ok)) {
          written <- c(written, html = fhtml)
          .ic2_step(rep, paste0("  Wrote: ", basename(fhtml)))
        }
      }
    }
  } else {
    .ic2_step(rep, "[2/2] Skipping plot (my.out.put = 'data', no format)")
  }

  if (my.out.put %in% c("data", "both")) {
    x@modalities[[modality]]$gene.model <- my.clust.genes
    x@misc$gene.model <- my.clust.genes

    model_stem <- file.path(
      x@storage$root_dir, "metadata", paste0(modality, "_gene_model")
    )
    dir.create(dirname(model_stem), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      data.frame(feature_id = my.clust.genes, stringsAsFactors = FALSE),
      file = paste0(model_stem, ".tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )
    model_path <- paste0(model_stem, ".tsv")
    .ic2_soda_unlink_rds(paste0(model_stem, ".rds"))
    if (!is.null(x@storage$metadata_paths)) {
      x@storage$metadata_paths$gene_model <- model_path
    }

    x@params$feature_selection[[modality]]$gene.model <- list(
      n_genes = length(my.clust.genes),
      gate.mode = gate.mode,
      dispersion.limit = dispersion.limit,
      mean.limit = mean.limit,
      base.mean.rank = base.mean.rank,
      gene.num.max = gene.num.max,
      no.mito.model = no.mito.model,
      no.cell.cycle = no.cell.cycle,
      path = model_path
    )

    x <- .ic2_log_command(
      x, "make.gene.model",
      arguments = list(
        modality = modality,
        gate.mode = gate.mode,
        dispersion.limit = dispersion.limit,
        mean.limit = mean.limit,
        base.mean.rank = base.mean.rank,
        gene.num.max = gene.num.max,
        format = if (length(write_formats)) write_formats else NULL,
        out.name = out.name
      ),
      outputs = list(
        n_genes = length(my.clust.genes),
        path = model_path,
        plot_files = if (length(written)) unname(written) else NULL
      )
    )

    if (!is.null(x@storage$root_dir)) {
      rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
    }
    rep <- .ic2_add_output(rep, "Gene model", model_path)
    for (wf in written) rep <- .ic2_add_output(rep, "Plot", wf)
    .ic2_finish(rep, extra_lines = paste0(
      "Selected: ", .ic2_fmt_num(length(my.clust.genes)), " genes",
      if (length(written)) {
        paste0("   |   wrote: ", paste(basename(written), collapse = ", "))
      } else ""
    ))

    if (identical(my.out.put, "both")) {
      out <- list(object = x, plot = plot_obj)
      if (length(written)) attr(out, "files") <- written
      return(out)
    }
    if (length(written)) attr(x, "files") <- written
    return(x)
  }

  .ic2_finish(rep, extra_lines = if (length(written)) {
    paste0("Wrote: ", paste(basename(written), collapse = ", "))
  } else NULL)
  if (length(written)) attr(plot_obj, "files") <- written
  plot_obj
}

#' One gene-model panel (classic or mean view).
#' @keywords internal
.ic2_gene_model_panel <- function(
    pdata,
    view = c("classic", "mean"),
    dispersion.limit,
    mean.limit,
    cellCountLimit,
    top.rank.line,
    base.mean.rank,
    gene.num.max,
    title_suffix = "",
    selected.col,
    other.col,
    mito.col,
    disp.line.col,
    mean.line.col,
    rank.line.col,
    cell.size,
    cell.transparency,
    mark.mito,
    mito.genes
) {
  view <- match.arg(view)
  df <- pdata
  if (identical(view, "mean")) {
    df$xplot <- log2(pmax(as.numeric(df$meanExp), 0) + 1)
    x_lab <- "log2(mean expression + 1)"
    x_cut <- log2(pmax(mean.limit, 0) + 1)
    x_cut_col <- mean.line.col
    plot_title <- paste0("Mean  x  dispersion", title_suffix)
    plot_sub <- paste0(
      "dispersion >= ", dispersion.limit,
      "   |   mean >= ", signif(mean.limit, 3),
      "   |   max ", gene.num.max
    )
  } else {
    df$xplot <- log2(pmax(as.numeric(df$numberOfCells), 0) + 1)
    x_lab <- "log2(n cells + 1)"
    x_cut <- top.rank.line
    x_cut_col <- rank.line.col
    plot_title <- paste0("Classic iCellR", title_suffix)
    plot_sub <- paste0(
      "dispersion >= ", dispersion.limit,
      "   |   base.mean.rank ", base.mean.rank,
      "   |   cell-count cut ", round(cellCountLimit)
    )
  }
  df$yplot <- as.numeric(df$SDs)
  col_map <- c(selected = selected.col, other = other.col, mito = mito.col)
  df$color <- factor(as.character(df$color), levels = c("other", "mito", "selected"))
  df <- df[order(df$color), , drop = FALSE]

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$xplot,
      y = .data$yplot,
      text = .data$genes,
      colour = .data$color
    )
  ) +
    ggplot2::geom_point(size = cell.size, alpha = cell.transparency, stroke = 0) +
    ggplot2::scale_color_manual(values = col_map, drop = FALSE) +
    ggplot2::scale_y_continuous(trans = "log1p") +
    ggplot2::geom_hline(
      yintercept = dispersion.limit, colour = disp.line.col,
      linetype = "dotted", linewidth = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = x_cut, colour = x_cut_col,
      linetype = "dotted", linewidth = 0.8
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", colour = "#14242A", size = 12),
      plot.subtitle = ggplot2::element_text(colour = "#5A6B72", size = 9),
      axis.title = ggplot2::element_text(colour = "#5A6B72", size = 10),
      axis.text = ggplot2::element_text(colour = "#5A6B72", size = 8),
      axis.line = ggplot2::element_line(colour = "#14242A", linewidth = 0.55)
    ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_sub,
      x = x_lab,
      y = "log1p(dispersion)"
    )

  if (isTRUE(mark.mito) && length(mito.genes) > 0L &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    top_labelled <- df[df$genes %in% mito.genes, , drop = FALSE]
    if (nrow(top_labelled) > 0L) {
      p <- p + ggrepel::geom_text_repel(
        data = top_labelled,
        mapping = ggplot2::aes(label = .data$genes),
        size = 2.4,
        fontface = "bold",
        color = mito.col,
        inherit.aes = TRUE,
        max.overlaps = 25
      )
    }
  }
  p
}

#' Normalize a plot output stem (create parent dir; strip extension).
#' @keywords internal
.ic2_plot_out_stem <- function(out.name, default_stem = "plot") {
  out.name <- path.expand(as.character(out.name)[1L])
  if (!nzchar(out.name)) out.name <- default_stem
  if (dir.exists(out.name) || grepl("[/\\\\]$", out.name)) {
    dir.create(out.name, recursive = TRUE, showWarnings = FALSE)
    return(file.path(sub("[/\\\\]$", "", out.name), default_stem))
  }
  parent <- dirname(out.name)
  if (!identical(parent, ".") && !dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  sub("\\.(html|pdf|png|svg)$", "", out.name, ignore.case = TRUE)
}

#' Save an htmlwidget; selfcontained when pandoc is available.
#' @return \code{TRUE} on success, \code{FALSE} otherwise.
#' @keywords internal
.ic2_save_htmlwidget <- function(widget, file, selfcontained = NULL) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    warning("Package 'htmlwidgets' is required to write HTML; skipping.")
    return(FALSE)
  }
  file <- path.expand(as.character(file)[1L])
  parent <- dirname(file)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  has_pandoc <- nzchar(Sys.which("pandoc"))
  if (!has_pandoc && requireNamespace("rmarkdown", quietly = TRUE)) {
    has_pandoc <- isTRUE(tryCatch(
      rmarkdown::pandoc_available(),
      error = function(e) FALSE
    ))
  }
  if (is.null(selfcontained)) selfcontained <- has_pandoc
  if (isTRUE(selfcontained) && !has_pandoc) {
    warning(
      "Pandoc not found; writing non-selfcontained HTML (",
      basename(file), "). Install pandoc for a single portable file."
    )
    selfcontained <- FALSE
  }
  ok <- tryCatch({
    htmlwidgets::saveWidget(widget, file, selfcontained = selfcontained)
    TRUE
  }, error = function(e) {
    if (isTRUE(selfcontained)) {
      warning(
        "selfcontained HTML failed (", conditionMessage(e),
        "); retrying with selfcontained = FALSE."
      )
      tryCatch({
        htmlwidgets::saveWidget(widget, file, selfcontained = FALSE)
        TRUE
      }, error = function(e2) {
        warning("Failed to write HTML: ", conditionMessage(e2))
        FALSE
      })
    } else {
      warning("Failed to write HTML: ", conditionMessage(e))
      FALSE
    }
  })
  isTRUE(ok)
}
