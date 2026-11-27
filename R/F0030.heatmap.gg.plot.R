#' Heatmap of genes (classic iCellR look + scalable means)
#'
#' Default matches classic iCellR: **genes  x  cells**, ordered / faceted by
#' cluster (\code{group.by = "cells"}, \code{cluster.by = "clusters"}).
#' For large objects, use \code{group.by = "clusters"} or
#' \code{"conditions"} for streamed group means (one column per group).
#'
#' @param x An \code{iCellR2} object.
#' @param gene Character vector of gene ids to plot.
#' @param group.by \code{"cells"} (default; classic per-cell tiles),
#'   \code{"clusters"} (means by \code{cluster_id}), or
#'   \code{"conditions"} (means by \code{condition_id}).
#' @param cluster.by Classic iCellR argument for cell-level heatmaps:
#'   \code{"clusters"} (default; facet / order by cluster),
#'   \code{"conditions"} (facet by condition), or \code{"none"}
#'   (no facet). Ignored when \code{group.by} is means mode.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Assay to use; overrides \code{data.type} when set.
#' @param data.type Which matrix to plot. Default \code{"main.data"}
#'   (normalized \code{"data"}). Also accepts \code{"raw.data"},
#'   \code{"imputed.data"}, \code{"scaled.data"}, and short aliases
#'   (\code{"main"}, \code{"imputed"}, \code{"scale.data"}, ...).
#' @param transform \code{"auto"} (default): \code{log2p1} unless assay is
#'   \code{scale.data}; else \code{"none"}, \code{"log2p1"}, \code{"log1p"}.
#' @param conds.to.plot Optional condition ids to include (vector OK).
#' @param clusts.to.plot Optional cluster ids to include (vector OK).
#' @param scale.rows Z-score each gene across columns. Default \code{TRUE}.
#' @param min.scale,max.scale Clip scaled values. Defaults \code{-2.5}, \code{2.5}.
#' @param heat.colors Length-3 color vector low/mid/high.
#'   Default \code{c("blue","white","red")}.
#' @param max.cells Max cells when \code{group.by = "cells"}. Default \code{500000}.
#' @param cell.sort Order cells by PC1 within each facet group
#'   (\code{group.by = "cells"} only). Default \code{FALSE} (classic code
#'   default).
#' @param chunk.size Cells per HDF5 read. Default \code{200000}.
#' @param cex.row Gene label size. Default \code{10}.
#' @param cex.col Cell / column label size (means mode). Default \code{10}.
#' @param no.key Hide color legend. Default \code{FALSE}.
#' @param interactive If \code{TRUE}, return plotly / write HTML.
#'   Default \code{FALSE} (scripts / SODA-friendly; classic iCellR defaulted
#'   to \code{TRUE}).
#' @param out.name If interactive and non-\code{NULL}, also write HTML.
#' @param seed RNG seed for cell downsampling.
#' @param verbose Progress style.
#'
#' @return A \code{ggplot} or plotly widget.
#'
#' @examples
#' \dontrun{
#' ## Classic look (default): genes  x  cells, faceted by cluster
#' heatmap.gg.plot(obj, gene = top.markers(markers, topde = 5))
#' ## Scalable summary: one column per cluster mean
#' heatmap.gg.plot(obj, gene = top.markers(markers, topde = 5),
#'                 group.by = "clusters")
#' heatmap.gg.plot(obj, gene = c("CD3D", "MS4A1"),
#'                 cluster.by = "conditions")
#' }
#'
#' @importFrom rlang .data
#' @export
heatmap.gg.plot <- function(
    x = NULL,
    gene = NULL,
    group.by = c("cells", "clusters", "conditions"),
    cluster.by = c("clusters", "conditions", "none"),
    modality = "rna",
    assay.name = NULL,
    data.type = "main.data",
    transform = c("auto", "none", "log2p1", "log1p"),
    conds.to.plot = NULL,
    clusts.to.plot = NULL,
    scale.rows = TRUE,
    min.scale = -2.5,
    max.scale = 2.5,
    heat.colors = c("blue", "white", "red"),
    max.cells = 500000L,
    cell.sort = FALSE,
    chunk.size = 200000L,
    cex.row = 10,
    cex.col = 10,
    no.key = FALSE,
    interactive = FALSE,
    out.name = NULL,
    seed = 1L,
    verbose = TRUE
) {
  .ic2_require_object(x)
  group.by <- match.arg(group.by)
  cluster.by <- match.arg(cluster.by)
  modality <- tolower(modality)
  transform <- match.arg(transform)
  chunk.size <- as.integer(chunk.size)
  max.cells <- as.integer(max.cells)

  if (is.null(gene) || length(gene) < 1L) {
    stop("Provide at least one gene in 'gene'.")
  }
  gene_req <- unique(as.character(gene))
  gene <- gene_req
  if (length(heat.colors) != 3L) {
    stop("'heat.colors' must be length 3 (low, mid, high).")
  }

  ## Classic data.type  ->  assay.name (assay.name wins when set)
  assay.name <- .ic2_resolve_plot_assay(
    x, modality, assay.name = assay.name, data.type = data.type, prefer.data = TRUE
  )
  if (identical(transform, "auto")) {
    transform <- if (identical(assay.name, "scale.data")) "none" else "log2p1"
  }

  rep <- .ic2_new_reporter("heatmap.gg.plot", verbose, soda_root = x@storage$root_dir)
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

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  size_factors <- .ic2_assay_size_factors(x, modality, assay.name)

  if (identical(group.by, "cells")) {
    .ic2_step(rep, "[1/2] Fetching genes for cells (classic layout)")
    cell_ids <- rownames(meta)
    n_all <- length(cell_ids)
    if (n_all > max.cells) {
      set.seed(as.integer(seed))
      if ("cluster_id" %in% colnames(meta)) {
        parts <- lapply(split(cell_ids, meta$cluster_id), function(ids) {
          ## Use doubles  -  max.cells * length(ids) overflows 32-bit int at atlas scale
          n_take <- max(
            1L,
            as.integer(round(
              as.numeric(max.cells) * as.numeric(length(ids)) / as.numeric(n_all)
            ))
          )
          n_take <- min(n_take, length(ids))
          if (length(ids) <= n_take) ids else sample(ids, n_take)
        })
        cell_ids <- unique(unlist(parts, use.names = FALSE))
        if (length(cell_ids) > max.cells) {
          cell_ids <- sample(cell_ids, max.cells)
        }
      } else {
        cell_ids <- cell_ids[sample.int(n_all, max.cells)]
      }
      meta <- meta[cell_ids, , drop = FALSE]
      .ic2_step(rep, paste0(
        "  Downsampled ", .ic2_fmt_num(n_all), "  ->  ",
        .ic2_fmt_num(length(cell_ids)), " cells"
      ))
    }

    expr <- .ic2_load_features_for_cells(
      handles = handles,
      features = gene,
      cell_ids = cell_ids,
      chunk.size = chunk.size,
      transform = transform,
      size_factors = size_factors
    )
    found <- intersect(gene_req, colnames(expr))
    missing <- setdiff(gene_req, found)
    if (length(missing)) {
      warning(
        "Gene(s) not available / not plotted: ",
        paste(utils::head(missing, 20L), collapse = ", "),
        if (length(missing) > 20L) ", ..." else ""
      )
    }
    gene <- found
    if (length(gene) < 1L) stop("None of the requested genes were found.")
    mat <- t(expr[, gene, drop = FALSE]) ## genes x cells

    if (identical(cluster.by, "clusters")) {
      if (!"cluster_id" %in% colnames(meta)) {
        stop("cluster.by = 'clusters' needs cluster_id. Run clustering first.")
      }
      ## Index by rownames  -  meta$cluster_id[cell_ids] is all-NA (unnamed vector)
      col_group <- as.character(meta[cell_ids, "cluster_id"])
      ord <- order(.ic2_factor_sorted(col_group))
      cell_ids <- cell_ids[ord]
      mat <- mat[, cell_ids, drop = FALSE]
      col_group <- col_group[ord]
    } else if (identical(cluster.by, "conditions")) {
      if (!"condition_id" %in% colnames(meta)) {
        stop("cluster.by = 'conditions' needs condition_id in cell_data.")
      }
      col_group <- as.character(meta[cell_ids, "condition_id"])
      ord <- order(.ic2_factor_sorted(col_group))
      cell_ids <- cell_ids[ord]
      mat <- mat[, cell_ids, drop = FALSE]
      col_group <- col_group[ord]
    } else {
      col_group <- NULL
    }

    if (isTRUE(cell.sort) && ncol(mat) > 2L && nrow(mat) > 1L) {
      groups <- if (is.null(col_group)) rep("all", length(cell_ids)) else col_group
      new_ord <- character(0)
      for (g in unique(groups)) {
        idx <- which(groups == g)
        sub <- mat[, idx, drop = FALSE]
        if (ncol(sub) < 3L || nrow(sub) < 2L) {
          new_ord <- c(new_ord, colnames(sub))
          next
        }
        ## sub is genes  x  cells. prcomp rows = observations  ->  cells must be rows.
        ## (Old code used prcomp(t(scale(t(sub)))) = genes as rows  ->  kept only
        ##  n_genes cells per cluster.)
        pc <- tryCatch(
          {
            zs <- t(scale(t(sub))) ## genes  x  cells, z-scored genes
            zs[!is.finite(zs)] <- 0
            stats::prcomp(t(zs), center = FALSE, scale. = FALSE)$x[, 1L]
          },
          error = function(e) seq_len(ncol(sub))
        )
        if (length(pc) != ncol(sub)) {
          new_ord <- c(new_ord, colnames(sub))
        } else {
          new_ord <- c(new_ord, colnames(sub)[order(pc)])
        }
      }
      mat <- mat[, new_ord, drop = FALSE]
      if (!is.null(col_group)) {
        col_group <- col_group[match(new_ord, cell_ids)]
      }
      cell_ids <- new_ord
    }

    col_labs <- cell_ids
    facet_var <- col_group
  } else {
    if (identical(group.by, "clusters")) {
      if (!"cluster_id" %in% colnames(meta)) {
        stop("No 'cluster_id' in cell_data. Run run.clustering() first.")
      }
      cell_to_grp <- setNames(
        as.character(meta$cluster_id),
        rownames(meta)
      )
    } else {
      if (!"condition_id" %in% colnames(meta)) {
        stop("No 'condition_id' in cell_data.")
      }
      cell_to_grp <- setNames(
        as.character(meta$condition_id),
        rownames(meta)
      )
    }
    .ic2_step(rep, paste0(
      "[1/2] Streaming means for ", length(gene), " gene(s) by ", group.by
    ))
    mat <- .ic2_feature_group_means(
      handles = handles,
      features = gene,
      cell_to_group = cell_to_grp,
      transform = transform,
      chunk.size = chunk.size,
      progress = rep,
      size_factors = size_factors
    )
    gene <- rownames(mat)
    if (length(gene) < 1L) stop("None of the requested genes were found.")
    cn <- colnames(mat)
    num <- suppressWarnings(as.numeric(cn))
    if (all(is.finite(num))) {
      mat <- mat[, cn[order(num)], drop = FALSE]
    } else {
      mat <- mat[, sort(cn), drop = FALSE]
    }
    col_labs <- colnames(mat)
    facet_var <- NULL
  }

  if (isTRUE(scale.rows) && ncol(mat) > 1L) {
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
  }
  mat[mat > max.scale] <- max.scale
  mat[mat < min.scale] <- min.scale

  mode_lab <- if (identical(group.by, "cells")) {
    paste0("cells   |   cluster.by=", cluster.by)
  } else {
    paste0("means   |   ", group.by)
  }
  rep <- .ic2_banner_start(rep, c(
    paste0("Genes: ", length(gene), "   |   ", mode_lab),
    paste0(
      "Assay: ", assay.name, "   |   transform: ", transform,
      "   |   scale.rows: ", scale.rows
    )
  ))
  .ic2_step(rep, "[2/2] Building heatmap")

  df <- data.frame(
    gene = factor(rep(rownames(mat), times = ncol(mat)), levels = rev(gene)),
    column = factor(
      rep(col_labs, each = nrow(mat)),
      levels = col_labs
    ),
    expression = as.numeric(mat),
    stringsAsFactors = FALSE
  )
  if (!is.null(facet_var) && length(unique(facet_var)) > 1L) {
    fl <- unique(facet_var)
    num <- suppressWarnings(as.numeric(fl))
    if (all(is.finite(num))) fl <- fl[order(num)] else fl <- sort(fl)
    df$facet <- factor(rep(facet_var, each = nrow(mat)), levels = fl)
  } else {
    facet_var <- NULL
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$column, y = .data$gene, fill = .data$expression)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = heat.colors[1], mid = heat.colors[2], high = heat.colors[3],
      midpoint = 0, name = if (scale.rows) "z" else "expr",
      guide = if (isTRUE(no.key)) {
        "none"
      } else {
        ggplot2::guide_colorbar(direction = "vertical", title.position = "left")
      }
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = cex.row),
      axis.text.x = if (identical(group.by, "cells")) {
        ggplot2::element_blank()
      } else {
        ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = cex.col)
      },
      axis.ticks.x = if (identical(group.by, "cells")) {
        ggplot2::element_blank()
      } else {
        ggplot2::element_line()
      },
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  if (!is.null(facet_var)) {
    p <- p +
      ggplot2::facet_wrap(~facet, nrow = 1L, scales = "free_x") +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = NA))
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
    "Heatmap: ", length(gene), " genes x ", length(col_labs),
    if (identical(group.by, "cells")) " cells" else paste0(" ", group.by)
  ))
  invisible(p)
}

#' Mean expression of selected features per group (streamed)
#'
#' @keywords internal
.ic2_feature_group_means <- function(
    handles,
    features,
    cell_to_group,
    transform = "log2p1",
    chunk.size = 200000L,
    progress = NULL,
    return.pct = FALSE,
    size_factors = NULL
) {
  features <- unique(as.character(features))
  groups <- sort(unique(as.character(unname(cell_to_group))))
  n_g <- length(groups)
  n_f <- length(features)
  if (n_g < 1L || n_f < 1L) {
    stop("Need features and groups for means.")
  }

  sum_mat <- matrix(0, nrow = n_f, ncol = n_g,
                    dimnames = list(features, groups))
  npos_mat <- if (isTRUE(return.pct)) {
    matrix(0, nrow = n_f, ncol = n_g, dimnames = list(features, groups))
  } else {
    NULL
  }
  count_vec <- setNames(integer(n_g), groups)
  found_feat <- logical(n_f)
  names(found_feat) <- features

  total <- length(cell_to_group)
  done <- 0L

  for (h in handles) {
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    row_idx <- match(features, h$feature_ids)
    present <- which(!is.na(row_idx))
    if (length(present)) found_feat[present] <- TRUE

    idx <- which(h$cell_ids %in% names(cell_to_group))
    if (length(idx) < 1L) next
    runs <- .ic2_contig_runs(sort(idx))

    for (run in runs) {
      starts <- seq.int(1L, length(run), by = as.integer(chunk.size))
      for (s in starts) {
        e <- min(s + as.integer(chunk.size) - 1L, length(run))
        cols <- run[s:e]
        raw <- .ic2_read_cols(h, cols)
        ## features x cells (requested order)
        block <- matrix(0, nrow = n_f, ncol = length(cols),
                        dimnames = list(features, h$cell_ids[cols]))
        if (length(present)) {
          sub <- raw[row_idx[present], , drop = FALSE]
          if (inherits(sub, "Matrix")) {
            sub <- as.matrix(.ic2_to_dgCMatrix(sub))
          } else {
            sub <- as.matrix(sub)
          }
          storage.mode(sub) <- "double"
          block[present, ] <- sub
        }
        block[!is.finite(block)] <- 0
        if (!is.null(size_factors)) {
          block <- .ic2_scale_cols_by_sf(block, size_factors)
        }
        block <- .ic2_assay_transform_block(block, transform = transform)

        cids <- h$cell_ids[cols]
        grp <- as.character(cell_to_group[cids])
        for (g in unique(grp)) {
          hit <- which(grp == g)
          if (!length(hit)) next
          sum_mat[, g] <- sum_mat[, g] + rowSums(block[, hit, drop = FALSE])
          if (!is.null(npos_mat)) {
            npos_mat[, g] <- npos_mat[, g] +
              rowSums(block[, hit, drop = FALSE] > 0)
          }
          count_vec[[g]] <- count_vec[[g]] + length(hit)
        }
        done <- done + length(cols)
        if (!is.null(progress)) {
          progress <- .ic2_progress_bar(progress, done, total, prefix = "  heat")
        }
      }
    }
  }

  missing <- features[!found_feat]
  if (length(missing)) {
    warning(
      "Gene(s) not found in assay: ",
      paste(utils::head(missing, 20L), collapse = ", "),
      if (length(missing) > 20L) ", ..." else ""
    )
  }
  keep <- found_feat
  if (!any(keep)) stop("None of the requested genes were found.")
  sum_mat <- sum_mat[keep, , drop = FALSE]
  mean_mat <- sweep(sum_mat, 2L, pmax(as.numeric(count_vec), 1), `/`)
  ## Drop empty groups
  keep_g <- count_vec > 0L
  mean_mat <- mean_mat[, keep_g, drop = FALSE]
  if (!isTRUE(return.pct)) {
    return(mean_mat)
  }
  npos_mat <- npos_mat[keep, keep_g, drop = FALSE]
  pct_mat <- sweep(npos_mat, 2L, pmax(as.numeric(count_vec[keep_g]), 1), `/`) * 100
  list(mean = mean_mat, pct = pct_mat)
}
