#' Differential expression between two groups (chunked, on-disk)
#'
#' Two-group DE without loading the full assay. Streams once, accumulates
#' sum / sum-of-squares for each side, then computes means, fold-change
#' (\code{group.2 / group.1}), and Welch \code{t.test} p-values from
#' aggregates. \code{wilcox.test} uses a second pass with optional cell
#' downsampling.
#'
#' \code{group.1}, \code{group.2}, and \code{within} all accept **multiple**
#' ids, e.g. \code{c(1, 3, 2)} or \code{c("KO", "base")}.
#'
#' @param x An \code{iCellR2} object.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Assay to use; default active / \code{"data"}.
#' @param transform \code{"auto"} (default): \code{log2p1} unless assay is
#'   \code{scale.data}; else \code{"none"}, \code{"log2p1"}, \code{"log1p"}.
#' @param compare Comparison mode:
#'   \itemize{
#'     \item \code{"clusters"}  -  \code{group.1}/\code{group.2} are cluster ids
#'     \item \code{"conditions"}  -  condition ids from \code{condition_id}
#'     \item \code{"conditions_in_cluster"}  -  conditions within \code{within} cluster(s)
#'     \item \code{"clusters_in_condition"}  -  clusters within \code{within} condition(s)
#'   }
#' @param group.1 First group id(s). Vector OK.
#' @param group.2 Second group id(s). Vector OK. FC is mean(group.2)/mean(group.1).
#' @param within Restricting cluster(s) or condition(s) for the
#'   \code{*_in_*} modes. Vector OK. Required there.
#' @param pval.test \code{"t.test"} (default) or \code{"wilcox.test"}.
#' @param p.adjust.method Passed to \code{p.adjust}. Default \code{"hochberg"}.
#' @param max.cells Max cells for \code{wilcox.test}. Default \code{5000000}.
#'   Ignored for \code{t.test}.
#' @param chunk.size Cells per HDF5 read. Default \code{200000}.
#' @param verbose Progress style.
#' @param de.by,cond.1,cond.2,base.cond Deprecated aliases for
#'   \code{compare}, \code{group.1}, \code{group.2}, \code{within}.
#'
#' @return A DE \code{data.frame} with columns \code{gene}, \code{baseMean},
#'   mean in each group (named), \code{foldChange}, \code{log2FoldChange},
#'   \code{pval}, \code{padj}.
#'
#' @examples
#' \dontrun{
#' de <- run.diff.exp(obj, compare = "clusters",
#'                    group.1 = 1, group.2 = 2)
#' de <- run.diff.exp(obj, compare = "clusters",
#'                    group.1 = c(1, 3), group.2 = c(2, 5))
#' de <- run.diff.exp(obj, compare = "conditions",
#'                    group.1 = "KO", group.2 = c("base", "WT"))
#' de <- run.diff.exp(obj, compare = "conditions_in_cluster",
#'                    within = 3, group.1 = "KO", group.2 = "base")
#' volcano.ma.plot(de)
#' }
#'
#' @export
run.diff.exp <- function(
    x = NULL,
    modality = "rna",
    assay.name = NULL,
    transform = c("auto", "none", "log2p1", "log1p"),
    compare = c(
      "clusters", "conditions",
      "conditions_in_cluster", "clusters_in_condition"
    ),
    group.1 = NULL,
    group.2 = NULL,
    within = NULL,
    pval.test = c("t.test", "wilcox.test"),
    p.adjust.method = "hochberg",
    max.cells = 5000000L,
    chunk.size = 200000L,
    verbose = TRUE,
    ## deprecated classic aliases
    de.by = NULL,
    cond.1 = NULL,
    cond.2 = NULL,
    base.cond = NULL
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  transform <- match.arg(transform)
  pval.test <- match.arg(pval.test)
  chunk.size <- as.integer(chunk.size)
  max.cells <- as.integer(max.cells)

  ## ---- resolve aliases -------------------------------------------------
  compare <- .ic2_de_resolve_compare(compare, de.by)
  if (!is.null(cond.1)) {
    if (!is.null(group.1)) {
      warning("Both 'group.1' and deprecated 'cond.1' set; using 'group.1'.")
    } else {
      warning("'cond.1' is deprecated; use 'group.1'.")
      group.1 <- cond.1
    }
  }
  if (!is.null(cond.2)) {
    if (!is.null(group.2)) {
      warning("Both 'group.2' and deprecated 'cond.2' set; using 'group.2'.")
    } else {
      warning("'cond.2' is deprecated; use 'group.2'.")
      group.2 <- cond.2
    }
  }
  if (!is.null(base.cond)) {
    if (!is.null(within)) {
      warning("Both 'within' and deprecated 'base.cond' set; using 'within'.")
    } else {
      warning("'base.cond' is deprecated; use 'within'.")
      within <- base.cond
    }
  }
  ## classic used base.cond = 0 to mean "unset"
  if (!is.null(within) && length(within) == 1L && is.numeric(within) &&
      isTRUE(within == 0)) {
    within <- NULL
  }

  if (is.null(group.1) || is.null(group.2)) {
    stop("Provide both 'group.1' and 'group.2' (vectors OK, e.g. c(1,3) or c('KO','base')).")
  }

  name_a <- paste(group.1, collapse = "_")
  name_b <- paste(group.2, collapse = "_")
  if (identical(compare, "conditions_in_cluster")) {
    if (is.null(within)) {
      stop("Set 'within' to the cluster id(s) for compare = 'conditions_in_cluster'.")
    }
    name_a <- paste0(name_a, ".inClust", paste(within, collapse = "_"))
    name_b <- paste0(name_b, ".inClust", paste(within, collapse = "_"))
  }
  if (identical(compare, "clusters_in_condition")) {
    if (is.null(within)) {
      stop("Set 'within' to the condition id(s) for compare = 'clusters_in_condition'.")
    }
    name_a <- paste0(name_a, ".inCond", paste(within, collapse = "_"))
    name_b <- paste0(name_b, ".inCond", paste(within, collapse = "_"))
  }

  assay.names <- names(x@modalities[[modality]]$assays %||% list())
  if (is.null(assay.name)) {
    assay.name <- if ("data" %in% assay.names) "data" else NULL
  }
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  if (identical(transform, "auto")) {
    transform <- if (identical(assay.name, "scale.data")) "none" else "log2p1"
  }

  rep <- .ic2_new_reporter("run.diff.exp", verbose, soda_root = x@storage$root_dir)
  meta <- x@cell_data
  if (nrow(meta) < 2L) stop("cell_data has fewer than 2 cells.")

  groups <- .ic2_de_resolve_groups(
    meta = meta,
    compare = compare,
    group.1 = group.1,
    group.2 = group.2,
    within = within
  )
  cells_a <- groups$cells_a
  cells_b <- groups$cells_b
  n_a <- length(cells_a)
  n_b <- length(cells_b)
  if (n_a < 2L || n_b < 2L) {
    stop(
      "Need >=2 cells per group. Got group.1=", n_a,
      "  group.2=", n_b, "."
    )
  }
  overlap <- intersect(cells_a, cells_b)
  if (length(overlap)) {
    stop(
      length(overlap), " cells are in both groups; ",
      "group.1 and group.2 must be disjoint."
    )
  }

  cell_to_grp <- c(
    setNames(rep(1L, n_a), cells_a),
    setNames(rep(2L, n_b), cells_b)
  )

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  h0 <- .ic2_enrich_handle(handles[[1]], load_dimnames = TRUE)
  genes <- h0$feature_ids
  n_feat <- length(genes)
  if (n_feat < 1L) stop("Assay has no features.")

  sum_a <- numeric(n_feat)
  sum_b <- numeric(n_feat)
  sumsq_a <- numeric(n_feat)
  sumsq_b <- numeric(n_feat)
  count_a <- 0L
  count_b <- 0L

  rep <- .ic2_banner_start(rep, c(
    paste0("Assay: ", assay.name, "   |   transform: ", transform),
    paste0("compare: ", compare, "   |   test: ", pval.test),
    paste0(
      name_a, " (n=", .ic2_fmt_num(n_a), ")  vs  ",
      name_b, " (n=", .ic2_fmt_num(n_b), ")"
    )
  ))

  .ic2_step(rep, "[1/2] Streaming assay in file order (group sum / sumsq)")
  total_plan <- n_a + n_b
  done_cells <- 0L

  for (hi in seq_along(handles)) {
    handle <- .ic2_enrich_handle(handles[[hi]], load_dimnames = TRUE)
    if (!identical(handle$feature_ids, genes)) {
      if (!all(genes %in% handle$feature_ids)) {
        stop("Feature ids differ across assay handles.")
      }
      row_map <- match(genes, handle$feature_ids)
    } else {
      row_map <- NULL
    }

    n_file <- length(handle$cell_ids)
    if (n_file < 1L) next
    ## One sequential read per chunk. Cluster members are interleaved in the
    ## file, so picking only those columns opens HDF5 once per short run.
    starts <- seq.int(1L, n_file, by = chunk.size)
    for (s in starts) {
      e <- min(s + chunk.size - 1L, n_file)
      cols <- s:e
      grp <- cell_to_grp[handle$cell_ids[cols]]
      keep <- which(!is.na(grp))
      if (length(keep) == 0L) next
      block <- .ic2_read_cols(handle, cols)
      if (!is.null(row_map)) block <- block[row_map, , drop = FALSE]
      block <- block[, keep, drop = FALSE]
      block <- .ic2_assay_transform_block(block, transform = transform)
      grp <- grp[keep]

      for (g in c(1L, 2L)) {
        hit <- which(grp == g)
        if (length(hit) == 0L) next
        sub <- block[, hit, drop = FALSE]
        if (inherits(sub, "Matrix")) {
          rs <- as.numeric(Matrix::rowSums(sub))
          rssq <- as.numeric(Matrix::rowSums(sub * sub))
        } else {
          sub <- as.matrix(sub)
          storage.mode(sub) <- "double"
          rs <- rowSums(sub)
          rssq <- rowSums(sub * sub)
        }
        if (g == 1L) {
          sum_a <- sum_a + rs
          sumsq_a <- sumsq_a + rssq
          count_a <- count_a + length(hit)
        } else {
          sum_b <- sum_b + rs
          sumsq_b <- sumsq_b + rssq
          count_b <- count_b + length(hit)
        }
      }
      done_cells <- done_cells + length(keep)
      rep <- .ic2_progress_bar(rep, done_cells, total_plan, prefix = "  de")
    }
  }

  if (count_a < 2L || count_b < 2L) {
    stop(
      "Accumulated <2 cells in a group (check assay cell ids). ",
      "A=", count_a, " B=", count_b
    )
  }

  .ic2_step(rep, "[2/2] Computing means, FC, and p-values")
  n1 <- as.numeric(count_a)
  n2 <- as.numeric(count_b)
  m1 <- sum_a / n1
  m2 <- sum_b / n2
  v1 <- (sumsq_a - (sum_a^2) / n1) / max(n1 - 1, 1)
  v2 <- (sumsq_b - (sum_b^2) / n2) / max(n2 - 1, 1)
  v1[!is.finite(v1) | v1 < 0] <- 0
  v2[!is.finite(v2) | v2 < 0] <- 0

  total_n <- n1 + n2
  base_mean <- (sum_a + sum_b) / total_n

  ## FC = mean(group.2) / mean(group.1)
  fc <- m2 / m1
  fc[!is.finite(fc) & m2 > 0 & m1 == 0] <- Inf
  fc[!is.finite(fc) & m2 == 0 & m1 == 0] <- NaN
  lfc <- log2(fc)
  ok <- !is.nan(fc)

  if (identical(pval.test, "t.test")) {
    se2 <- v1 / n1 + v2 / n2
    se2[se2 <= 0 | !is.finite(se2)] <- .Machine$double.eps
    tstat <- (m1 - m2) / sqrt(se2)
    a <- v1 / n1
    b <- v2 / n2
    dfree <- (a + b)^2 / (a^2 / max(n1 - 1, 1) + b^2 / max(n2 - 1, 1))
    dfree[!is.finite(dfree) | dfree < 1] <- 1
    pval <- 2 * stats::pt(-abs(tstat), df = dfree)
  } else {
    all_ids <- names(cell_to_grp)
    if (length(all_ids) > max.cells) {
      set.seed(1L)
      n_take_a <- max(
        2L,
        as.integer(round(
          as.numeric(max.cells) * as.numeric(n_a) / as.numeric(n_a + n_b)
        ))
      )
      n_take_a <- min(n_take_a, n_a)
      n_take_b <- max(2L, as.integer(max.cells) - n_take_a)
      n_take_b <- min(n_take_b, n_b)
      wa <- if (n_a <= n_take_a) cells_a else sample(cells_a, n_take_a)
      wb <- if (n_b <= n_take_b) cells_b else sample(cells_b, n_take_b)
      wilcox_cells <- c(wa, wb)
      .ic2_step(rep, paste0(
        "  wilcox downsample: ", .ic2_fmt_num(length(wilcox_cells)),
        " / ", .ic2_fmt_num(length(all_ids)), " cells"
      ))
    } else {
      wilcox_cells <- all_ids
    }
    use_genes <- genes[ok]
    if (length(use_genes) < 1L) stop("No genes with defined fold-change.")
    .ic2_step(rep, paste0(
      "  wilcox on ", .ic2_fmt_num(length(use_genes)), " genes (batched)"
    ))
    grp_s <- cell_to_grp[wilcox_cells]
    in_a <- grp_s == 1L
    pval <- rep(NA_real_, n_feat)
    names(pval) <- genes
    gene_chunk <- 500L
    gstarts <- seq.int(1L, length(use_genes), by = gene_chunk)
    for (gs in gstarts) {
      ge <- min(gs + gene_chunk - 1L, length(use_genes))
      gbatch <- use_genes[gs:ge]
      expr <- .ic2_load_features_for_cells(
        handles = handles,
        features = gbatch,
        cell_ids = wilcox_cells,
        chunk.size = chunk.size,
        transform = transform
      )
      for (g in gbatch) {
        if (!g %in% colnames(expr)) next
        x1 <- expr[wilcox_cells[in_a], g]
        x2 <- expr[wilcox_cells[!in_a], g]
        pval[[g]] <- tryCatch(
          stats::wilcox.test(x = x1, y = x2)$p.value,
          error = function(e) NA_real_
        )
      }
      rep <- .ic2_progress_bar(rep, ge, length(use_genes), prefix = "  wilcox")
    }
    pval <- unname(pval)
  }

  pval[!ok] <- NA_real_
  padj <- stats::p.adjust(pval, method = p.adjust.method)

  out <- data.frame(
    gene = genes,
    baseMean = base_mean,
    mean_a = m1,
    mean_b = m2,
    foldChange = fc,
    log2FoldChange = lfc,
    pval = pval,
    padj = padj,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out)[names(out) == "mean_a"] <- name_a
  names(out)[names(out) == "mean_b"] <- name_b
  out <- out[ok, , drop = FALSE]
  rownames(out) <- out$gene

  if (!is.null(x@storage$root_dir) && nzchar(x@storage$root_dir)) {
    out_dir <- file.path(x@storage$root_dir, "metadata", "diffexp")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    safe <- gsub("[^A-Za-z0-9._-]+", "_", paste(compare, name_a, name_b, sep = "__"))
    out_path <- .ic2_soda_write_table(
      out, file.path(out_dir, paste0(modality, "_", safe))
    )
    invisible(.ic2_log_command(
      x, "run.diff.exp",
      arguments = list(
        modality = modality, assay.name = assay.name, transform = transform,
        compare = compare, group.1 = group.1, group.2 = group.2, within = within,
        pval.test = pval.test, p.adjust.method = p.adjust.method,
        max.cells = max.cells
      ),
      outputs = list(n_genes = nrow(out), path = out_path)
    ))
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
    rep <- .ic2_add_output(rep, "DE table", out_path)
  }

  .ic2_finish(rep, extra_lines = paste0(
    "DE genes: ", .ic2_fmt_num(nrow(out)),
    "   |   ", name_a, " vs ", name_b
  ))
  out
}

#' Map compare / legacy de.by to canonical compare value
#' @keywords internal
.ic2_de_resolve_compare <- function(compare, de.by = NULL) {
  legacy_map <- c(
    clusters = "clusters",
    conditions = "conditions",
    clustBase.condComp = "conditions_in_cluster",
    condBase.clustComp = "clusters_in_condition",
    conditions_in_cluster = "conditions_in_cluster",
    clusters_in_condition = "clusters_in_condition"
  )
  if (!is.null(de.by)) {
    if (!de.by %in% names(legacy_map)) {
      stop("Unknown deprecated 'de.by': ", de.by)
    }
    mapped <- unname(legacy_map[[de.by]])
    ## If user also passed compare explicitly as default vector, prefer de.by
    if (length(compare) > 1L) {
      warning("'de.by' is deprecated; use 'compare = \"", mapped, "\"'.")
      return(mapped)
    }
    if (!identical(compare, mapped)) {
      warning(
        "Both 'compare' and deprecated 'de.by' set; using 'compare = \"",
        compare, "\"'."
      )
    }
  }
  allowed <- c(
    "clusters", "conditions",
    "conditions_in_cluster", "clusters_in_condition"
  )
  if (length(compare) > 1L) compare <- compare[[1L]]
  if (!compare %in% allowed) {
    ## allow legacy string in compare itself
    if (compare %in% names(legacy_map)) {
      compare <- unname(legacy_map[[compare]])
    } else {
      stop(
        "'compare' must be one of: ",
        paste(allowed, collapse = ", ")
      )
    }
  }
  compare
}

#' Resolve two DE cell groups from cell_data
#' @keywords internal
.ic2_de_resolve_groups <- function(meta, compare, group.1, group.2, within) {
  ids <- rownames(meta)
  if (is.null(ids) || anyDuplicated(ids)) {
    stop("cell_data needs unique rownames (cell ids).")
  }

  if (identical(compare, "clusters")) {
    if (!"cluster_id" %in% colnames(meta)) {
      stop("No 'cluster_id' in cell_data. Run run.clustering() first.")
    }
    cl <- meta$cluster_id
    return(list(
      cells_a = ids[cl %in% group.1],
      cells_b = ids[cl %in% group.2]
    ))
  }

  if (identical(compare, "conditions")) {
    if (!"condition_id" %in% colnames(meta)) {
      stop("No 'condition_id' in cell_data.")
    }
    cd <- as.character(meta$condition_id)
    return(list(
      cells_a = ids[cd %in% as.character(group.1)],
      cells_b = ids[cd %in% as.character(group.2)]
    ))
  }

  if (identical(compare, "conditions_in_cluster")) {
    if (is.null(within) || length(within) < 1L) {
      stop("Set 'within' to the cluster id(s) for compare = 'conditions_in_cluster'.")
    }
    if (!"cluster_id" %in% colnames(meta) || !"condition_id" %in% colnames(meta)) {
      stop("Need both 'cluster_id' and 'condition_id' in cell_data.")
    }
    keep <- meta$cluster_id %in% within
    cd <- as.character(meta$condition_id)
    return(list(
      cells_a = ids[keep & cd %in% as.character(group.1)],
      cells_b = ids[keep & cd %in% as.character(group.2)]
    ))
  }

  if (identical(compare, "clusters_in_condition")) {
    if (is.null(within) || length(within) < 1L) {
      stop("Set 'within' to the condition id(s) for compare = 'clusters_in_condition'.")
    }
    if (!"cluster_id" %in% colnames(meta) || !"condition_id" %in% colnames(meta)) {
      stop("Need both 'cluster_id' and 'condition_id' in cell_data.")
    }
    keep <- as.character(meta$condition_id) %in% as.character(within)
    cl <- meta$cluster_id
    return(list(
      cells_a = ids[keep & cl %in% group.1],
      cells_b = ids[keep & cl %in% group.2]
    ))
  }

  stop("Unknown compare: ", compare)
}
