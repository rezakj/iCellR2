#' iCellR batch alignment (IBA): CPCA or CCCA
#'
#' Same idea as classic \code{iba}: each cell is replaced by the mean of its
#' \code{k} nearest neighbors \emph{inside every condition}, then PCA is run
#' again on that smoothed matrix.
#'
#' The classic \eqn{n \times n} distance matrix is not built. Neighbors are
#' queried per condition with \pkg{FNN} (\code{get.knnx}). \code{CPCA}
#' (default) averages PCA coordinates only  -  no expression I/O. \code{CCCA}
#' streams the expression mean with the same single-pass writer as
#' \code{\link{run.impute}}.
#'
#' Aligned scores are written to \code{reductions/aligned_pca}.
#' \code{reductions/pca} is left unchanged.
#'
#' Conditions come from \code{cell_data[[condition.by]]} (default
#' \code{condition_id}), not from splitting cell barcodes on \code{"_"}.
#'
#' @param x An \code{iCellR2} object.
#' @param dims PCA columns used for neighbor search. Default \code{1:30}.
#' @param k Neighbors taken inside each condition. Default \code{10}.
#' @param ba.method \code{"CPCA"} (default) or \code{"CCCA"}.
#' @param method,top.rank,features Passed to \code{\link{run.pca}} when
#'   \code{refit = TRUE} (and again after CCCA).
#' @param scale.data If \code{TRUE} (default), the PCA fit centers and scales.
#'   Mapped onto \code{run.pca(center, scale.)}; \code{plus.log.value} is
#'   ignored (iCellR2 uses \code{log2p1}).
#' @param gene.list Deprecated. Use \code{features}.
#' @param plus.log.value Ignored. Kept so classic calls do not error.
#' @param refit Recompute PCA before alignment. Default \code{FALSE} when a
#'   PCA reduction is already on disk.
#' @param condition.by Column in \code{cell_data}. Default \code{"condition_id"}.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Source assay for a refit / CCCA. Default: active \code{data}.
#' @param chunk.cells CCCA gather size. Default \code{30000}. Ignored by CPCA.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} object with aligned PCA embeddings.
#'
#' @examples
#' \dontrun{
#' obj <- iba(obj)                              # CPCA, k = 10, dims 1:30
#' obj <- iba(obj, ba.method = "CCCA", k = 10)
#' }
#'
#' @export
iba <- function(
    x = NULL,
    dims = 1:30,
    k = 10L,
    ba.method = "CPCA",
    method = c("base.mean.rank", "gene.model", "features"),
    top.rank = 500L,
    plus.log.value = 0.1,
    scale.data = TRUE,
    gene.list = NULL,
    refit = NULL,
    condition.by = "condition_id",
    modality = "rna",
    assay.name = NULL,
    features = NULL,
    chunk.cells = 30000L,
    verbose = TRUE
) {
  .ic2_require_object(x)
  ba.method <- toupper(as.character(ba.method)[1L])
  ba.method <- match.arg(ba.method, c("CPCA", "CCCA"))
  method <- match.arg(method)
  modality <- tolower(modality)
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) stop("'k' must be >= 1.")
  chunk.cells <- as.integer(chunk.cells)
  if (!is.finite(chunk.cells) || chunk.cells < 1L) chunk.cells <- 30000L
  if (!is.null(gene.list) && !identical(gene.list, "character")) {
    features <- gene.list
    method <- "features"
  }

  has_pca <- !is.null(x@modalities[[modality]]$reductions$pca)
  if (is.null(refit)) refit <- !isTRUE(has_pca)

  rep <- .ic2_new_reporter("iba", verbose, soda_root = x@storage$root_dir)
  rep <- .ic2_banner_start(rep, c(
    paste0("Method: ", ba.method, "   |   k = ", k, " / condition"),
    paste0("Conditions: cell_data$", condition.by)
  ))

  if (isTRUE(refit) || !isTRUE(has_pca)) {
    .ic2_step(rep, "[1/4] PCA (streaming)")
    x <- run.pca(
      x,
      modality = modality,
      assay.name = assay.name,
      method = method,
      top.rank = as.integer(top.rank),
      features = features,
      npcs = max(as.integer(dims), na.rm = TRUE),
      center = isTRUE(scale.data),
      scale. = isTRUE(scale.data),
      verbose = verbose
    )
  } else {
    .ic2_step(rep, "[1/4] Using PCA already on disk")
  }

  .ic2_step(rep, "[2/4] Per-condition kNN (no full distance matrix)")
  emb <- .ic2_get_plot_embedding(x, modality, "pca")
  dims <- as.integer(dims)
  dims <- dims[is.finite(dims) & dims >= 1L & dims <= ncol(emb)]
  if (length(dims) < 1L) {
    stop("No valid 'dims' (PCA has ", ncol(emb), " columns).")
  }
  emb <- emb[, dims, drop = FALSE]
  cell_ids <- rownames(emb)
  cond <- .ic2_iba_conditions(x, cell_ids, condition.by)
  knn <- .ic2_iba_knn_by_condition(emb, cond, k = k)
  n_cond <- length(unique(cond))

  if (identical(ba.method, "CPCA")) {
    .ic2_step(rep, "[3/4] CPCA: mean of neighbor PCs")
    acc <- matrix(0, nrow = nrow(emb), ncol = ncol(emb))
    for (j in seq_len(ncol(knn))) {
      acc <- acc + emb[knn[, j], , drop = FALSE]
    }
    acc <- acc / ncol(knn)
    rownames(acc) <- cell_ids

    .ic2_step(rep, "[4/4] PCA of smoothed PCs")
    fit <- stats::prcomp(acc, center = TRUE, scale. = FALSE, rank. = ncol(acc))
    aligned <- fit$x
    rownames(aligned) <- cell_ids
    colnames(aligned) <- paste0("PC", seq_len(ncol(aligned)))
    sdev <- as.numeric(fit$sdev)
    rotation <- fit$rotation
    rownames(rotation) <- colnames(emb)
    x <- .ic2_iba_store_pca(
      x, modality,
      emb = aligned,
      rotation = rotation,
      sdev = sdev,
      ba.method = ba.method,
      k = k,
      dims = dims,
      n_cond = n_cond,
      condition.by = condition.by
    )
  } else {
    ## CCCA: same per-condition knn, stream neighbor-mean expression, then PCA
    .ic2_step(rep, "[3/4] CCCA: stream neighbor-mean expression")
    assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
    handles <- .ic2_get_assay_handles(x, modality, assay.name)
    feature_universe <- character(0)
    for (h in handles) {
      h2 <- .ic2_enrich_handle(h, load_dimnames = TRUE)
      feature_universe <- union(feature_universe, as.character(h2$feature_ids))
    }
    feature_universe <- as.character(feature_universe)
    if (length(feature_universe) < 1L) stop("No features found for CCCA.")

    dest <- .ic2_soda_assay_path(
      x, modality, "iba_ccca", sample_id = NULL, ext = ".h5"
    )
    handle <- .ic2_impute_write_csc(
      handles = handles,
      cell_ids = cell_ids,
      knn_idx = knn,
      feature_universe = feature_universe,
      dest_path = dest,
      chunk.cells = chunk.cells,
      round.num = TRUE,
      rounding.digits = 3L,
      drop.zero.genes = TRUE,
      progress = function(cur, tot) {
        rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  ccca")
      }
    )
    handle$modality <- modality
    handle$sample_id <- "iba_ccca"
    x <- .ic2_assign_assay(
      x,
      modality = modality,
      assay.name = "iba_ccca",
      handles = list(handle),
      virtual = FALSE,
      feature_universe = NULL,
      set.active = FALSE
    )

    ## Keep original reductions/pca untouched: snapshot, PCA into it, copy to
    ## aligned_pca, then restore the snapshot.
    .ic2_step(rep, "[4/4] PCA on CCCA-smoothed assay \u2192 aligned_pca")
    root <- x@storage$root_dir
    pca_dir <- file.path(root, "reductions", "pca")
    pca_bak <- file.path(root, "reductions", "pca_pre_ccca")
    pca_slot <- x@modalities[[modality]]$reductions$pca
    if (dir.exists(pca_dir)) {
      if (dir.exists(pca_bak)) unlink(pca_bak, recursive = TRUE)
      dir.create(pca_bak, recursive = TRUE, showWarnings = FALSE)
      file.copy(list.files(pca_dir, full.names = TRUE), pca_bak, recursive = TRUE)
    }

    x <- run.pca(
      x,
      modality = modality,
      assay.name = "iba_ccca",
      method = method,
      top.rank = as.integer(top.rank),
      features = features,
      npcs = max(as.integer(dims), na.rm = TRUE),
      center = isTRUE(scale.data),
      scale. = isTRUE(scale.data),
      verbose = verbose
    )
    pca_emb <- .ic2_get_plot_embedding(x, modality, "pca")
    red <- x@modalities[[modality]]$reductions$pca
    rot_path <- red$feature_loadings
    rotation <- if (!is.null(rot_path) && file.exists(rot_path)) {
      .ic2_h5_read_dense(rot_path)
    } else {
      matrix(numeric(0), nrow = 0L, ncol = ncol(pca_emb))
    }
    sdev <- as.numeric(red$stdev %||% rep(1, ncol(pca_emb)))
    x <- .ic2_iba_store_pca(
      x, modality,
      emb = pca_emb,
      rotation = rotation,
      sdev = sdev,
      ba.method = ba.method,
      k = k,
      dims = dims,
      n_cond = n_cond,
      condition.by = condition.by
    )

    ## Restore original pca from snapshot
    if (dir.exists(pca_bak)) {
      unlink(pca_dir, recursive = TRUE)
      dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(list.files(pca_bak, full.names = TRUE), pca_dir, recursive = TRUE)
      if (!is.null(pca_slot)) {
        x@modalities[[modality]]$reductions$pca <- pca_slot
      }
    }
  }

  x <- .ic2_log_command(
    x, "iba",
    arguments = list(
      ba.method = ba.method,
      k = as.integer(k),
      dims = as.integer(dims),
      condition.by = condition.by,
      refit = isTRUE(refit)
    ),
    outputs = list(
      aligned_pca = file.path(x@storage$root_dir, "reductions", "aligned_pca")
    )
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  .ic2_finish(rep, extra_lines = paste0(
    ba.method, "   |   ", .ic2_fmt_num(nrow(emb)), " cells   |   ",
    n_cond, " conditions   |   k = ", k
  ))
  x
}


#' Condition labels aligned to embedding cell ids
#' @keywords internal
.ic2_iba_conditions <- function(x, cell_ids, condition.by = "condition_id") {
  cd <- x@cell_data
  if (is.null(cd) || nrow(cd) < 1L) {
    stop("cell_data is empty; iba needs a condition column.")
  }
  if (!condition.by %in% colnames(cd)) {
    if (identical(condition.by, "condition_id") && "sample_id" %in% colnames(cd)) {
      condition.by <- "sample_id"
    } else {
      stop("cell_data has no column '", condition.by, "'.")
    }
  }
  if (!is.null(rownames(cd)) && all(cell_ids %in% rownames(cd))) {
    cond <- as.character(cd[cell_ids, condition.by])
  } else if ("cell_id" %in% colnames(cd)) {
    m <- match(cell_ids, as.character(cd$cell_id))
    if (anyNA(m)) {
      stop(sum(is.na(m)), " embedding cells are missing from cell_data.")
    }
    cond <- as.character(cd[[condition.by]][m])
  } else {
    stop("Cannot align cell_data to PCA cell ids.")
  }
  if (anyNA(cond) || any(!nzchar(cond))) {
    stop("Some cells have an empty ", condition.by, ".")
  }
  if (length(unique(cond)) < 2L) {
    warning("Only one condition in '", condition.by, "'; alignment is a within-batch smooth.")
  }
  cond
}


#' k nearest cells inside each condition (positions into \code{emb})
#' @keywords internal
.ic2_iba_knn_by_condition <- function(emb, cond, k = 10L) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required for iba.")
  }
  emb <- as.matrix(emb)
  cond <- as.character(cond)
  levels <- unique(cond)
  blocks <- vector("list", length(levels))
  for (i in seq_along(levels)) {
    sel <- which(cond == levels[[i]])
    kk <- min(as.integer(k), length(sel))
    data <- emb[sel, , drop = FALSE]
    if (kk == 1L) {
      idx <- matrix(sel[[1L]], nrow = nrow(emb), ncol = 1L)
    } else {
      kn <- FNN::get.knnx(data = data, query = emb, k = kk)
      idx <- matrix(sel[kn$nn.index], nrow = nrow(emb), ncol = kk)
    }
    if (kk < k) {
      pad <- matrix(idx[, kk], nrow = nrow(emb), ncol = k - kk)
      idx <- cbind(idx, pad)
    }
    blocks[[i]] <- idx
  }
  out <- do.call(cbind, blocks)
  storage.mode(out) <- "integer"
  out
}


#' Write CPCA scores to reductions/aligned_pca. Does not touch reductions/pca.
#' @keywords internal
.ic2_iba_store_pca <- function(
    x, modality, emb, rotation, sdev, ba.method, k, dims, n_cond, condition.by
) {
  root <- x@storage$root_dir
  name <- "aligned_pca"
  red_dir <- file.path(root, "reductions", name)
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  load_path <- file.path(red_dir, "feature_loadings.h5")
  stdev_path <- file.path(red_dir, "stdev.tsv")
  .ic2_h5_write_dense(emb, emb_path)
  .ic2_h5_write_dense(rotation, load_path)
  utils::write.table(
    data.frame(pc = seq_along(sdev), stdev = as.numeric(sdev), stringsAsFactors = FALSE),
    file = stdev_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  meta <- list(
    name = name,
    mark = "Aligned",
    aligned = TRUE,
    ba.method = ba.method,
    k = as.integer(k),
    dims = as.integer(dims),
    n_conditions = as.integer(n_cond),
    condition.by = condition.by,
    npcs = ncol(emb),
    n_cells = nrow(emb),
    created = as.character(Sys.time()),
    paths = list(embeddings = emb_path, loadings = load_path, stdev = stdev_path)
  )
  .ic2_soda_write_json(meta, file.path(red_dir, "reduction.json"))
  x@modalities[[modality]]$reductions[[name]] <- list(
    cell_embeddings = emb_path,
    feature_loadings = load_path,
    stdev = as.numeric(sdev),
    npcs = ncol(emb),
    meta_path = file.path(red_dir, "reduction.json"),
    aligned = TRUE,
    mark = "Aligned"
  )
  x@storage$reduction_paths[[name]] <- red_dir
  x@params$batch_alignment <- list(
    ba.method = ba.method,
    reduction = name,
    mark = "Aligned",
    k = as.integer(k),
    dims = as.integer(dims),
    n_conditions = as.integer(n_cond),
    condition.by = condition.by
  )
  x
}


#' Pick aligned_<family> when it exists, otherwise <family>.
#'
#' \code{use = "auto"} (default) prefers \code{aligned_pca} / \code{aligned_umap}.
#' \code{use = "pca"} or \code{use = "reductions/pca"} forces the original.
#' \code{use = "aligned_pca"} or \code{use = "reductions/aligned_pca"} forces aligned.
#' @keywords internal
.ic2_resolve_use <- function(x, modality = "rna", family = c("pca", "umap"), use = "auto") {
  family <- match.arg(family)
  modality <- tolower(modality)
  u <- tolower(as.character(use %||% "auto")[1L])
  u <- sub("^reductions/", "", u)
  aligned_name <- paste0("aligned_", family)
  reds <- x@modalities[[modality]]$reductions
  has <- function(nm) !is.null(reds[[nm]])
  if (u %in% c("auto", "")) {
    if (has(aligned_name)) {
      name <- aligned_name
    } else if (has(family)) {
      name <- family
    } else {
      stop(
        "No '", aligned_name, "' or '", family, "' reduction for modality '",
        modality, "'. Run ", if (identical(family, "pca")) "run.pca()" else "run.umap()",
        " first."
      )
    }
  } else if (u %in% c(family, aligned_name)) {
    name <- u
    if (!has(name)) {
      stop(
        "No '", name, "' reduction. ",
        if (identical(name, aligned_name)) "Run iba() first." else
          paste0("Run ", if (identical(family, "pca")) "run.pca()" else "run.umap()", " first.")
      )
    }
  } else {
    stop(
      "'use' must be 'auto', '", family, "', '", aligned_name,
      "', or 'reductions/", aligned_name, "'. Got '", use, "'."
    )
  }
  aligned <- identical(name, aligned_name)
  list(
    name = name,
    aligned = aligned,
    mark = if (aligned) "Aligned" else NULL,
    family = family
  )
}


#' Output reduction name. Aligned inputs are stored as aligned_<base>.
#' @keywords internal
.ic2_aligned_out_name <- function(base, aligned) {
  base <- as.character(base)[1L]
  if (isTRUE(aligned)) {
    if (startsWith(base, "aligned_")) base else paste0("aligned_", base)
  } else {
    base
  }
}


#' Default scatter reduction for \code{cluster.plot} / \code{gene.plot}.
#'
#' Prefers \code{active_reduction} when it is plottable, otherwise the first
#' available among aligned UMAP, UMAP, then other common embeddings. This keeps
#' scatter plots on \code{aligned_umap} after \code{iba()} + \code{run.umap()}
#' instead of hard-coding \code{"umap"}.
#' @keywords internal
.ic2_default_scatter_reduction <- function(
    x,
    modality = "rna",
    allowed = c(
      "aligned_umap", "umap", "aligned_knetl", "knetl",
      "aligned_tsne", "tsne", "cgumap", "ptumap", "cgmap",
      "aligned_pca", "pca", "cgpca"
    )
) {
  modality <- tolower(as.character(modality)[1L])
  reds <- x@modalities[[modality]]$reductions
  if (is.null(reds) || !length(reds)) {
    stop(
      "No reductions for modality '", modality,
      "'. Run run.pca() / run.umap() first."
    )
  }
  active <- x@modalities[[modality]]$active_reduction
  if (!is.null(active) && nzchar(as.character(active)[1L])) {
    a <- as.character(active)[1L]
    if (a %in% allowed && !is.null(reds[[a]])) {
      return(a)
    }
  }
  prefer <- c(
    "aligned_umap", "umap", "aligned_knetl", "knetl",
    "aligned_tsne", "tsne", "cgumap", "ptumap", "cgmap",
    "aligned_pca", "pca", "cgpca"
  )
  for (cand in prefer) {
    if (cand %in% allowed && !is.null(reds[[cand]])) {
      return(cand)
    }
  }
  stop(
    "No plottable reduction for modality '", modality,
    "'. Run run.umap() (or run.pca()) first."
  )
}
