#' Impute an assay via kNN smoothing (streamed, on-disk)
#'
#' Classic iCellR \code{run.impute} (\code{iCellR.imp}) ported to SODA /
#' \code{iCellR2}: neighbors are found in a reduction (PCA/UMAP), then each
#' cell's expression is replaced by the mean of its \code{nn} neighbors
#' (including self, matching classic behaviour). The full expression matrix
#' is never densified  -  columns are gathered in chunks and written as a
#' sparse CSC assay under SODA.
#'
#' Requires \code{\link{run.pca}} (or \code{\link{run.umap}}) first. Source
#' assay defaults to normalized \code{"data"} when present. Prefer a
#' \strong{materialized} \code{norm.data} assay for speed; lazy norm works
#' via the same handles but pays an extra size-factor path when applicable.
#'
#' @param x An \code{iCellR2} object.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Source assay (default: active / \code{"data"}).
#' @param new.assay.name Output assay. Default \code{"imputed"}
#'   (\code{"atac.imputed"} when \code{modality = "atac"}).
#' @param nn Number of neighbors including self (classic default \code{10}).
#' @param dims Dimensions of the reduction to use. Default \code{1:10}
#'   (classic iCellR). Values outside available columns are dropped.
#' @param data.type Alias for \code{reduction}.
#' @param reduction Embedding for kNN: \code{"pca"} (default) or \code{"umap"}.
#' @param use Which stored reduction to read. \code{"auto"} uses an aligned
#'   embedding when present, otherwise \code{reduction}.
#' @param knn.dist.method Distance for kNN when FNN is unavailable.
#'   Default \code{"euclidean"}.
#' @param rounding.digits Digits when \code{round.num = TRUE}.
#' @param round.num Round nonzero imputed values (default \code{TRUE}).
#' @param chunk.cells Target cells per write/gather chunk. Default
#'   \code{30000}. Larger chunks amortize HDF5 I/O (more neighbor overlap);
#'   lower on memory-constrained machines.
#' @param drop.zero.genes Drop features with no coverage in the written
#'   assay (post-impute). Default \code{TRUE} (classic behaviour).
#' @param set.active Set the imputed assay active. Default \code{FALSE}.
#' @param overwrite Overwrite existing imputed assay file.
#' @param verbose Progress style.
#' @param ATAC.data Deprecated; use \code{modality = "atac"}.
#'
#' @return Updated \code{iCellR2} object with a new on-disk imputed assay.
#'
#' @examples
#' \dontrun{
#' obj <- run.impute(obj)                       # PCA dims 1:10, nn = 10
#' obj <- run.impute(obj, nn = 20, dims = 1:15)
#' gene.plot(obj, gene = "CD3D", data.type = "imputed.data")
#' }
#'
#' @export
run.impute <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    new.assay.name = NULL,
    nn = 10L,
    dims = 1:10,
    data.type = NULL,
    reduction = c("pca", "umap"),
    use = "auto",
    knn.dist.method = "euclidean",
    rounding.digits = 3L,
    round.num = TRUE,
    chunk.cells = 30000L,
    drop.zero.genes = TRUE,
    set.active = FALSE,
    overwrite = TRUE,
    verbose = TRUE,
    ATAC.data = FALSE
) {
  .ic2_require_object(x)
  if (isTRUE(ATAC.data)) modality <- "atac"
  modality <- tolower(modality)

  if (!is.null(data.type)) {
    reduction <- match.arg(data.type, c("pca", "umap", "tsne", "knetl", "diffusion"))
    if (!reduction %in% c("pca", "umap")) {
      stop(
        "reduction / data.type = '", reduction,
        "' is not supported in iCellR2 yet; use 'pca' or 'umap'."
      )
    }
  } else {
    reduction <- match.arg(reduction)
  }
  src <- .ic2_resolve_use(x, modality, reduction, use)
  reduction <- src$name

  assay_name_set <- !is.null(new.assay.name) && nzchar(as.character(new.assay.name)[1L])
  if (!assay_name_set) {
    new.assay.name <- if (isTRUE(src$aligned)) {
      "aligned.imputed"
    } else if (identical(modality, "atac")) {
      "atac.imputed"
    } else {
      "imputed"
    }
  }
  new.assay.name <- as.character(new.assay.name)[1L]
  nn <- as.integer(nn)
  if (!is.finite(nn) || nn < 1L) stop("'nn' must be >= 1.")
  chunk.cells <- as.integer(chunk.cells)
  if (!is.finite(chunk.cells) || chunk.cells < 1L) chunk.cells <- 30000L

  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  rep <- .ic2_new_reporter("run.impute", verbose, soda_root = x@storage$root_dir)

  ## ---- iCellR.imp (streamed) ----

  .ic2_step(rep, paste0("[1/4] Loading ", reduction, " embeddings"))
  emb <- .ic2_get_plot_embedding(x, modality, reduction)
  cell_ids <- rownames(emb)
  if (is.null(cell_ids) || length(cell_ids) < 1L) {
    stop("Embedding has no cell ids (rownames).")
  }

  ## Restrict to cells still in cell_data when available
  cd_ids <- tryCatch(rownames(x@cell_data), error = function(e) NULL)
  if (!is.null(cd_ids) && length(cd_ids) > 0L) {
    keep <- cell_ids %in% cd_ids
    if (!all(keep)) {
      emb <- emb[keep, , drop = FALSE]
      cell_ids <- rownames(emb)
    }
  }
  n_cells <- nrow(emb)
  if (n_cells < 2L) stop("Need at least 2 cells to impute.")
  if (nn > n_cells) {
    .ic2_step(rep, paste0("nn = ", nn, " > n_cells; using nn = ", n_cells))
    nn <- n_cells
  }

  dims <- as.integer(dims)
  dims <- dims[is.finite(dims) & dims >= 1L & dims <= ncol(emb)]
  if (length(dims) < 1L) {
    stop(
      "No valid 'dims' for embedding (have ", ncol(emb),
      " columns). Classic default is dims = 1:10 on PCA."
    )
  }
  emb <- emb[, dims, drop = FALSE]

  ## Order cells by first embedding axis so consecutive chunks share more
  ## kNN neighbors (fewer HDF5 gathers). Assay colnames follow this order.
  if (n_cells > 1L) {
    ord <- order(emb[, 1L], method = "radix")
    if (!isTRUE(all(ord == seq_len(n_cells)))) {
      emb <- emb[ord, , drop = FALSE]
      cell_ids <- rownames(emb)
    }
  }

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  rep <- .ic2_register_assay_inputs(
    rep, x, modality, assay.name, label = paste0("Assay[", assay.name, "]")
  )

  feature_universe <- character(0)
  for (h in handles) {
    h2 <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    feature_universe <- union(feature_universe, as.character(h2$feature_ids))
  }
  feature_universe <- as.character(feature_universe)
  if (length(feature_universe) < 1L) stop("No features found in source assay.")

  ## Ensure every embedding cell is present in the assay
  assay_cells <- character(0)
  for (h in handles) {
    h2 <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    assay_cells <- union(assay_cells, as.character(h2$cell_ids))
  }
  missing_cells <- setdiff(cell_ids, assay_cells)
  if (length(missing_cells) > 0L) {
    stop(
      length(missing_cells), " embedding cells are missing from assay '",
      assay.name, "' (e.g. ", missing_cells[[1L]], ")."
    )
  }

  rep <- .ic2_banner_start(rep, c(
    paste0("Method: iCellR.imp   |   nn = ", nn, " (includes self)"),
    paste0("Reduction: ", reduction, "   |   dims: ", paste(range(dims), collapse = ":")),
    paste0(
      "Assay: ", assay.name, "  ->  ", new.assay.name,
      "   |   cells: ", .ic2_fmt_num(n_cells),
      "   |   features: ", .ic2_fmt_num(length(feature_universe))
    )
  ))

  .ic2_step(rep, paste0("[2/4] Building kNN (FNN, k = ", nn - 1L, " + self)"))
  knn_idx <- .ic2_impute_knn_index(emb, nn = nn, dist.method = knn.dist.method)

  ## Persist knn for audit / future lazy path
  graph_dir <- file.path(x@storage$root_dir, "graphs", paste0("impute_", new.assay.name))
  dir.create(graph_dir, recursive = TRUE, showWarnings = FALSE)
  knn_path <- file.path(graph_dir, "knn_index.h5")
  .ic2_h5_write_dense(knn_idx, knn_path, group = "data")
  .ic2_soda_write_json(
    list(
      method = "iCellR.imp",
      reduction = reduction,
      dims = as.integer(dims),
      nn = as.integer(nn),
      n_cells = as.integer(n_cells),
      cell_ids_path = "cell_ids.tsv",
      knn_path = basename(knn_path),
      assay.name = assay.name,
      new.assay.name = new.assay.name
    ),
    file.path(graph_dir, "impute.json")
  )
  utils::write.table(
    data.frame(cell_id = cell_ids, stringsAsFactors = FALSE),
    file = file.path(graph_dir, "cell_ids.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  dest <- .ic2_soda_assay_path(
    x, modality, new.assay.name, sample_id = NULL, ext = ".h5"
  )
  if (file.exists(dest) && !overwrite) {
    stop("Imputed assay exists at '", dest, "'. Use overwrite = TRUE.")
  }

  .ic2_step(rep, paste0(
    "[3/4] Streaming impute  ->  CSC (single pass, ",
    .ic2_fmt_num(chunk.cells), " cells/chunk)"
  ))
  handle <- .ic2_impute_write_csc(
    handles = handles,
    cell_ids = cell_ids,
    knn_idx = knn_idx,
    feature_universe = feature_universe,
    dest_path = dest,
    chunk.cells = chunk.cells,
    round.num = round.num,
    rounding.digits = as.integer(rounding.digits),
    drop.zero.genes = drop.zero.genes,
    progress = function(cur, tot) {
      rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  impute")
    }
  )
  handle$modality <- modality
  handle$sample_id <- "imputed"

  .ic2_step(rep, "[4/4] Registering imputed assay in SODA")
  feat_out <- handle$feature_ids
  if (is.null(feat_out) || length(feat_out) < 1L) {
    h2 <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
    feat_out <- h2$feature_ids
  }
  x <- .ic2_assign_assay(
    x,
    modality = modality,
    assay.name = new.assay.name,
    handles = list(handle),
    virtual = FALSE,
    feature_universe = NULL,
    set.active = set.active
  )

  x@params$imputation <- x@params$imputation %||% list()
  x@params$imputation[[modality]] <- list(
    assay.name = assay.name,
    new.assay.name = new.assay.name,
    imp.method = "iCellR.imp",
    reduction = reduction,
    based_on = reduction,
    dims = as.integer(dims),
    nn = as.integer(nn),
    aligned = isTRUE(src$aligned),
    mark = src$mark,
    n_cells = as.integer(n_cells),
    n_features = length(feat_out),
    knn_dir = graph_dir,
    chunk.cells = as.integer(chunk.cells)
  )

  x <- .ic2_log_command(
    x,
    "run.impute",
    arguments = list(
      modality = modality,
      assay.name = assay.name,
      new.assay.name = new.assay.name,
      imp.method = "iCellR.imp",
      reduction = reduction,
      dims = as.integer(dims),
      nn = as.integer(nn),
      chunk.cells = as.integer(chunk.cells)
    ),
    outputs = list(
      assay = new.assay.name,
      path = dest,
      knn_dir = graph_dir,
      n_cells = as.integer(n_cells)
    )
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  rep <- .ic2_add_output(rep, "Imputed assay", dest)
  rep <- .ic2_add_output(rep, "kNN graph", graph_dir)
  rep <- .ic2_add_scale(rep, c(
    paste0("Method:   iCellR.imp (streamed kNN mean, nn = ", nn, ")"),
    paste0("Cells:    ", .ic2_fmt_num(n_cells)),
    paste0("Features: ", .ic2_fmt_num(length(feat_out))),
    paste0("Assay:    ", assay.name, "  ->  ", new.assay.name),
    "Memory:   single-pass gather  x  sparse kNN mean (no full densify)"
  ))
  .ic2_finish(rep, extra_lines = paste0(
    "Imputed assay: ", .ic2_fmt_num(length(feat_out)), "  x  ", .ic2_fmt_num(n_cells),
    "   |   nn = ", nn, "   |   ", new.assay.name
  ))
  x
}


#' kNN index including self as neighbor 1 (classic iCellR parity)
#' @keywords internal
.ic2_impute_knn_index <- function(emb, nn = 10L, dist.method = "euclidean") {
  emb <- as.matrix(emb)
  n <- nrow(emb)
  nn <- as.integer(min(as.integer(nn), n))
  if (nn < 1L) stop("'nn' must be >= 1.")

  if (nn == 1L) {
    return(matrix(seq_len(n), nrow = n, ncol = 1L))
  }

  k_other <- nn - 1L
  if (requireNamespace("FNN", quietly = TRUE) && identical(dist.method, "euclidean")) {
    kn <- FNN::get.knn(emb, k = k_other)
    idx <- cbind(seq_len(n), kn$nn.index)
  } else {
    D <- as.matrix(stats::dist(emb, method = dist.method))
    idx <- matrix(0L, nrow = n, ncol = nn)
    for (i in seq_len(n)) {
      idx[i, ] <- as.integer(order(D[i, ])[seq_len(nn)])
    }
  }
  storage.mode(idx) <- "integer"
  idx
}


#' Gather full assay columns for selected cell ids (feature-aligned)
#' @keywords internal
.ic2_gather_assay_cols <- function(handles, cell_ids, feature_universe, chunk.size = 5000L) {
  cell_ids <- as.character(cell_ids)
  feature_universe <- as.character(feature_universe)
  if (length(cell_ids) < 1L) stop("No cells to gather.")

  blocks <- vector("list", 0L)
  remaining <- cell_ids

  for (h in handles) {
    if (length(remaining) < 1L) break
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    local <- intersect(remaining, h$cell_ids)
    if (length(local) < 1L) next
    col_idx <- match(local, h$cell_ids)
    ord <- order(col_idx)
    local <- local[ord]
    col_idx <- col_idx[ord]
    starts <- seq.int(1L, length(col_idx), by = as.integer(chunk.size))
    same_features <- !is.null(h$feature_ids) &&
      length(h$feature_ids) == length(feature_universe) &&
      identical(as.character(h$feature_ids), feature_universe)
    for (s in starts) {
      e <- min(s + as.integer(chunk.size) - 1L, length(col_idx))
      raw <- .ic2_read_cols(h, col_idx[s:e])
      mat <- .ic2_to_dgCMatrix(raw)
      if (!same_features) {
        ## Fast path when rownames already match (common for single-sample SODA)
        rn <- rownames(mat)
        if (is.null(rn) || length(rn) != length(feature_universe) ||
            !identical(as.character(rn), feature_universe)) {
          mat <- .ic2_expand_sparse_to_universe(mat, feature_universe)
        }
      }
      colnames(mat) <- local[s:e]
      blocks[[length(blocks) + 1L]] <- mat
    }
    remaining <- setdiff(remaining, local)
  }
  if (length(remaining) > 0L) {
    stop(
      "Could not gather ", length(remaining), " cells from assay handles ",
      "(e.g. ", remaining[[1L]], ")."
    )
  }
  if (length(blocks) == 1L) {
    out <- blocks[[1L]]
  } else {
    out <- do.call(cbind, blocks)
  }
  ## Reorder to requested cell_ids
  out <- out[, cell_ids, drop = FALSE]
  out
}


#' Stream kNN-mean impute into a CSC HDF5 assay (single pass).
#'
#' Uses growable HDF5 datasets so we do not need a counting pass (the old
#' two-pass path roughly doubled wall time on large atlases).
#' @keywords internal
.ic2_impute_write_csc <- function(
    handles,
    cell_ids,
    knn_idx,
    feature_universe,
    dest_path,
    group = "counts",
    chunk.cells = 30000L,
    round.num = TRUE,
    rounding.digits = 3L,
    drop.zero.genes = TRUE,
    progress = NULL
) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  cell_ids <- as.character(cell_ids)
  feature_universe <- as.character(feature_universe)
  n_cells <- length(cell_ids)
  n_feat <- length(feature_universe)
  if (nrow(knn_idx) != n_cells) {
    stop("knn_idx rows (", nrow(knn_idx), ") != n_cells (", n_cells, ").")
  }
  nn <- ncol(knn_idx)
  chunk.cells <- as.integer(chunk.cells)
  if (!is.finite(chunk.cells) || chunk.cells < 1L) chunk.cells <- 30000L
  w <- 1 / as.numeric(nn)

  if (file.exists(dest_path)) unlink(dest_path)
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)

  row_keep <- if (isTRUE(drop.zero.genes)) {
    rep(FALSE, n_feat)
  } else {
    NULL
  }

  ## Growable CSC sinks  -  append each target chunk once (no counting pass).
  h5 <- hdf5r::H5File$new(dest_path, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5$create_group(group)
  g[["format"]] <- "iCellR2_csc_v1"
  ## dims / rownames finalized after optional drop.zero.genes
  g[["colnames"]] <- cell_ids

  chunk_nx <- as.numeric(5e6)
  space_i <- hdf5r::H5S$new(dims = 0, maxdims = Inf)
  space_x <- hdf5r::H5S$new(dims = 0, maxdims = Inf)
  ds_i <- g$create_dataset(
    "i",
    robj = NULL,
    dtype = hdf5r::h5types$H5T_NATIVE_INT,
    space = space_i,
    chunk_dims = chunk_nx
  )
  ds_x <- g$create_dataset(
    "x",
    robj = NULL,
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    space = space_x,
    chunk_dims = chunk_nx
  )
  ds_p <- g$create_dataset(
    "p",
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = n_cells + 1L,
    chunk_dims = as.integer(min(n_cells + 1L, 1e6L))
  )
  .ic2_h5_write_1d(ds_p, start1 = 1, values = 0)

  out_nnz <- 0
  out_col <- 0L
  starts <- seq.int(1L, n_cells, by = chunk.cells)
  ## Reuse neighbor columns across consecutive chunks (PC1-ordered cells
  ## share many kNN neighbors  ->  large gather savings).
  cache_ids <- character(0)
  cache_M <- NULL

  for (s in starts) {
    e <- min(s + chunk.cells - 1L, n_cells)
    tgt <- s:e
    n_tgt <- length(tgt)
    neigh_i <- unique(as.integer(knn_idx[tgt, , drop = FALSE]))
    neigh_ids <- cell_ids[neigh_i]
    hit <- intersect(neigh_ids, cache_ids)
    miss <- setdiff(neigh_ids, cache_ids)
    parts <- list()
    if (length(hit) > 0L && !is.null(cache_M)) {
      parts[[length(parts) + 1L]] <- cache_M[, hit, drop = FALSE]
    }
    if (length(miss) > 0L) {
      parts[[length(parts) + 1L]] <- .ic2_gather_assay_cols(
        handles, miss, feature_universe,
        chunk.size = max(as.integer(chunk.cells), 5000L)
      )
    }
    if (length(parts) == 1L) {
      M <- parts[[1L]]
    } else if (length(parts) > 1L) {
      M <- do.call(cbind, parts)
    } else {
      stop("Internal error: empty neighbor gather.")
    }
    M <- M[, neigh_ids, drop = FALSE]
    cache_M <- M
    cache_ids <- neigh_ids

    ## Build sparse weight matrix neighbors  x  targets.
    ## knn rows are targets; R stores matrices column-major, so flatten
    ## with t(nb_mat) (or use times=nn with as.integer(nb_mat)) so each
    ## target is paired with its own nn neighbors \u2014 matching classic
    ## rowMeans(my.data[, knn[i, ]]). Using as.integer(nb_mat) with
    ## each=nn scrambled neighbors across targets (MS4A1 "bleed").
    nb_mat <- knn_idx[tgt, , drop = FALSE]
    loc <- match(cell_ids[as.integer(t(nb_mat))], colnames(M))
    if (anyNA(loc)) {
      stop("Internal error: neighbor cell missing from gathered assay block.")
    }
    W <- Matrix::sparseMatrix(
      i = as.integer(loc),
      j = as.integer(rep(seq_len(n_tgt), each = nn)),
      x = rep(w, length.out = length(loc)),
      dims = c(ncol(M), n_tgt)
    )
    imp <- .ic2_to_dgCMatrix(M %*% W)
    rm(W)
    if (isTRUE(round.num) && length(imp@x)) {
      imp@x <- round(imp@x, digits = rounding.digits)
      imp <- Matrix::drop0(imp)
    }
    if (!is.null(row_keep) && length(imp@x)) {
      row_keep[imp@i + 1L] <- TRUE
    }

    bi <- as.integer(imp@i)
    bx <- as.numeric(imp@x)
    bp <- as.numeric(imp@p)
    nadd <- length(bx)
    if (nadd > 0L) {
      new_nnz <- out_nnz + nadd
      ds_i$set_extent(new_nnz)
      ds_x$set_extent(new_nnz)
      .ic2_h5_write_1d(ds_i, start1 = out_nnz + 1, values = bi)
      .ic2_h5_write_1d(ds_x, start1 = out_nnz + 1, values = bx)
    }
    .ic2_h5_write_1d(
      ds_p,
      start1 = out_col + 2L,
      values = out_nnz + bp[-1L]
    )
    out_nnz <- out_nnz + nadd
    out_col <- out_col + n_tgt
    if (is.function(progress)) progress(out_col, n_cells)
    rm(imp, bi, bx, bp)
  }
  rm(cache_M)

  if (out_col != n_cells) {
    stop("Internal error: wrote ", out_col, " columns, expected ", n_cells, ".")
  }

  ## Optional gene filter: remap row indices in-place when some features never
  ## appear (avoids a second full impute pass).
  feature_write <- feature_universe
  if (isTRUE(drop.zero.genes) && !is.null(row_keep) && !all(row_keep)) {
    feature_write <- feature_universe[row_keep]
    if (length(feature_write) < 1L) {
      stop("All features were zero after imputation; check input assay.")
    }
    if (length(feature_write) < n_feat && out_nnz > 0) {
      old_to_new <- rep(NA_integer_, n_feat)
      old_to_new[row_keep] <- seq_along(feature_write) - 1L ## 0-based for CSC @i
      ## Rewrite i in chunks
      pos <- 0
      while (pos < out_nnz) {
        take <- min(as.integer(chunk_nx), out_nnz - pos)
        ii <- as.integer(.ic2_h5_read_1d(ds_i, start1 = pos + 1, count = take))
        mapped <- old_to_new[ii + 1L]
        if (anyNA(mapped)) {
          stop("Internal error: imputed row index not in kept feature set.")
        }
        .ic2_h5_write_1d(ds_i, start1 = pos + 1, values = as.integer(mapped))
        pos <- pos + take
      }
    }
  }

  g[["dims"]] <- as.integer(c(length(feature_write), n_cells))
  g[["rownames"]] <- feature_write
  g[["nnz"]] <- as.numeric(out_nnz)
  h5$flush()
  h5$close_all()
  on.exit(NULL)

  .ic2_assay_handle(
    path = dest_path,
    backend = "hdf5_sparse",
    group = group,
    n_features = length(feature_write),
    n_cells = n_cells,
    feature_ids = feature_write,
    cell_ids = cell_ids
  )
}

