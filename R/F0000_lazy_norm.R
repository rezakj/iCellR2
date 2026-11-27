#' Lazy normalization helpers
#'
#' Size-factor-only assays (\code{aggregation_type = "lazy_norm"}) keep raw
#' counts on disk and apply scaling on read.
#'
#' @name lazy_norm
#' @keywords internal
NULL

#' Read a size-factor TSV (cell_id, size_factor) to a named numeric vector
#' @keywords internal
.ic2_read_size_factors <- function(path) {
  if (is.null(path) || !nzchar(as.character(path)[1L]) || !file.exists(path)) {
    return(NULL)
  }
  df <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(df)) return(NULL)
  id_col <- if ("cell_id" %in% names(df)) "cell_id" else names(df)[[1L]]
  sf_col <- if ("size_factor" %in% names(df)) {
    "size_factor"
  } else if (ncol(df) >= 2L) {
    names(df)[[2L]]
  } else {
    stop("size_factors table needs a size_factor column: ", path)
  }
  sf <- as.numeric(df[[sf_col]])
  names(sf) <- as.character(df[[id_col]])
  sf
}

#' Size factors for a lazy_norm assay (NULL if already materialized)
#' @keywords internal
.ic2_assay_size_factors <- function(x, modality, assay.name) {
  modality <- tolower(modality)
  ref <- x@modalities[[modality]]$assays[[assay.name]]
  if (!is.list(ref) || !identical(ref$aggregation_type, "lazy_norm")) {
    return(NULL)
  }
  path <- ref$size_factors_path %||%
    x@params$normalization[[modality]]$size_factors_path
  .ic2_read_size_factors(path)
}

#' Is this assay a size-factor-only lazy norm?
#' @keywords internal
.ic2_assay_is_lazy_norm <- function(x, modality, assay.name) {
  ref <- x@modalities[[modality]]$assays[[assay.name]]
  is.list(ref) && identical(ref$aggregation_type, "lazy_norm")
}

#' Divide columns of a features x cells block by size factors
#' @keywords internal
.ic2_scale_cols_by_sf <- function(block, size_factors) {
  if (is.null(size_factors) || length(size_factors) < 1L) return(block)
  cn <- colnames(block)
  if (is.null(cn)) stop("Block must have colnames to apply size factors.")
  sf <- as.numeric(size_factors[cn])
  if (anyNA(sf) || any(!is.finite(sf) | sf <= 0)) {
    stop("Missing or non-positive size factors for ", sum(is.na(sf) | !is.finite(sf) | sf <= 0), " cells.")
  }
  if (inherits(block, "Matrix")) {
    block <- block %*% Matrix::Diagonal(x = 1 / sf)
    return(.ic2_to_dgCMatrix(block))
  }
  sweep(as.matrix(block), 2L, sf, "/")
}

#' Register a lazy_norm assay (no normalized HDF5)
#' @keywords internal
.ic2_assign_lazy_norm_assay <- function(
    x,
    modality,
    assay.name,
    source_assay,
    size_factors_path,
    feature_universe,
    n_cells,
    set.active = TRUE
) {
  assay_ref <- list(
    aggregation_type = "lazy_norm",
    source_assay = source_assay,
    size_factors_path = normalizePath(size_factors_path, winslash = "/", mustWork = FALSE),
    feature_universe = as.character(feature_universe),
    n_cells = as.integer(n_cells),
    n_features = length(feature_universe)
  )
  x@modalities[[modality]]$assays[[assay.name]] <- assay_ref
  x@storage$assay_paths <- list(
    modality = modality,
    assay = assay.name,
    aggregation = "lazy_norm",
    source_assay = source_assay,
    size_factors_path = assay_ref$size_factors_path
  )
  ## Persist enough for open.soda to rebuild the lazy assay (SF-only dir)
  assay_dir <- dirname(assay_ref$size_factors_path)
  tryCatch(
    .ic2_soda_write_json(
      list(
        aggregation_type = "lazy_norm",
        source_assay = source_assay,
        size_factors_path = basename(assay_ref$size_factors_path),
        n_cells = as.integer(n_cells),
        n_features = length(feature_universe),
        ## Keep universe for PCA / markers without re-scanning counts
        feature_universe = as.character(feature_universe)
      ),
      file.path(assay_dir, "assay.json")
    ),
    error = function(e) invisible(NULL)
  )
  if (isTRUE(set.active)) {
    x@modalities[[modality]]$active_assay <- assay.name
  }
  x
}

#' Can we materialize by copy + in-place scale of a native 10x H5?
#' @keywords internal
.ic2_tenx_inplace_norm_eligible <- function(handle, norm.method, feature_universe = NULL) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  if (!identical(handle$backend, "tenx_h5")) return(FALSE)
  if (!norm.method %in% c("global.glsf", "rpm", "LogNormalize", "no.norm")) {
    return(FALSE)
  }
  if (is.null(handle$meta_path) || !file.exists(handle$meta_path)) return(FALSE)
  if (is.null(handle$path) || !file.exists(handle$path)) return(FALSE)
  meta <- .ic2_soda_read_tenx_meta(handle$meta_path)
  if (is.null(meta)) return(FALSE)
  if (!is.null(meta$row_keep) && length(meta$row_keep)) return(FALSE)
  grp <- as.character(meta$group %||% "matrix")[1L]
  if (!identical(grp, "matrix")) return(FALSE)
  fu <- as.character(feature_universe %||% character(0))
  if (length(fu) > 0L) {
    feats <- as.character(handle$feature_ids %||% meta$feature_ids)
    if (length(fu) != length(feats) || any(is.na(match(feats, fu)))) {
      return(FALSE)
    }
  }
  TRUE
}

#' Fast materialize: copy 10x H5 and divide matrix/data by size factors in place
#'
#' Mirrors rawdata/materialize/normalize_10x_h5_chunked.R using hdf5r.
#' Cells without a size factor are left unchanged (divisor 1). Optional
#' \code{log1p} for LogNormalize materialize.
#'
#' @keywords internal
.ic2_materialize_tenx_inplace <- function(
    handle,
    dest_path,
    size_factors,
    chunk.size = 100000L,
    transform = c("none", "log1p", "log2p1"),
    round.num = FALSE,
    rounding.digits = 3,
    progress = NULL,
    overwrite = TRUE
) {
  transform <- match.arg(transform)
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  if (!identical(handle$backend, "tenx_h5")) {
    stop("Expected tenx_h5 handle for in-place normalize.")
  }
  src <- normalizePath(handle$path, winslash = "/", mustWork = TRUE)
  dest_path <- as.character(dest_path)[1L]
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest_path)) {
    if (!isTRUE(overwrite)) stop("Destination exists: ", dest_path)
    unlink(dest_path)
  }
  if (!identical(src, normalizePath(dest_path, winslash = "/", mustWork = FALSE))) {
    if (!file.copy(src, dest_path, overwrite = TRUE)) {
      stop("Failed to copy '", src, "'  ->  '", dest_path, "'.")
    }
  }
  dest_path <- normalizePath(dest_path, winslash = "/", mustWork = TRUE)

  meta <- .ic2_soda_read_tenx_meta(handle$meta_path)
  if (is.null(meta)) stop("Could not read tenx meta: ", handle$meta_path)
  group <- meta$group %||% "matrix"
  cell_ids <- as.character(handle$cell_ids %||% meta$cell_ids)
  n_cells <- length(cell_ids)
  if (n_cells < 1L) stop("No cells in tenx handle.")

  sf <- rep(1.0, n_cells)
  if (!identical(transform, "none") || length(size_factors) > 0L) {
    if (!is.null(names(size_factors))) {
      m <- match(cell_ids, names(size_factors))
      ok <- !is.na(m)
      sf[ok] <- as.numeric(size_factors[m[ok]])
    } else if (length(size_factors) == n_cells) {
      sf <- as.numeric(size_factors)
    }
  }
  if (any(!is.finite(sf) | sf <= 0)) {
    stop(
      "Non-positive size factors for ",
      sum(!is.finite(sf) | sf <= 0),
      " cells in tenx in-place normalize."
    )
  }

  chunk.size <- as.integer(chunk.size)
  if (!is.finite(chunk.size) || chunk.size < 1L) chunk.size <- 100000L

  h5 <- hdf5r::H5File$new(dest_path, mode = "r+")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  if (!group %in% names(h5)) stop("HDF5 group '", group, "' missing in ", dest_path)
  g <- h5[[group]]
  if (!all(c("data", "indptr") %in% names(g))) {
    stop("Not a modern 10x matrix group (need data + indptr): ", dest_path)
  }
  ds_data <- g[["data"]]
  ds_indptr <- g[["indptr"]]
  ## Full indptr is small relative to nnz (n_cells+1)
  indptr <- as.numeric(ds_indptr[])
  if (length(indptr) != n_cells + 1L) {
    stop(
      "indptr length (", length(indptr), ") != n_cells+1 (", n_cells + 1L,
      "). Cell id / matrix mismatch."
    )
  }

  ## no.norm + no transform: copy only
  need_scale <- !all(abs(sf - 1) < 1e-15)
  need_tf <- !identical(transform, "none")
  need_round <- isTRUE(round.num)
  if (need_scale || need_tf || need_round) {
    starts <- seq.int(1L, n_cells, by = chunk.size)
    for (si in seq_along(starts)) {
      c0 <- starts[si]
      c1 <- min(c0 + chunk.size - 1L, n_cells)
      x_start0 <- indptr[c0]
      x_end0 <- indptr[c1 + 1L]
      n_x <- x_end0 - x_start0
      if (n_x >= 1) {
        block <- as.numeric(.ic2_h5_read_1d(ds_data, start1 = x_start0 + 1, count = n_x))
        local_ptr <- indptr[c0:(c1 + 1L)] - x_start0
        for (j in seq_len(c1 - c0 + 1L)) {
          s <- local_ptr[j] + 1L
          e <- local_ptr[j + 1L]
          if (e >= s) {
            if (need_scale) {
              block[s:e] <- block[s:e] / sf[c0 + j - 1L]
            }
            if (identical(transform, "log1p")) {
              block[s:e] <- log1p(block[s:e])
            } else if (identical(transform, "log2p1")) {
              block[s:e] <- log2(block[s:e] + 1)
            }
            if (need_round) {
              block[s:e] <- round(block[s:e], digits = rounding.digits)
            }
          }
        }
        .ic2_h5_write_1d(ds_data, start1 = x_start0 + 1, values = block)
      }
      if (is.function(progress)) progress(c1, n_cells)
    }
  } else if (is.function(progress)) {
    progress(n_cells, n_cells)
  }
  h5$flush()
  h5$close_all()
  on.exit(NULL)

  ## Write SODA tenx meta (preserve logical cell_ids from source handle)
  meta$source <- dest_path
  meta$cell_ids <- cell_ids
  meta$feature_ids <- as.character(handle$feature_ids %||% meta$feature_ids)
  meta$n_cells <- n_cells
  meta$n_features <- length(meta$feature_ids)
  meta_dest <- file.path(dirname(dest_path), "tenx_meta.json")
  meta_dest <- .ic2_soda_write_tenx_meta(meta, meta_dest)
  .ic2_soda_unlink_rds(file.path(dirname(dest_path), "tenx_meta.rds"))

  h_out <- .ic2_assay_handle(
    path = dest_path,
    backend = "tenx_h5",
    group = group,
    n_features = as.integer(meta$n_features),
    n_cells = as.integer(n_cells),
    feature_ids = NULL,
    cell_ids = NULL,
    sample_id = handle$sample_id,
    modality = handle$modality
  )
  h_out$meta_path <- normalizePath(meta_dest, winslash = "/", mustWork = TRUE)
  h_out
}
