#' Scale / transform an assay for dimensionality reduction
#'
#' Default matches \code{\link{run.pca}} feature prep:
#' \code{log2(x+1)}  ->  center  ->  sample SD \eqn{/ (n-1)}  ->  no clip, on
#' \strong{all assay features}, written as \code{scale.data}. Then
#' \code{run.pca} skips re-scaling (\code{transform = "none"},
#' \code{center = FALSE}, \code{scale. = FALSE}) when this assay was
#' z-scored. For large datasets prefer
#' \code{features = gene.model} (PCA still uses the gene model).
#'
#' Use \code{method = "log2"} / \code{"log1p"} for log-only (z-score in PCA).
#' Pass \code{clip = 10} for Seurat-like clipping after z-score.
#'
#' @param x An \code{iCellR2} object.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Source assay; default active (usually normalized
#'   \code{"data"}).
#' @param new.assay.name Output assay name. Default \code{"scale.data"}.
#' @param method \code{"log2_z"} (default: log2 + center + sample SD),
#'   \code{"log1p_z"}, or log-only \code{"log2"} / \code{"log1p"}.
#' @param center,scale Center / scale features (z-score). Only used for
#'   \code{*_z} methods. Ignored for \code{log2} / \code{log1p}.
#' @param features Feature ids to transform. If \code{NULL} (default),
#'   scales \strong{all} features in the assay. Pass a gene model
#'   (e.g. \code{features = x@misc$gene.model}) to restrict to HVGs.
#' @param clip Absolute clip after z-scoring (e.g. \code{10}). \code{NULL}
#'   (default) disables clipping. Set \code{10} for Seurat-like clip.
#' @param chunk.size Cells per chunk. Default \code{200000}.
#' @param set.active Set the output assay active. Default \code{FALSE}.
#' @param overwrite Overwrite existing scaled assay files.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} object.
#'
#' @examples
#' \dontrun{
#' obj <- data.scale(obj)   # all genes; log2 + center + sample SD
#' obj <- data.scale(obj, features = obj@misc$gene.model)  # HVGs only
#' obj <- run.pca(obj)      # uses scale.data; skips re-scale if z-scored
#' }
#'
#' @export
data.scale <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    new.assay.name = "scale.data",
    method = c("log2_z", "log1p_z", "log2", "log1p"),
    center = TRUE,
    scale = TRUE,
    features = NULL,
    clip = NULL,
    chunk.size = 200000,
    set.active = FALSE,
    overwrite = TRUE,
    verbose = TRUE
) {
  .ic2_require_object(x)
  method <- match.arg(method)
  modality <- tolower(modality)
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  chunk.size <- as.integer(chunk.size)
  rep <- .ic2_new_reporter("data.scale", verbose, soda_root = x@storage$root_dir)

  log_only <- method %in% c("log2", "log1p")
  if (log_only) {
    center <- FALSE
    scale <- FALSE
    clip <- NULL
  }
  ## Default: no clip (match run.pca). Pass clip = 10 for Seurat-like.

  ## Resolve features (default = all assay features)
  features <- if (!is.null(features)) unique(as.character(features)) else NULL

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  rep <- .ic2_register_assay_inputs(rep, x, modality, assay.name, label = paste0("Assay[", assay.name, "]"))
  is_virtual <- .ic2_assay_is_virtual(x, modality, assay.name)

  ## Peek dims / feature universe from first handle
  h0 <- .ic2_enrich_handle(handles[[1L]], load_dimnames = TRUE)
  if (is.null(features)) {
    features <- as.character(h0$feature_ids)
    .ic2_step(rep, paste0(
      "Scaling all ", .ic2_fmt_num(length(features)),
      " features (pass features = gene.model to restrict)."
    ))
  } else {
    missing <- setdiff(features, h0$feature_ids)
    features <- intersect(features, h0$feature_ids)
    if (length(features) < 1L) stop("None of the requested features were found in the assay.")
    if (length(missing) > 0L) {
      .ic2_step(rep, paste0("  Dropping ", length(missing), " features not in assay."))
    }
  }

  n_cells_guess <- sum(vapply(handles, function(h) {
    as.integer(.ic2_enrich_handle(h, FALSE)$n_cells %||% 0L)
  }, integer(1)))

  ## Dense scale.data for all genes  x  large n can be huge  -  warn early
  approx_gb <- (as.numeric(length(features)) * as.numeric(n_cells_guess) * 8) / (1024^3)
  if (length(features) > 5000L && approx_gb > 8) {
    warning(
      "data.scale will write a dense matrix of ~",
      .ic2_fmt_num(length(features)), " features  x  ",
      .ic2_fmt_num(n_cells_guess), " cells (~",
      sprintf("%.1f", approx_gb), " GB). ",
      "For large SODAs prefer features = gene.model (PCA only needs HVGs).",
      call. = FALSE
    )
  }

  rep <- .ic2_banner_start(rep, c(
    paste0("Method: ", method, " | center=", center, " | scale=", scale),
    paste0("Assay: ", assay.name, "  ->  ", new.assay.name),
    paste0(
      "Features: ", .ic2_fmt_num(length(features)),
      " | Cells: ", .ic2_fmt_num(n_cells_guess)
    )
  ))

  log_fun <- if (method %in% c("log2", "log2_z")) {
    function(m) log2(m + 1)
  } else {
    function(m) log1p(m)
  }

  ## Pass 1: feature means / vars across all cells (for z-score)
  feat_mean <- setNames(numeric(length(features)), features)
  feat_m2 <- setNames(numeric(length(features)), features) ## sum of squares of log values
  n_tot <- 0L

  if (isTRUE(center) || isTRUE(scale)) {
    .ic2_step(rep, "[1/2] Estimating feature mean / variance (streaming)")
    for (hi in seq_along(handles)) {
      h <- .ic2_enrich_handle(handles[[hi]], load_dimnames = TRUE)
      row_idx <- match(features, h$feature_ids)
      keep_f <- !is.na(row_idx)
      row_idx <- row_idx[keep_f]
      feat_local <- features[keep_f]
      n_cells <- h$n_cells
      starts <- seq.int(1L, n_cells, by = chunk.size)
      for (si in seq_along(starts)) {
        start <- starts[si]
        end <- min(start + chunk.size - 1L, n_cells)
        block <- .ic2_read_cols(h, start:end)
        block <- block[row_idx, , drop = FALSE]
        block <- log_fun(block)
        ## densify small feature block only (HVGs  x  chunk)  -  intentional
        mat <- as.matrix(block)
        feat_mean[feat_local] <- feat_mean[feat_local] + rowSums(mat)
        feat_m2[feat_local] <- feat_m2[feat_local] + rowSums(mat * mat)
        n_tot <- n_tot + ncol(mat)
        rep <- .ic2_progress_bar(
          rep,
          (hi - 1L) * n_cells + end,
          n_cells_guess,
          prefix = "  stats"
        )
      }
    }
    n <- as.numeric(n_tot)
    feat_mean <- feat_mean / pmax(n, 1)
    ## Sample SD (n-1), same convention as run.pca / prcomp
    feat_var <- pmax(0, (feat_m2 - n * feat_mean * feat_mean) / pmax(n - 1, 1))
    feat_sd <- sqrt(feat_var)
    feat_sd[!is.finite(feat_sd) | feat_sd < 1e-8] <- 1
  } else {
    .ic2_step(rep, "[1/2] No centering/scaling  -  log transform only")
    feat_mean <- setNames(rep(0, length(features)), features)
    feat_sd <- setNames(rep(1, length(features)), features)
  }

  ## Pass 2: write scaled assay per handle
  .ic2_step(rep, "[2/2] Writing scaled assay to HDF5")
  new_handles <- vector("list", length(handles))
  dest_paths <- character(0)

  for (hi in seq_along(handles)) {
    h <- .ic2_enrich_handle(handles[[hi]], load_dimnames = TRUE)
    row_idx <- match(features, h$feature_ids)
    keep_f <- !is.na(row_idx)
    row_idx <- row_idx[keep_f]
    feat_local <- features[keep_f]
    mu <- feat_mean[feat_local]
    sdv <- feat_sd[feat_local]

    sid <- if (!is.null(h$sample_id) && nzchar(as.character(h$sample_id))) {
      h$sample_id
    } else if (is_virtual) {
      paste0("sample", hi)
    } else NULL

    dest <- .ic2_soda_assay_path(x, modality, new.assay.name, sample_id = sid, ext = ".h5")
    if (file.exists(dest) && !overwrite) {
      stop("Scaled assay exists at '", dest, "'. Use overwrite = TRUE.")
    }

    blocks <- list()
    n_cells <- h$n_cells
    starts <- seq.int(1L, n_cells, by = chunk.size)
    for (si in seq_along(starts)) {
      start <- starts[si]
      end <- min(start + chunk.size - 1L, n_cells)
      block <- .ic2_read_cols(h, start:end)[row_idx, , drop = FALSE]
      block <- log_fun(block)
      mat <- as.matrix(block)
      if (isTRUE(center)) {
        mat <- sweep(mat, 1L, mu, "-")
      }
      if (isTRUE(scale)) {
        mat <- sweep(mat, 1L, sdv, "/")
      }
      if (!is.null(clip) && is.finite(clip)) {
        mat[mat > clip] <- clip
        mat[mat < -clip] <- -clip
      }
      rownames(mat) <- feat_local
      colnames(mat) <- colnames(block)
      ## Dense on disk  -  scaled values are nearly dense; CSC drop0 zeros them out
      blocks[[length(blocks) + 1L]] <- mat
      rep <- .ic2_progress_bar(
        rep,
        (hi - 1L) * n_cells + end,
        n_cells_guess,
        prefix = "  scale"
      )
    }
    mat_all <- do.call(cbind, blocks)
    storage.mode(mat_all) <- "double"
    rm(blocks)
    invisible(gc)
    new_handles[[hi]] <- write_dense_hdf5(mat_all, dest, group = "counts", overwrite = TRUE)
    new_handles[[hi]]$sample_id <- sid
    new_handles[[hi]]$modality <- modality
    new_handles[[hi]]$feature_ids <- NULL
    new_handles[[hi]]$cell_ids <- NULL

    ## Verify round-trip (catches hdf5r flatten bugs immediately)
    probe <- read_dense_hdf5_cols(dest, group = "counts", cols = seq_len(min(8L, ncol(mat_all))))
    if (max(abs(probe), na.rm = TRUE) == 0 && max(abs(mat_all), na.rm = TRUE) > 0) {
      stop(
        "Dense scale write verification failed for '", dest,
        "' (on-disk block is all zeros). Please report this hdf5r issue."
      )
    }
    dest_paths <- c(dest_paths, dest)
  }

  x <- .ic2_assign_assay(
    x,
    modality = modality,
    assay.name = new.assay.name,
    handles = new_handles,
    virtual = is_virtual,
    feature_universe = features,
    set.active = set.active
  )

  ## store scaling params (compact JSON  -  no .rds)
  scale_path <- file.path(
    x@storage$root_dir, "modalities", modality, "assays", new.assay.name, "scale_params.json"
  )
  dir.create(dirname(scale_path), recursive = TRUE, showWarnings = FALSE)
  .ic2_soda_write_json(
    list(
      features = features,
      mean = feat_mean,
      sd = feat_sd,
      method = method,
      center = isTRUE(center),
      scale = isTRUE(scale),
      sd_type = if (isTRUE(scale)) "sample" else "none",
      clip = clip
    ),
    scale_path
  )

  x@params$normalization[[paste0(modality, "_scale")]] <- list(
    assay.name = assay.name,
    new.assay.name = new.assay.name,
    method = method,
    center = center,
    scale = scale,
    n_features = length(features),
    clip = clip,
    params_path = scale_path
  )

  x <- .ic2_log_command(
    x, "data.scale",
    arguments = list(
      modality = modality, assay.name = assay.name,
      new.assay.name = new.assay.name, method = method,
      n_features = length(features)
    ),
    outputs = list(paths = dest_paths, params_path = scale_path)
  )

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  for (p in dest_paths) rep <- .ic2_add_output(rep, "Scaled assay", p)
  rep <- .ic2_add_output(rep, "Scale params", scale_path)

  wall <- as.numeric(difftime(Sys.time(), rep$wall0, units = "secs"))
  if (is.finite(wall) && wall > 0) {
    rep$stats$cells_per_sec <- n_cells_guess / wall
  }
  .ic2_finish(rep, extra_lines = paste0(
    if (log_only) "Log-transformed " else "Scaled ",
    .ic2_fmt_num(length(features)), " features  x  ",
    .ic2_fmt_num(n_cells_guess), " cells (", method, ")",
    if (log_only) {
      "  -  z-score in run.pca()"
    } else {
      "  -  log2/log1p + center + sample SD; run.pca skips re-scale"
    }
  ))
  x
}
