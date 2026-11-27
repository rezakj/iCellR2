#' Run PCA on scaled / normalized features (scalable)
#'
#' Computes principal components from an on-disk assay. Prefer normalized
#' \code{"data"} with on-the-fly \code{log2(x+1)} on the selected features
#' (typically ~2k from \code{gene.model}); storing a full \code{scale.data}
#' is optional.
#'
#' When the feature count is modest (default \eqn{\le} 5000  -  typical
#' gene models), \strong{exact} PCA is computed out-of-core by streaming
#' cells from HDF5 and accumulating a \eqn{p \times p} covariance
#' (\code{pca.method = "cov"} / \code{"auto"}). For wider feature sets,
#' falls back to sketch + truncated SVD (\pkg{irlba}).
#'
#' @param x An \code{iCellR2} object.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Source assay. Default: \code{"scale.data"} if present
#'   (already gene.model-sized), else \code{"data"} (with on-the-fly
#'   \code{log2p1}), else active.
#' @param method Feature selection: \code{"gene.model"} (default),
#'   \code{"base.mean.rank"}, or \code{"features"}.
#' @param top.rank Top features by mean for \code{base.mean.rank}.
#' @param features Explicit feature ids when \code{method = "features"}.
#' @param npcs Number of principal components. Default \code{50}.
#' @param transform Expression transform applied to assay values for the
#'   selected features only: \code{"auto"} (default) uses \code{"log2p1"}
#'   when the assay is \code{"data"} and \code{"none"} for
#'   \code{"scale.data"}; or set \code{"log2p1"}, \code{"log1p"}, \code{"none"}.
#' @param pca.method \code{"auto"} (default): covariance stream when
#'   \code{length(features) <= max.features.cov}, else sketch/direct.
#'   \code{"cov"}: exact out-of-core covariance PCA on all cells.
#'   \code{"sketch"}: sketch + irlba + project.
#'   \code{"direct"}: load all cells and irlba/prcomp.
#' @param max.features.cov Feature-count cutoff for \code{pca.method = "auto"}
#'   to choose covariance PCA. Default \code{5000}.
#' @param use.sketch If \code{TRUE} and sketch path is selected, fit on a
#'   sketch when available (or cells exceed \code{max.cells.direct}).
#' @param sketch.name Sketch name from \code{\link{run.sketch}}.
#' @param max.cells.direct Max cells for direct / sketch-fallback PCA.
#'   Default \code{100000}.
#' @param center,scale. Feature centering / sample-SD scaling. Defaults
#'   \code{TRUE}. When assay is \code{scale.data} already z-scored by
#'   \code{\link{data.scale}} (\code{log2_z} / \code{log1p_z}), both
#'   auto-flip to \code{FALSE} unless you pass them explicitly.
#' @param chunk.size Cells per HDF5 / projection chunk.
#' @param seed RNG seed for auto-sketch.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} with PCA embeddings / loadings on disk.
#'
#' @examples
#' \dontrun{
#' obj <- norm.data(obj)
#' obj <- run.pca(obj, npcs = 50)         # scale.data if present, else data+log2p1
#' }
#'
#' @export
run.pca <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    method = c("gene.model", "base.mean.rank", "features"),
    top.rank = 500,
    features = NULL,
    npcs = 50,
    transform = c("auto", "none", "log2p1", "log1p"),
    pca.method = c("auto", "cov", "sketch", "direct"),
    max.features.cov = 5000L,
    use.sketch = TRUE,
    sketch.name = "sketch1",
    max.cells.direct = 100000L,
    center = TRUE,
    scale. = TRUE,
    chunk.size = 200000L,
    seed = 1,
    verbose = TRUE
) {
  .ic2_require_object(x)
  method <- match.arg(method)
  pca.method <- match.arg(pca.method)
  transform <- match.arg(transform)
  modality <- tolower(modality)
  chunk.size <- as.integer(chunk.size)
  npcs <- as.integer(npcs)
  max.cells.direct <- as.integer(max.cells.direct)
  max.features.cov <- as.integer(max.features.cov)
  rep <- .ic2_new_reporter("run.pca", verbose, soda_root = x@storage$root_dir)

  ## Resolve assay  -  prefer scale.data (skinny gene.model) over full data
  assay.names <- names(x@modalities[[modality]]$assays %||% list())
  if (is.null(assay.name)) {
    if ("scale.data" %in% assay.names) {
      assay.name <- "scale.data"
    } else if ("data" %in% assay.names) {
      assay.name <- "data"
    } else {
      assay.name <- .ic2_resolve_assay_name(x, modality, NULL)
    }
  } else {
    assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  }
  if (identical(transform, "auto")) {
    transform <- if (assay.name %in% c("scale.data", "tfidf")) "none" else "log2p1"
  }
  ## TF-IDF is already weighted. Do not z-score peaks; rare peaks would dominate.
  if (identical(assay.name, "tfidf") && missing(scale.)) scale. <- FALSE

  ## If scale.data was written with center+sample SD, skip re-scaling
  ## unless the caller passed center / scale. explicitly.
  if (identical(assay.name, "scale.data")) {
    scale_path <- file.path(
      x@storage$root_dir, "modalities", modality, "assays",
      assay.name, "scale_params.json"
    )
    sp <- NULL
    if (is.character(scale_path) && file.exists(scale_path)) {
      sp <- tryCatch(.ic2_soda_read_meta(scale_path), error = function(e) NULL)
    }
    if (is.null(sp)) {
      sp <- x@params$normalization[[paste0(modality, "_scale")]]
    }
    already_z <- is.list(sp) && isTRUE(sp$center) && isTRUE(sp$scale) &&
      (
        identical(as.character(sp$method %||% ""), "log2_z") ||
          identical(as.character(sp$method %||% ""), "log1p_z") ||
          identical(as.character(sp$sd_type %||% ""), "sample")
      )
    if (isTRUE(already_z)) {
      if (missing(center)) center <- FALSE
      if (missing(scale.)) scale. <- FALSE
      .ic2_step(rep, "scale.data already z-scored  -  PCA center/scale. = FALSE")
    }
  }

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  size_factors <- .ic2_assay_size_factors(x, modality, assay.name)
  if (!is.null(size_factors)) {
    .ic2_step(rep, paste0(
      "Lazy norm: applying ", .ic2_fmt_num(length(size_factors)),
      " size factors on the fly from '", assay.name, "'"
    ))
  }
  rep <- .ic2_register_assay_inputs(rep, x, modality, assay.name, label = paste0("Assay[", assay.name, "]"))
  assay_ref <- x@modalities[[modality]]$assays[[assay.name]]
  if (is.list(assay_ref) && !is.null(assay_ref$feature_universe) &&
      length(assay_ref$feature_universe) > 0L) {
    all_features <- as.character(assay_ref$feature_universe)
  } else {
    all_features <- character(0)
    for (hh in handles) {
      hh <- .ic2_enrich_handle(hh, load_dimnames = TRUE)
      all_features <- union(all_features, as.character(hh$feature_ids))
    }
  }

  ## Resolve features
  if (identical(method, "gene.model")) {
    features <- x@modalities[[modality]]$gene.model %||% x@misc$gene.model
    if (is.null(features) || length(features) < 2L) {
      ## Portable SODA: metadata/<mod>_gene_model.tsv (Python / new R)
      gm_tsv <- file.path(
        x@storage$root_dir, "metadata", paste0(modality, "_gene_model.tsv")
      )
      gm_rds <- file.path(
        x@storage$root_dir, "metadata", paste0(modality, "_gene_model.rds")
      )
      if (file.exists(gm_tsv)) {
        tab <- utils::read.delim(gm_tsv, sep = "\t", stringsAsFactors = FALSE)
        col <- if ("feature_id" %in% names(tab)) "feature_id" else names(tab)[[1L]]
        features <- as.character(tab[[col]])
      } else if (file.exists(gm_rds)) {
        features <- as.character(readRDS(gm_rds))
      }
    }
    if (is.null(features) || length(features) < 2L) {
      if (identical(modality, "atac")) {
        features <- all_features
        .ic2_step(rep, "ATAC has no gene.model  -  PCA uses every peak in the assay")
      } else {
        stop("No gene.model found. Run make.gene.model() or set method = 'base.mean.rank'.")
      }
    }
    x@modalities[[modality]]$gene.model <- features
    x@misc$gene.model <- features
  } else if (identical(method, "base.mean.rank")) {
    gd <- x@modalities[[modality]]$gene.data
    if (is.null(gd)) {
      sp <- x@storage$metadata_paths$gene_stats
      cands <- c(
        sp,
        file.path(x@storage$root_dir, "metadata", paste0(modality, "_gene_stats.tsv")),
        file.path(x@storage$root_dir, "metadata", paste0(modality, "_gene_stats.rds"))
      )
      for (p in unique(cands[!vapply(cands, function(z) is.null(z) || !nzchar(as.character(z)[1L]), logical(1))])) {
        if (!file.exists(p)) next
        gd <- tryCatch({
          if (grepl("\\.rds$", p, ignore.case = TRUE)) readRDS(p) else .ic2_soda_read_table(p)
        }, error = function(e) NULL)
        if (is.data.frame(gd) && nrow(gd) > 0L) break
      }
    }
    if (is.null(gd)) {
      stop("No gene.stats found for base.mean.rank. Run gene.stats() first.")
    }
    d <- gd[gd$condition == "all", , drop = FALSE]
    features <- as.character(d$genes[order(d$meanExp, decreasing = TRUE)])
    features <- utils::head(features, top.rank)
  } else {
    if (is.null(features) || length(features) < 2L) {
      stop("Provide 'features' when method = 'features'.")
    }
    features <- unique(as.character(features))
  }
  features <- intersect(features, all_features)
  if (length(features) < 2L) stop("Fewer than 2 features available for PCA.")

  cell_ids <- rownames(x@cell_data)
  n_cells <- length(cell_ids)
  npcs <- min(npcs, length(features) - 1L, max(2L, n_cells - 1L))

  ## Resolve engine
  if (identical(pca.method, "auto")) {
    pca.method <- if (length(features) <= max.features.cov) {
      "cov"
    } else if (isTRUE(use.sketch) && n_cells > max.cells.direct) {
      "sketch"
    } else {
      "direct"
    }
  }

  sketch_ids <- NULL
  sketch_mode <- FALSE
  if (identical(pca.method, "sketch")) {
    if (isTRUE(use.sketch)) {
      sk <- x@modalities[[modality]]$sketches$sketch_models[[sketch.name]]
      if (!is.null(sk$path) && file.exists(sk$path)) {
        meta <- readRDS(sk$path)
        sketch_ids <- intersect(meta$cell_ids, cell_ids)
      } else if (!is.null(x@modalities[[modality]]$sketches$sketch_ids)) {
        sketch_ids <- intersect(x@modalities[[modality]]$sketches$sketch_ids, cell_ids)
      }
      if (length(sketch_ids) >= 10L && length(sketch_ids) < n_cells) {
        sketch_mode <- TRUE
      } else if (length(sketch_ids) >= n_cells) {
        sketch_mode <- FALSE
        sketch_ids <- NULL
        pca.method <- "direct"
      } else if (n_cells > max.cells.direct) {
        set.seed(seed)
        nsk <- min(max.cells.direct, n_cells)
        sketch_ids <- cell_ids[sort(sample.int(n_cells, nsk))]
        sketch_mode <- TRUE
        .ic2_step(rep, paste0(
          "No sketch found  -  auto-sketching ", .ic2_fmt_num(nsk),
          " / ", .ic2_fmt_num(n_cells), " cells for PCA fit"
        ))
      } else {
        pca.method <- "direct"
      }
    } else if (n_cells > max.cells.direct) {
      warning(
        "n_cells (", n_cells, ") > max.cells.direct (", max.cells.direct,
        "). Consider pca.method = 'cov' or use.sketch = TRUE."
      )
      pca.method <- "direct"
    } else {
      pca.method <- "direct"
    }
  }

  fit_ids <- if (sketch_mode) sketch_ids else cell_ids
  engine_label <- switch(
    pca.method,
    cov = "cov stream (exact, all cells)",
    sketch = "sketch  ->  irlba  ->  project",
    direct = "direct irlba/prcomp",
    pca.method
  )

  rep <- .ic2_banner_start(rep, c(
    paste0("Assay: ", assay.name, "   |   transform: ", transform, "   |   features: ", method),
    paste0(
      "Features: ", .ic2_fmt_num(length(features)),
      "   |   Cells: ", .ic2_fmt_num(n_cells),
      "   |   npcs: ", npcs
    ),
    paste0("Engine: ", engine_label)
  ))

  rotation <- NULL
  sdev <- NULL
  center_vec <- NULL
  scale_vec <- NULL
  pca_fit <- NULL

  ## Stream PCA from assay handles (SF + transform on the fly).
  ## When the assay has far more rows than gene.model / selected features,
  ## materialize a skinny CSC cache once so cov/project do not re-read all genes.
  pca_handles <- handles
  pca_sf <- size_factors
  pca_transform <- transform

  max_assay_rows <- 0L
  for (hh in handles) {
    hh <- .ic2_enrich_handle(hh, load_dimnames = FALSE)
    max_assay_rows <- max(max_assay_rows, as.integer(hh$n_features %||% 0L))
  }
  if (max_assay_rows > max(2L * length(features), length(features) + 50L)) {
    cache_dir <- file.path(x@storage$root_dir, "cache")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    key <- paste(
      modality, assay.name, length(features), length(cell_ids),
      paste(utils::head(features, 3L), collapse = ","),
      paste(utils::tail(features, 3L), collapse = ","),
      paste(utils::head(cell_ids, 2L), collapse = ","),
      paste(utils::tail(cell_ids, 2L), collapse = ","),
      sep = "|"
    )
    digest <- sprintf("%08x", sum(utf8ToInt(key)) %% 4294967296)
    cache_path <- file.path(
      cache_dir,
      paste0("pca_features_", modality, "_", assay.name, "_", digest, ".h5")
    )
    meta_path <- paste0(cache_path, ".meta.rds")
    reuse <- FALSE
    if (file.exists(cache_path) && file.exists(meta_path)) {
      meta <- tryCatch(readRDS(meta_path), error = function(e) NULL)
      if (is.list(meta) &&
          identical(as.character(meta$features), as.character(features)) &&
          identical(as.character(meta$cell_ids), as.character(cell_ids)) &&
          identical(as.character(meta$assay), as.character(assay.name))) {
        reuse <- TRUE
      }
    }
    if (!isTRUE(reuse)) {
      .ic2_step(rep, paste0(
        "Building gene.model feature cache (",
        .ic2_fmt_num(length(features)), " features  x  ",
        .ic2_fmt_num(length(cell_ids)), " cells)  ->  ",
        basename(cache_path)
      ))
      ch <- .ic2_materialize_feature_subset_hdf5(
        handles = handles,
        features = features,
        cell_ids = cell_ids,
        dest_path = cache_path,
        chunk.size = chunk.size,
        size_factors = size_factors,
        transform = transform,
        overwrite = TRUE,
        progress = function(cur, tot) {
          rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  hvg cache")
        }
      )
      saveRDS(
        list(
          features = features,
          cell_ids = cell_ids,
          assay = assay.name,
          path = cache_path,
          n_features = length(features),
          n_cells = length(cell_ids)
        ),
        meta_path
      )
    } else {
      .ic2_step(rep, paste0("Reusing gene.model feature cache: ", basename(cache_path)))
      ch <- .ic2_assay_handle(
        path = cache_path,
        backend = "hdf5_sparse",
        group = "counts",
        n_features = length(features),
        n_cells = length(cell_ids),
        feature_ids = NULL,
        cell_ids = NULL
      )
    }
    pca_handles <- list(ch)
    pca_sf <- NULL
    pca_transform <- "none"
  }

  if (identical(pca.method, "cov")) {
    .ic2_step(rep, "[1/2] Out-of-core covariance PCA (stream all cells)")
    pca_fit <- .ic2_pca_cov_stream(
      handles = pca_handles,
      features = features,
      cell_ids = cell_ids,
      npcs = npcs,
      center = center,
      scale. = scale.,
      transform = pca_transform,
      chunk.size = chunk.size,
      size_factors = pca_sf,
      progress = function(cur, tot, phase = "") {
        rep <<- .ic2_progress_bar(
          rep, cur, tot,
          prefix = if (nzchar(phase)) paste0("  ", phase) else "  cov"
        )
      }
    )
    rotation <- pca_fit$rotation
    sdev <- pca_fit$sdev
    center_vec <- pca_fit$center
    scale_vec <- pca_fit$scale
    features <- rownames(rotation)
    npcs <- ncol(rotation)
    fit_ids <- cell_ids
    sketch_mode <- FALSE

    .ic2_step(rep, "[2/2] Projecting all cells onto loadings")
    emb <- .ic2_pca_project(
      handles = pca_handles,
      features = features,
      cell_ids = cell_ids,
      rotation = rotation,
      center = center_vec,
      scale = scale_vec,
      transform = pca_transform,
      chunk.size = chunk.size,
      size_factors = pca_sf,
      progress = function(cur, tot) {
        rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  project")
      }
    )
  } else {
    ## sketch or direct: load cells  x  features, irlba/prcomp, project rest
    .ic2_step(rep, "[1/3] Loading features for PCA fit (sparse preferred)")
    X_fit <- .ic2_load_cells_by_features_sparse(
      handles = pca_handles,
      features = features,
      cell_ids = fit_ids,
      chunk.size = chunk.size,
      transform = pca_transform,
      size_factors = pca_sf,
      progress = function(cur, tot) {
        rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  load")
      }
    )

    .ic2_step(rep, "[2/3] Computing PCA (irlba / truncated SVD)")
    pca_fit <- .ic2_pca_fit(
      X_fit,
      npcs = npcs,
      center = center,
      scale. = scale.
    )
    rep <- .ic2_progress_bar(rep, 1L, 1L, prefix = "  svd")

    rotation <- pca_fit$rotation
    sdev <- pca_fit$sdev
    center_vec <- pca_fit$center
    scale_vec <- pca_fit$scale
    features <- rownames(rotation)
    npcs <- ncol(rotation)

    .ic2_step(rep, "[3/3] Writing embeddings (project all cells)")
    emb <- matrix(NA_real_, nrow = n_cells, ncol = npcs,
                  dimnames = list(cell_ids, paste0("PC", seq_len(npcs))))

    if (!is.null(pca_fit$x) && sketch_mode) {
      common <- intersect(rownames(pca_fit$x), cell_ids)
      emb[common, ] <- pca_fit$x[common, , drop = FALSE]
      project_ids <- setdiff(cell_ids, common)
    } else if (!is.null(pca_fit$x) && !sketch_mode) {
      emb[rownames(pca_fit$x), ] <- pca_fit$x
      project_ids <- character(0)
    } else {
      project_ids <- cell_ids
    }

    if (length(project_ids) > 0L) {
      emb_proj <- .ic2_pca_project(
        handles = pca_handles,
        features = features,
        cell_ids = project_ids,
        rotation = rotation,
        center = center_vec,
        scale = scale_vec,
        transform = pca_transform,
        chunk.size = chunk.size,
        size_factors = pca_sf,
        progress = function(cur, tot) {
          rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  project")
        }
      )
      emb[rownames(emb_proj), ] <- emb_proj
    }
  }

  ## Suggest opt.pcs (knee on variance; typically ~10 for scRNA elbows)
  OPTpcs <- .ic2_suggest_opt_pcs(sdev, method = "knee")

  ## Guard: never leave NA/NaN/Inf in PCA embeddings (breaks UMAP/neighbors)
  bad <- which(rowSums(!is.finite(emb)) > 0L)
  if (length(bad) > 0L) {
    stop(
      length(bad), " cells have non-finite PCA embeddings after fit/project. ",
      "First ids: ", paste(utils::head(rownames(emb)[bad], 5L), collapse = ", "),
      ". Check assay handles / feature overlap, then re-run run.pca()."
    )
  }

  ## ---- Persist (HDF5 dense embeddings + loadings; JSON meta) ----
  ## RNA stays at reductions/pca. Other modalities must not overwrite it.
  red_dir <- if (identical(modality, "rna")) {
    file.path(x@storage$root_dir, "reductions", "pca")
  } else {
    file.path(x@storage$root_dir, "modalities", modality, "reductions", "pca")
  }
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  load_path <- file.path(red_dir, "feature_loadings.h5")
  stdev_path <- file.path(red_dir, "stdev.tsv")
  meta_path <- file.path(red_dir, "reduction.json")

  .ic2_h5_write_dense(emb, emb_path)
  .ic2_h5_write_dense(rotation, load_path)
  utils::write.table(
    data.frame(pc = seq_along(sdev), stdev = as.numeric(sdev), stringsAsFactors = FALSE),
    file = stdev_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  ## Remove legacy .rds siblings if present
  unlink(c(
    file.path(red_dir, "cell_embeddings.rds"),
    file.path(red_dir, "feature_loadings.rds"),
    file.path(red_dir, "stdev.rds"),
    file.path(red_dir, "pca_meta.rds")
  ), force = TRUE)

  pca_meta <- list(
    name = "pca",
    assay.name = assay.name,
    assay = assay.name,
    method = method,
    pca.method = pca.method,
    transform = transform,
    features = features,
    npcs = npcs,
    opt.pcs = OPTpcs,
    sketch_mode = sketch_mode,
    n_fit_cells = length(fit_ids),
    n_cells = n_cells,
    center = center,
    scale. = scale.,
    created = as.character(Sys.time()),
    paths = list(
      embeddings = emb_path,
      loadings = load_path,
      stdev = stdev_path
    )
  )
  .ic2_soda_write_json(pca_meta, meta_path)

  pca_ref <- list(
    cell_embeddings = emb_path,
    feature_loadings = load_path,
    stdev = sdev,
    opt.pcs = OPTpcs,
    npcs = npcs,
    features = features,
    assay = assay.name,
    transform = transform,
    meta_path = meta_path
  )

  x@modalities[[modality]]$reductions$pca <- pca_ref
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x@modalities[[modality]]$active_reduction <- "pca"
  if (identical(modality, "rna")) {
    x@reductions$integrated_pca <- NULL
    x@misc$opt.pcs <- OPTpcs
    x@storage$reduction_paths$pca <- red_dir
  }

  x@params$dimensionality_reduction[[modality]] <- list(
    pca = list(
      assay.name = assay.name,
      method = method,
      pca.method = pca.method,
      transform = transform,
      npcs = npcs,
      opt.pcs = OPTpcs,
      sketch_mode = sketch_mode,
      n_features = length(features)
    )
  )

  x <- .ic2_log_command(
    x, "run.pca",
    arguments = list(
      modality = modality, assay.name = assay.name, method = method,
      pca.method = pca.method, transform = transform, npcs = npcs,
      use.sketch = use.sketch, sketch_mode = sketch_mode
    ),
    outputs = list(
      embeddings = emb_path, loadings = load_path, opt.pcs = OPTpcs
    )
  )

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  rep <- .ic2_add_output(rep, "Embeddings", emb_path)
  rep <- .ic2_add_output(rep, "Loadings", load_path)
  rep <- .ic2_add_output(rep, "Stdev", stdev_path)
  wall <- as.numeric(difftime(Sys.time(), rep$wall0, units = "secs"))
  if (is.finite(wall) && wall > 0) {
    rep$stats$cells_per_sec <- n_cells / wall
  }
  .ic2_finish(rep, extra_lines = paste0(
    "PCs: ", npcs, "   |   suggested opt.pcs: ", OPTpcs,
    "   |   ", pca.method
  ))
  x
}

#' Iterate assay columns as feature-subset batches (features  x  cells)
#' @keywords internal
.ic2_for_feature_batches <- function(
    handles,
    features,
    cell_ids,
    chunk.size = 200000L,
    FUN,
    progress = NULL,
    progress.phase = "",
    transform = c("none", "log2p1", "log1p"),
    size_factors = NULL
) {
  transform <- match.arg(transform)
  features <- as.character(features)
  cell_ids <- as.character(cell_ids)
  n_cells <- length(cell_ids)
  remaining <- cell_ids
  done <- 0L

  for (h in handles) {
    if (length(remaining) < 1L) break
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    local <- intersect(remaining, h$cell_ids)
    if (length(local) < 1L) next
    row_idx <- match(features, h$feature_ids)
    col_idx <- match(local, h$cell_ids)
    ord <- order(col_idx)
    local <- local[ord]
    col_idx <- col_idx[ord]
    ## Chunk by chunk.size on ordered columns (gaps OK). H5 reader keeps one
    ## open and splits contig runs internally  -  not one open per tiny run.
    starts <- seq.int(1L, length(col_idx), by = as.integer(chunk.size))
    for (s in starts) {
      e <- min(s + as.integer(chunk.size) - 1L, length(col_idx))
      cols <- col_idx[s:e]
      raw <- .ic2_read_cols(h, cols)
      cn <- colnames(raw)
      if (is.null(cn)) cn <- h$cell_ids[cols]
      ## Missing genes in this sample  ->  zeros (virtual multi-sample)
      block <- matrix(0, nrow = length(features), ncol = length(cols),
                      dimnames = list(features, cn))
      present <- which(!is.na(row_idx))
      if (length(present) > 0L) {
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
      block <- .ic2_scale_cols_by_sf(block, size_factors)
      block <- .ic2_assay_transform_block(block, transform = transform)
      FUN(block)
      done <- done + length(cn)
      if (is.function(progress)) {
        if (length(formals(progress)) >= 3L) {
          progress(done, n_cells, progress.phase)
        } else {
          progress(done, n_cells)
        }
      }
    }
    remaining <- setdiff(remaining, local)
  }
  if (length(remaining) > 0L) {
    stop(
      "Could not stream assay data for ", length(remaining),
      " cells (first: ", paste(utils::head(remaining, 5L), collapse = ", "), ")."
    )
  }
  if (is.function(progress)) {
    if (length(formals(progress)) >= 3L) {
      progress(n_cells, n_cells, progress.phase)
    } else {
      progress(n_cells, n_cells)
    }
  }
  invisible(done)
}

#' Exact out-of-core PCA via streamed p x p covariance (p << n)
#'
#' Two HDF5 passes: (1) feature moments, (2) Gram accumulation, then
#' symmetric eigen on the \eqn{p \times p} covariance. Returns loadings
#' compatible with \code{.ic2_pca_project} (no cell embeddings).
#'
#' @keywords internal
.ic2_pca_cov_stream <- function(
    handles,
    features,
    cell_ids,
    npcs = 50L,
    center = TRUE,
    scale. = FALSE,
    transform = c("none", "log2p1", "log1p"),
    chunk.size = 200000L,
    progress = NULL,
    size_factors = NULL
) {
  transform <- match.arg(transform)
  features <- as.character(features)
  cell_ids <- as.character(cell_ids)
  p <- length(features)
  n_cells <- length(cell_ids)
  if (p < 2L || n_cells < 2L) {
    stop("Covariance PCA needs >= 2 features and >= 2 cells.")
  }
  npcs <- as.integer(npcs)
  npcs <- min(npcs, p - 1L, n_cells - 1L)
  if (npcs < 1L) stop("'npcs' must be >= 1 after rank checks.")

  n <- as.numeric(n_cells)
  cm <- rep(0, p)
  sds <- rep(1, p)

  ## ---- Pass 1: moments (skip when already centered/scaled) ----
  if (isTRUE(center) || isTRUE(scale.)) {
    sum_g <- numeric(p)
    sumsq_g <- numeric(p)
    n_seen <- 0L
    .ic2_for_feature_batches(
      handles, features, cell_ids,
      chunk.size = chunk.size,
      progress = progress,
      progress.phase = "moments",
      transform = transform,
      size_factors = size_factors,
      FUN = function(block) {
        ## block: features  x  cells
        sum_g <<- sum_g + rowSums(block)
        sumsq_g <<- sumsq_g + rowSums(block * block)
        n_seen <<- n_seen + ncol(block)
      }
    )
    if (n_seen != n_cells) {
      stop("Moment pass saw ", n_seen, " cells but expected ", n_cells, ".")
    }
    n <- as.numeric(n_seen)
    cm <- if (isTRUE(center)) sum_g / n else rep(0, p)
    if (isTRUE(scale.)) {
      ## sample sd with explicit zeros (same convention as dense prcomp)
      sds <- sqrt(pmax(0, (sumsq_g - n * cm * cm) / (n - 1)))
    } else {
      sds <- rep(1, p)
    }
    sds[!is.finite(sds) | sds < 1e-8] <- 1

    ## Drop near-zero-variance features when scaling (or always if var ~ 0 and center)
    var0 <- (sumsq_g - n * (sum_g / n) * (sum_g / n)) / (n - 1)
    keep <- which(is.finite(var0) & var0 > 1e-8)
    if (length(keep) < 2L) {
      stop(
        "Fewer than 2 variable features for covariance PCA ",
        "(n_features=", p, "). Check assay / gene.model."
      )
    }
    if (length(keep) < p) {
      features <- features[keep]
      cm <- cm[keep]
      sds <- sds[keep]
      p <- length(features)
      npcs <- min(npcs, p - 1L)
    }
  }

  ## ---- Pass 2: Gram of centered/scaled features (p  x  p) ----
  S <- matrix(0, nrow = p, ncol = p)
  .ic2_for_feature_batches(
    handles, features, cell_ids,
    chunk.size = chunk.size,
    progress = progress,
    progress.phase = "gram",
    transform = transform,
    size_factors = size_factors,
    FUN = function(block) {
      if (isTRUE(center)) {
        block <- block - cm
      }
      if (isTRUE(scale.)) {
        block <- block / sds
      }
      S <<- S + tcrossprod(block)
    }
  )
  C <- S / (n - 1)

  ## ---- Eigen ----
  eg <- eigen(C, symmetric = TRUE)
  k <- min(npcs, sum(eg$values > 1e-12), p)
  if (k < 1L) stop("Covariance PCA produced no positive eigenvalues.")
  rotation <- eg$vectors[, seq_len(k), drop = FALSE]
  sdev <- sqrt(pmax(eg$values[seq_len(k)], 0))
  colnames(rotation) <- paste0("PC", seq_len(k))
  rownames(rotation) <- features
  names(sdev) <- colnames(rotation)
  names(cm) <- features
  names(sds) <- features

  if (any(!is.finite(rotation)) || any(!is.finite(sdev))) {
    stop("Covariance PCA produced non-finite loadings/sdev.")
  }

  list(
    x = NULL,
    rotation = rotation,
    sdev = as.numeric(sdev),
    center = as.numeric(cm),
    scale = as.numeric(sds)
  )
}

#' Load cells  x  features as a sparse dgCMatrix from assay handles
#' @keywords internal
.ic2_load_cells_by_features_sparse <- function(
    handles,
    features,
    cell_ids,
    chunk.size = 200000L,
    transform = c("none", "log2p1", "log1p"),
    progress = NULL,
    size_factors = NULL
) {
  transform <- match.arg(transform)
  features <- as.character(features)
  cell_ids <- as.character(cell_ids)
  n_cells <- length(cell_ids)
  n_feat <- length(features)
  if (n_cells < 1L || n_feat < 1L) {
    stop("Need at least one cell and one feature for PCA load.")
  }

  ## Collect feature x cell blocks, then transpose to cells x features
  blocks <- list()
  remaining <- cell_ids
  done <- 0L

  for (h in handles) {
    if (length(remaining) < 1L) break
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    local <- intersect(remaining, h$cell_ids)
    if (length(local) < 1L) next
    row_idx <- match(features, h$feature_ids)
    col_idx <- match(local, h$cell_ids)
    ord <- order(col_idx)
    local <- local[ord]
    col_idx <- col_idx[ord]
    starts <- seq.int(1L, length(col_idx), by = as.integer(chunk.size))
    for (s in starts) {
      e <- min(s + as.integer(chunk.size) - 1L, length(col_idx))
      cols <- col_idx[s:e]
      raw <- .ic2_read_cols(h, cols)
      cn <- colnames(raw)
      if (is.null(cn)) cn <- h$cell_ids[cols]
      present <- which(!is.na(row_idx))
      if (length(present) == length(features) && !anyNA(row_idx)) {
        block <- raw[row_idx, , drop = FALSE]
        block <- .ic2_to_dgCMatrix(block)
        rownames(block) <- features
        colnames(block) <- cn
      } else {
        ## Missing genes in this sample  ->  explicit zero rows
        block <- Matrix::Matrix(
          0, nrow = length(features), ncol = length(cols),
          sparse = TRUE, dimnames = list(features, cn)
        )
        if (length(present) > 0L) {
          sub <- .ic2_to_dgCMatrix(raw[row_idx[present], , drop = FALSE])
          block[present, ] <- sub
        }
      }
      block <- .ic2_scale_cols_by_sf(block, size_factors)
      block <- .ic2_assay_transform_block(block, transform = transform)
      block <- .ic2_to_dgCMatrix(block)
      rownames(block) <- features
      colnames(block) <- cn
      blocks[[length(blocks) + 1L]] <- block
      done <- done + length(cn)
      if (is.function(progress)) progress(done, n_cells)
    }
    remaining <- setdiff(remaining, local)
  }
  if (length(remaining) > 0L) {
    stop("Could not load assay data for ", length(remaining), " cells.")
  }
  mat <- if (length(blocks) == 1L) {
    blocks[[1L]]
  } else {
    do.call(cbind, blocks)
  }
  ## cells  x  features, ordered as cell_ids
  mat <- Matrix::t(mat)
  if (!identical(rownames(mat), cell_ids)) {
    mat <- mat[cell_ids, , drop = FALSE]
  }
  if (is.function(progress)) progress(n_cells, n_cells)
  mat
}

#' Dense loader (compat); prefer \code{.ic2_load_cells_by_features_sparse}
#' @keywords internal
.ic2_load_cells_by_features <- function(
    handles,
    features,
    cell_ids,
    chunk.size = 200000L,
    transform = c("none", "log2p1", "log1p"),
    progress = NULL,
    size_factors = NULL
) {
  sp <- .ic2_load_cells_by_features_sparse(
    handles, features, cell_ids,
    chunk.size = chunk.size,
    transform = transform,
    progress = progress,
    size_factors = size_factors
  )
  as.matrix(sp)
}

#' Column SDs for a cells x features matrix (sparse-friendly)
#' @keywords internal
.ic2_col_sds <- function(X) {
  n <- nrow(X)
  if (n < 2L) return(rep(1, ncol(X)))
  if (inherits(X, "dgCMatrix")) {
    cm <- as.numeric(Matrix::colMeans(X))
    ## column sum of squares via CSC
    ss <- numeric(ncol(X))
    p <- X@p
    xx <- X@x
    for (j in seq_len(ncol(X))) {
      a <- p[j] + 1L
      b <- p[j + 1L]
      if (b >= a) ss[j] <- sum(xx[a:b]^2)
    }
    ## sample variance treating missing entries as explicit zeros
    sds <- sqrt(pmax(0, (ss - n * cm * cm) / (n - 1)))
    sds[!is.finite(sds) | sds < 1e-8] <- 1
    return(sds)
  }
  X <- as.matrix(X)
  sds <- apply(X, 2L, stats::sd)
  sds[!is.finite(sds) | sds < 1e-8] <- 1
  sds
}

#' Fit PCA on cells  x  features (sparse + irlba preferred)
#' @keywords internal
.ic2_pca_fit <- function(X, npcs, center = TRUE, scale. = FALSE) {
  if (inherits(X, "Matrix")) {
    X <- methods::as(X, "dgCMatrix")
    X@x[!is.finite(X@x)] <- 0
  } else {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    X[!is.finite(X)] <- 0
  }

  ## Drop near-zero-variance features
  sds0 <- .ic2_col_sds(X)
  keep <- which(is.finite(sds0) & sds0 > 1e-8)
  if (length(keep) < 2L) {
    mx <- if (inherits(X, "Matrix")) {
      max(abs(X@x), na.rm = TRUE)
    } else {
      max(abs(X), na.rm = TRUE)
    }
    stop(
      "Fewer than 2 variable features available for PCA ",
      "(max|value|=", signif(mx, 3), ", n_features=", ncol(X), "). ",
      "If you are using scale.data, re-run data.scale(..., overwrite = TRUE)."
    )
  }
  if (length(keep) < ncol(X)) {
    X <- X[, keep, drop = FALSE]
  }

  npcs <- as.integer(npcs)
  npcs <- min(npcs, ncol(X), max(1L, nrow(X) - 1L))
  if (npcs < 1L) stop("'npcs' must be >= 1 after rank checks.")

  cell_names <- rownames(X)
  feat_names <- colnames(X)

  pack <- function(emb, rotation, sdev, center_vec, scale_vec) {
    colnames(emb) <- paste0("PC", seq_len(ncol(emb)))
    colnames(rotation) <- colnames(emb)
    rownames(emb) <- cell_names
    rownames(rotation) <- feat_names
    if (any(!is.finite(emb)) || any(!is.finite(rotation)) || any(!is.finite(sdev))) {
      return(NULL)
    }
    list(
      x = emb,
      rotation = rotation,
      sdev = as.numeric(sdev),
      center = as.numeric(center_vec),
      scale = as.numeric(scale_vec)
    )
  }

  out <- NULL
  use_irlba <- requireNamespace("irlba", quietly = TRUE) &&
    nrow(X) > 100L &&
    npcs < min(nrow(X), ncol(X))

  if (use_irlba) {
    cm <- if (isTRUE(center)) {
      if (inherits(X, "Matrix")) {
        as.numeric(Matrix::colMeans(X))
      } else {
        as.numeric(colMeans(X))
      }
    } else {
      rep(0, ncol(X))
    }
    sds <- if (isTRUE(scale.)) {
      .ic2_col_sds(X)
    } else {
      rep(1, ncol(X))
    }
    ## irlba can center/scale without forming a dense copy
    fit <- tryCatch(
      suppressWarnings(
        irlba::irlba(
          X,
          nv = npcs,
          center = if (isTRUE(center)) cm else FALSE,
          scale = if (isTRUE(scale.)) sds else FALSE,
          maxit = 3000L,
          work = max(npcs + 20L, as.integer(ceiling(npcs * 2)))
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(fit) && length(fit$d) >= 1L) {
      k <- min(npcs, length(fit$d), ncol(fit$u), ncol(fit$v))
      emb <- fit$u[, seq_len(k), drop = FALSE] %*%
        diag(fit$d[seq_len(k)], nrow = k)
      rotation <- fit$v[, seq_len(k), drop = FALSE]
      sdev <- fit$d[seq_len(k)] / sqrt(max(1, nrow(X) - 1))
      out <- pack(emb, rotation, sdev, cm, sds)
    }
  }

  if (is.null(out)) {
    if (use_irlba) {
      message("  irlba PCA unstable  -  falling back to prcomp().")
    }
    Xd <- as.matrix(X)
    storage.mode(Xd) <- "double"
    fit <- stats::prcomp(Xd, center = center, scale. = scale., rank. = npcs)
    k <- min(npcs, ncol(fit$x), ncol(fit$rotation), length(fit$sdev))
    emb <- fit$x[, seq_len(k), drop = FALSE]
    rotation <- fit$rotation[, seq_len(k), drop = FALSE]
    sdev <- fit$sdev[seq_len(k)]
    cm <- if (is.null(fit$center)) rep(0, ncol(Xd)) else fit$center
    sds <- if (is.null(fit$scale) || isTRUE(all(fit$scale == FALSE))) {
      rep(1, ncol(Xd))
    } else {
      as.numeric(fit$scale)
    }
    sds[!is.finite(sds) | sds < 1e-8] <- 1
    out <- pack(emb, rotation, sdev, cm, sds)
  }

  if (is.null(out)) {
    stop("PCA produced non-finite embeddings (irlba and prcomp). Check scale.data / gene.model.")
  }
  out
}

#' Project cells onto PCA rotation (chunked; sparse multiply when possible)
#' @keywords internal
.ic2_pca_project <- function(
    handles,
    features,
    cell_ids,
    rotation,
    center,
    scale,
    transform = c("none", "log2p1", "log1p"),
    chunk.size = 200000L,
    progress = NULL,
    size_factors = NULL
) {
  transform <- match.arg(transform)
  features <- as.character(features)
  n_cells <- length(cell_ids)
  npcs <- ncol(rotation)
  out <- matrix(0, nrow = n_cells, ncol = npcs,
                dimnames = list(cell_ids, colnames(rotation)))
  remaining <- cell_ids
  done <- 0L

  center <- as.numeric(center)
  scale <- as.numeric(scale)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  ## emb = t(block) %*% (rotation/scale) - (center/scale)%*%rotation
  R2 <- rotation / scale
  offset <- as.numeric(crossprod(center / scale, rotation))

  for (h in handles) {
    if (length(remaining) < 1L) break
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    local <- intersect(remaining, h$cell_ids)
    if (length(local) < 1L) next
    row_idx <- match(features, h$feature_ids)
    col_idx <- match(local, h$cell_ids)
    ord <- order(col_idx)
    local <- local[ord]
    col_idx <- col_idx[ord]
    starts <- seq.int(1L, length(col_idx), by = as.integer(chunk.size))
    for (s in starts) {
      e <- min(s + as.integer(chunk.size) - 1L, length(col_idx))
      cols <- col_idx[s:e]
      raw <- .ic2_read_cols(h, cols)
      cn <- colnames(raw)
      if (is.null(cn)) cn <- h$cell_ids[cols]
      present <- which(!is.na(row_idx))
      block <- matrix(0, nrow = length(features), ncol = length(cols),
                      dimnames = list(features, cn))
      if (length(present) > 0L) {
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
      block <- .ic2_scale_cols_by_sf(block, size_factors)
      block <- .ic2_assay_transform_block(block, transform = transform)
      emb_block <- crossprod(block, R2)
      emb_block <- sweep(emb_block, 2L, offset, "-")
      rownames(emb_block) <- cn
      out[cn, ] <- emb_block
      done <- done + length(cn)
      if (is.function(progress)) progress(done, n_cells)
    }
    remaining <- setdiff(remaining, local)
  }
  if (length(remaining) > 0L) {
    stop(
      "PCA project could not find assay data for ", length(remaining),
      " cells (first: ", paste(utils::head(remaining, 5L), collapse = ", "), ")."
    )
  }
  if (any(!is.finite(out))) {
    out[!is.finite(out)] <- 0
  }
  if (is.function(progress)) progress(n_cells, n_cells)
  out
}

#' Suggest an optimal number of PCs from standard deviations
#'
#' @param sdev Numeric PC standard deviations (length = npcs).
#' @param method \code{"knee"} (default): maximum distance from the
#'   chord joining the first and last variance points  -  usually lands
#'   near the visual elbow (~10 for many scRNA runs).
#'   \code{"legacy"}: old iCellR rule \code{sum(sdev > 2*mean(sdev)) + 1}.
#'   \code{"variance"}: fewest PCs explaining at least \code{var.explained}
#'   of the computed PC variance.
#' @param var.explained Cumulative variance target for \code{"variance"}.
#' @param min.pcs,max.pcs Bounds on the suggestion.
#' @return Integer scalar.
#' @keywords internal
.ic2_suggest_opt_pcs <- function(
    sdev,
    method = c("knee", "legacy", "variance"),
    var.explained = 0.8,
    min.pcs = 2L,
    max.pcs = NULL
) {
  method <- match.arg(method)
  sdev <- as.numeric(sdev)
  sdev <- sdev[is.finite(sdev) & sdev > 0]
  n <- length(sdev)
  if (n < 2L) return(max(1L, n))
  if (is.null(max.pcs)) max.pcs <- n
  max.pcs <- as.integer(min(as.integer(max.pcs), n))
  min.pcs <- as.integer(max(1L, min(as.integer(min.pcs), max.pcs)))

  opt <- switch(
    method,
    legacy = {
      as.integer(sum(sdev > mean(sdev) * 2) + 1L)
    },
    variance = {
      ve <- sdev^2
      ve <- ve / sum(ve)
      as.integer(which(cumsum(ve) >= var.explained)[1L])
    },
    knee = {
      ## Distance from each point to the chord (1, y1)-(n, yn) on variance
      y <- sdev^2
      y <- y / max(y[1L], .Machine$double.eps)
      x <- seq_len(n)
      x1 <- 1
      y1 <- y[1L]
      x2 <- n
      y2 <- y[n]
      ## perpendicular distance
      num <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
      den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
      if (!is.finite(den) || den < .Machine$double.eps) {
        as.integer(max(min.pcs, round(n / 5)))
      } else {
        as.integer(which.max(num / den))
      }
    }
  )

  if (!is.finite(opt) || is.na(opt)) opt <- min.pcs
  opt <- max(min.pcs, min(as.integer(opt), max.pcs))
  opt
}
