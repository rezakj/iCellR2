#' Run UMAP on PCA embeddings (scalable)
#'
#' Fits UMAP with \pkg{uwot} on PCA cell embeddings. By default runs
#' \strong{full} UMAP when \code{n_cells <= max.cells.direct} (5M). Above
#' that (or if \code{use.sketch = TRUE}), fits on a cell sketch
#' (\code{ret_model = TRUE}) and projects the rest with
#' \code{umap_transform()} in chunks.
#'
#' @param x An \code{iCellR2} object after \code{\link{run.pca}}.
#' @param modality Modality. Default \code{"rna"}.
#' @param use Which PCA to embed. \code{"auto"} (default) uses
#'   \code{reductions/aligned_pca} when that file exists, otherwise
#'   \code{reductions/pca}. To force one of them, pass
#'   \code{"aligned_pca"} or \code{"pca"}
#'   (\code{"reductions/aligned_pca"} and \code{"reductions/pca"}
#'   are accepted too).
#' @param dims PC indices to use. Default \code{1:10}. Pass \code{NULL}
#'   to use \code{1:opt.pcs} when available.
#' @param my.seed RNG seed. Default \code{1}.
#' @param use.sketch If \code{TRUE}, prefer sketch -> project even when
#'   \code{n_cells <= max.cells.direct}. If \code{FALSE} (default), run
#'   full UMAP up to \code{max.cells.direct}; above that, sketch is
#'   forced automatically.
#' @param sketch.name Sketch name from \code{\link{run.sketch}}.
#' @param max.cells.direct Max cells for direct (full) UMAP.
#'   Default \code{5000000}. Larger objects auto-switch to sketch -> project.
#' @param chunk.size Cells per \code{umap_transform} chunk.
#'   Default \code{200000}.
#' @param n_neighbors,n_components,metric,n_epochs,learning_rate,scale,init,init_sdev,spread,min_dist,set_op_mix_ratio,local_connectivity,bandwidth,repulsion_strength,negative_sample_rate,a,b,nn_method,n_trees,search_k,approx_pow,y,target_n_neighbors,target_metric,target_weight,pca,pca_center,pcg_rand,fast_sgd,n_threads,n_sgd_threads,grain_size,tmpdir
#'   Passed to \code{uwot::umap} / \code{uwot::umap_transform}.
#' @param verbose Progress style for the **iCellR2** reporter:
#'   \code{TRUE}/\code{"full"}, \code{"compact"}, or \code{FALSE}.
#'   Underlying \pkg{uwot} console logging is off by default (set
#'   \code{options(iCellR2.uwot.verbose = TRUE)} only for package debugging).
#' @param save.model If \code{TRUE}, also write an R-only \code{.uwot} model
#'   for later \code{umap_transform}. Default \code{FALSE}: portable archive
#'   keeps embeddings (\code{cell_embeddings.h5}) + params in
#'   \code{reduction.json} only.
#'
#' @return Updated \code{iCellR2} with UMAP embeddings on disk
#'   (columns \code{V1}, \code{V2}, ...  -  classic iCellR naming).
#'
#' @examples
#' \dontrun{
#' obj <- run.pca(obj, npcs = 50)
#' obj <- run.umap(obj)                    # auto: aligned_pca if present
#' obj <- run.umap(obj, use = "aligned_pca")  # reductions/aligned_pca
#' obj <- run.umap(obj, use = "pca")          # reductions/pca
#' }
#'
#' @importFrom uwot umap umap_transform
#' @export
run.umap <- function(
    x,
    modality = "rna",
    use = "auto",
    dims = 1:10,
    my.seed = 1,
    use.sketch = FALSE,
    sketch.name = "sketch1",
    max.cells.direct = 5000000L,
    chunk.size = 200000L,
    n_neighbors = 15,
    n_components = 2,
    metric = "euclidean",
    n_epochs = NULL,
    learning_rate = 1,
    scale = FALSE,
    init = "spectral",
    init_sdev = NULL,
    spread = 1,
    min_dist = 0.01,
    set_op_mix_ratio = 1,
    local_connectivity = 1,
    bandwidth = 1,
    repulsion_strength = 1,
    negative_sample_rate = 5,
    a = NULL,
    b = NULL,
    nn_method = NULL,
    n_trees = 50,
    search_k = NULL,
    approx_pow = FALSE,
    y = NULL,
    target_n_neighbors = n_neighbors,
    target_metric = "euclidean",
    target_weight = 0.5,
    pca = NULL,
    pca_center = TRUE,
    pcg_rand = TRUE,
    fast_sgd = FALSE,
    n_threads = 4,
    n_sgd_threads = 0,
    grain_size = 1,
    tmpdir = tempdir(),
    verbose = TRUE,
    save.model = FALSE
) {
  .ic2_require_object(x)
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for run.umap(). Install it from CRAN.")
  }
  modality <- tolower(modality)
  max.cells.direct <- as.integer(max.cells.direct)
  chunk.size <- as.integer(chunk.size)
  n_components <- as.integer(n_components)
  n_neighbors <- as.integer(n_neighbors)
  if (is.null(search_k)) search_k <- 2L * n_neighbors * as.integer(n_trees)

  rep <- .ic2_new_reporter("run.umap", verbose, soda_root = x@storage$root_dir)
  ## Keep uwot quiet  -  Seurat surfaces uwot's Annoy/epoch bars; we use iCellR2 progress only
  uwot_verbose <- isTRUE(getOption("iCellR2.uwot.verbose", FALSE))

  src <- .ic2_resolve_use(x, modality, "pca", use)
  out_name <- .ic2_aligned_out_name("umap", src$aligned)
  pca_ref <- x@modalities[[modality]]$reductions[[src$name]]
  if (!is.null(pca_ref$cell_embeddings) && nzchar(as.character(pca_ref$cell_embeddings)[1L])) {
    rep <- .ic2_add_input(rep, src$name, pca_ref$cell_embeddings)
  }
  emb_pca <- .ic2_get_plot_embedding(x, modality, src$name)
  cell_ids <- rownames(x@cell_data)
  common <- intersect(cell_ids, rownames(emb_pca))
  if (length(common) < 2L) {
    stop("PCA embeddings do not match object cell ids. Re-run run.pca().")
  }
  if (length(common) < length(cell_ids)) {
    warning(
      length(cell_ids) - length(common),
      " cells lack PCA embeddings and will be skipped."
    )
  }
  cell_ids <- common
  emb_pca <- emb_pca[cell_ids, , drop = FALSE]
  n_cells <- length(cell_ids)
  npcs_avail <- ncol(emb_pca)

  ## Resolve dims
  if (is.null(dims)) {
    opt <- x@misc$opt.pcs %||% pca_ref$opt.pcs %||% min(10L, npcs_avail)
    dims <- seq_len(min(as.integer(opt), npcs_avail))
  } else {
    dims <- as.integer(dims)
    dims <- dims[dims >= 1L & dims <= npcs_avail]
  }
  if (length(dims) < 2L) stop("'dims' must select at least 2 PCs.")
  X_all <- emb_pca[, dims, drop = FALSE]

  ## uwot rejects NA/NaN; drop incomplete PCA rows (e.g. failed project)
  finite_row <- rowSums(!is.finite(X_all)) == 0L
  n_bad <- sum(!finite_row)
  if (n_bad > 0L) {
    warning(
      n_bad, " / ", nrow(X_all),
      " cells have non-finite PCA values and will be skipped. ",
      "Re-run run.pca() if this is unexpected."
    )
    X_all <- X_all[finite_row, , drop = FALSE]
    cell_ids <- rownames(X_all)
  }
  if (nrow(X_all) < 2L) {
    stop(
      "Fewer than 2 cells with finite PCA embeddings for dims ",
      paste(range(dims), collapse = "-"),
      ". Re-run run.pca() then run.umap()."
    )
  }
  if (any(!is.finite(X_all))) {
    stop("PCA embeddings still contain non-finite values after filtering.")
  }
  n_cells <- length(cell_ids)

  ## Sketch decision:
  ## - n > max.cells.direct  ->  always sketch -> project (even if use.sketch = FALSE)
  ## - else if use.sketch = TRUE  ->  sketch when a proper subset exists / auto-sample
  ## - else  ->  full UMAP on all cells
  sketch_ids <- NULL
  sketch_mode <- FALSE
  force_sketch <- n_cells > max.cells.direct
  want_sketch <- isTRUE(use.sketch) || force_sketch
  if (force_sketch && !isTRUE(use.sketch)) {
    .ic2_step(rep, paste0(
      "n_cells (", .ic2_fmt_num(n_cells), ") > max.cells.direct (",
      .ic2_fmt_num(max.cells.direct), ")  -  forcing sketch -> project"
    ))
  }
  if (want_sketch) {
    sk <- x@modalities[[modality]]$sketches$sketch_models[[sketch.name]]
    if (!is.null(sk$path) && file.exists(sk$path)) {
      meta <- readRDS(sk$path)
      sketch_ids <- intersect(meta$cell_ids, cell_ids)
    } else if (!is.null(x@modalities[[modality]]$sketches$sketch_ids)) {
      sketch_ids <- intersect(x@modalities[[modality]]$sketches$sketch_ids, cell_ids)
    }
    ## Keep only sketch cells that have finite PCA rows
    sketch_ids <- intersect(sketch_ids, cell_ids)
    if (length(sketch_ids) >= max(n_neighbors + 1L, 50L) &&
        length(sketch_ids) < n_cells) {
      sketch_mode <- TRUE
    } else if (length(sketch_ids) >= n_cells) {
      sketch_mode <- FALSE
      sketch_ids <- NULL
    } else if (force_sketch || n_cells > max.cells.direct) {
      set.seed(my.seed)
      ## Fit sketch stays modest (full UMAP on 5M+ is the problem we avoid)
      nsk <- min(100000L, n_cells - 1L)
      nsk <- max(nsk, max(n_neighbors + 2L, 50L))
      nsk <- min(nsk, n_cells - 1L)
      sketch_ids <- cell_ids[sort(sample.int(n_cells, nsk))]
      sketch_mode <- TRUE
      .ic2_step(rep, paste0(
        "No sketch found  -  auto-sketching ", .ic2_fmt_num(nsk),
        " / ", .ic2_fmt_num(n_cells), " cells for UMAP fit"
      ))
    } else if (isTRUE(use.sketch) && n_cells <= max.cells.direct) {
      ## Optional sketch on smaller data only if a stored sketch exists
      sketch_mode <- FALSE
      sketch_ids <- NULL
    }
  }

  fit_ids <- if (sketch_mode) sketch_ids else cell_ids
  ## If sketch cells were dropped as non-finite, fall back to direct / auto-sketch
  if (sketch_mode && length(fit_ids) <= n_neighbors) {
    warning(
      "Too few finite sketch cells for UMAP; fitting on all finite cells instead."
    )
    sketch_mode <- FALSE
    fit_ids <- cell_ids
  }
  if (length(fit_ids) <= n_neighbors) {
    stop(
      "Need more than n_neighbors (", n_neighbors,
      ") cells to fit UMAP; have ", length(fit_ids), "."
    )
  }

  rep <- .ic2_banner_start(rep, c(
    paste0("Modality: ", modality, "   |   input: PCA dims ", paste(range(dims), collapse = "-")),
    paste0(
      "Cells: ", .ic2_fmt_num(n_cells),
      "   |   components: ", n_components,
      "   |   neighbors: ", n_neighbors
    ),
    paste0(
      "Engine: sketch -> project on SODA" ,
      if (!sketch_mode) " (direct fit; all cells)" else paste0(
        "   |   fit cells: ", .ic2_fmt_num(length(fit_ids))
      )
    )
  ))

  set.seed(my.seed)
  umap_args <- list(
    n_neighbors = n_neighbors,
    n_components = n_components,
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    scale = scale,
    init = init,
    init_sdev = init_sdev,
    spread = spread,
    min_dist = min_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    nn_method = nn_method,
    n_trees = n_trees,
    search_k = search_k,
    approx_pow = approx_pow,
    y = if (!sketch_mode) y else NULL,
    target_n_neighbors = target_n_neighbors,
    target_metric = target_metric,
    target_weight = target_weight,
    pca = pca,
    pca_center = pca_center,
    pcg_rand = pcg_rand,
    fast_sgd = fast_sgd,
    n_threads = n_threads,
    n_sgd_threads = n_sgd_threads,
    grain_size = grain_size,
    tmpdir = tmpdir,
    verbose = uwot_verbose
  )

  ## Classic iCellR embedding names (not Seurat-style UMAP_1)
  colnames_umap <- paste0("V", seq_len(n_components))
  emb_umap <- matrix(
    NA_real_,
    nrow = n_cells,
    ncol = n_components,
    dimnames = list(cell_ids, colnames_umap)
  )

  model <- NULL
  ## Overall stages: prepare  ->  fit  ->  project/write
  n_stages <- if (sketch_mode) 3L else 2L
  rep <- .ic2_progress_bar(rep, 0L, n_stages, prefix = "  embed")

  if (sketch_mode) {
    .ic2_step(rep, paste0(
      "[1/2] Building manifold model on sketch (", .ic2_fmt_num(length(fit_ids)), " cells)"
    ))
    rep <- .ic2_progress_bar(rep, 1L, n_stages, prefix = "  embed")
    X_fit <- X_all[fit_ids, , drop = FALSE]
    model <- do.call(
      uwot::umap,
      c(list(X = X_fit, ret_model = TRUE), umap_args)
    )
    emb_fit <- model$embedding
    rownames(emb_fit) <- fit_ids
    colnames(emb_fit) <- colnames_umap
    emb_umap[fit_ids, ] <- emb_fit

    project_ids <- setdiff(cell_ids, fit_ids)
    .ic2_step(rep, paste0(
      "[2/2] Projecting ", .ic2_fmt_num(length(project_ids)),
      " cells onto sketch model (chunks of ", .ic2_fmt_num(chunk.size), ")"
    ))
    rep <- .ic2_progress_bar(rep, 2L, n_stages, prefix = "  embed")
    if (length(project_ids) > 0L) {
      starts <- seq.int(1L, length(project_ids), by = chunk.size)
      n_chunk <- length(starts)
      for (i in seq_along(starts)) {
        s <- starts[[i]]
        e <- min(s + chunk.size - 1L, length(project_ids))
        ids <- project_ids[s:e]
        proj <- uwot::umap_transform(
          X = X_all[ids, , drop = FALSE],
          model = model,
          n_epochs = n_epochs,
          learning_rate = learning_rate,
          verbose = FALSE,
          n_threads = n_threads,
          n_sgd_threads = n_sgd_threads,
          grain_size = grain_size
        )
        rownames(proj) <- ids
        colnames(proj) <- colnames_umap
        emb_umap[ids, ] <- proj
        rep <- .ic2_progress_bar(rep, i, n_chunk, prefix = "  project")
      }
    } else {
      .ic2_step(rep, "Sketch covers all cells  -  no projection pass.")
    }
    rep <- .ic2_progress_bar(rep, 3L, n_stages, prefix = "  embed")
  } else {
    .ic2_step(rep, paste0(
      "[1/1] Embedding all ", .ic2_fmt_num(n_cells), " cells"
    ))
    rep <- .ic2_progress_bar(rep, 1L, n_stages, prefix = "  embed")
    emb_fit <- do.call(uwot::umap, c(list(X = X_all, ret_model = FALSE), umap_args))
    if (is.list(emb_fit) && !is.null(emb_fit$embedding)) {
      emb_fit <- emb_fit$embedding
    }
    emb_fit <- as.matrix(emb_fit)
    rownames(emb_fit) <- cell_ids
    colnames(emb_fit) <- colnames_umap
    emb_umap[, ] <- emb_fit
    rep <- .ic2_progress_bar(rep, 2L, n_stages, prefix = "  embed")
  }

  ## Persist: embeddings + JSON params (portable). Optional R-only .uwot.
  red_dir <- if (identical(modality, "rna")) {
    file.path(x@storage$root_dir, "reductions", out_name)
  } else {
    file.path(x@storage$root_dir, "modalities", modality, "reductions", out_name)
  }
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  model_path <- file.path(red_dir, "umap_model.uwot")
  meta_path <- file.path(red_dir, "reduction.json")

  .ic2_h5_write_dense(emb_umap, emb_path)
  unlink(c(
    file.path(red_dir, "cell_embeddings.rds"),
    file.path(red_dir, "umap_meta.rds"),
    file.path(red_dir, "umap_model.rds")
  ), force = TRUE)

  saved_model_path <- NULL
  if (isTRUE(save.model) && !is.null(model)) {
    .ic2_soda_write_uwot_model(model, model_path)
    saved_model_path <- model_path
  } else {
    unlink(c(model_path, file.path(red_dir, "umap_model.rds")), force = TRUE)
    if (!is.null(model) && is.list(model) && !is.null(model$mod_dir)) {
      try(uwot::unload_uwot(model), silent = TRUE)
    }
  }

  umap_meta <- list(
    name = out_name,
    mark = src$mark,
    aligned = isTRUE(src$aligned),
    source = src$name,
    based_on = src$name,
    modality = modality,
    dims = as.integer(dims),
    n_components = as.integer(n_components),
    n_neighbors = as.integer(n_neighbors),
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    scale = isTRUE(scale),
    init = init,
    spread = spread,
    min_dist = min_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    nn_method = nn_method,
    n_trees = as.integer(n_trees),
    search_k = search_k,
    sketch_mode = sketch_mode,
    n_fit_cells = length(fit_ids),
    n_cells = n_cells,
    seed = my.seed,
    save.model = isTRUE(save.model),
    created = as.character(Sys.time()),
    paths = list(
      embeddings = emb_path,
      model = saved_model_path
    )
  )
  .ic2_soda_write_json(umap_meta, meta_path)

  umap_ref <- list(
    cell_embeddings = emb_path,
    model = saved_model_path,
    dims = dims,
    n_components = n_components,
    n_neighbors = n_neighbors,
    metric = metric,
    min_dist = min_dist,
    sketch_mode = sketch_mode,
    meta_path = meta_path
  )
  umap_ref$aligned <- isTRUE(src$aligned)
  umap_ref$mark <- src$mark
  umap_ref$source <- src$name
  umap_ref$based_on <- src$name
  x@modalities[[modality]]$reductions[[out_name]] <- umap_ref
  x@modalities[[modality]]$active_reduction <- out_name
  if (identical(modality, "rna")) {
    x@storage$reduction_paths[[out_name]] <- red_dir
  }
  try(.ic2_soda_write_manifest(x), silent = TRUE)

  if (is.null(x@params$dimensionality_reduction)) {
    x@params$dimensionality_reduction <- list()
  }
  x@params$dimensionality_reduction[[modality]][[out_name]] <- list(
    source = src$name,
    based_on = src$name,
    aligned = isTRUE(src$aligned),
    mark = src$mark,
    dims = dims,
    n_components = n_components,
    n_neighbors = n_neighbors,
    metric = metric,
    min_dist = min_dist,
    sketch_mode = sketch_mode,
    n_fit_cells = length(fit_ids),
    save.model = isTRUE(save.model)
  )

  x <- .ic2_log_command(
    x, "run.umap",
    arguments = list(
      modality = modality, dims = dims, n_components = n_components,
      n_neighbors = n_neighbors, use.sketch = use.sketch,
      sketch_mode = sketch_mode, seed = my.seed,
      save.model = isTRUE(save.model)
    ),
    outputs = list(
      embeddings = emb_path,
      model = saved_model_path %||% NA_character_,
      n_cells = n_cells
    )
  )

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
    rep <- .ic2_add_output(rep, "Log", file.path(x@storage$root_dir, "iCellR2.log"))
  }
  rep <- .ic2_add_output(rep, "UMAP", emb_path)
  if (!is.null(saved_model_path)) rep <- .ic2_add_output(rep, "Model", saved_model_path)

  wall <- as.numeric(difftime(Sys.time(), rep$wall0, units = "secs"))
  if (is.finite(wall) && wall > 0) {
    rep$stats$cells_per_sec <- n_cells / wall
  }

  .ic2_finish(rep, extra_lines = paste0(
    "Stored V1..V", n_components, "   |   cells: ", .ic2_fmt_num(n_cells),
    if (sketch_mode) "   |   sketch -> project" else "   |   direct",
    "   |   embeddings only"
  ))
  x
}
