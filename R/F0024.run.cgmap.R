#' Internal cluster-guided MAP engine
#'
#' Used by \code{\link{run.cgpca}} and \code{\link{run.cgumap}}.
#' Not exported  -  prefer those wrappers.
#'
#' Cluster-guided MAP (cgMAP)
#'
#' Build a cluster-aware layout from an existing reduction plus
#' \code{cluster_id}:
#' \enumerate{
#'   \item \strong{input}: PCA dims (default \code{1:3}) or UMAP dims
#'     (default \code{1:2})
#'   \item Append \code{cluster_id} (optionally scaled / weighted)
#'   \item \strong{method}: run PCA (fast, up to 3D) or UMAP (nonlinear)
#'     on that guided matrix
#' }
#'
#' Internal engine for \code{\link{run.cgpca}} and \code{\link{run.cgumap}}.
#'
#' @param x An \code{iCellR2} object after clustering.
#' @param modality Modality. Default \code{"rna"}.
#' @param input Source embeddings: \code{"pca"} or \code{"umap"}.
#' @param dims Dimensions from \code{input}. Default \code{NULL}  -> 
#'   \code{1:3} for PCA, \code{1:2} for UMAP.
#' @param method Embedding engine on the guided matrix:
#'   \code{"pca"} (default, fast) or \code{"umap"} (slower, nonlinear).
#' @param npcs Output components when \code{method = "pca"}. Default \code{3}.
#' @param scale If \code{TRUE} (default), z-score input dims and cluster
#'   before embedding.
#' @param cluster.weight Multiplier on the cluster column after optional
#'   scaling. Default \code{1}.
#' @param cluster.col Metadata column. Default \code{"cluster_id"}.
#' @param center Passed to \code{prcomp} when \code{method = "pca"} and
#'   \code{scale = FALSE}.
#' @param my.seed RNG seed for the UMAP path.
#' @param use.sketch,sketch.name,max.cells.direct,chunk.size
#'   Sketch controls for \code{method = "umap"} (full until 5M).
#' @param n_neighbors,n_components,metric,n_epochs,learning_rate,init,spread,min_dist,n_threads,n_sgd_threads,grain_size,tmpdir
#'   Passed to \code{uwot} when \code{method = "umap"}.
#' @param reduction.name SODA / slot name for the reduction.
#'   Default \code{"cgmap"}. Wrappers may use \code{"cgpca"} / \code{"cgumap"}.
#' @param verbose Progress style.
#'
#' @return Updated \code{iCellR2} with a \code{cgmap} (or named) reduction.
#'   Plot with \code{cluster.plot(x, plot.type = "cgmap")}.
#'
#' @examples
#' \dontrun{
#' # Prefer public wrappers:
#' obj <- run.cgpca(obj)
#' obj <- run.cgumap(obj)
#' cluster.plot(obj, plot.type = "cgpca")
#' cluster.plot(obj, plot.type = "cgumap")
#' }
#'
#' @keywords internal
.ic2_run_cgmap <- function(
    x,
    modality = "rna",
    input = c("pca", "umap"),
    dims = NULL,
    method = c("pca", "umap"),
    npcs = 3L,
    scale = TRUE,
    cluster.weight = 1,
    cluster.col = "cluster_id",
    center = TRUE,
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
    init = "spectral",
    spread = 1,
    min_dist = 0.01,
    n_threads = 1,
    n_sgd_threads = 0,
    grain_size = 1,
    tmpdir = tempdir(),
    reduction.name = "cgmap",
    use = "auto",
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  input <- match.arg(input)
  method <- match.arg(method)
  src <- .ic2_resolve_use(x, modality, input, use)
  input <- src$name
  if (isTRUE(src$aligned)) {
    reduction.name <- .ic2_aligned_out_name(reduction.name, TRUE)
  }
  reduction.name <- as.character(reduction.name)[1L]
  if (!nzchar(reduction.name)) stop("'reduction.name' must be non-empty.")
  npcs <- as.integer(npcs)
  max.cells.direct <- as.integer(max.cells.direct)
  chunk.size <- as.integer(chunk.size)
  n_components <- as.integer(n_components)
  n_neighbors <- as.integer(n_neighbors)
  cluster.weight <- as.numeric(cluster.weight)[1L]
  if (!is.finite(cluster.weight) || cluster.weight == 0) {
    stop("'cluster.weight' must be a non-zero finite number.")
  }
  if (identical(method, "umap") && !requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for method = 'umap'.")
  }

  cmd <- paste0("run.", reduction.name)
  if (!startsWith(cmd, "run.")) cmd <- "run.cgmap"
  ## Prefer a stable command name for the unified API
  if (identical(reduction.name, "cgmap")) cmd <- "run.cgmap"
  if (identical(reduction.name, "cgpca")) cmd <- "run.cgpca"
  if (identical(reduction.name, "cgumap")) cmd <- "run.cgumap"

  rep <- .ic2_new_reporter(cmd, verbose, soda_root = x@storage$root_dir)
  uwot_verbose <- isTRUE(getOption("iCellR2.uwot.verbose", FALSE))

  if (!cluster.col %in% colnames(x@cell_data)) {
    stop(
      "No '", cluster.col, "' in cell_data. Run iclust() / run.clustering() first."
    )
  }

  .ic2_step(rep, paste0("[1/3] Loading ", toupper(input), " embeddings + cluster"))
  if (!grepl("umap", input)) {
    emb_src <- .ic2_get_plot_embedding(x, modality, input)
    default_dims <- 1:3
    dim_lab <- "PC"
  } else {
    emb_src <- .ic2_get_plot_embedding(x, modality, input)
    default_dims <- 1:2
    dim_lab <- "UMAP"
  }

  cell_ids <- rownames(x@cell_data)
  common <- intersect(cell_ids, rownames(emb_src))
  if (length(common) < 3L) {
    stop(toupper(input), " embeddings do not match cell_data.")
  }
  if (length(common) < length(cell_ids)) {
    warning(
      length(cell_ids) - length(common),
      " cells lack ", input, " embeddings and will be skipped."
    )
  }
  cell_ids <- common
  emb_src <- emb_src[cell_ids, , drop = FALSE]
  n_cells <- length(cell_ids)
  nd_avail <- ncol(emb_src)

  if (is.null(dims)) dims <- default_dims
  dims <- as.integer(dims)
  dims <- dims[dims >= 1L & dims <= nd_avail]
  if (length(dims) < 1L) {
    stop("No valid '", input, "' dims selected (available: ", nd_avail, ").")
  }
  X_in <- emb_src[, dims, drop = FALSE]
  storage.mode(X_in) <- "double"

  cl <- x@cell_data[cell_ids, cluster.col]
  cl <- as.numeric(as.character(cl))
  if (anyNA(cl)) {
    stop(
      sum(is.na(cl)), " cells have missing '", cluster.col,
      "'. Re-run clustering or drop those cells."
    )
  }

  colnames_in <- c(
    if (!is.null(colnames(X_in)) && all(nzchar(colnames(X_in)))) {
      colnames(X_in)
    } else {
      paste0(dim_lab, dims)
    },
    cluster.col
  )

  if (isTRUE(scale)) {
    X_in <- scale(X_in, center = TRUE, scale = TRUE)
    X_in[!is.finite(X_in)] <- 0
    cl_z <- as.numeric(scale(cl, center = TRUE, scale = TRUE))
    cl_z[!is.finite(cl_z)] <- 0
    cl_col <- cl_z * cluster.weight
  } else {
    cl_col <- cl * cluster.weight
  }
  X_all <- cbind(X_in, cluster = cl_col)
  colnames(X_all) <- colnames_in
  if (any(!is.finite(X_all))) X_all[!is.finite(X_all)] <- 0

  model <- NULL
  sketch_mode <- FALSE
  fit_ids <- cell_ids
  rotation <- NULL
  sdev <- NULL

  if (identical(method, "pca")) {
    n_out <- min(as.integer(npcs), ncol(X_all), max(1L, n_cells - 1L))
    rep <- .ic2_banner_start(rep, c(
      paste0(
        "Input: ", toupper(input), " dims ",
        paste(range(dims), collapse = "-"), " + ", cluster.col
      ),
      paste0(
        "Cells: ", .ic2_fmt_num(n_cells),
        "   |   features: ", ncol(X_all),
        "   |   npcs: ", n_out
      ),
      paste0(
        "method=pca   |   scale=", scale,
        "   |   cluster.weight=", cluster.weight
      )
    ))
    .ic2_step(rep, "[2/3] PCA on guided matrix")
    if (isTRUE(scale)) {
      fit <- stats::prcomp(X_all, center = FALSE, scale. = FALSE, rank. = n_out)
    } else {
      fit <- stats::prcomp(
        X_all,
        center = isTRUE(center),
        scale. = FALSE,
        rank. = n_out
      )
    }
    k <- min(n_out, ncol(fit$x), length(fit$sdev))
    if (k < 1L) stop("cgMAP PCA produced no components.")
    emb_out <- fit$x[, seq_len(k), drop = FALSE]
    colnames(emb_out) <- paste0("V", seq_len(k))
    rownames(emb_out) <- cell_ids
    sdev <- as.numeric(fit$sdev[seq_len(k)])
    rotation <- fit$rotation[, seq_len(k), drop = FALSE]
    colnames(rotation) <- colnames(emb_out)
    rownames(rotation) <- colnames_in
    n_components <- k
  } else {
    sketch_ids <- NULL
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
      sketch_ids <- intersect(sketch_ids, cell_ids)
      if (length(sketch_ids) >= max(n_neighbors + 1L, 50L) &&
          length(sketch_ids) < n_cells) {
        sketch_mode <- TRUE
      } else if (length(sketch_ids) >= n_cells) {
        sketch_mode <- FALSE
        sketch_ids <- NULL
      } else if (force_sketch || n_cells > max.cells.direct) {
        set.seed(my.seed)
        nsk <- min(100000L, n_cells - 1L)
        nsk <- max(nsk, max(n_neighbors + 2L, 50L))
        nsk <- min(nsk, n_cells - 1L)
        sketch_ids <- cell_ids[sort(sample.int(n_cells, nsk))]
        sketch_mode <- TRUE
        .ic2_step(rep, paste0(
          "No sketch found  -  auto-sketching ", .ic2_fmt_num(nsk),
          " / ", .ic2_fmt_num(n_cells), " cells for cgMAP fit"
        ))
      } else if (isTRUE(use.sketch) && n_cells <= max.cells.direct) {
        sketch_mode <- FALSE
        sketch_ids <- NULL
      }
    }

    fit_ids <- if (sketch_mode) sketch_ids else cell_ids
    if (sketch_mode && length(fit_ids) <= n_neighbors) {
      warning("Too few sketch cells for cgMAP; fitting on all cells instead.")
      sketch_mode <- FALSE
      fit_ids <- cell_ids
    }
    if (length(fit_ids) <= n_neighbors) {
      stop(
        "Need more than n_neighbors (", n_neighbors,
        ") cells to fit cgMAP; have ", length(fit_ids), "."
      )
    }

    rep <- .ic2_banner_start(rep, c(
      paste0(
        "Input: ", toupper(input), " dims ",
        paste(range(dims), collapse = "-"), " + ", cluster.col
      ),
      paste0(
        "Cells: ", .ic2_fmt_num(n_cells),
        "   |   features: ", ncol(X_all),
        "   |   components: ", n_components
      ),
      paste0(
        "method=umap   |   scale=", scale,
        "   |   cluster.weight=", cluster.weight,
        if (!sketch_mode) "   |   direct" else paste0(
          "   |   sketch fit: ", .ic2_fmt_num(length(fit_ids))
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
      scale = FALSE,
      init = init,
      spread = spread,
      min_dist = min_dist,
      n_threads = n_threads,
      n_sgd_threads = n_sgd_threads,
      grain_size = grain_size,
      tmpdir = tmpdir,
      verbose = uwot_verbose
    )
    colnames_out <- paste0("V", seq_len(n_components))
    emb_out <- matrix(
      NA_real_, nrow = n_cells, ncol = n_components,
      dimnames = list(cell_ids, colnames_out)
    )
    n_stages <- if (sketch_mode) 3L else 2L
    rep <- .ic2_progress_bar(rep, 0L, n_stages, prefix = "  embed")
    .ic2_step(rep, "[2/3] UMAP on guided matrix")
    if (sketch_mode) {
      rep <- .ic2_progress_bar(rep, 1L, n_stages, prefix = "  embed")
      model <- do.call(
        uwot::umap,
        c(list(X = X_all[fit_ids, , drop = FALSE], ret_model = TRUE), umap_args)
      )
      emb_fit <- model$embedding
      rownames(emb_fit) <- fit_ids
      colnames(emb_fit) <- colnames_out
      emb_out[fit_ids, ] <- emb_fit
      project_ids <- setdiff(cell_ids, fit_ids)
      .ic2_step(rep, paste0(
        "  Projecting ", .ic2_fmt_num(length(project_ids)),
        " cells (chunks of ", .ic2_fmt_num(chunk.size), ")"
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
          colnames(proj) <- colnames_out
          emb_out[ids, ] <- proj
          rep <- .ic2_progress_bar(rep, i, n_chunk, prefix = "  project")
        }
      }
      rep <- .ic2_progress_bar(rep, 3L, n_stages, prefix = "  embed")
    } else {
      .ic2_step(rep, paste0(
        "  Embedding all ", .ic2_fmt_num(n_cells), " cells"
      ))
      rep <- .ic2_progress_bar(rep, 1L, n_stages, prefix = "  embed")
      emb_fit <- do.call(uwot::umap, c(list(X = X_all, ret_model = FALSE), umap_args))
      if (is.list(emb_fit) && !is.null(emb_fit$embedding)) {
        emb_fit <- emb_fit$embedding
      }
      emb_fit <- as.matrix(emb_fit)
      rownames(emb_fit) <- cell_ids
      colnames(emb_fit) <- colnames_out
      emb_out[, ] <- emb_fit
      rep <- .ic2_progress_bar(rep, 2L, n_stages, prefix = "  embed")
    }
  }

  if (any(!is.finite(emb_out))) emb_out[!is.finite(emb_out)] <- 0

  .ic2_step(rep, paste0("[3/3] Writing ", reduction.name, " reduction"))
  red_dir <- file.path(x@storage$root_dir, "reductions", reduction.name)
  dir.create(red_dir, recursive = TRUE, showWarnings = FALSE)
  emb_path <- file.path(red_dir, "cell_embeddings.h5")
  load_path <- file.path(red_dir, "feature_loadings.h5")
  stdev_path <- file.path(red_dir, "stdev.tsv")
  model_path <- file.path(red_dir, paste0(reduction.name, "_model.uwot"))
  meta_path <- file.path(red_dir, "reduction.json")

  .ic2_h5_write_dense(emb_out, emb_path)
  unlink(c(
    file.path(red_dir, "cell_embeddings.rds"),
    file.path(red_dir, paste0(reduction.name, "_meta.rds")),
    file.path(red_dir, paste0(reduction.name, "_model.rds"))
  ), force = TRUE)
  if (!is.null(rotation)) {
    .ic2_h5_write_dense(rotation, load_path)
    unlink(file.path(red_dir, "feature_loadings.rds"), force = TRUE)
  } else if (file.exists(load_path)) unlink(load_path)
  if (!is.null(sdev)) {
    utils::write.table(
      data.frame(pc = seq_along(sdev), stdev = as.numeric(sdev), stringsAsFactors = FALSE),
      file = stdev_path, sep = "\t", quote = FALSE, row.names = FALSE
    )
    unlink(file.path(red_dir, "stdev.rds"), force = TRUE)
  } else if (file.exists(stdev_path)) unlink(stdev_path)
  if (!is.null(model)) {
    ## In-memory model used for projection only; do not archive .uwot by default
    unlink(c(
      model_path,
      file.path(red_dir, paste0(reduction.name, "_model.rds")),
      file.path(red_dir, paste0(reduction.name, "_model.uwot"))
    ), force = TRUE)
    if (is.list(model) && !is.null(model$mod_dir)) {
      try(uwot::unload_uwot(model), silent = TRUE)
    }
    model_path <- NULL
  } else if (file.exists(model_path)) {
    unlink(model_path)
    model_path <- NULL
  }

  meta <- list(
    name = reduction.name,
    modality = modality,
    input = input,
    based_on = input,
    method = method,
    dims = dims,
    n_components = ncol(emb_out),
    scale = isTRUE(scale),
    cluster.weight = cluster.weight,
    cluster.col = cluster.col,
    sketch_mode = sketch_mode,
    n_fit_cells = length(fit_ids),
    n_cells = n_cells,
    seed = my.seed,
    created = as.character(Sys.time()),
    paths = list(
      embeddings = emb_path,
      loadings = if (!is.null(rotation)) load_path else NULL,
      stdev = if (!is.null(sdev)) stdev_path else NULL,
      model = NULL
    )
  )
  .ic2_soda_write_json(meta, meta_path)

  red_ref <- list(
    cell_embeddings = emb_path,
    feature_loadings = if (!is.null(rotation)) load_path else NULL,
    stdev = sdev,
    model = NULL,
    input = input,
    method = method,
    dims = dims,
    n_components = ncol(emb_out),
    scale = isTRUE(scale),
    cluster.weight = cluster.weight,
    sketch_mode = sketch_mode,
    meta_path = meta_path
  )
  ## Legacy field names for older plot / meta readers
  if (identical(input, "pca")) red_ref$pca_dims <- dims
  if (identical(input, "umap")) red_ref$umap_dims <- dims
  if (identical(method, "pca")) red_ref$npcs <- ncol(emb_out)

  x@modalities[[modality]]$reductions[[reduction.name]] <- red_ref
  x@modalities[[modality]]$active_reduction <- reduction.name
  if (is.null(x@storage$reduction_paths)) x@storage$reduction_paths <- list()
  x@storage$reduction_paths[[reduction.name]] <- red_dir

  if (is.null(x@params$dimensionality_reduction)) {
    x@params$dimensionality_reduction <- list()
  }
  x@params$dimensionality_reduction[[modality]][[reduction.name]] <- list(
    input = input,
    method = method,
    dims = dims,
    n_components = ncol(emb_out),
    scale = isTRUE(scale),
    cluster.weight = cluster.weight,
    cluster.col = cluster.col,
    sketch_mode = sketch_mode
  )

  x <- .ic2_log_command(
    x, cmd,
    arguments = list(
      modality = modality, input = input, method = method, dims = dims,
      n_components = ncol(emb_out), scale = isTRUE(scale),
      cluster.weight = cluster.weight, cluster.col = cluster.col,
      reduction.name = reduction.name, sketch_mode = sketch_mode
    ),
    outputs = list(
      embeddings = emb_path,
      n_cells = n_cells,
      n_components = ncol(emb_out)
    )
  )

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  rep <- .ic2_add_output(rep, "Embeddings", emb_path)
  if (!is.null(rotation)) rep <- .ic2_add_output(rep, "Loadings", load_path)
  if (!is.null(rotation)) rep <- .ic2_add_output(rep, "Loadings", load_path)

  .ic2_finish(rep, extra_lines = paste0(
    "V1..V", ncol(emb_out),
    "   |   cells: ", .ic2_fmt_num(n_cells),
    "   |   input=", input,
    "   |   method=", method,
    "   |   scale=", scale
  ))
  x
}
