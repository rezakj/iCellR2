#' Sketch cells for scalable dimensionality reduction
#'
#' Draws a cell sketch (uniform or density-biased placeholder) and stores
#' sketch ids on disk / in the object. Downstream PCA/UMAP should run on
#' the sketch and project remaining cells  -  not naive full UMAP.
#'
#' The default path is base-R sampling over on-disk cell ids.
#'
#' @param x An \code{iCellR2} object.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Assay name; default active assay.
#' @param n.sketch Number of cells to sketch. Default \code{5000}.
#' @param method \code{"uniform"} (default). Reserved: \code{"geosketch"}
#'   when an optional backend provides it.
#' @param seed Random seed.
#' @param sketch.name Name for this sketch. Default \code{"sketch1"}.
#' @param verbose Progress style: \code{TRUE}/\code{"full"}, \code{"compact"},
#'   or \code{FALSE}/\code{"silent"}.
#'
#' @return Updated \code{iCellR2} with sketch metadata (no full embedding yet).
#'
#' @export
run.sketch <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    n.sketch = 5000,
    method = c("uniform", "geosketch"),
    seed = 1,
    sketch.name = "sketch1",
    verbose = TRUE
) {
  method <- match.arg(method)
  .ic2_require_object(x)
  modality <- tolower(modality)
  if (is.null(assay.name)) {
    assay.name <- x@modalities[[modality]]$active_assay
  }
  if (is.null(assay.name)) stop("No active assay for modality '", modality, "'.")

  rep <- .ic2_new_reporter("run.sketch", verbose, soda_root = x@storage$root_dir)

  cell_ids <- rownames(x@cell_data)
  n_cells <- length(cell_ids)
  n.sketch <- as.integer(min(as.integer(n.sketch), n_cells))
  if (n.sketch < 2L) stop("'n.sketch' must be >= 2.")

  rep <- .ic2_banner_start(rep, c(
    paste0("Modality: ", modality, "   |   assay: ", assay.name),
    paste0(
      "Cells: ", .ic2_fmt_num(n_cells),
      "   |   sketch: ", .ic2_fmt_num(n.sketch),
      "   |   method: ", method
    )
  ))

  .ic2_step(rep, "[1/2] Sampling sketch cells")
  if (identical(method, "geosketch")) {
    message("geosketch is not available in this build; using uniform sampling.")
    method_used <- "uniform"
  } else {
    method_used <- "uniform"
  }

  set.seed(seed)
  sketch_idx <- sort(sample.int(n_cells, size = n.sketch, replace = FALSE))
  sketch_ids <- cell_ids[sketch_idx]
  rep <- .ic2_progress_bar(rep, 1L, 2L, prefix = "  sketch")

  .ic2_step(rep, "[2/2] Writing sketch metadata to SODA")
  sketch_dir <- file.path(
    x@storage$root_dir, "modalities", modality, "sketches", sketch.name
  )
  dir.create(sketch_dir, recursive = TRUE, showWarnings = FALSE)
  sketch_tsv <- file.path(sketch_dir, "sketch_ids.tsv")
  utils::write.table(
    data.frame(cell_id = as.character(sketch_ids), stringsAsFactors = FALSE),
    file = sketch_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  sketch_meta <- list(
    sketch_name = sketch.name,
    method = method_used,
    n_sketch = n.sketch,
    n_cells = n_cells,
    seed = seed,
    modality = modality,
    assay.name = assay.name,
    created = as.character(Sys.time())
  )
  sketch_path <- .ic2_soda_write_json(sketch_meta, file.path(sketch_dir, "sketch_meta.json"))
  .ic2_soda_unlink_rds(file.path(sketch_dir, "sketch_ids.rds"))

  if (is.null(x@modalities[[modality]]$sketches)) {
    x@modalities[[modality]]$sketches <- list(
      sketch_ids = NULL, sketch_models = list(), projections = list()
    )
  }
  x@modalities[[modality]]$sketches$sketch_ids <- sketch_ids
  x@modalities[[modality]]$sketches$sketch_models[[sketch.name]] <- list(
    path = sketch_path,
    method = method_used,
    n_sketch = n.sketch,
    seed = seed
  )

  x@params$sketch[[modality]] <- list(
    sketch.name = sketch.name,
    method = method_used,
    n.sketch = n.sketch,
    seed = seed
  )

  x <- .ic2_log_command(
    x,
    "run.sketch",
    arguments = list(
      modality = modality, assay.name = assay.name,
      n.sketch = n.sketch, method = method_used, seed = seed,
      sketch.name = sketch.name
    ),
    outputs = list(n_sketch = n.sketch, path = sketch_path)
  )

  rep <- .ic2_progress_bar(rep, 2L, 2L, prefix = "  sketch")
  rep <- .ic2_add_output(rep, "sketch", sketch_path)
  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Log", file.path(x@storage$root_dir, "iCellR2.log"))
  }
  .ic2_finish(rep, extra_lines = paste0(
    "Sketched ", .ic2_fmt_num(n.sketch), " / ", .ic2_fmt_num(n_cells),
    " cells (", method_used, ")"
  ))

  x
}

#' Placeholder: project full dataset from a sketch model
#'
#' Scaffold for sketch-based projection (PCA/UMAP). Full implementations
#' will load sketch models and project remaining cells without running
#' naive full UMAP on all cells.
#'
#' @param x An \code{iCellR2} object.
#' @param sketch.name Sketch name created by \code{run.sketch()}.
#' @param modality Modality. Default \code{"rna"}.
#' @param verbose Progress style.
#' @export
project.from.sketch <- function(
    x,
    sketch.name = "sketch1",
    modality = "rna",
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  sk <- x@modalities[[modality]]$sketches$sketch_models[[sketch.name]]
  if (is.null(sk)) {
    stop(
      "Sketch '", sketch.name, "' not found. Run run.sketch() first."
    )
  }
  rep <- .ic2_new_reporter("project.from.sketch", verbose, soda_root = x@storage$root_dir)
  rep <- .ic2_banner_start(rep, c(
    paste0("Sketch: ", sketch.name, "   |   modality: ", modality),
    paste0("Registered cells: ", .ic2_fmt_num(sk$n_sketch))
  ))
  .ic2_finish(rep, extra_lines = paste0(
    "Scaffold only  -  PCA/UMAP projection from sketch comes next; ",
    "use sketch cell ids for scalable reduction."
  ))
  invisible(x)
}
