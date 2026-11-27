#' Aggregate multiple samples with flexible output formats
#'
#' For large projects use \code{output = "virtual"} (default): returns a
#' logical aggregation of **assay handles / paths and maps only**  -  no
#' merged mega-matrix is built in RAM.
#'
#' Inputs may be sparse matrices, assay handles, or \code{iCellR2_sample}
#' objects from \code{\link{load10x}}. Prefer setting a unique
#' \code{sample.id} in each \code{load10x()} call; then aggregation can
#' omit \code{condition.names} and use those ids for barcode prefixes /
#' conditions.
#'
#' @param samples List of matrices, handles, or virtual samples.
#' @param condition.names Optional condition labels (same length as
#'   \code{samples}). If \code{NULL} (default), uses each sample's
#'   \code{sample.id} from \code{load10x} / \code{sample.ids}.
#' @param modality Modality string. Default \code{"rna"}.
#' @param sample.ids Optional sample ids. Default \code{NULL} uses each
#'   sample's \code{sample_id} (required from \code{load10x}).
#' @param output One of \code{"virtual"}, \code{"sparse"}, \code{"matrix"},
#'   \code{"data.frame"}, \code{"hdf5"}. Default \code{"virtual"}.
#' @param h5.file Required for physical \code{"hdf5"} merge; for
#'   \code{"virtual"}, optional directory prefix for per-sample HDF5 files.
#' @param dataset.name HDF5 group name. Default \code{"counts"}.
#' @param overwrite Logical.
#' @param out.dir Directory for per-sample \strong{10x HDF5} files
#'   (\code{*_tenx.h5} + \code{*_tenx_meta.json}) when materializing
#'   \strong{matrix} inputs. Unused when every sample is already a virtual
#'   handle from \code{\link{load10x}}. Default \code{"iCellR2_agg_tmp"}.
#' @param verbose Progress style.
#'
#' @return Aggregation descriptor or merged matrix / handle.
#'
#' @examples
#' \dontrun{
#' s1 <- load10x(path1, sample.id = "my1.3m")
#' s2 <- load10x(path2, sample.id = "my340k")
#' agg <- data.aggregation(samples = list(s1, s2), output = "virtual")
#' }
#'
#' @importFrom methods as
#' @importFrom Matrix Matrix
#' @export
data.aggregation <- function(
    samples = NULL,
    condition.names = NULL,
    modality = "rna",
    sample.ids = NULL,
    output = c("virtual", "sparse", "matrix", "data.frame", "hdf5"),
    h5.file = NULL,
    dataset.name = "counts",
    overwrite = FALSE,
    out.dir = "iCellR2_agg_tmp",
    verbose = TRUE
) {
  output <- match.arg(output)
  rep <- .ic2_new_reporter("data.aggregation", verbose)

  if (is.null(samples) || !is.list(samples) || length(samples) < 2L) {
    stop("'samples' must be a list with at least two elements.")
  }
  if (!is.character(modality) || length(modality) != 1L || !nzchar(modality)) {
    stop("'modality' must be a non-empty character string.")
  }
  modality <- tolower(modality)

  if (is.null(sample.ids)) {
    nm <- names(samples)
    sample.ids <- vapply(seq_along(samples), function(i) {
      s <- samples[[i]]
      sid <- NULL
      if (inherits(s, "iCellR2_sample")) {
        sid <- s$sample_id %||% s$handle$sample_id
      } else if (.ic2_is_assay_handle(s)) {
        sid <- s$sample_id
      }
      ## Named list: list(WT = wt, KO = ko, ...)  ->  sample.ids from names
      if ((is.null(sid) || !nzchar(as.character(sid)[1L])) &&
          !is.null(nm) && length(nm) >= i && nzchar(nm[[i]])) {
        sid <- nm[[i]]
      }
      if (is.null(sid) || !nzchar(as.character(sid)[1L])) {
        stop(
          "Sample ", i, " has no sample.id. ",
          "Pass sample.ids=, use a named list, or set sample.id in load10x()."
        )
      }
      as.character(sid)[1L]
    }, character(1))
  }
  if (!is.character(sample.ids) || length(sample.ids) != length(samples)) {
    stop("'sample.ids' must match length of 'samples'.")
  }
  if (any(is.na(sample.ids)) || any(!nzchar(sample.ids))) {
    stop("All 'sample.ids' must be non-empty.")
  }
  sample.ids <- .ic2_clean_names(sample.ids)
  if (anyDuplicated(sample.ids)) {
    stop(
      "Duplicate sample.ids after cleaning: ",
      paste(unique(sample.ids[duplicated(sample.ids)]), collapse = ", "),
      ". Choose unique sample.id values in load10x()."
    )
  }

  if (is.null(condition.names)) {
    condition.names <- sample.ids
  }
  if (!is.character(condition.names)) {
    stop("'condition.names' must be a character vector or NULL.")
  }
  if (length(samples) != length(condition.names)) {
    stop("'samples' and 'condition.names' must have the same length.")
  }
  if (any(is.na(condition.names)) || any(!nzchar(condition.names))) {
    stop("All entries in 'condition.names' must be non-empty strings.")
  }
  condition.names <- .ic2_clean_names(condition.names)

  if (output == "hdf5") {
    if (is.null(h5.file) || !is.character(h5.file) || length(h5.file) != 1L || !nzchar(h5.file)) {
      stop("'h5.file' must be provided when output = 'hdf5'.")
    }
  }

  ## Physical merge is not atlas-safe  -  refuse tenx_h5 / huge samples
  if (!identical(output, "virtual")) {
    n_guess <- 0
    has_tenx <- FALSE
    for (s in samples) {
      h <- NULL
      if (inherits(s, "iCellR2_sample") && !is.null(s$handle)) h <- s$handle
      else if (.ic2_is_assay_handle(s)) h <- s
      if (!is.null(h)) {
        h <- .ic2_enrich_handle(h, load_dimnames = FALSE)
        if (identical(h$backend, "tenx_h5")) has_tenx <- TRUE
        n_guess <- n_guess + as.integer(h$n_cells %||% 0L)
      } else if (!is.null(ncol(s))) {
        n_guess <- n_guess + as.integer(ncol(s))
      }
    }
    if (isTRUE(has_tenx) || n_guess > 200000L) {
      stop(
        "output = '", output, "' would merge assays in memory and is not safe for ",
        "native tenx_h5 and/or large cell counts (~", .ic2_fmt_num(n_guess), "). ",
        "Use output = \"virtual\" (default)."
      )
    }
  }

  clean_names <- function(x) {
    x <- gsub("_", ".", x, fixed = TRUE)
    x <- gsub("-", ".", x, fixed = TRUE)
    x <- gsub("/", ".or.", x, fixed = TRUE)
    x <- gsub(" ", ".", x, fixed = TRUE)
    make.names(x, unique = TRUE)
  }

  resolve_one <- function(obj, sample_id, condition_id, dest_h5) {
    if (inherits(obj, "iCellR2_sample")) {
      if (!is.null(obj$handle) && .ic2_is_assay_handle(obj$handle)) {
        h <- .ic2_enrich_handle(obj$handle, load_dimnames = TRUE)
        return(list(handle = h, feature_ids = h$feature_ids, cell_ids = h$cell_ids))
      }
      if (!is.null(obj$assay)) obj <- obj$assay
    }
    if (.ic2_is_assay_handle(obj)) {
      h <- .ic2_enrich_handle(obj, load_dimnames = TRUE)
      return(list(handle = h, feature_ids = h$feature_ids, cell_ids = h$cell_ids))
    }

    mat <- .ic2_to_dgCMatrix(obj)
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop("Each sample must have rownames and colnames.")
    }
    rownames(mat) <- clean_names(rownames(mat))
    colnames(mat) <- clean_names(colnames(mat))
    ## Raw / agg temp: Cell Ranger-style 10x H5 (not legacy CSC)
    h <- .ic2_write_matrix_tenx_h5(
      mat,
      filepath = dest_h5,
      sample_id = sample_id,
      modality = modality,
      overwrite = TRUE
    )
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    rm(mat)
    invisible(gc)
    list(handle = h, feature_ids = h$feature_ids, cell_ids = h$cell_ids)
  }

  needs_materialize <- vapply(samples, function(s) {
    !(inherits(s, "iCellR2_sample") || .ic2_is_assay_handle(s))
  }, logical(1))
  if (any(needs_materialize)) {
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    rep <- .ic2_add_output(
      rep, "Agg dir",
      normalizePath(out.dir, winslash = "/", mustWork = FALSE)
    )
  }
  rep <- .ic2_set_io_dirs(
    rep,
    input_dir = getwd(),
    output_dir = if (!is.null(h5.file)) {
      h5.file
    } else if (any(needs_materialize)) {
      out.dir
    } else {
      getwd()
    }
  )

  handles <- vector("list", length(samples))
  feature_maps <- vector("list", length(samples))
  cell_maps <- vector("list", length(samples))
  input_types <- character(length(samples))
  all.features <- character(0)

  rep <- .ic2_banner_start(rep, c(
    paste0("Samples: ", length(samples), "   |   modality: ", modality, "   |   output: ", output),
    paste0("Conditions: ", paste(unique(condition.names), collapse = ", ")),
    if (!any(needs_materialize)) {
      "Virtual handles only  -  no agg temp files (data copied at make.obj)"
    } else {
      paste0("Materializing ", sum(needs_materialize), " matrix sample(s) under ", out.dir)
    }
  ))
  .ic2_step(rep, paste0("[1/2] Resolving ", length(samples), " sample handle(s)"))

  for (i in seq_along(samples)) {
    dest <- if (needs_materialize[i]) {
      file.path(out.dir, paste0(sample.ids[i], "_tenx.h5"))
    } else {
      NA_character_
    }
    if (.ic2_is_assay_handle(samples[[i]]) || inherits(samples[[i]], "iCellR2_sample")) {
      input_types[i] <- "handle"
    } else if (inherits(samples[[i]], "dgCMatrix") || inherits(samples[[i]], "Matrix")) {
      input_types[i] <- "sparse"
    } else if (is.data.frame(samples[[i]])) {
      input_types[i] <- "data.frame"
    } else if (is.matrix(samples[[i]])) {
      input_types[i] <- "matrix"
    } else {
      input_types[i] <- "matrix-like"
    }

    resolved <- resolve_one(samples[[i]], sample.ids[i], condition.names[i], dest)
    feature_names <- clean_names(resolved$feature_ids)
    cell_names <- clean_names(resolved$cell_ids)
    prefixed_cells <- paste(condition.names[i], cell_names, sep = "_")

    ## Update colnames for consistency (iCellR2 CSC only; tenx_h5 is read-only source)
    if (identical(resolved$handle$backend, "hdf5_sparse")) {
      h5 <- hdf5r::H5File$new(resolved$handle$path, mode = "r+")
      g <- h5[[resolved$handle$group]]
      if ("colnames" %in% names(g)) g$link_delete("colnames")
      g[["colnames"]] <- as.character(prefixed_cells)
      h5$close_all()
    } else if (identical(resolved$handle$backend, "tenx_h5")) {
      ## Do not mutate the original 10x H5  -  update sidecar meta only (JSON)
      mp <- resolved$handle$meta_path
      if (!is.null(mp) && file.exists(mp)) {
        meta <- .ic2_soda_read_tenx_meta(mp)
        if (!is.null(meta)) {
          meta$cell_ids <- as.character(prefixed_cells)
          resolved$handle$meta_path <- .ic2_soda_write_tenx_meta(meta, mp)
        }
      }
      resolved$handle$cell_ids <- as.character(prefixed_cells)
    }

    resolved$handle$sample_id <- sample.ids[i]
    resolved$handle$modality <- modality
    resolved$handle$n_features <- length(feature_names)
    resolved$handle$n_cells <- length(prefixed_cells)
    ## keep handles lightweight: dimnames live on disk / in maps
    resolved$handle$feature_ids <- NULL
    resolved$handle$cell_ids <- NULL
    handles[[i]] <- resolved$handle

    all.features <- union(all.features, feature_names)
    feature_maps[[i]] <- data.frame(
      sample_id = sample.ids[i],
      condition_id = condition.names[i],
      local_feature = feature_names,
      global_feature = feature_names,
      stringsAsFactors = FALSE
    )
    cell_maps[[i]] <- data.frame(
      sample_id = sample.ids[i],
      condition_id = condition.names[i],
      local_cell = cell_names,
      global_cell = prefixed_cells,
      stringsAsFactors = FALSE
    )
    rep <- .ic2_progress_bar(rep, i, length(samples), prefix = "  samples")
  }

  sample_data <- data.frame(
    sample_id = sample.ids,
    condition_id = condition.names,
    modality = modality,
    input_type = input_types,
    assay_path = vapply(handles, function(h) h$path, character(1)),
    assay_backend = vapply(handles, function(h) h$backend, character(1)),
    stringsAsFactors = FALSE
  )

  .ic2_step(rep, "[2/2] Building aggregation descriptor")
  n_cells_total <- sum(vapply(handles, function(h) as.integer(h$n_cells), integer(1)))

  if (output == "virtual") {
    agg <- list(
      type = "iCellR2_aggregation",
      modality = modality,
      created = as.character(Sys.time()),
      output = "virtual",
      handles = handles,
      sample_data = sample_data,
      condition.names = condition.names,
      sample.ids = sample.ids,
      feature_universe = all.features,
      feature_maps = feature_maps,
      cell_maps = cell_maps,
      n_samples = length(samples)
    )
    class(agg) <- c("iCellR2_aggregation", "list")
    .ic2_finish(rep, extra_lines = paste0(
      "Virtual agg  |  features: ", .ic2_fmt_num(length(all.features)),
      "   |   cells: ", .ic2_fmt_num(n_cells_total),
      "   |   samples: ", length(samples)
    ))
    return(agg)
  }

  ## Physical merge (small / medium data only)
  .ic2_step(rep, "Merging samples into one matrix (non-virtual)")
  aligned_samples <- lapply(seq_along(handles), function(i) {
    dn <- .ic2_handle_dimnames(handles[[i]])
    ## read full sample  -  intentional for non-virtual outputs
    if (identical(handles[[i]]$backend, "hdf5_sparse")) {
      x <- read_sparse_hdf5_cols(
        handles[[i]]$path,
        group = handles[[i]]$group,
        cols = seq_len(handles[[i]]$n_cells)
      )
    } else {
      x <- readRDS(handles[[i]]$path)
    }
    rownames(x) <- clean_names(rownames(x))
    missing.features <- setdiff(all.features, rownames(x))
    if (length(missing.features) > 0L) {
      zero.block <- Matrix::Matrix(
        0, nrow = length(missing.features), ncol = ncol(x), sparse = TRUE
      )
      rownames(zero.block) <- missing.features
      colnames(zero.block) <- colnames(x)
      x <- rbind(x, zero.block)
    }
    x[all.features, , drop = FALSE]
  })

  merged <- do.call(cbind, aligned_samples)
  merged <- .ic2_to_dgCMatrix(merged)

  if (output == "sparse") {
    .ic2_finish(rep, extra_lines = paste0(
      "Merged sparse: ", .ic2_fmt_num(nrow(merged)), "  x  ", .ic2_fmt_num(ncol(merged))
    ))
    return(merged)
  }
  if (output == "matrix") {
    .ic2_finish(rep)
    return(as.matrix(merged))
  }
  if (output == "data.frame") {
    .ic2_finish(rep)
    return(as.data.frame(as.matrix(merged), stringsAsFactors = FALSE))
  }

  out_h <- write_sparse_hdf5(merged, filepath = h5.file, group = dataset.name, overwrite = overwrite)
  rep <- .ic2_add_output(rep, "HDF5", normalizePath(h5.file, winslash = "/", mustWork = FALSE))
  .ic2_finish(rep, extra_lines = paste0(
    "Merged HDF5  |  features: ", .ic2_fmt_num(nrow(merged)),
    "   |   cells: ", .ic2_fmt_num(ncol(merged))
  ))
  out_h
}
