#' Load 10x Genomics count data with flexible input and output formats
#'
#' Reads either:
#' \itemize{
#'   \item a 10x Genomics directory containing `matrix.mtx`, `barcodes.tsv`,
#'   and either `genes.tsv` or `features.tsv` (optionally gzipped), or
#'   \item a 10x Genomics HDF5 file such as `filtered_feature_bc_matrix.h5`.
#' }
#'
#' Default \code{output = "virtual"}:
#' \itemize{
#'   \item \strong{10x HDF5 input:} native \code{tenx_h5} handle pointing at the
#'   original file (no matrix copy).
#'   \item \strong{10x directory (MTX):} converts to a modern \strong{10x HDF5}
#'   file, then returns a virtual sample wrapping a \code{tenx_h5} handle.
#' }
#'
#' Use \code{output = "10x_hdf5"} to always produce an owned 10x HDF5 file
#' (copy for H5 input; MTX -> 10x H5 convert for directories).
#'
#' @param input Path to a 10x input source (directory or \code{.h5}/\code{.hdf5}).
#'   A Visium HD \code{outs} folder is sent to \code{binned_outputs/square_*um}.
#' @param bin.size Visium HD bin in microns: 8 (default), 16, or 2.
#'   Classic Visium directories are unchanged.
#' @param input.format One of \code{"auto"}, \code{"10x_dir"}, \code{"10x_h5"},
#'   or \code{"icellr2_h5"} (legacy iCellR2 CSC reload). Default \code{"auto"}.
#' @param feature.column Integer column in the feature file for directory input.
#'   Default \code{2}.
#' @param feature.field Field for HDF5 feature names: \code{"name"}, \code{"id"},
#'   or \code{"feature_type"}. Default \code{"name"}.
#' @param output One of \code{"virtual"}, \code{"10x_hdf5"}, \code{"sparse"},
#'   \code{"matrix"}, or \code{"data.frame"}. Default \code{"virtual"}.
#' @param h5.file Output path. If \code{NULL}:
#'   \itemize{
#'     \item \code{virtual} + 10x H5  -> 
#'       \code{iCellR2_import/<sample>_tenx_meta.json}, or
#'       \code{iCellR2_import/<sample>_<feature.type>_tenx_meta.json}
#'       when \code{feature.type} is set, so Gene Expression and Peaks
#'       from the same file do not overwrite each other
#'     \item otherwise  -> 
#'       \code{iCellR2_import/<sample>_tenx.h5}
#'   }
#' @param dataset.name Unused for 10x outputs (kept for legacy CSC reload).
#' @param overwrite Overwrite existing outputs. Default \code{FALSE}.
#' @param sample.id \strong{Required} sample id. Used for default import
#'   paths and as the SODA / aggregation sample name (barcode prefix).
#'   Choose a stable, unique id (e.g. \code{"my1.3m"}, \code{"donorA"}).
#' @param modality Modality label. Default \code{"rna"}.
#' @param feature.type Optional feature-type filter. Use a label from the
#'   file (\code{"Gene Expression"}, \code{"Peaks"}, \code{"Antibody Capture"},
#'   \code{"Multiplexing Capture"}, \code{"CRISPR Guide Capture"}, \code{"Custom"},
#'   or any other type that file contains). Short names \code{"rna"},
#'   \code{"atac"}, \code{"adt"}, and \code{"hto"} resolve to the matching
#'   type when it is present.
#' @param chunk.cells Reserved (chunk size for legacy helpers).
#' @param stream.mtx Ignored (kept for compatibility).
#' @param materialize Reserved; errors if \code{TRUE}.
#' @param verbose Progress style.
#'
#' @return Depending on \code{output}: sparse/dense matrix, data frame,
#'   \code{tenx_h5} assay handle (\code{"10x_hdf5"}), or \code{iCellR2_sample}
#'   (\code{"virtual"}).
#'
#' @importFrom Matrix readMM
#' @importFrom methods as new
#' @export
load10x <- function(
    input,
    input.format = c("auto", "10x_dir", "10x_h5", "icellr2_h5"),
    feature.column = 2,
    feature.field = c("name", "id", "feature_type"),
    output = c("virtual", "10x_hdf5", "sparse", "matrix", "data.frame"),
    h5.file = NULL,
    dataset.name = "counts",
    overwrite = FALSE,
    sample.id,
    modality = "rna",
    feature.type = NULL,
    chunk.cells = 50000L,
    stream.mtx = FALSE,
    materialize = FALSE,
    verbose = TRUE,
    bin.size = 8
) {
  input.format <- match.arg(input.format)
  feature.field <- match.arg(feature.field)
  output <- match.arg(output)
  rep <- .ic2_new_reporter("load10x", verbose)

  if (isTRUE(materialize)) {
    stop(
      "'materialize = TRUE' is reserved. ",
      "Use output = \"10x_hdf5\" for an owned 10x HDF5 copy, ",
      "or output = \"virtual\" for a native tenx_h5 handle."
    )
  }
  if (isTRUE(stream.mtx)) {
    warning("'stream.mtx' is ignored; MTX -> 10x HDF5 uses the Rcpp path.", call. = FALSE)
  }

  if (missing(input) || is.null(input) || !is.character(input) ||
      length(input) != 1L || !nzchar(input)) {
    stop("'input' must be a non-empty character string.")
  }

  if (dir.exists(input)) {
    alt <- .ic2_spatial_resolve_matrix(input, bin.size)
    if (!identical(normalizePath(alt, winslash = "/", mustWork = FALSE),
                   normalizePath(input, winslash = "/", mustWork = FALSE))) {
      .ic2_msg(rep, "Visium HD ", bin.size, " um bins: ", alt)
      input <- alt
    }
  }

  is_dir <- dir.exists(input)
  is_h5_file <- file.exists(input) && grepl("\\.(h5|hdf5)$", input, ignore.case = TRUE)

  if (input.format == "auto") {
    if (is_dir) {
      input.format <- "10x_dir"
    } else if (is_h5_file) {
      kind <- .ic2_detect_h5_kind(input)
      if (identical(kind, "icellr2_csc")) {
        input.format <- "icellr2_h5"
      } else if (identical(kind, "tenx_h5")) {
        input.format <- "10x_h5"
      } else {
        stop(
          "Could not determine HDF5 type for '", input, "'. ",
          "Expected 10x Genomics HDF5 or legacy iCellR2 CSC HDF5."
        )
      }
    } else {
      stop(
        "Could not determine input type. 'input' must be a 10x directory ",
        "or a .h5/.hdf5 file."
      )
    }
  }

  if (input.format == "10x_dir" && !is_dir) {
    stop("'input.format' = '10x_dir' but 'input' is not an existing directory.")
  }
  if (input.format %in% c("10x_h5", "icellr2_h5") && !is_h5_file) {
    stop("'input.format' = '", input.format, "' but 'input' is not an existing .h5/.hdf5 file.")
  }
  if (identical(input.format, "10x_h5") && is_h5_file) {
    kind <- .ic2_detect_h5_kind(input)
    if (identical(kind, "icellr2_csc")) input.format <- "icellr2_h5"
  }

  if (!is.numeric(feature.column) || length(feature.column) != 1L || is.na(feature.column)) {
    stop("'feature.column' must be a single numeric value.")
  }
  feature.column <- as.integer(feature.column)
  if (feature.column < 1L) stop("'feature.column' must be >= 1.")

  if (!is.null(feature.type)) {
    if (!is.character(feature.type) || !length(feature.type) || anyNA(feature.type)) {
      stop("'feature.type' must be NULL or a non-empty character vector.")
    }
  }

  if (!is.character(modality) || length(modality) != 1L || !nzchar(modality)) {
    stop("'modality' must be a non-empty character string.")
  }
  modality <- tolower(modality)

  if (missing(sample.id) || is.null(sample.id)) {
    stop(
      "'sample.id' is required. Choose a unique id for this sample ",
      "(e.g. sample.id = \"donorA\"). It names import files and SODA handles."
    )
  }
  if (!is.character(sample.id) || length(sample.id) != 1L || !nzchar(sample.id)) {
    stop("'sample.id' must be a non-empty character string.")
  }
  sample.id <- .ic2_clean_names(sample.id)[1L]
  if (!nzchar(sample.id)) {
    stop("'sample.id' is empty after cleaning; use letters/numbers (e.g. \"sampleA\").")
  }

  if (output %in% c("virtual", "10x_hdf5") ||
      identical(input.format, "10x_h5") ||
      identical(input.format, "icellr2_h5")) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
      stop("Package 'hdf5r' is required for HDF5 input/output.")
    }
  }

  ## Default output paths
  meta_path <- NULL
  tenx_path <- NULL
  if (identical(input.format, "icellr2_h5") && identical(output, "virtual")) {
    h5.file <- NULL
  } else if (output %in% c("virtual", "10x_hdf5")) {
    if (is.null(h5.file)) {
      if (identical(input.format, "10x_h5") && identical(output, "virtual")) {
        tag <- sample.id
        if (!is.null(feature.type)) {
          ft <- .ic2_clean_names(paste(feature.type, collapse = "_"))[1L]
          if (nzchar(ft)) tag <- paste0(sample.id, "_", ft)
        }
        h5.file <- file.path("iCellR2_import", paste0(tag, "_tenx_meta.json"))
      } else {
        h5.file <- file.path("iCellR2_import", paste0(sample.id, "_tenx.h5"))
      }
    }
    if (!is.character(h5.file) || length(h5.file) != 1L || !nzchar(h5.file)) {
      stop("'h5.file' must be a non-empty character string.")
    }
    h5.file <- path.expand(h5.file)
    if (grepl("\\.(rds|json)$", h5.file, ignore.case = TRUE)) {
      meta_path <- sub("\\.rds$", ".json", h5.file, ignore.case = TRUE)
      tenx_path <- paste0(tools::file_path_sans_ext(h5.file), ".h5")
      if (identical(input.format, "10x_h5") && identical(output, "virtual")) {
        tenx_path <- NULL ## pointer  -  no owned H5
      }
    } else {
      tenx_path <- h5.file
      meta_path <- paste0(tools::file_path_sans_ext(h5.file), "_tenx_meta.json")
    }
    check_exists <- unique(c(meta_path, tenx_path))
    check_exists <- check_exists[!vapply(check_exists, is.null, logical(1))]
    for (p in check_exists) {
      if (file.exists(p) && !overwrite) {
        stop("File already exists at '", p, "'. Use overwrite = TRUE.")
      }
    }
  }

  out_for_log <- if (!is.null(tenx_path)) tenx_path else if (!is.null(meta_path)) meta_path else input
  rep <- .ic2_set_io_dirs(rep, input_dir = input, output_dir = out_for_log)
  rep <- .ic2_add_input(rep, "Input", normalizePath(input, winslash = "/", mustWork = FALSE))

  make_sample <- function(handle, source_path) {
    sample_desc <- list(
      type = "iCellR2_sample",
      modality = modality,
      sample_id = sample.id,
      source = normalizePath(source_path, winslash = "/", mustWork = FALSE),
      handle = handle,
      feature_ids = NULL,
      cell_ids = NULL,
      n_features = handle$n_features,
      n_cells = handle$n_cells,
      created = as.character(Sys.time())
    )
    class(sample_desc) <- c("iCellR2_sample", "list")
    sample_desc
  }

  ## ---- Legacy iCellR2 CSC reload ----
  if (identical(input.format, "icellr2_h5")) {
    if (identical(output, "10x_hdf5")) {
      stop(
        "output = \"10x_hdf5\" is not supported for legacy iCellR2 CSC input. ",
        "Use output = \"virtual\" to attach the CSC, or re-import from 10x sources."
      )
    }
    rep <- .ic2_banner_start(rep, c(
      paste0("Source: ", input),
      paste0("Input: icellr2_h5 (legacy CSC)   |   sample: ", sample.id),
      paste0("Output: ", output)
    ))
    group <- dataset.name
    kind_check <- .ic2_detect_h5_kind(input)
    if (!identical(kind_check, "icellr2_csc")) {
      stop("'", input, "' is not an iCellR2 CSC HDF5 file. Detected: ", kind_check, ".")
    }
    h5tmp <- hdf5r::H5File$new(input, mode = "r")
    on.exit(try(h5tmp$close_all(), silent = TRUE), add = TRUE)
    gnames <- names(h5tmp)
    if (!group %in% gnames) {
      group <- if ("counts" %in% gnames) "counts" else gnames[1L]
    }
    dims <- .ic2_hdf5_dims(input, group = group)
    h5tmp$close_all()
    on.exit(NULL)

    if (identical(output, "virtual")) {
      .ic2_step(rep, "[1/1] Attaching legacy iCellR2 CSC HDF5")
      handle <- .ic2_assay_handle(
        path = input,
        backend = "hdf5_sparse",
        group = group,
        n_features = dims$n_features,
        n_cells = dims$n_cells,
        sample_id = sample.id,
        modality = modality
      )
      rep <- .ic2_add_output(rep, "HDF5", normalizePath(input, winslash = "/", mustWork = FALSE))
      .ic2_finish(rep, extra_lines = paste0(
        "Features: ", .ic2_fmt_num(handle$n_features),
        "   |   Cells: ", .ic2_fmt_num(handle$n_cells),
        "   |   legacy CSC"
      ))
      return(make_sample(handle, input))
    }

    .ic2_step(rep, "[1/1] Materializing legacy CSC into memory")
    handle <- .ic2_assay_handle(
      path = input, backend = "hdf5_sparse", group = group,
      n_features = dims$n_features, n_cells = dims$n_cells
    )
    handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
    mtx <- .ic2_read_cols(handle, seq_len(handle$n_cells))
    if (identical(output, "sparse")) { .ic2_finish(rep); return(mtx) }
    if (identical(output, "matrix")) { .ic2_finish(rep); return(as.matrix(mtx)) }
    .ic2_finish(rep)
    return(as.data.frame(as.matrix(mtx), stringsAsFactors = FALSE))
  }

  open_file <- function(x) {
    if (grepl("\\.gz$", x, ignore.case = TRUE)) gzfile(x) else x
  }
  to_num <- function(x) {
    if (inherits(x, "integer64")) return(as.numeric(x))
    x
  }

  read_10x_dir <- function(path, feature.column) {
    files <- .ic2_find_10x_dir_files(path)
    mtx <- Matrix::readMM(open_file(files$matrix))
    mtx <- .ic2_to_dgCMatrix(mtx)
    cell.barcodes <- readLines(open_file(files$barcodes), warn = FALSE)
    gene.info <- utils::read.table(
      file = open_file(files$features), header = FALSE, sep = "\t",
      stringsAsFactors = FALSE, quote = "", comment.char = ""
    )
    if (feature.column > ncol(gene.info)) {
      stop("'feature.column' = ", feature.column, " out of range (", ncol(gene.info), " columns).")
    }
    list(mtx = mtx, barcodes = cell.barcodes, features = as.character(gene.info[[feature.column]]))
  }

  read_10x_h5_memory <- function(path, feature.field, feature.type = NULL) {
    h5 <- hdf5r::H5File$new(path, mode = "r")
    on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
    if (h5$exists("matrix")) {
      g <- h5[["matrix"]]
      needed <- c("data", "indices", "indptr", "shape", "barcodes", "features")
      missing_needed <- needed[!vapply(needed, g$exists, logical(1))]
      if (length(missing_needed)) {
        stop("Invalid modern 10x HDF5 file: missing ", paste(missing_needed, collapse = ", "))
      }
      data <- g[["data"]]$read()
      indices <- to_num(g[["indices"]]$read())
      indptr <- to_num(g[["indptr"]]$read())
      shape <- to_num(g[["shape"]]$read())
      barcodes <- as.character(g[["barcodes"]]$read())
      fg <- g[["features"]]
      if (!fg$exists(feature.field)) {
        stop("Requested feature field '", feature.field, "' not found in HDF5 file.")
      }
      features <- as.character(fg[[feature.field]]$read())
      feature_type_all <- if (fg$exists("feature_type")) {
        as.character(fg[["feature_type"]]$read())
      } else NULL
    } else {
      top <- names(h5)
      if (!length(top)) stop("HDF5 file appears empty.")
      g <- h5[[top[1]]]
      needed <- c("data", "indices", "indptr", "shape", "barcodes")
      missing_needed <- needed[!vapply(needed, g$exists, logical(1))]
      if (length(missing_needed)) {
        stop("Invalid legacy 10x HDF5 file: missing ", paste(missing_needed, collapse = ", "))
      }
      data <- g[["data"]]$read()
      indices <- to_num(g[["indices"]]$read())
      indptr <- to_num(g[["indptr"]]$read())
      shape <- to_num(g[["shape"]]$read())
      barcodes <- as.character(g[["barcodes"]]$read())
      legacy_map <- c(id = "genes", name = "gene_names")
      selected_field <- legacy_map[[feature.field]]
      if (is.null(selected_field)) {
        stop("Unsupported feature field for legacy 10x HDF5 input: ", feature.field)
      }
      if (!g$exists(selected_field)) {
        stop("Requested feature field '", selected_field, "' not found in legacy HDF5 file.")
      }
      features <- as.character(g[[selected_field]]$read())
      feature_type_all <- rep("Gene Expression", length(features))
    }

    nr <- as.integer(shape[1]); nc <- as.integer(shape[2])
    if (length(indices) > .Machine$integer.max || max(indptr, na.rm = TRUE) > .Machine$integer.max) {
      stop(
        "Matrix too large for in-memory dgCMatrix. ",
        "Use output = \"virtual\" or output = \"10x_hdf5\"."
      )
    }
    if (!is.null(feature.type)) {
      matched <- .ic2_match_feature_types(feature_type_all, feature.type)
      keep <- matched$keep
      mtx <- methods::new("dgCMatrix", i = as.integer(indices), p = as.integer(indptr),
                          x = as.numeric(data), Dim = c(nr, nc))
      mtx <- mtx[keep, , drop = FALSE]
      features <- features[keep]
      return(list(mtx = .ic2_to_dgCMatrix(mtx), barcodes = barcodes, features = features))
    }
    mtx <- methods::new("dgCMatrix", i = as.integer(indices), p = as.integer(indptr),
                        x = as.numeric(data), Dim = c(nr, nc))
    list(mtx = mtx, barcodes = barcodes, features = features)
  }

  rep <- .ic2_banner_start(rep, c(
    paste0("Source: ", input),
    paste0("Input: ", input.format),
    paste0("Output: ", output, "   |   sample: ", sample.id, "   |   modality: ", modality),
    if (!is.null(tenx_path)) paste0("Out: ", tenx_path)
    else if (!is.null(meta_path)) paste0("Out: ", meta_path)
    else NULL
  ))

  ## ---- 10x H5 + virtual: native pointer ----
  if (identical(input.format, "10x_h5") && identical(output, "virtual")) {
    .ic2_step(rep, "[1/2] Inspecting 10x HDF5 (native handle, no convert)")
    handle <- .ic2_tenx_h5_handle(
      src_path = input,
      meta_path = meta_path,
      feature.field = feature.field,
      feature.type = feature.type,
      sample_id = sample.id,
      modality = modality,
      overwrite = overwrite
    )
    rep <- .ic2_progress_bar(rep, 1L, 2L, prefix = "  inspect")
    .ic2_step(rep, "[2/2] Building sample descriptor")
    rep <- .ic2_progress_bar(rep, 2L, 2L, prefix = "  inspect")
    rep <- .ic2_add_output(rep, "Meta", normalizePath(meta_path, winslash = "/", mustWork = FALSE))
    .ic2_finish(rep, extra_lines = paste0(
      "Features: ", .ic2_fmt_num(handle$n_features),
      "   |   Cells: ", .ic2_fmt_num(handle$n_cells),
      "   |   native tenx_h5 (no copy)"
    ))
    return(make_sample(handle, input))
  }

  ## ---- 10x H5 + 10x_hdf5: owned copy ----
  if (identical(input.format, "10x_h5") && identical(output, "10x_hdf5")) {
    .ic2_step(rep, "[1/2] Copying 10x HDF5 into owned file")
    dir.create(dirname(tenx_path), recursive = TRUE, showWarnings = FALSE)
    if (!identical(normalizePath(input, mustWork = FALSE),
                   normalizePath(tenx_path, mustWork = FALSE))) {
      if (file.exists(tenx_path)) unlink(tenx_path)
      if (!file.copy(input, tenx_path, overwrite = TRUE)) {
        stop("Failed to copy '", input, "'  ->  '", tenx_path, "'.")
      }
    }
    .ic2_step(rep, "[2/2] Building tenx_h5 handle")
    handle <- .ic2_tenx_h5_handle(
      src_path = tenx_path,
      meta_path = meta_path,
      feature.field = feature.field,
      feature.type = feature.type,
      sample_id = sample.id,
      modality = modality,
      overwrite = TRUE
    )
    rep <- .ic2_add_output(rep, "10x H5", normalizePath(tenx_path, winslash = "/", mustWork = FALSE))
    .ic2_finish(rep, extra_lines = paste0(
      "Features: ", .ic2_fmt_num(handle$n_features),
      "   |   Cells: ", .ic2_fmt_num(handle$n_cells),
      "   |   owned 10x HDF5 copy"
    ))
    return(handle)
  }

  ## ---- MTX dir  ->  10x HDF5 ----
  if (identical(input.format, "10x_dir") && output %in% c("virtual", "10x_hdf5")) {
    pfx <- new.env(parent = emptyenv())
    pfx$v <- "  mtx parse"
    tenx_out <- .ic2_convert_mtx_dir_to_tenx_hdf5(
      dir_path = input,
      dest_path = tenx_path,
      feature.column = feature.column,
      overwrite = overwrite,
      progress = function(cur, tot) {
        rep <<- .ic2_progress_bar(rep, cur, tot, prefix = pfx$v)
      },
      on_phase = function(phase) {
        if (identical(phase, "parse")) {
          rep <<- .ic2_step(rep, "[1/3] Rcpp MTX parse (in memory)")
          pfx$v <- "  mtx parse"
        } else if (identical(phase, "write")) {
          rep <<- .ic2_step(rep, "[2/3] Writing 10x HDF5")
          pfx$v <- "  tenx write"
        }
      }
    )
    .ic2_step(rep, "[3/3] Building tenx_h5 handle")
    handle <- .ic2_tenx_h5_handle(
      src_path = tenx_out,
      meta_path = meta_path,
      feature.field = "name",
      feature.type = feature.type,
      sample_id = sample.id,
      modality = modality,
      overwrite = TRUE
    )
    rep <- .ic2_set_io_dirs(rep, input_dir = input, output_dir = dirname(tenx_out))
    rep <- .ic2_add_output(rep, "10x H5", normalizePath(tenx_out, winslash = "/", mustWork = FALSE))
    .ic2_finish(rep, extra_lines = paste0(
      "Features: ", .ic2_fmt_num(handle$n_features),
      "   |   Cells: ", .ic2_fmt_num(handle$n_cells),
      "   |   MTX -> 10x HDF5"
    ))
    if (identical(output, "10x_hdf5")) return(handle)
    return(make_sample(handle, input))
  }

  ## ---- In-memory outputs ----
  .ic2_step(rep, paste0(
    "[1/1] Reading 10x input into memory",
    if (identical(input.format, "10x_dir")) " (Matrix::readMM)" else ""
  ))
  dat <- switch(
    input.format,
    "10x_dir" = read_10x_dir(input, feature.column),
    "10x_h5"  = read_10x_h5_memory(input, feature.field, feature.type)
  )
  mtx <- dat$mtx
  cell.barcodes <- .ic2_clean_names(dat$barcodes)
  gene.names.ids <- .ic2_clean_names(dat$features)
  colnames(mtx) <- cell.barcodes
  rownames(mtx) <- gene.names.ids

  if (output == "sparse") {
    .ic2_finish(rep, extra_lines = paste0(
      "Returned sparse matrix: ", .ic2_fmt_num(nrow(mtx)), "  x  ", .ic2_fmt_num(ncol(mtx))
    ))
    return(mtx)
  }
  if (output == "matrix") {
    .ic2_finish(rep)
    return(as.matrix(mtx))
  }
  if (output == "data.frame") {
    .ic2_finish(rep)
    return(as.data.frame(as.matrix(mtx), stringsAsFactors = FALSE))
  }

  stop("Unhandled output = '", output, "'.")
}
