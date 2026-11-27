#' Internal: CRAN-friendly sparse HDF5 storage and assay handles
#'
#' Default on-disk backend for iCellR2. Matrices are stored as CSC
#' (`i`, `p`, `x`) via \pkg{hdf5r}.
#'
#' @keywords internal
#' @name iCellR2-storage
NULL

.ic2_clean_names <- function(x) {
  x <- gsub("-", ".", x, fixed = TRUE)
  make.names(x, unique = TRUE)
}

#' Default sample.id from an input path (avoids multi-sample collisions)
#' @keywords internal
.ic2_default_sample_id <- function(input) {
  input <- as.character(input)[1L]
  if (!nzchar(input)) return("sample1")
  bn <- basename(path.expand(input))
  bn <- sub("\\.(h5|hdf5)$", "", bn, ignore.case = TRUE)
  bn <- sub("_filtered_feature_bc_matrix$", "", bn, ignore.case = TRUE)
  bn <- sub("\\.mtx(\\.gz)?$", "", bn, ignore.case = TRUE)
  if (!nzchar(bn) || bn %in% c(".", "..")) bn <- "sample1"
  ## If input was the standard 10x folder name, use parent directory name
  if (identical(tolower(bn), "filtered_feature_bc_matrix") ||
      identical(tolower(bn), "raw_feature_bc_matrix")) {
    parent <- basename(dirname(path.expand(input)))
    if (nzchar(parent) && !parent %in% c(".", "..")) bn <- parent
  }
  sid <- .ic2_clean_names(bn)[1L]
  if (!nzchar(sid)) sid <- "sample1"
  sid
}

.ic2_to_dgCMatrix <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.matrix(x)) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  if (inherits(x, "Matrix") || !is.null(dim(x))) {
    ## Matrix >= 1.4-2: prefer CsparseMatrix over deprecated as(., "dgCMatrix")
    x <- methods::as(x, "CsparseMatrix")
    if (!inherits(x, "dgCMatrix")) {
      x <- methods::as(x, "generalMatrix")
      x <- methods::as(x, "CsparseMatrix")
    }
    return(x)
  }
  stop("Unsupported matrix type for sparse conversion.")
}

.ic2_is_assay_handle <- function(x) {
  inherits(x, "iCellR2_assay_handle") ||
    (is.list(x) && isTRUE(x$type == "iCellR2_assay_handle"))
}

.ic2_backend_available <- function(backend = c("hdf5", "hdf5array")) {
  backend <- match.arg(backend)
  if (backend == "hdf5") {
    return(requireNamespace("hdf5r", quietly = TRUE))
  }
  if (backend == "hdf5array") {
    return(requireNamespace("HDF5Array", quietly = TRUE))
  }
  FALSE
}

#' Create a lightweight assay handle (no matrix payload)
#'
#' @param path Path to on-disk assay file.
#' @param backend Storage backend: `"hdf5_sparse"` or `"rds"`.
#' @param group HDF5 group name (hdf5_sparse only).
#' @param n_features,n_cells Dimensions.
#' @param feature_ids,cell_ids Optional dimnames (prefer on-disk for huge data).
#' @param sample_id Optional sample id.
#' @param modality Optional modality.
#' @keywords internal
.ic2_assay_handle <- function(
    path,
    backend = "hdf5_sparse",
    group = "counts",
    n_features = NA_integer_,
    n_cells = NA_integer_,
    feature_ids = NULL,
    cell_ids = NULL,
    sample_id = NULL,
    modality = NULL
) {
  h <- list(
    type = "iCellR2_assay_handle",
    backend = backend,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    group = group,
    n_features = as.integer(n_features),
    n_cells = as.integer(n_cells),
    feature_ids = feature_ids,
    cell_ids = cell_ids,
    sample_id = sample_id,
    modality = modality,
    created = as.character(Sys.time())
  )
  class(h) <- c("iCellR2_assay_handle", "list")
  h
}

#' Match a requested feature.type against the types in this file.
#'
#' Exact labels are compared without case or punctuation. Short names such as
#' \code{rna}, \code{atac}, \code{adt}, and \code{hto} are used only when the
#' corresponding type is present. Anything else must match a type in the file.
#' @keywords internal
.ic2_feature_type_key <- function(x) {
  gsub("[^a-z0-9]", "", tolower(as.character(x)))
}

#' @keywords internal
.ic2_resolve_feature_types <- function(requested, available) {
  requested <- as.character(requested)
  requested <- requested[nzchar(requested) & !is.na(requested)]
  available <- unique(as.character(available))
  available <- available[nzchar(available) & !is.na(available)]
  if (!length(requested)) stop("'feature.type' must contain at least one name.")
  if (!length(available)) {
    stop("'feature.type' was supplied, but this input has no feature types.")
  }
  aliases <- list(
    rna = "Gene Expression", gene = "Gene Expression", genes = "Gene Expression",
    gex = "Gene Expression", scrna = "Gene Expression",
    atac = "Peaks", peak = "Peaks", peaks = "Peaks", scatac = "Peaks",
    adt = "Antibody Capture", cite = "Antibody Capture", citeseq = "Antibody Capture",
    antibody = "Antibody Capture",
    hto = c("Multiplexing Capture", "Antibody Capture"),
    hashtag = c("Multiplexing Capture", "Antibody Capture"),
    hashtags = c("Multiplexing Capture", "Antibody Capture"),
    hashing = c("Multiplexing Capture", "Antibody Capture"),
    cmo = "Multiplexing Capture", cellplex = "Multiplexing Capture",
    crispr = "CRISPR Guide Capture", guide = "CRISPR Guide Capture",
    custom = "Custom", antigen = "Antigen Capture"
  )
  avail_key <- .ic2_feature_type_key(available)
  out <- character(0)
  for (req in requested) {
    key <- .ic2_feature_type_key(req)
    exact <- available[avail_key == key]
    if (length(exact)) {
      out <- c(out, exact)
      next
    }
    pref <- aliases[[key]]
    if (!is.null(pref)) {
      pref_key <- .ic2_feature_type_key(pref)
      hit <- pref_key %in% avail_key
      if (any(hit)) {
        out <- c(out, available[match(pref_key[hit][1L], avail_key)])
        next
      }
    }
    sub <- available[nzchar(key) & vapply(avail_key, function(k) grepl(key, k, fixed = TRUE), logical(1))]
    if (length(sub) == 1L) {
      out <- c(out, sub)
      next
    }
    stop(
      "No features matched feature.type '", req, "'. This input has: ",
      paste(available, collapse = ", "), "."
    )
  }
  unique(out)
}

#' @keywords internal
.ic2_match_feature_types <- function(feature_type_all, feature.type) {
  if (is.null(feature.type)) {
    return(list(keep = rep(TRUE, length(feature_type_all)), feature.type = NULL))
  }
  if (is.null(feature_type_all)) {
    stop("'feature.type' was supplied, but no feature_type information is available.")
  }
  resolved <- .ic2_resolve_feature_types(feature.type, feature_type_all)
  keep <- feature_type_all %in% resolved
  if (!any(keep)) {
    stop(
      "No features matched 'feature.type'. This input has: ",
      paste(unique(feature_type_all), collapse = ", "), "."
    )
  }
  list(keep = keep, feature.type = resolved)
}

#' Inspect a 10x HDF5 file (modern or legacy) without loading the matrix
#' @keywords internal
.ic2_tenx_h5_inspect <- function(path, feature.field = "name", feature.type = NULL) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  to_num <- function(x) {
    if (inherits(x, "integer64")) return(as.numeric(x))
    as.numeric(x)
  }
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  if (h5$exists("matrix")) {
    layout <- "modern"
    group <- "matrix"
    g <- h5[["matrix"]]
    needed <- c("data", "indices", "indptr", "shape", "barcodes", "features")
    missing_needed <- needed[!vapply(needed, g$exists, logical(1))]
    if (length(missing_needed)) {
      stop("Invalid modern 10x HDF5 file: missing ", paste(missing_needed, collapse = ", "))
    }
    shape <- to_num(g[["shape"]]$read())
    barcodes <- as.character(g[["barcodes"]]$read())
    fg <- g[["features"]]
    if (!fg$exists(feature.field)) {
      stop("Requested feature field '", feature.field, "' not found in HDF5 file.")
    }
    features <- as.character(fg[[feature.field]]$read())
    feature_type_all <- if (fg$exists("feature_type")) {
      as.character(fg[["feature_type"]]$read())
    } else {
      NULL
    }
  } else {
    layout <- "legacy"
    top <- names(h5)
    if (!length(top)) stop("HDF5 file appears empty.")
    group <- top[1]
    g <- h5[[group]]
    needed <- c("data", "indices", "indptr", "shape", "barcodes")
    missing_needed <- needed[!vapply(needed, g$exists, logical(1))]
    if (length(missing_needed)) {
      stop("Invalid legacy 10x HDF5 file: missing ", paste(missing_needed, collapse = ", "))
    }
    shape <- to_num(g[["shape"]]$read())
    barcodes <- as.character(g[["barcodes"]]$read())
    legacy_map <- c(id = "genes", name = "gene_names")
    selected_field <- legacy_map[[feature.field]]
    if (is.null(selected_field)) {
      if (identical(feature.field, "feature_type")) {
        stop("'feature.field = \"feature_type\"' is not available for legacy 10x HDF5 input.")
      }
      stop("Unsupported feature field for legacy 10x HDF5 input: ", feature.field)
    }
    if (!g$exists(selected_field)) {
      stop("Requested feature field '", selected_field, "' not found in legacy HDF5 file.")
    }
    features <- as.character(g[[selected_field]]$read())
    feature_type_all <- rep("Gene Expression", length(features))
  }

  nr0 <- as.integer(shape[1])
  nc <- as.integer(shape[2])
  if (length(features) != nr0) stop("Feature count does not match matrix rows.")
  if (length(barcodes) != nc) stop("Barcode count does not match matrix columns.")

  row_keep <- NULL
  if (!is.null(feature.type)) {
    matched <- .ic2_match_feature_types(feature_type_all, feature.type)
    feature.type <- matched$feature.type
    row_keep <- which(matched$keep)
    features <- features[row_keep]
  }

  list(
    layout = layout,
    group = group,
    n_features_src = nr0,
    n_features = length(features),
    n_cells = nc,
    features = .ic2_clean_names(features),
    barcodes = .ic2_clean_names(barcodes),
    row_keep = row_keep,
    feature.field = feature.field,
    feature.type = feature.type
  )
}

#' Build a native tenx_h5 assay handle (no matrix copy)
#' @keywords internal
.ic2_tenx_h5_handle <- function(
    src_path,
    meta_path,
    feature.field = "name",
    feature.type = NULL,
    sample_id = NULL,
    modality = "rna",
    overwrite = FALSE
) {
  info <- .ic2_tenx_h5_inspect(src_path, feature.field = feature.field, feature.type = feature.type)
  src_path <- normalizePath(src_path, winslash = "/", mustWork = TRUE)
  if (file.exists(meta_path) && !isTRUE(overwrite)) {
    stop("Meta file already exists at '", meta_path, "'. Use overwrite = TRUE.")
  }
  dir.create(dirname(meta_path), recursive = TRUE, showWarnings = FALSE)

  meta <- list(
    type = "iCellR2_tenx_h5_meta",
    source = src_path,
    layout = info$layout,
    group = info$group,
    feature.field = feature.field,
    feature.type = feature.type,
    row_keep = info$row_keep,
    feature_ids = info$features,
    cell_ids = info$barcodes,
    n_features = info$n_features,
    n_features_src = info$n_features_src,
    n_cells = info$n_cells,
    created = as.character(Sys.time())
  )
  meta_path <- .ic2_soda_write_tenx_meta(meta, meta_path)

  h <- .ic2_assay_handle(
    path = src_path,
    backend = "tenx_h5",
    group = info$group,
    n_features = info$n_features,
    n_cells = info$n_cells,
    feature_ids = NULL,
    cell_ids = NULL,
    sample_id = sample_id,
    modality = modality
  )
  h$meta_path <- normalizePath(meta_path, winslash = "/", mustWork = FALSE)
  h
}

#' Write a dgCMatrix as modern 10x HDF5 + tenx_h5 handle
#'
#' Preferred on-disk format for raw counts (agg temp + SODA).
#' @keywords internal
.ic2_write_matrix_tenx_h5 <- function(
    mat,
    filepath,
    sample_id = NULL,
    modality = "rna",
    feature_types = "Gene Expression",
    genome = "unknown",
    overwrite = TRUE
) {
  mat <- .ic2_to_dgCMatrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("'mat' must have rownames and colnames.")
  }
  if (!is.character(filepath) || length(filepath) != 1L || !nzchar(filepath)) {
    stop("'filepath' must be a non-empty character string.")
  }
  filepath <- path.expand(filepath)
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  stem <- tools::file_path_sans_ext(basename(filepath))
  meta_path <- file.path(dirname(filepath), paste0(stem, "_tenx_meta.json"))

  .ic2_write_tenx_hdf5(
    filepath = filepath,
    i = as.integer(mat@i),
    p = as.integer(mat@p),
    x = as.numeric(mat@x),
    barcodes = as.character(colnames(mat)),
    feature_ids = as.character(rownames(mat)),
    feature_names = as.character(rownames(mat)),
    feature_types = feature_types,
    genome = genome,
    overwrite = overwrite
  )
  .ic2_tenx_h5_handle(
    src_path = filepath,
    meta_path = meta_path,
    feature.field = "name",
    sample_id = sample_id,
    modality = modality,
    overwrite = TRUE
  )
}

#' Convert an on-disk CSC (hdf5_sparse) handle to tenx_h5 at dest_path
#' @keywords internal
.ic2_convert_csc_handle_to_tenx <- function(handle, dest_path, overwrite = TRUE) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  if (!identical(handle$backend, "hdf5_sparse")) {
    stop("Expected hdf5_sparse handle for CSC -> tenx conversion.")
  }
  nc <- as.integer(handle$n_cells)
  if (!is.finite(nc) || nc < 1L) {
    stop("CSC handle has invalid n_cells.")
  }
  mat <- read_sparse_hdf5_cols(
    handle$path,
    group = handle$group %||% "counts",
    cols = seq_len(nc)
  )
  if (is.null(rownames(mat)) && !is.null(handle$feature_ids)) {
    rownames(mat) <- as.character(handle$feature_ids)
  }
  if (is.null(colnames(mat)) && !is.null(handle$cell_ids)) {
    colnames(mat) <- as.character(handle$cell_ids)
  }
  .ic2_write_matrix_tenx_h5(
    mat,
    filepath = dest_path,
    sample_id = handle$sample_id,
    modality = handle$modality %||% "rna",
    overwrite = overwrite
  )
}

#' Write a modern 10x Genomics HDF5 file (CSC under /matrix)
#' @keywords internal
.ic2_write_tenx_hdf5 <- function(
    filepath,
    i,
    p,
    x,
    barcodes,
    feature_ids,
    feature_names = NULL,
    feature_types = NULL,
    genome = "unknown",
    overwrite = FALSE,
    progress = NULL
) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  nr <- length(feature_ids)
  nc <- length(barcodes)
  nnz <- length(i)
  if (length(x) != nnz) stop("length(i) != length(x).")
  if (length(p) != nc + 1L) stop("length(p) != ncol + 1.")
  if (is.null(feature_names)) feature_names <- feature_ids
  if (length(feature_names) != nr) stop("feature_names length != n_features.")
  if (is.null(feature_types)) feature_types <- rep("Gene Expression", nr)
  if (length(feature_types) == 1L) feature_types <- rep(as.character(feature_types), nr)
  if (length(feature_types) != nr) stop("feature_types length != n_features.")
  if (length(genome) == 1L) genome <- rep(as.character(genome), nr)
  if (length(genome) != nr) stop("genome length != n_features.")

  if (file.exists(filepath)) {
    if (!isTRUE(overwrite)) {
      stop("File already exists at '", filepath, "'. Use overwrite = TRUE.")
    }
    unlink(filepath)
  }
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)

  h5 <- hdf5r::H5File$new(filepath, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5$create_group("matrix")
  fg <- g$create_group("features")

  ## Cell Ranger-style datasets
  g[["barcodes"]] <- as.character(barcodes)
  g[["shape"]] <- as.integer(c(nr, nc))
  ## indptr as integer64-friendly double stored; hdf5r writes numeric as float64.
  ## 10x often uses int64; float64 is fine for nnz < 2^53.
  chunk <- as.integer(min(max(nnz, 1), 5000000L))
  ds_data <- g$create_dataset(
    "data",
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = max(nnz, 1),
    chunk_dims = chunk
  )
  ds_idx <- g$create_dataset(
    "indices",
    dtype = hdf5r::h5types$H5T_NATIVE_INT,
    dims = max(nnz, 1),
    chunk_dims = chunk
  )
  if (nnz > 0) {
    if (is.function(progress)) progress(0, nnz)
    pos <- 0
    while (pos < nnz) {
      take <- min(chunk, nnz - pos)
      idx <- (pos + 1):(pos + take)
      .ic2_h5_write_1d(ds_data, start1 = pos + 1, values = as.numeric(x[idx]))
      .ic2_h5_write_1d(ds_idx, start1 = pos + 1, values = as.integer(i[idx]))
      pos <- pos + take
      if (is.function(progress)) progress(pos, nnz)
    }
  }
  g[["indptr"]] <- as.numeric(p)

  fg[["id"]] <- as.character(feature_ids)
  fg[["name"]] <- as.character(feature_names)
  fg[["feature_type"]] <- as.character(feature_types)
  fg[["genome"]] <- as.character(genome)

  h5$flush()
  invisible(normalizePath(filepath, winslash = "/", mustWork = TRUE))
}

#' Convert 10x MTX directory  ->  modern 10x HDF5 via Rcpp
#' @keywords internal
.ic2_convert_mtx_dir_to_tenx_hdf5 <- function(
    dir_path,
    dest_path,
    feature.column = 2L,
    overwrite = FALSE,
    progress = NULL,
    on_phase = NULL,
    force = c("auto", "memory", "stream")
) {
  force <- match.arg(force)
  if (!exists("cpp_mtx_to_csc", mode = "function") ||
      !exists("cpp_mtx_info", mode = "function")) {
    stop("Rcpp MTX parser is not available; rebuild the package.")
  }
  files <- .ic2_find_10x_dir_files(dir_path)
  feature.column <- as.integer(feature.column)[1L]

  bc_con <- .ic2_open_text(files$barcodes)
  on.exit(try(close(bc_con), silent = TRUE), add = TRUE)
  barcodes <- readLines(bc_con, warn = FALSE)
  close(bc_con)

  ft_con <- .ic2_open_text(files$features)
  on.exit(try(close(ft_con), silent = TRUE), add = TRUE)
  gene.info <- utils::read.table(
    file = ft_con,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  close(ft_con)
  on.exit(NULL)

  if (ncol(gene.info) < 1L) stop("Empty features/genes file.")
  feature_ids <- as.character(gene.info[[1L]])
  name_col <- if (ncol(gene.info) >= 2L) min(as.integer(feature.column), ncol(gene.info)) else 1L
  feature_names <- as.character(gene.info[[name_col]])
  feature_types <- if (ncol(gene.info) >= 3L) {
    as.character(gene.info[[3L]])
  } else {
    rep("Gene Expression", nrow(gene.info))
  }
  barcodes <- .ic2_clean_names(barcodes)
  feature_ids <- as.character(feature_ids)
  feature_names <- as.character(feature_names)

  info <- cpp_mtx_info(files$matrix)
  nr <- as.integer(info$nrow)
  nc <- as.integer(info$ncol)
  nnz <- as.numeric(info$nnz)
  if (length(feature_ids) != nr) {
    stop("Features (", length(feature_ids), ") != MTX rows (", nr, ").")
  }
  if (length(barcodes) != nc) {
    stop("Barcodes (", length(barcodes), ") != MTX columns (", nc, ").")
  }

  need <- .ic2_mtx_csc_ram_bytes(nnz, nc)
  avail <- .ic2_available_ram_bytes()
  use_memory <- switch(
    force,
    memory = TRUE,
    stream = FALSE,
    auto = is.finite(need) && is.finite(avail) && need < avail * 0.5
  )
  call_phase <- function(phase) if (is.function(on_phase)) on_phase(phase)

  ## Stream path still needs full CSC for tenx write for now  ->  require memory
  ## or fall back to memory when stream requested but unsorted.
  if (!use_memory) {
    ## Prefer memory if possible; else parse memory-style (same as CSC stream
    ## then write). For true low-RAM we'd stream into tenx datasets  -  use memory
    ## parser when estimate fits 85% avail; else error with guidance.
    if (need < avail * 0.85) {
      use_memory <- TRUE
    } else {
      stop(
        "MTX -> 10x HDF5 needs ~", sprintf("%.1f GB", need / 1024^3),
        " RAM to assemble CSC before writing. Free memory or convert on a larger machine."
      )
    }
  }

  call_phase("parse")
  csc <- cpp_mtx_to_csc(files$matrix, progress = progress)
  call_phase("write")
  path <- .ic2_write_tenx_hdf5(
    filepath = dest_path,
    i = csc$i,
    p = csc$p,
    x = csc$x,
    barcodes = barcodes,
    feature_ids = feature_ids,
    feature_names = feature_names,
    feature_types = feature_types,
    overwrite = overwrite,
    progress = progress
  )
  attr(path, "mtx_mode") <- "rcpp-tenx"
  path
}

#' Read column block from a native 10x HDF5 assay
#' @keywords internal
read_tenx_hdf5_cols <- function(filepath, meta, cols) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  cols <- as.integer(cols)
  if (length(cols) < 1L || any(is.na(cols)) || any(cols < 1L)) {
    stop("'cols' must be positive 1-based column indices.")
  }
  if (is.character(meta) && length(meta) == 1L && file.exists(meta)) {
    meta <- .ic2_soda_read_tenx_meta(meta)
  }
  if (!is.list(meta)) stop("Invalid tenx_h5 meta.")

  ncol_m <- as.integer(meta$n_cells)
  if (max(cols) > ncol_m) {
    stop("Requested columns exceed matrix ncol (", ncol_m, ").")
  }

  ## One open per call; contig runs share the same handle
  h5 <- hdf5r::H5File$new(filepath, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5[[meta$group]]

  rn <- as.character(meta$feature_ids)
  nrow_out <- length(rn)
  row_keep <- meta$row_keep
  remap <- NULL
  if (!is.null(row_keep)) {
    remap <- rep(NA_integer_, as.integer(meta$n_features_src))
    remap[as.integer(row_keep)] <- seq_len(length(row_keep)) - 1L
  }

  .read_run <- function(run) {
    indptr <- as.numeric(.ic2_h5_read_1d(g[["indptr"]], start1 = run[1L], count = length(run) + 1L))
    x_start0 <- indptr[1L]
    n_x <- indptr[length(indptr)] - x_start0
    if (n_x < 1) {
      i_vals <- integer(0)
      x_vals <- numeric(0)
    } else {
      i_vals <- as.integer(.ic2_h5_read_1d(g[["indices"]], start1 = x_start0 + 1, count = n_x))
      x_vals <- as.numeric(.ic2_h5_read_1d(g[["data"]], start1 = x_start0 + 1, count = n_x))
    }
    p_local <- as.integer(indptr - indptr[1L])
    cn <- as.character(meta$cell_ids[run])

    if (!is.null(remap)) {
      new_i <- integer(0)
      new_x <- numeric(0)
      new_p <- integer(length(run) + 1L)
      new_p[1L] <- 0L
      for (j in seq_along(run)) {
        a <- p_local[j] + 1L
        b <- p_local[j + 1L]
        if (b >= a) {
          ii <- i_vals[a:b]
          xx <- x_vals[a:b]
          mapped <- remap[ii + 1L]
          keep_e <- !is.na(mapped)
          new_i <- c(new_i, as.integer(mapped[keep_e]))
          new_x <- c(new_x, xx[keep_e])
        }
        new_p[j + 1L] <- length(new_i)
      }
      i_vals <- new_i
      x_vals <- new_x
      p_local <- new_p
    }

    Matrix::sparseMatrix(
      i = i_vals + 1L,
      p = p_local,
      x = x_vals,
      dims = c(nrow_out, length(run)),
      dimnames = list(rn, cn),
      index1 = TRUE
    )
  }

  if (length(cols) > 1L && (is.unsorted(cols) || any(diff(cols) != 1L))) {
    runs <- .ic2_contig_runs(cols)
    mats <- lapply(runs, .read_run)
    return(.ic2_to_dgCMatrix(do.call(cbind, mats)))
  }
  .read_run(cols)
}

#' Read a 1D HDF5 dataset hyperslab (1-based start; supports large offsets)
#' @keywords internal
.ic2_h5_read_1d <- function(ds, start1, count) {
  count <- as.numeric(count)[1L]
  start1 <- as.numeric(start1)[1L]
  if (!is.finite(count) || count < 1) {
    return(NULL)
  }
  if (!is.finite(start1) || start1 < 1) {
    stop("Invalid hyperslab start: ", start1)
  }
  end1 <- start1 + count - 1
  ## Bracket indexing is reliable with current hdf5r. read_low_level hyperslabs
  ## were observed to treat start=0 and start=1 as the same first element,
  ## shifting all later column chunks (wrong colnames / CSC slices).
  if (end1 <= .Machine$integer.max && start1 <= .Machine$integer.max) {
    return(ds[as.integer(start1):as.integer(end1)])
  }
  fs <- ds$get_space()
  fs$select_hyperslab(start = start1 - 1, count = count)
  ms <- hdf5r::H5S$new(dims = count, maxdims = count)
  ds$read_low_level(file_space = fs, mem_space = ms)
}

#' Write a 1D HDF5 dataset hyperslab (1-based start; supports large offsets)
#' @keywords internal
.ic2_h5_write_1d <- function(ds, start1, values) {
  n <- length(values)
  if (n < 1L) return(invisible(NULL))
  start1 <- as.numeric(start1)[1L]
  if (!is.finite(start1) || start1 < 1) {
    stop("Invalid hyperslab write start: ", start1)
  }
  end1 <- start1 + as.numeric(n) - 1
  ## Prefer bracket assignment  -  write_low_level hyperslabs are unreliable here
  if (end1 <= .Machine$integer.max && start1 <= .Machine$integer.max) {
    ds[as.integer(start1):as.integer(end1)] <- values
  } else {
    fs <- ds$get_space()
    fs$select_hyperslab(start = start1 - 1, count = as.numeric(n))
    ms <- hdf5r::H5S$new(dims = as.numeric(n), maxdims = as.numeric(n))
    ds$write_low_level(values, file_space = fs, mem_space = ms)
  }
  invisible(NULL)
}

#' Convert a 10x HDF5 matrix to iCellR2 CSC HDF5 without building a dgCMatrix
#'
#' Streams column chunks from a modern (\code{/matrix}) or legacy 10x H5 file
#' into an on-disk CSC assay. Supports matrices whose nnz exceeds R's 32-bit
#' sparse index limit by storing column pointers \code{p} as double.
#'
#' @param src_path Path to 10x \code{.h5}.
#' @param dest_path Output iCellR2 CSC \code{.h5}.
#' @param group Destination group name.
#' @param feature.field Feature name field for modern layout.
#' @param feature.type Optional feature types to keep (modern layout).
#' @param overwrite Overwrite destination.
#' @param chunk.cells Columns per streaming chunk.
#' @param progress Optional \code{function(current, total)} callback.
#' @return An \code{iCellR2_assay_handle}.
#' @keywords internal
.ic2_convert_10x_h5_to_csc_hdf5 <- function(
    src_path,
    dest_path,
    group = "counts",
    feature.field = "name",
    feature.type = NULL,
    overwrite = FALSE,
    chunk.cells = 2000L,
    progress = NULL,
    row.map = NULL,
    dest.features = NULL,
    cell.ids = NULL
) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  if (file.exists(dest_path)) {
    if (!isTRUE(overwrite)) {
      stop("File already exists at '", dest_path, "'. Use overwrite = TRUE.")
    }
    unlink(dest_path)
  }
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  chunk.cells <- as.integer(max(1L, as.integer(chunk.cells)))

  to_num <- function(x) {
    if (inherits(x, "integer64")) return(as.numeric(x))
    as.numeric(x)
  }

  h_src <- hdf5r::H5File$new(src_path, mode = "r")
  on.exit(try(h_src$close_all(), silent = TRUE), add = TRUE)

  legacy <- FALSE
  if (h_src$exists("matrix")) {
    g <- h_src[["matrix"]]
    needed <- c("data", "indices", "indptr", "shape", "barcodes", "features")
    missing_needed <- needed[!vapply(needed, g$exists, logical(1))]
    if (length(missing_needed)) {
      stop("Invalid modern 10x HDF5 file: missing ", paste(missing_needed, collapse = ", "))
    }
    shape <- to_num(g[["shape"]]$read())
    barcodes <- as.character(g[["barcodes"]]$read())
    indptr <- to_num(g[["indptr"]]$read())
    fg <- g[["features"]]
    if (!fg$exists(feature.field)) {
      stop("Requested feature field '", feature.field, "' not found in HDF5 file.")
    }
    features <- as.character(fg[[feature.field]]$read())
    feature_type_all <- if (fg$exists("feature_type")) {
      as.character(fg[["feature_type"]]$read())
    } else {
      NULL
    }
    ds_i_src <- g[["indices"]]
    ds_x_src <- g[["data"]]
  } else {
    legacy <- TRUE
    top <- names(h_src)
    if (!length(top)) stop("HDF5 file appears empty.")
    g <- h_src[[top[1]]]
    needed <- c("data", "indices", "indptr", "shape", "barcodes")
    missing_needed <- needed[!vapply(needed, g$exists, logical(1))]
    if (length(missing_needed)) {
      stop("Invalid legacy 10x HDF5 file: missing ", paste(missing_needed, collapse = ", "))
    }
    shape <- to_num(g[["shape"]]$read())
    barcodes <- as.character(g[["barcodes"]]$read())
    indptr <- to_num(g[["indptr"]]$read())
    legacy_map <- c(id = "genes", name = "gene_names")
    selected_field <- legacy_map[[feature.field]]
    if (is.null(selected_field)) {
      if (identical(feature.field, "feature_type")) {
        stop("'feature.field = \"feature_type\"' is not available for legacy 10x HDF5 input.")
      }
      stop("Unsupported feature field for legacy 10x HDF5 input: ", feature.field)
    }
    if (!g$exists(selected_field)) {
      stop("Requested feature field '", selected_field, "' not found in legacy HDF5 file.")
    }
    features <- as.character(g[[selected_field]]$read())
    feature_type_all <- rep("Gene Expression", length(features))
    ds_i_src <- g[["indices"]]
    ds_x_src <- g[["data"]]
  }

  nr0 <- as.integer(shape[1])
  nc <- as.integer(shape[2])
  if (length(indptr) != nc + 1L) {
    stop("Invalid CSC structure: length(indptr) != ncol + 1.")
  }
  if (length(features) != nr0) stop("Feature count does not match matrix rows.")
  if (length(barcodes) != nc) stop("Barcode count does not match matrix columns.")

  ## Optional feature-type filter  ->  row keep map (0-based source  ->  0-based dest)
  row_map <- NULL
  if (!is.null(feature.type)) {
    matched <- .ic2_match_feature_types(feature_type_all, feature.type)
    feature.type <- matched$feature.type
    keep <- matched$keep
    features <- features[keep]
    row_map <- rep(NA_integer_, nr0)
    row_map[keep] <- as.integer(seq_len(sum(keep)) - 1L)
  }
  ## Optional full-matrix remap (0-based destination row, NA = drop).
  ## Used to place sample peaks on a consensus set without loading the matrix.
  if (!is.null(row.map)) {
    if (length(row.map) != nr0) {
      stop("'row.map' must have one entry per source feature.")
    }
    if (is.null(dest.features) || !length(dest.features)) {
      stop("'dest.features' is required when 'row.map' is set.")
    }
    row_map <- as.integer(row.map)
    features <- as.character(dest.features)
  }
  nr <- length(features)

  barcodes <- .ic2_clean_names(barcodes)
  if (!is.null(cell.ids)) {
    if (length(cell.ids) != length(barcodes)) {
      stop("'cell.ids' must have one id per source barcode.")
    }
    barcodes <- as.character(cell.ids)
  }
  if (is.null(row.map)) features <- .ic2_clean_names(features)

  nnz_src <- indptr[length(indptr)]
  ## Destination nnz: same if no filter; unknown upper bound = nnz_src if filtering
  nnz_est <- nnz_src

  h_out <- hdf5r::H5File$new(dest_path, mode = "w")
  on.exit({
    try(h_src$close_all(), silent = TRUE)
    try(h_out$close_all(), silent = TRUE)
  }, add = FALSE)

  g_out <- h_out$create_group(group)
  g_out[["format"]] <- "iCellR2_csc_v1"
  g_out[["dims"]] <- as.integer(c(nr, nc))
  g_out[["rownames"]] <- as.character(features)
  g_out[["colnames"]] <- as.character(barcodes)

  ## p always as double so atlas-scale nnz is representable
  chunk_p <- as.integer(min(length(indptr), max(1000L, min(1e6L, length(indptr)))))
  ds_p <- g_out$create_dataset(
    "p",
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = length(indptr),
    chunk_dims = chunk_p
  )
  chunk_nx <- as.numeric(min(nnz_est, 1e7))
  if (chunk_nx < 1) chunk_nx <- 1
  ## Pre-allocate to source nnz; shrink via attribute of true nnz if filtered
  ds_i <- g_out$create_dataset(
    "i",
    dtype = hdf5r::h5types$H5T_NATIVE_INT,
    dims = max(nnz_est, 1),
    chunk_dims = chunk_nx
  )
  ds_x <- g_out$create_dataset(
    "x",
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = max(nnz_est, 1),
    chunk_dims = chunk_nx
  )

  out_nnz <- 0
  p_out <- numeric(nc + 1L)
  p_out[1L] <- 0

  n_chunks <- as.integer(ceiling(nc / chunk.cells))
  for (ci in seq_len(n_chunks)) {
    c0 <- (ci - 1L) * chunk.cells + 1L
    c1 <- min(ci * chunk.cells, nc)
    ## Source nnz range for columns c0..c1 (0-based indptr)
    src_off <- indptr[c0]
    n_src <- indptr[c1 + 1L] - src_off
    s0 <- src_off + 1

    if (is.null(row_map)) {
      ## Fast path: copy CSC values for this column block as-is
      if (n_src > 0) {
        idx <- as.integer(.ic2_h5_read_1d(ds_i_src, start1 = s0, count = n_src))
        xx <- as.numeric(.ic2_h5_read_1d(ds_x_src, start1 = s0, count = n_src))
        .ic2_h5_write_1d(ds_i, start1 = out_nnz + 1, values = idx)
        .ic2_h5_write_1d(ds_x, start1 = out_nnz + 1, values = xx)
      }
      p_out[(c0 + 1L):(c1 + 1L)] <- indptr[(c0 + 1L):(c1 + 1L)]
      out_nnz <- indptr[c1 + 1L]
    } else {
      if (n_src > 0) {
        idx <- as.integer(.ic2_h5_read_1d(ds_i_src, start1 = s0, count = n_src))
        xx <- as.numeric(.ic2_h5_read_1d(ds_x_src, start1 = s0, count = n_src))
      } else {
        idx <- integer(0)
        xx <- numeric(0)
      }
      ip_chunk <- indptr[c0:(c1 + 1L)] - src_off
      for (j in seq_len(c1 - c0 + 1L)) {
        a <- ip_chunk[j] + 1L
        b <- ip_chunk[j + 1L]
        if (b < a) {
          p_out[c0 + j] <- out_nnz
          next
        }
        ii <- idx[a:b]
        xv <- xx[a:b]
        new_i <- row_map[ii + 1L]
        keep_e <- !is.na(new_i)
        ii <- as.integer(new_i[keep_e])
        xv <- xv[keep_e]
        if (length(ii) > 1L) {
          ord <- order(ii)
          ii <- ii[ord]
          xv <- xv[ord]
          if (anyDuplicated(ii)) {
            change <- which(c(TRUE, ii[-1L] != ii[-length(ii)]))
            ends <- c(change[-1L] - 1L, length(ii))
            xv <- vapply(seq_along(change), function(k) sum(xv[change[k]:ends[k]]), numeric(1))
            ii <- ii[change]
          }
        }
        n_add <- length(ii)
        if (n_add > 0L) {
          .ic2_h5_write_1d(ds_i, start1 = out_nnz + 1, values = ii)
          .ic2_h5_write_1d(ds_x, start1 = out_nnz + 1, values = xv)
          out_nnz <- out_nnz + n_add
        }
        p_out[c0 + j] <- out_nnz
      }
    }

    if (is.function(progress)) {
      progress(c1, nc)
    }
  }

  .ic2_h5_write_1d(ds_p, start1 = 1, values = as.numeric(p_out))
  g_out[["nnz"]] <- as.numeric(out_nnz)

  h_out$flush()
  h_out$close_all()
  h_src$close_all()
  on.exit(NULL)

  .ic2_assay_handle(
    path = dest_path,
    backend = "hdf5_sparse",
    group = group,
    n_features = nr,
    n_cells = nc,
    feature_ids = NULL,
    cell_ids = NULL
  )
}

#' Write a peaks-only or genes-only 10x HDF5 from a mixed multiome file.
#'
#' Streams. The SODA assay file contains one feature type and no row filter,
#' so RNA \code{main.data} and ATAC normalization never see each other's rows.
#' @keywords internal
.ic2_split_10x_feature_type <- function(
    src_path,
    dest_path,
    feature.type,
    cell.ids = NULL,
    sample_id = NULL,
    modality = "rna"
) {
  if (file.exists(dest_path)) unlink(dest_path)
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = "ic2_split_", fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  .ic2_convert_10x_h5_to_csc_hdf5(
    src_path = src_path,
    dest_path = tmp,
    feature.type = feature.type,
    overwrite = TRUE,
    cell.ids = cell.ids
  )
  label <- as.character(feature.type)[1L]
  if (!nzchar(label)) label <- if (identical(modality, "atac")) "Peaks" else "Gene Expression"
  .ic2_stream_csc_to_tenx_h5(
    csc_path = tmp,
    dest_path = dest_path,
    feature_type = label
  )
  meta_path <- file.path(dirname(dest_path), "tenx_meta.json")
  .ic2_tenx_h5_handle(
    src_path = dest_path,
    meta_path = meta_path,
    feature.field = "name",
    feature.type = NULL,
    sample_id = sample_id,
    modality = modality,
    overwrite = TRUE
  )
}

#' Copy an iCellR2 CSC file into a 10x HDF5 without building a dgCMatrix.
#' @keywords internal
.ic2_stream_csc_to_tenx_h5 <- function(csc_path, dest_path, feature_type, group = "counts") {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  src <- hdf5r::H5File$new(csc_path, mode = "r")
  on.exit(try(src$close_all(), silent = TRUE), add = TRUE)
  g <- src[[group]]
  features <- as.character(g[["rownames"]]$read())
  barcodes <- as.character(g[["colnames"]]$read())
  p <- as.numeric(g[["p"]]$read())
  nnz <- if (g$exists("nnz")) as.numeric(g[["nnz"]]$read())[1L] else p[length(p)]
  if (!is.finite(nnz) || nnz < 0) nnz <- p[length(p)]
  nr <- length(features)
  nc <- length(barcodes)

  if (file.exists(dest_path)) unlink(dest_path)
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  h5 <- hdf5r::H5File$new(dest_path, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  mg <- h5$create_group("matrix")
  fg <- mg$create_group("features")
  mg[["barcodes"]] <- barcodes
  mg[["shape"]] <- as.integer(c(nr, nc))
  mg[["indptr"]] <- p
  chunk <- as.integer(min(max(nnz, 1), 5000000))
  ds_data <- mg$create_dataset(
    "data", dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = max(nnz, 1), chunk_dims = chunk
  )
  ds_idx <- mg$create_dataset(
    "indices", dtype = hdf5r::h5types$H5T_NATIVE_INT,
    dims = max(nnz, 1), chunk_dims = chunk
  )
  if (nnz > 0) {
    pos <- 0
    while (pos < nnz) {
      take <- min(as.numeric(chunk), nnz - pos)
      .ic2_h5_write_1d(ds_idx, pos + 1, as.integer(.ic2_h5_read_1d(g[["i"]], pos + 1, take)))
      .ic2_h5_write_1d(ds_data, pos + 1, as.numeric(.ic2_h5_read_1d(g[["x"]], pos + 1, take)))
      pos <- pos + take
    }
  }
  fg[["id"]] <- features
  fg[["name"]] <- features
  fg[["feature_type"]] <- rep(as.character(feature_type)[1L], nr)
  fg[["genome"]] <- rep("unknown", nr)
  h5$flush()
  invisible(normalizePath(dest_path, winslash = "/", mustWork = TRUE))
}

#' Locate matrix.mtx / barcodes / features in a 10x directory
#' @keywords internal
.ic2_find_10x_dir_files <- function(path) {
  matrix.file <- list.files(path = path, pattern = "^matrix\\.mtx(\\.gz)?$", full.names = TRUE)
  barcode.file <- list.files(path = path, pattern = "^barcodes\\.tsv(\\.gz)?$", full.names = TRUE)
  feature.file <- list.files(path = path, pattern = "^(genes|features)\\.tsv(\\.gz)?$", full.names = TRUE)
  if (length(matrix.file) != 1L) stop("Expected exactly one matrix.mtx(.gz) in '", path, "'.")
  if (length(barcode.file) != 1L) stop("Expected exactly one barcodes.tsv(.gz) in '", path, "'.")
  if (length(feature.file) != 1L) {
    stop("Expected exactly one genes.tsv(.gz) or features.tsv(.gz) in '", path, "'.")
  }
  list(matrix = matrix.file, barcodes = barcode.file, features = feature.file)
}

#' Open a text or gzipped file for reading
#' @keywords internal
.ic2_open_text <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
}

#' Read Matrix Market header: returns list(nr, nc, nnz, con) with con positioned at first triplet
#' @keywords internal
.ic2_mtx_open_header <- function(mtx_path) {
  con <- .ic2_open_text(mtx_path)
  banner <- readLines(con, n = 1L, warn = FALSE)
  if (!length(banner) || !grepl("MatrixMarket", banner, fixed = TRUE)) {
    close(con)
    stop("Not a Matrix Market file: '", mtx_path, "'.")
  }
  if (!grepl("coordinate", banner, ignore.case = TRUE)) {
    close(con)
    stop("Only coordinate Matrix Market files are supported: '", mtx_path, "'.")
  }
  ## Skip comment lines
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) {
      close(con)
      stop("Unexpected end of MTX header in '", mtx_path, "'.")
    }
    if (!startsWith(ln, "%")) {
      dims <- scan(text = ln, what = double(), quiet = TRUE)
      if (length(dims) < 3L) {
        close(con)
        stop("Invalid MTX size line: '", ln, "'.")
      }
      return(list(
        nr = as.numeric(dims[1L]),
        nc = as.numeric(dims[2L]),
        nnz = as.numeric(dims[3L]),
        con = con,
        path = mtx_path
      ))
    }
  }
}

#' Estimate available RAM in bytes (best-effort)
#' @keywords internal
.ic2_available_ram_bytes <- function() {
  ## Linux
  if (file.exists("/proc/meminfo")) {
    lines <- tryCatch(readLines("/proc/meminfo", n = 40L, warn = FALSE), error = function(e) character())
    avail <- lines[grepl("^MemAvailable:", lines)]
    if (length(avail) >= 1L) {
      kb <- suppressWarnings(as.numeric(sub(".*?([0-9]+).*", "\\1", avail[1L])))
      if (is.finite(kb) && kb > 0) return(kb * 1024)
    }
    total <- lines[grepl("^MemTotal:", lines)]
    if (length(total) >= 1L) {
      kb <- suppressWarnings(as.numeric(sub(".*?([0-9]+).*", "\\1", total[1L])))
      if (is.finite(kb) && kb > 0) return(kb * 1024 * 0.5)
    }
  }
  ## macOS: physical memory; treat half as usable
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    sysctl <- Sys.which("sysctl")
    if (nzchar(sysctl)) {
      out <- suppressWarnings(tryCatch(
        system2(sysctl, c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE),
        error = function(e) character()
      ))
      bytes <- suppressWarnings(as.numeric(out[1L]))
      if (is.finite(bytes) && bytes > 0) return(bytes * 0.5)
    }
  }
  ## Conservative fallback
  4 * 1024^3
}

#' Estimated bytes for a full in-memory CSC (i int + x double + p double)
#' @keywords internal
.ic2_mtx_csc_ram_bytes <- function(nnz, ncol) {
  nnz <- as.numeric(nnz)
  ncol <- as.numeric(ncol)
  ## i(4) + x(8) + p(8)*(ncol+1) + ~20% overhead for R SEXPs / temps
  (nnz * 12 + (ncol + 1) * 8) * 1.2
}

#' Write CSC arrays (i, p, x) to iCellR2 HDF5
#' @keywords internal
.ic2_write_csc_arrays_hdf5 <- function(
    filepath,
    i,
    p,
    x,
    rownames,
    colnames,
    group = "counts",
    overwrite = FALSE,
    progress = NULL
) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  nr <- length(rownames)
  nc <- length(colnames)
  nnz <- length(i)
  if (length(x) != nnz) stop("length(i) != length(x).")
  if (length(p) != nc + 1L) stop("length(p) != ncol + 1.")
  if (file.exists(filepath)) {
    if (!isTRUE(overwrite)) {
      stop("File already exists at '", filepath, "'. Use overwrite = TRUE.")
    }
    unlink(filepath)
  }
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)

  h5 <- hdf5r::H5File$new(filepath, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5$create_group(group)
  g[["format"]] <- "iCellR2_csc_v1"
  g[["dims"]] <- as.integer(c(nr, nc))
  g[["rownames"]] <- as.character(rownames)
  g[["colnames"]] <- as.character(colnames)

  ## Chunked writes keep peak temp lower and allow progress ticks
  chunk <- as.integer(min(max(nnz, 1), 5000000L))
  ds_i <- g$create_dataset(
    "i",
    dtype = hdf5r::h5types$H5T_NATIVE_INT,
    dims = max(nnz, 1),
    chunk_dims = chunk
  )
  ds_x <- g$create_dataset(
    "x",
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = max(nnz, 1),
    chunk_dims = chunk
  )
  if (nnz > 0) {
    if (is.function(progress)) progress(0, nnz)
    pos <- 0
    while (pos < nnz) {
      take <- min(chunk, nnz - pos)
      idx <- (pos + 1):(pos + take)
      .ic2_h5_write_1d(ds_i, start1 = pos + 1, values = as.integer(i[idx]))
      .ic2_h5_write_1d(ds_x, start1 = pos + 1, values = as.numeric(x[idx]))
      pos <- pos + take
      if (is.function(progress)) progress(pos, nnz)
    }
  } else if (is.function(progress)) {
    progress(0, 1)
    progress(1, 1)
  }
  g[["p"]] <- as.numeric(p)
  g[["nnz"]] <- as.numeric(nnz)
  h5$flush()
  invisible(.ic2_assay_handle(
    path = filepath,
    backend = "hdf5_sparse",
    group = group,
    n_features = nr,
    n_cells = nc,
    feature_ids = NULL,
    cell_ids = NULL
  ))
}

#' Convert 10x MTX directory  ->  iCellR2 CSC HDF5 via Rcpp (auto memory/stream)
#'
#' @param on_phase Optional \code{function(phase)} called with
#'   \code{"parse"}, \code{"write"}, or \code{"stream"} before that phase.
#' @keywords internal
.ic2_convert_mtx_dir_rcpp <- function(
    dir_path,
    dest_path,
    group = "counts",
    feature.column = 2L,
    overwrite = FALSE,
    progress = NULL,
    on_phase = NULL,
    force = c("auto", "memory", "stream")
) {
  force <- match.arg(force)
  if (!exists("cpp_mtx_to_csc", mode = "function") ||
      !exists("cpp_mtx_info", mode = "function") ||
      !exists("cpp_mtx_stream_sorted", mode = "function")) {
    stop("Rcpp MTX parsers are not available; rebuild the package.")
  }
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }

  files <- .ic2_find_10x_dir_files(dir_path)
  feature.column <- as.integer(feature.column)[1L]

  bc_con <- .ic2_open_text(files$barcodes)
  on.exit(try(close(bc_con), silent = TRUE), add = TRUE)
  barcodes <- readLines(bc_con, warn = FALSE)
  close(bc_con)

  ft_con <- .ic2_open_text(files$features)
  on.exit(try(close(ft_con), silent = TRUE), add = TRUE)
  gene.info <- utils::read.table(
    file = ft_con,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  close(ft_con)
  on.exit(NULL)

  if (feature.column < 1L || feature.column > ncol(gene.info)) {
    stop("'feature.column' = ", feature.column, " out of range (", ncol(gene.info), " columns).")
  }
  features <- .ic2_clean_names(as.character(gene.info[[feature.column]]))
  barcodes <- .ic2_clean_names(barcodes)

  info <- cpp_mtx_info(files$matrix)
  nr <- as.integer(info$nrow)
  nc <- as.integer(info$ncol)
  nnz <- as.numeric(info$nnz)
  if (length(features) != nr) {
    stop("Features (", length(features), ") != MTX rows (", nr, ").")
  }
  if (length(barcodes) != nc) {
    stop("Barcodes (", length(barcodes), ") != MTX columns (", nc, ").")
  }

  need <- .ic2_mtx_csc_ram_bytes(nnz, nc)
  avail <- .ic2_available_ram_bytes()
  use_memory <- switch(
    force,
    memory = TRUE,
    stream = FALSE,
    auto = is.finite(need) && is.finite(avail) && need < avail * 0.5
  )

  prog <- progress
  mode_label <- if (use_memory) "memory" else "stream"
  call_phase <- function(phase) {
    if (is.function(on_phase)) on_phase(phase)
  }

  if (use_memory) {
    call_phase("parse")
    csc <- cpp_mtx_to_csc(files$matrix, progress = prog)
    call_phase("write")
    handle <- .ic2_write_csc_arrays_hdf5(
      filepath = dest_path,
      i = csc$i,
      p = csc$p,
      x = csc$x,
      rownames = features,
      colnames = barcodes,
      group = group,
      overwrite = overwrite,
      progress = prog
    )
    attr(handle, "mtx_mode") <- paste0("rcpp-", mode_label)
    return(handle)
  }

  ## Low-RAM: stream sorted MTX chunks into preallocated HDF5
  call_phase("stream")
  if (file.exists(dest_path)) {
    if (!isTRUE(overwrite)) {
      stop("File already exists at '", dest_path, "'. Use overwrite = TRUE.")
    }
    unlink(dest_path)
  }
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)

  h5 <- hdf5r::H5File$new(dest_path, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5$create_group(group)
  g[["format"]] <- "iCellR2_csc_v1"
  g[["dims"]] <- as.integer(c(nr, nc))
  g[["rownames"]] <- as.character(features)
  g[["colnames"]] <- as.character(barcodes)

  chunk <- as.integer(min(max(as.numeric(nnz), 1), 5000000))
  ds_i <- g$create_dataset(
    "i",
    dtype = hdf5r::h5types$H5T_NATIVE_INT,
    dims = max(nnz, 1),
    chunk_dims = chunk
  )
  ds_x <- g$create_dataset(
    "x",
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = max(nnz, 1),
    chunk_dims = chunk
  )
  env <- new.env(parent = emptyenv())
  env$out_nnz <- 0

  flush_chunk <- function(ii, xx) {
    n <- length(ii)
    if (n < 1L) return(invisible(NULL))
    .ic2_h5_write_1d(ds_i, start1 = env$out_nnz + 1, values = as.integer(ii))
    .ic2_h5_write_1d(ds_x, start1 = env$out_nnz + 1, values = as.numeric(xx))
    env$out_nnz <- env$out_nnz + n
    invisible(NULL)
  }

  streamed <- cpp_mtx_stream_sorted(
    files$matrix,
    flush_chunk = flush_chunk,
    progress = prog,
    chunk_nnz = 5000000
  )

  if (!isTRUE(streamed$sorted)) {
    ## Partial writes may exist; fall back to full-memory if it can fit
    try(h5$close_all(), silent = TRUE)
    on.exit(NULL)
    unlink(dest_path)
    if (need < avail * 0.85) {
      call_phase("parse")
      csc <- cpp_mtx_to_csc(files$matrix, progress = prog)
      call_phase("write")
      handle <- .ic2_write_csc_arrays_hdf5(
        filepath = dest_path,
        i = csc$i,
        p = csc$p,
        x = csc$x,
        rownames = features,
        colnames = barcodes,
        group = group,
        overwrite = TRUE,
        progress = prog
      )
      attr(handle, "mtx_mode") <- "rcpp-memory-unsorted-fallback"
      return(handle)
    }
    stop(
      "MTX is not column-sorted and estimated CSC RAM (",
      sprintf("%.1f GB", need / 1024^3),
      ") exceeds available memory. ",
      "Sort the MTX (Cell Ranger style) or use a machine with more RAM."
    )
  }

  g[["p"]] <- as.numeric(streamed$p)
  g[["nnz"]] <- as.numeric(env$out_nnz)
  h5$flush()
  h5$close_all()
  on.exit(NULL)

  handle <- .ic2_assay_handle(
    path = dest_path,
    backend = "hdf5_sparse",
    group = group,
    n_features = nr,
    n_cells = nc,
    feature_ids = NULL,
    cell_ids = NULL
  )
  attr(handle, "mtx_mode") <- paste0("rcpp-", mode_label)
  handle
}

#' Stream a 10x MTX directory to iCellR2 CSC HDF5 (no full dgCMatrix)
#'
#' Prefer column-sorted MTX (Cell Ranger default): streams triplets and writes
#' CSC hyperslabs. If the file is not sorted by column, falls back to a
#' coordinate load + order (still avoids a \code{dgCMatrix} S4 object).
#'
#' @keywords internal
.ic2_convert_mtx_dir_to_csc_hdf5 <- function(
    dir_path,
    dest_path,
    group = "counts",
    feature.column = 2L,
    overwrite = FALSE,
    chunk.lines = 2000000L,
    progress = NULL
) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  files <- .ic2_find_10x_dir_files(dir_path)
  feature.column <- as.integer(feature.column)[1L]

  bc_con <- .ic2_open_text(files$barcodes)
  on.exit(try(close(bc_con), silent = TRUE), add = TRUE)
  barcodes <- readLines(bc_con, warn = FALSE)
  close(bc_con)

  ft_con <- .ic2_open_text(files$features)
  on.exit(try(close(ft_con), silent = TRUE), add = TRUE)
  gene.info <- utils::read.table(
    file = ft_con,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  close(ft_con)
  on.exit(NULL)
  if (feature.column < 1L || feature.column > ncol(gene.info)) {
    stop("'feature.column' = ", feature.column, " out of range (", ncol(gene.info), " columns).")
  }
  features <- as.character(gene.info[[feature.column]])
  barcodes <- .ic2_clean_names(barcodes)
  features <- .ic2_clean_names(features)

  hdr <- .ic2_mtx_open_header(files$matrix)
  on.exit(try(close(hdr$con), silent = TRUE), add = TRUE)
  nr <- as.integer(hdr$nr)
  nc <- as.integer(hdr$nc)
  nnz <- as.numeric(hdr$nnz)
  if (!is.finite(nnz) || nnz < 0) stop("Invalid MTX nnz.")
  if (length(features) != nr) {
    stop("Features (", length(features), ") != MTX rows (", nr, ").")
  }
  if (length(barcodes) != nc) {
    stop("Barcodes (", length(barcodes), ") != MTX columns (", nc, ").")
  }

  if (file.exists(dest_path)) {
    if (!overwrite) stop("File already exists at '", dest_path, "'. Use overwrite = TRUE.")
    unlink(dest_path)
  }
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)

  h_out <- hdf5r::H5File$new(dest_path, mode = "w")
  on.exit({
    try(close(hdr$con), silent = TRUE)
    try(h_out$close_all(), silent = TRUE)
  }, add = FALSE)

  g_out <- h_out$create_group(group)
  g_out[["format"]] <- "iCellR2_csc_v1"
  g_out[["dims"]] <- as.integer(c(nr, nc))
  g_out[["rownames"]] <- as.character(features)
  g_out[["colnames"]] <- as.character(barcodes)

  chunk_p <- as.integer(min(nc + 1L, max(1000L, min(1e6L, nc + 1L))))
  ds_p <- g_out$create_dataset(
    "p",
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = nc + 1L,
    chunk_dims = chunk_p
  )
  nnz_store <- max(nnz, 1)
  chunk_nx <- as.numeric(min(nnz_store, 1e7))
  if (chunk_nx < 1) chunk_nx <- 1
  ds_i <- g_out$create_dataset(
    "i",
    dtype = hdf5r::h5types$H5T_NATIVE_INT,
    dims = nnz_store,
    chunk_dims = chunk_nx
  )
  ds_x <- g_out$create_dataset(
    "x",
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = nnz_store,
    chunk_dims = chunk_nx
  )

  ## Use an environment so buffer updates are reliable (subsetted <<- is unsafe)
  buf <- new.env(parent = emptyenv())
  buf$cap <- 500000L
  buf$i <- integer(buf$cap)
  buf$x <- numeric(buf$cap)
  buf$n <- 0L
  buf$out_nnz <- 0
  buf$completed_col <- 0L
  buf$p_out <- numeric(nc + 1L)
  buf$p_out[1L] <- 0

  flush_buf <- function() {
    if (buf$n < 1L) return(invisible(NULL))
    .ic2_h5_write_1d(ds_i, start1 = buf$out_nnz + 1, values = buf$i[seq_len(buf$n)])
    .ic2_h5_write_1d(ds_x, start1 = buf$out_nnz + 1, values = buf$x[seq_len(buf$n)])
    buf$out_nnz <- buf$out_nnz + buf$n
    buf$n <- 0L
    invisible(NULL)
  }
  append_buf <- function(ii, xx) {
    n_add <- length(ii)
    if (n_add < 1L) return(invisible(NULL))
    pos <- 0L
    while (pos < n_add) {
      space <- buf$cap - buf$n
      if (space < 1L) {
        flush_buf()
        space <- buf$cap
      }
      take <- min(space, n_add - pos)
      idx <- (buf$n + 1L):(buf$n + take)
      src <- (pos + 1L):(pos + take)
      buf$i[idx] <- as.integer(ii[src])
      buf$x[idx] <- as.numeric(xx[src])
      buf$n <- buf$n + take
      pos <- pos + take
      if (buf$n >= buf$cap) flush_buf()
    }
    invisible(NULL)
  }
  close_cols_through <- function(last_col) {
    nnz_now <- buf$out_nnz + buf$n
    while (buf$completed_col < last_col) {
      buf$completed_col <- buf$completed_col + 1L
      buf$p_out[buf$completed_col + 1L] <- nnz_now
    }
  }

  sorted_ok <- TRUE
  n_read <- 0
  while (n_read < nnz && sorted_ok) {
    take <- as.integer(min(chunk.lines, nnz - n_read))
    trip <- scan(
      hdr$con,
      what = list(i = integer(), j = integer(), x = double()),
      nmax = take,
      quiet = TRUE
    )
    got <- length(trip$i)
    if (got < 1L) break
    n_read <- n_read + got

    ii <- trip$i
    jj <- trip$j
    xx <- trip$x
    if (any(ii < 1L | ii > nr | jj < 1L | jj > nc)) {
      close(hdr$con)
      h_out$close_all()
      on.exit(NULL)
      stop("MTX index out of range while streaming '", files$matrix, "'.")
    }
    if (jj[1L] < buf$completed_col + 1L || (got > 1L && any(diff(jj) < 0L))) {
      sorted_ok <- FALSE
      break
    }

    starts <- c(1L, which(diff(jj) != 0L) + 1L)
    ends <- c(starts[-1L] - 1L, got)
    for (r in seq_along(starts)) {
      a <- starts[r]
      b <- ends[r]
      cc <- jj[a]
      if (cc > buf$completed_col + 1L) {
        flush_buf()
        close_cols_through(cc - 1L)
      }
      append_buf(ii[a:b] - 1L, xx[a:b])
      next_cc <- if (r < length(starts)) jj[starts[r + 1L]] else NA_integer_
      if (!is.na(next_cc) && next_cc != cc) {
        flush_buf()
        close_cols_through(cc)
      } else if (is.na(next_cc) && n_read >= nnz) {
        flush_buf()
        close_cols_through(cc)
      }
    }

    if (is.function(progress)) progress(min(n_read, nnz), nnz)
  }

  if (sorted_ok) {
    flush_buf()
    close_cols_through(nc)
    out_nnz <- buf$out_nnz
    p_out <- buf$p_out
    if (out_nnz > nnz) {
      stop("Wrote more entries than MTX nnz header (", out_nnz, " > ", nnz, ").")
    }
    .ic2_h5_write_1d(ds_p, start1 = 1, values = as.numeric(p_out))
    g_out[["nnz"]] <- as.numeric(out_nnz)
    h_out$flush()
    h_out$close_all()
    try(close(hdr$con), silent = TRUE)
    on.exit(NULL)
    if (is.function(progress)) progress(nnz, nnz)
    return(.ic2_assay_handle(
      path = dest_path,
      backend = "hdf5_sparse",
      group = group,
      n_features = nr,
      n_cells = nc,
      feature_ids = NULL,
      cell_ids = NULL
    ))
  }

  ## ---- Unsorted fallback: coordinate vectors  ->  order  ->  write CSC ----
  try(close(hdr$con), silent = TRUE)
  h_out$close_all()
  unlink(dest_path)
  on.exit(NULL)

  if (is.function(progress)) {
    ## Restart progress for fallback
  }
  hdr2 <- .ic2_mtx_open_header(files$matrix)
  on.exit(try(close(hdr2$con), silent = TRUE), add = TRUE)
  i_all <- integer(nnz)
  j_all <- integer(nnz)
  x_all <- numeric(nnz)
  filled <- 0
  while (filled < nnz) {
    take <- as.integer(min(chunk.lines, nnz - filled))
    trip <- scan(
      hdr2$con,
      what = list(i = integer(), j = integer(), x = double()),
      nmax = take,
      quiet = TRUE
    )
    got <- length(trip$i)
    if (got < 1L) break
    idx <- (filled + 1):(filled + got)
    i_all[idx] <- trip$i
    j_all[idx] <- trip$j
    x_all[idx] <- trip$x
    filled <- filled + got
    if (is.function(progress)) progress(filled, nnz * 2)
  }
  close(hdr2$con)
  on.exit(NULL)
  if (filled != nnz) {
    i_all <- i_all[seq_len(filled)]
    j_all <- j_all[seq_len(filled)]
    x_all <- x_all[seq_len(filled)]
    nnz <- filled
  }
  ord <- order(j_all, i_all)
  i_all <- as.integer(i_all[ord] - 1L)
  x_all <- as.numeric(x_all[ord])
  j_all <- j_all[ord]
  p_out <- numeric(nc + 1L)
  ## tabulate column counts
  counts <- tabulate(j_all, nbins = nc)
  p_out[1L] <- 0
  p_out[-1L] <- cumsum(counts)

  h_out <- hdf5r::H5File$new(dest_path, mode = "w")
  on.exit(try(h_out$close_all(), silent = TRUE), add = TRUE)
  g_out <- h_out$create_group(group)
  g_out[["format"]] <- "iCellR2_csc_v1"
  g_out[["dims"]] <- as.integer(c(nr, nc))
  g_out[["rownames"]] <- as.character(features)
  g_out[["colnames"]] <- as.character(barcodes)
  g_out[["i"]] <- i_all
  g_out[["x"]] <- x_all
  g_out[["p"]] <- as.numeric(p_out)
  g_out[["nnz"]] <- as.numeric(length(i_all))
  h_out$flush()
  h_out$close_all()
  on.exit(NULL)
  if (is.function(progress)) progress(nnz * 2, nnz * 2)

  .ic2_assay_handle(
    path = dest_path,
    backend = "hdf5_sparse",
    group = group,
    n_features = nr,
    n_cells = nc,
    feature_ids = NULL,
    cell_ids = NULL
  )
}

#' Write a dgCMatrix to CSC HDF5 (CRAN default backend)
#'
#' @param mat A \code{dgCMatrix}.
#' @param filepath Output \code{.h5} path.
#' @param group Group name inside the file.
#' @param overwrite Overwrite existing file.
#' @return An \code{iCellR2_assay_handle}.
#' @export
write_sparse_hdf5 <- function(mat, filepath, group = "counts", overwrite = FALSE) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required for the default HDF5 sparse backend.")
  }
  mat <- .ic2_to_dgCMatrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("'mat' must have rownames and colnames.")
  }
  if (!is.character(filepath) || length(filepath) != 1L || !nzchar(filepath)) {
    stop("'filepath' must be a non-empty character string.")
  }
  if (file.exists(filepath)) {
    if (!overwrite) {
      stop("File already exists at '", filepath, "'. Use overwrite = TRUE.")
    }
    unlink(filepath)
  }
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)

  h5 <- hdf5r::H5File$new(filepath, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  if (group %in% names(h5)) {
    h5$link_delete(group)
  }
  g <- h5$create_group(group)

  g[["format"]] <- "iCellR2_csc_v1"
  g[["dims"]] <- as.integer(c(nrow(mat), ncol(mat)))
  g[["i"]] <- as.integer(mat@i)
  g[["p"]] <- as.integer(mat@p)
  g[["x"]] <- as.numeric(mat@x)
  g[["rownames"]] <- as.character(rownames(mat))
  g[["colnames"]] <- as.character(colnames(mat))

  invisible(.ic2_assay_handle(
    path = filepath,
    backend = "hdf5_sparse",
    group = group,
    n_features = nrow(mat),
    n_cells = ncol(mat),
    feature_ids = NULL,
    cell_ids = NULL
  ))
}

#' Write a dense numeric matrix to HDF5 (for scaled / dense assays)
#'
#' Scaled data is nearly dense (many tiny non-zeros). Storing it as CSC
#' can effectively become all-zeros after drop0; use this for
#' \code{scale.data}.
#'
#' @param mat Numeric matrix (features  x  cells) with dimnames.
#' @param filepath Output \code{.h5} path.
#' @param group Group name inside the file.
#' @param overwrite Overwrite existing file.
#' @return An \code{iCellR2_assay_handle} with backend \code{hdf5_dense}.
#' @export
write_dense_hdf5 <- function(mat, filepath, group = "counts", overwrite = FALSE) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required for the dense HDF5 backend.")
  }
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("'mat' must have rownames and colnames.")
  }
  if (!is.character(filepath) || length(filepath) != 1L || !nzchar(filepath)) {
    stop("'filepath' must be a non-empty character string.")
  }
  if (file.exists(filepath)) {
    if (!overwrite) {
      stop("File already exists at '", filepath, "'. Use overwrite = TRUE.")
    }
    unlink(filepath)
  }
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)

  nr <- nrow(mat)
  nc <- ncol(mat)
  h5 <- hdf5r::H5File$new(filepath, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  if (group %in% names(h5)) h5$link_delete(group)
  g <- h5$create_group(group)
  g[["format"]] <- "iCellR2_dense_v1"
  g[["dims"]] <- as.integer(c(nr, nc))
  ## IMPORTANT: g[["x"]] <- mat flattens to 1D in hdf5r; create a true 2D dataset
  chunk_c <- as.integer(min(nc, max(1L, min(4096L, nc))))
  g$create_dataset(
    name = "x",
    robj = mat,
    chunk_dims = c(nr, chunk_c)
  )
  g[["rownames"]] <- as.character(rownames(mat))
  g[["colnames"]] <- as.character(colnames(mat))

  ## Sanity: refuse to write an all-zero scaled matrix
  if (!any(mat != 0)) {
    warning("write_dense_hdf5: matrix is all zeros  -  check upstream scaling.")
  }

  invisible(.ic2_assay_handle(
    path = filepath,
    backend = "hdf5_dense",
    group = group,
    n_features = nr,
    n_cells = nc,
    feature_ids = NULL,
    cell_ids = NULL
  ))
}

#' Read a contiguous column block from dense HDF5 as a matrix
#' @param filepath Path to the dense HDF5 file.
#' @param group HDF5 group. Default \code{"counts"}.
#' @param cols Contiguous 1-based column indices.
#' @export
read_dense_hdf5_cols <- function(filepath, group = "counts", cols) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  cols <- as.integer(cols)
  if (length(cols) < 1L || any(is.na(cols)) || any(cols < 1L)) {
    stop("'cols' must be positive 1-based column indices.")
  }
  h5 <- hdf5r::H5File$new(filepath, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5[[group]]
  dims <- as.integer(g[["dims"]][])
  nrow_m <- dims[1L]
  ncol_m <- dims[2L]
  if (max(cols) > ncol_m) {
    stop("Requested columns exceed matrix ncol (", ncol_m, ").")
  }

  ds <- g[["x"]]
  ds_dims <- as.integer(ds$dims)
  if (length(ds_dims) == 2L) {
    block <- ds[, cols, drop = FALSE]
  } else if (length(ds_dims) == 1L && ds_dims[1L] == nrow_m * ncol_m) {
    ## Legacy broken flat write  -  reshape column-major then subset
    flat <- ds[]
    full <- matrix(flat, nrow = nrow_m, ncol = ncol_m)
    block <- full[, cols, drop = FALSE]
  } else {
    stop(
      "Dense HDF5 dataset 'x' has unexpected dims: ",
      paste(ds_dims, collapse = "x"),
      " (expected ", nrow_m, "x", ncol_m, "). Re-run data.scale(overwrite = TRUE)."
    )
  }
  if (!is.matrix(block)) {
    block <- matrix(block, nrow = nrow_m, ncol = length(cols))
  }
  storage.mode(block) <- "double"
  rownames(block) <- as.character(g[["rownames"]][])
  colnames(block) <- as.character(g[["colnames"]][cols])
  block
}

#' Detect whether an HDF5 file is 10x Genomics or iCellR2 CSC
#' @return One of \code{"tenx_h5"}, \code{"icellr2_csc"}, \code{"unknown"}
#' @keywords internal
.ic2_detect_h5_kind <- function(filepath) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  if (!file.exists(filepath)) return("unknown")
  h5 <- tryCatch(hdf5r::H5File$new(filepath, mode = "r"), error = function(e) NULL)
  if (is.null(h5)) return("unknown")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  ## iCellR2 CSC: group "counts" (or any group) with format / i,p,x
  for (nm in names(h5)) {
    g <- tryCatch(h5[[nm]], error = function(e) NULL)
    if (is.null(g) || !inherits(g, "H5Group")) next
    nms <- names(g)
    if ("format" %in% nms) {
      fmt <- as.character(g[["format"]][1L])
      if (identical(fmt, "iCellR2_csc_v1")) return("icellr2_csc")
      if (identical(fmt, "iCellR2_dense_v1")) return("icellr2_dense")
    }
    if (all(c("i", "p", "x", "dims") %in% nms)) return("icellr2_csc")
  }

  ## Modern 10x
  if ("matrix" %in% names(h5)) {
    g <- h5[["matrix"]]
    if (all(c("data", "indices", "indptr", "shape", "barcodes") %in% names(g))) {
      return("tenx_h5")
    }
  }

  ## Legacy 10x: top-level genome group with data/indices/indptr/barcodes
  top <- names(h5)
  if (length(top) >= 1L) {
    g <- tryCatch(h5[[top[1L]]], error = function(e) NULL)
    if (!is.null(g) && inherits(g, "H5Group")) {
      if (all(c("data", "indices", "indptr", "shape", "barcodes") %in% names(g))) {
        return("tenx_h5")
      }
    }
  }
  "unknown"
}

#' Read format string from an iCellR2 HDF5 assay group
#' @keywords internal
.ic2_hdf5_format <- function(filepath, group = "counts") {
  h5 <- hdf5r::H5File$new(filepath, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5[[group]]
  if ("format" %in% names(g)) {
    return(as.character(g[["format"]][1L]))
  }
  if ("i" %in% names(g) && "p" %in% names(g) && "x" %in% names(g)) {
    return("iCellR2_csc_v1")
  }
  if ("x" %in% names(g)) return("iCellR2_dense_v1")
  "unknown"
}

#' Read dimnames / dims from a sparse HDF5 assay
#' @keywords internal
.ic2_hdf5_dims <- function(filepath, group = "counts") {
  h5 <- hdf5r::H5File$new(filepath, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5[[group]]
  dims <- g[["dims"]][]
  list(
    n_features = as.integer(dims[1L]),
    n_cells = as.integer(dims[2L]),
    feature_ids = as.character(g[["rownames"]][]),
    cell_ids = as.character(g[["colnames"]][])
  )
}

#' Read a contiguous column block from CSC HDF5 as dgCMatrix
#'
#' Split 1-based indices into contiguous runs (input should be ordered)
#' @keywords internal
.ic2_contig_runs <- function(idx) {
  idx <- as.integer(idx)
  if (length(idx) < 1L) return(list())
  if (length(idx) == 1L) return(list(idx))
  brk <- which(diff(idx) != 1L)
  starts <- c(1L, brk + 1L)
  ends <- c(brk, length(idx))
  lapply(seq_along(starts), function(k) idx[starts[k]:ends[k]])
}

#' Read selected columns from a CSC HDF5 file
#'
#' @param filepath HDF5 path.
#' @param group Group name.
#' @param cols 1-based column indices (preferably contiguous).
#' @return A \code{dgCMatrix} for the requested columns.
#' @export
read_sparse_hdf5_cols <- function(filepath, group = "counts", cols) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  cols <- as.integer(cols)
  if (length(cols) < 1L || any(is.na(cols)) || any(cols < 1L)) {
    stop("'cols' must be positive 1-based column indices.")
  }

  ## One open per call. Non-contiguous cols  ->  contig runs under the same file
  ## handle (never open/close once per filtered cell).
  h5 <- hdf5r::H5File$new(filepath, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5[[group]]

  dims <- as.integer(g[["dims"]][])
  nrow_m <- dims[1L]
  ncol_m <- dims[2L]
  if (max(cols) > ncol_m) {
    stop("Requested columns exceed matrix ncol (", ncol_m, ").")
  }
  rn <- as.character(g[["rownames"]][])

  .read_run <- function(run) {
    p_slice <- as.numeric(.ic2_h5_read_1d(
      g[["p"]],
      start1 = run[1L],
      count = length(run) + 1L
    ))
    x_start0 <- p_slice[1L]
    x_end0 <- p_slice[length(p_slice)]
    n_x <- x_end0 - x_start0
    if (n_x < 1) {
      i_vals <- integer(0)
      x_vals <- numeric(0)
    } else {
      i_vals <- as.integer(.ic2_h5_read_1d(g[["i"]], start1 = x_start0 + 1, count = n_x))
      x_vals <- as.numeric(.ic2_h5_read_1d(g[["x"]], start1 = x_start0 + 1, count = n_x))
    }
    p_local <- as.integer(p_slice - p_slice[1L])
    cn <- as.character(.ic2_h5_read_1d(g[["colnames"]], start1 = run[1L], count = length(run)))
    Matrix::sparseMatrix(
      i = i_vals + 1L,
      p = p_local,
      x = x_vals,
      dims = c(nrow_m, length(run)),
      dimnames = list(rn, cn),
      index1 = TRUE
    )
  }

  if (length(cols) > 1L && (is.unsorted(cols) || any(diff(cols) != 1L))) {
    runs <- .ic2_contig_runs(cols)
    mats <- lapply(runs, .read_run)
    return(.ic2_to_dgCMatrix(do.call(cbind, mats)))
  }
  .read_run(cols)
}

#' Read selected feature rows for given columns from CSC HDF5
#' @keywords internal
.ic2_hdf5_rowcol_block <- function(filepath, group, row_idx, cols) {
  block <- read_sparse_hdf5_cols(filepath, group = group, cols = cols)
  if (length(row_idx) == 0L) {
    return(NULL)
  }
  block[row_idx, , drop = FALSE]
}

#' Resolve an assay input to an on-disk handle
#'
#' Raw counts are written as modern \strong{10x HDF5} (\code{tenx_h5},
#' \code{/matrix}) by default when \code{backend = "hdf5_sparse"}.
#' Use \code{backend = "rds"} for an RDS matrix dump.
#'
#' @keywords internal
.ic2_materialize_handle <- function(
    x,
    filepath,
    group = "counts",
    overwrite = TRUE,
    backend = c("hdf5_sparse", "rds"),
    sample_id = NULL,
    modality = "rna"
) {
  backend <- match.arg(backend)
  if (.ic2_is_assay_handle(x)) {
    return(x)
  }
  if (inherits(x, "iCellR2_sample") && !is.null(x$handle)) {
    return(x$handle)
  }
  if (inherits(x, "iCellR2_sample") && !is.null(x$assay)) {
    x <- x$assay
  }

  mat <- .ic2_to_dgCMatrix(x)
  if (backend == "hdf5_sparse") {
    ## Prefer Cell Ranger-style 10x H5 for raw assays in SODA / agg temp
    h5_path <- if (grepl("\\.h5$", filepath, ignore.case = TRUE)) {
      filepath
    } else {
      paste0(tools::file_path_sans_ext(filepath), ".h5")
    }
    return(.ic2_write_matrix_tenx_h5(
      mat,
      filepath = h5_path,
      sample_id = sample_id,
      modality = modality,
      overwrite = overwrite
    ))
  }

  ## Never write matrix .rds into SODA  -  always materialize as 10x H5
  h5_path <- if (grepl("\\.h5$", filepath, ignore.case = TRUE)) {
    filepath
  } else {
    paste0(tools::file_path_sans_ext(filepath), ".h5")
  }
  .ic2_write_matrix_tenx_h5(
    mat,
    filepath = h5_path,
    sample_id = sample_id,
    modality = modality,
    overwrite = overwrite
  )
}

#' Ensure handle has dims; optionally load dimnames
#' @keywords internal
.ic2_enrich_handle <- function(handle, load_dimnames = FALSE) {
  if (!.ic2_is_assay_handle(handle)) {
    stop("Expected an iCellR2_assay_handle.")
  }
  if (identical(handle$backend, "tenx_h5")) {
    ## Recover tenx_meta.json next to the matrix when open.soda omitted meta_path
    if (is.null(handle$meta_path) || !nzchar(as.character(handle$meta_path)[1L]) ||
        !file.exists(handle$meta_path)) {
      if (!is.null(handle$path) && nzchar(as.character(handle$path)[1L])) {
        meta_cands <- c(
          file.path(dirname(handle$path), "tenx_meta.json"),
          file.path(dirname(handle$path), "tenx_meta.rds")
        )
        meta_hit <- meta_cands[file.exists(meta_cands)]
        if (length(meta_hit)) handle$meta_path <- meta_hit[[1L]]
      }
    }
    meta <- NULL
    if (!is.null(handle$meta_path) && file.exists(handle$meta_path)) {
      meta <- .ic2_soda_read_tenx_meta(handle$meta_path)
    }
    if (!is.null(meta)) {
      handle$n_features <- as.integer(meta$n_features)
      handle$n_cells <- as.integer(meta$n_cells)
      handle$group <- meta$group %||% handle$group
      if (load_dimnames) {
        ## Prefer in-memory overrides (e.g. prefixed barcodes after aggregation)
        handle$feature_ids <- handle$feature_ids %||% as.character(meta$feature_ids %||% character(0))
        handle$cell_ids <- handle$cell_ids %||% as.character(meta$cell_ids %||% character(0))
      }
    }
    return(handle)
  }
  if (handle$backend %in% c("hdf5_sparse", "hdf5_dense") && file.exists(handle$path)) {
    need_dims <- is.null(handle$n_features) || is.null(handle$n_cells) ||
      !is.finite(as.numeric(handle$n_features)[1L]) || !is.finite(as.numeric(handle$n_cells)[1L])
    need_names <- isTRUE(load_dimnames) &&
      (is.null(handle$feature_ids) || is.null(handle$cell_ids) ||
         length(handle$feature_ids) < 1L || length(handle$cell_ids) < 1L)
    ## Avoid reopening HDF5 on every tiny column batch once dims/ids are known
    if (need_dims || need_names) {
      ## Auto-detect dense vs CSC if format was upgraded on disk
      fmt <- tryCatch(
        .ic2_hdf5_format(handle$path, group = handle$group),
        error = function(e) NA_character_
      )
      if (identical(fmt, "iCellR2_dense_v1")) {
        handle$backend <- "hdf5_dense"
      } else if (identical(fmt, "iCellR2_csc_v1")) {
        handle$backend <- "hdf5_sparse"
      }
      info <- .ic2_hdf5_dims(handle$path, group = handle$group)
      handle$n_features <- info$n_features
      handle$n_cells <- info$n_cells
      if (need_names) {
        handle$feature_ids <- info$feature_ids
        handle$cell_ids <- info$cell_ids
      }
    }
  } else if (identical(handle$backend, "rds") && file.exists(handle$path)) {
    mat <- readRDS(handle$path)
    handle$n_features <- nrow(mat)
    handle$n_cells <- ncol(mat)
    if (load_dimnames) {
      handle$feature_ids <- rownames(mat)
      handle$cell_ids <- colnames(mat)
    }
  }
  handle
}

#' Read a column chunk from any supported assay handle
#' @keywords internal
.ic2_read_cols <- function(handle, cols) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = FALSE)
  if (identical(handle$backend, "hdf5_dense")) {
    return(read_dense_hdf5_cols(handle$path, group = handle$group, cols = cols))
  }
  if (identical(handle$backend, "hdf5_sparse")) {
    return(read_sparse_hdf5_cols(handle$path, group = handle$group, cols = cols))
  }
  if (identical(handle$backend, "tenx_h5")) {
    if (is.null(handle$meta_path) || !file.exists(handle$meta_path)) {
      stop("tenx_h5 handle is missing meta_path.")
    }
    meta <- .ic2_soda_read_tenx_meta(handle$meta_path)
    if (is.null(meta)) stop("Could not read tenx_h5 meta: ", handle$meta_path)
    ## Allow runtime overrides of dimnames (aggregation prefixes)
    if (!is.null(handle$feature_ids)) meta$feature_ids <- handle$feature_ids
    if (!is.null(handle$cell_ids)) meta$cell_ids <- handle$cell_ids
    return(read_tenx_hdf5_cols(handle$path, meta = meta, cols = cols))
  }
  if (identical(handle$backend, "rds")) {
    mat <- readRDS(handle$path)
    return(mat[, cols, drop = FALSE])
  }
  stop("Unsupported assay backend: ", handle$backend)
}

#' Get feature / cell ids for a handle
#' @keywords internal
.ic2_handle_dimnames <- function(handle) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  list(feature_ids = handle$feature_ids, cell_ids = handle$cell_ids)
}

print.iCellR2_assay_handle <- function(x, ...) {
  cat(
    "<iCellR2 assay handle>\n",
    "  backend: ", x$backend, "\n",
    "  path: ", x$path, "\n",
    "  dims: ", x$n_features, " features x ", x$n_cells, " cells\n",
    sep = ""
  )
  invisible(x)
}

## ---- Assay resolution / subset / streaming helpers --------------------

#' @keywords internal
.ic2_require_object <- function(x) {
  if (is.null(x)) {
    stop("'x' must be an object of class 'iCellR2'.")
  }
  ## Proper S4 required  -  a plain list/S3 with class "iCellR2" will fail on @
  if (!isS4(x) || !methods::is(x, "iCellR2")) {
    stop(
      "'x' must be a valid S4 object of class 'iCellR2'.\n",
      "  class(x) = ", paste(class(x), collapse = "/"), "\n",
      "  isS4(x)  = ", isS4(x), "\n",
      "If you reloaded the package with devtools::load_all() / reinstall, ",
      "recreate the object with make.obj() (old S4 objects can become invalid)."
    )
  }
  invisible(x)
}

#' @keywords internal
.ic2_resolve_assay_name <- function(x, modality, assay.name = NULL) {
  modality <- tolower(modality)
  if (is.null(x@modalities[[modality]])) {
    stop("Modality '", modality, "' was not found.")
  }
  if (is.null(assay.name)) {
    assay.name <- x@modalities[[modality]]$active_assay
  }
  if (is.null(assay.name) || !nzchar(assay.name)) {
    stop("No assay name provided and no active assay for modality '", modality, "'.")
  }
  if (is.null(x@modalities[[modality]]$assays[[assay.name]])) {
    stop("Assay '", assay.name, "' not found in modality '", modality, "'.")
  }
  assay.name
}

#' Map classic / friendly \code{data.type} labels to assay names.
#'
#' Accepts both slot-style names (\code{main.data}, \code{raw.data}, ...)
#' and short aliases (\code{main}, \code{raw}, \code{imputed}, ...).
#' Unknown labels return \code{NULL} (caller may warn / fall back).
#'
#' @keywords internal
.ic2_data_type_to_assay <- function(data.type) {
  if (is.null(data.type) || !nzchar(as.character(data.type)[1L])) {
    return(NULL)
  }
  dt <- tolower(trimws(as.character(data.type)[1L]))
  # allow dots / underscores interchangeably
  dt <- gsub("_", ".", dt, fixed = TRUE)
  switch(
    dt,
    # raw counts
    "raw.data" = "counts",
    "raw" = "counts",
    "counts" = "counts",
    # normalized main / data
    "main.data" = "data",
    "main" = "data",
    "data" = "data",
    "norm.data" = "data",
    "normalized" = "data",
    "norm" = "data",
    # imputed
    "imputed.data" = "imputed",
    "imputed" = "imputed",
    # scaled
    "scaled.data" = "scale.data",
    "scale.data" = "scale.data",
    "scaled" = "scale.data",
    "scale" = "scale.data",
    # ATAC
    "atac" = "atac",
    "atac.imputed" = "atac.imputed",
    "atac.imputed.data" = "atac.imputed",
    NULL
  )
}

#' Resolve assay from optional \code{assay.name} and/or \code{data.type}.
#' \code{assay.name} wins when set. Default \code{data.type} is
#' \code{"main.data"} (normalized \code{"data"} assay when present).
#'
#' @keywords internal
.ic2_resolve_plot_assay <- function(
    x,
    modality,
    assay.name = NULL,
    data.type = "main.data",
    prefer.data = TRUE
) {
  modality <- tolower(modality)
  if (!is.null(assay.name) && nzchar(as.character(assay.name)[1L])) {
    return(.ic2_resolve_assay_name(x, modality, as.character(assay.name)[1L]))
  }
  mapped <- .ic2_data_type_to_assay(data.type)
  if (is.null(mapped) && !is.null(data.type) && nzchar(as.character(data.type)[1L])) {
    warning(
      "Unknown data.type '", as.character(data.type)[1L],
      "'; falling back to default assay. ",
      "Try: raw.data, main.data, imputed.data, scaled.data."
    )
  }
  assay.names <- names(x@modalities[[modality]]$assays %||% list())
  if (is.null(mapped)) {
    if (isTRUE(prefer.data) && "data" %in% assay.names) {
      mapped <- "data"
    } else {
      mapped <- NULL
    }
  }
  .ic2_resolve_assay_name(x, modality, mapped)
}

#' Resolve assay for gene.stats: raw counts by default, or normalized "data"
#' @keywords internal
.ic2_resolve_gene_stats_assay <- function(x, modality, assay.name = NULL, use.raw = TRUE) {
  modality <- tolower(modality)
  assays <- x@modalities[[modality]]$assays %||% list()
  assay_names <- names(assays)

  if (!is.null(assay.name) && nzchar(as.character(assay.name)[1L])) {
    return(.ic2_resolve_assay_name(x, modality, assay.name))
  }

  if (isTRUE(use.raw)) {
    if ("counts" %in% assay_names) return("counts")
    active <- x@modalities[[modality]]$active_assay
    if (!is.null(active) && active %in% assay_names) return(active)
    stop(
      "use.raw = TRUE but no 'counts' (or active) assay found in modality '",
      modality, "'."
    )
  }

  if ("data" %in% assay_names) return("data")
  stop(
    "Normalized assay 'data' not found. Run norm.data() first, ",
    "or keep use.raw = TRUE / assay.name = \"counts\" for raw counts."
  )
}

#' Infer backend + restore tenx_h5 meta_path from a filepath on disk
#' @keywords internal
.ic2_handle_from_path <- function(path, group = "counts", backend = NULL) {
  path <- as.character(path)[1L]
  if (!nzchar(path)) stop("Empty assay path.")

  ## Prefer explicit tenx pointer sidecar (written by make.obj)
  ptr_cands <- c(
    file.path(dirname(path), "tenx_pointer.json"),
    file.path(dirname(path), "tenx_pointer.rds") ## legacy
  )
  for (ptr in ptr_cands) {
    if (!file.exists(ptr)) next
    h <- if (grepl("\\.json$", ptr, ignore.case = TRUE)) {
      meta <- jsonlite::fromJSON(ptr, simplifyVector = TRUE, simplifyDataFrame = FALSE)
      sh <- list(
        type = "iCellR2_assay_handle",
        backend = "tenx_h5",
        path = meta$source %||% meta$path %||% path,
        group = meta$group %||% "matrix",
        n_features = meta$n_features,
        n_cells = meta$n_cells,
        meta_path = meta$meta_path
      )
      class(sh) <- c("iCellR2_assay_handle", "list")
      sh
    } else {
      readRDS(ptr)
    }
    if (.ic2_is_assay_handle(h) && identical(h$backend, "tenx_h5")) return(h)
  }

  ## Sidecar meta next to import path: <stem>_tenx_meta.json (.rds legacy)
  stem <- sub("\\.[^.]+$", "", basename(path))
  cand_meta <- c(
    file.path(dirname(path), "tenx_meta.json"),
    file.path(dirname(path), paste0(stem, "_tenx_meta.json")),
    paste0(sub("\\.[^.]+$", "", path), "_tenx_meta.json"),
    file.path(dirname(path), "tenx_meta.rds"),
    file.path(dirname(path), paste0(stem, "_tenx_meta.rds")),
    paste0(sub("\\.[^.]+$", "", path), "_tenx_meta.rds")
  )
  for (mp in cand_meta) {
    if (file.exists(mp)) {
      h <- .ic2_assay_handle(path = path, backend = "tenx_h5", group = group)
      h$meta_path <- normalizePath(mp, winslash = "/", mustWork = FALSE)
      return(h)
    }
  }

  if (is.null(backend)) {
    backend <- if (grepl("\\.h5$", path, ignore.case = TRUE)) "hdf5_sparse" else "hdf5_sparse"
  }
  .ic2_assay_handle(path = path, backend = backend, group = group)
}

#' @keywords internal
.ic2_get_assay_handles <- function(x, modality, assay.name) {
  assay_ref <- x@modalities[[modality]]$assays[[assay.name]]
  if (is.list(assay_ref) && identical(assay_ref$aggregation_type, "lazy_norm")) {
    src <- assay_ref$source_assay %||% "counts"
    return(.ic2_get_assay_handles(x, modality, src))
  }
  if (is.list(assay_ref) && identical(assay_ref$aggregation_type, "virtual")) {
    if (!is.null(assay_ref$sample_handles)) {
      return(assay_ref$sample_handles)
    }
    if (!is.null(assay_ref$sample_assay_paths)) {
      return(lapply(assay_ref$sample_assay_paths, function(p) {
        .ic2_handle_from_path(p, group = "counts")
      }))
    }
    stop("Virtual assay is missing sample handles/paths.")
  }
  if (is.list(assay_ref) && !is.null(assay_ref$handle)) {
    return(list(assay_ref$handle))
  }
  if (.ic2_is_assay_handle(assay_ref)) {
    return(list(assay_ref))
  }
  if (inherits(assay_ref, "Matrix") || is.matrix(assay_ref) || is.data.frame(assay_ref)) {
    tmp <- tempfile(fileext = ".h5")
    return(list(write_sparse_hdf5(.ic2_to_dgCMatrix(assay_ref), tmp, overwrite = TRUE)))
  }
  ap <- x@storage$assay_paths
  if (!is.null(ap$filepath) && file.exists(ap$filepath)) {
    backend <- if (!is.null(ap$backend)) ap$backend else NULL
    h <- .ic2_handle_from_path(
      ap$filepath,
      group = if (!is.null(ap$group)) ap$group else "counts",
      backend = backend
    )
    if (!is.null(ap$meta_path) && identical(h$backend, "tenx_h5")) {
      h$meta_path <- ap$meta_path
    }
    return(list(h))
  }
  stop("Could not resolve on-disk assay handles.")
}

#' @keywords internal
.ic2_assay_is_virtual <- function(x, modality, assay.name) {
  assay_ref <- x@modalities[[modality]]$assays[[assay.name]]
  is.list(assay_ref) && identical(assay_ref$aggregation_type, "virtual")
}

#' @keywords internal
.ic2_soda_assay_path <- function(x, modality, assay.name, sample_id = NULL, ext = ".h5") {
  base <- file.path(x@storage$root_dir, "modalities", modality, "assays", assay.name)
  if (is.null(sample_id) || !nzchar(as.character(sample_id)[1L])) {
    return(file.path(base, paste0("assay_source", ext)))
  }
  sid <- as.character(sample_id)[1L]
  file.path(base, sid, paste0(sid, ext))
}

#' Subset an assay handle to selected cell ids; write a new on-disk assay
#' @keywords internal
.ic2_subset_handle_cells <- function(
    handle,
    keep_cell_ids,
    dest_path,
    group = "counts",
    chunk.size = 200000L,
    progress_cb = NULL
) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  cell_ids <- handle$cell_ids
  keep_idx <- match(keep_cell_ids, cell_ids)
  keep_idx <- keep_idx[!is.na(keep_idx)]
  if (length(keep_idx) < 1L) {
    stop("No overlapping cells between keep list and assay handle '", handle$path, "'.")
  }
  ## Prefer original assay column order for contiguous HDF5 runs
  keep_idx <- sort(unique(keep_idx))

  ## Fast path: dropped few cells  -  read covering ranges as contiguous runs
  runs <- .ic2_contig_runs(keep_idx)
  n_runs <- length(runs)
  blocks <- vector("list", n_runs)
  for (i in seq_along(runs)) {
    run <- runs[[i]]
    ## Further split huge runs by chunk.size
    if (length(run) <= chunk.size) {
      blocks[[i]] <- .ic2_read_cols(handle, run)
    } else {
      sub <- list()
      starts <- seq.int(1L, length(run), by = as.integer(chunk.size))
      for (s in starts) {
        e <- min(s + as.integer(chunk.size) - 1L, length(run))
        sub[[length(sub) + 1L]] <- .ic2_read_cols(handle, run[s:e])
      }
      blocks[[i]] <- .ic2_to_dgCMatrix(do.call(cbind, sub))
    }
    if (is.function(progress_cb)) progress_cb(i, n_runs)
  }
  mat <- .ic2_to_dgCMatrix(do.call(cbind, blocks))
  ## Restore keep_cell_ids order if it differed from assay order
  want <- intersect(keep_cell_ids, colnames(mat))
  if (!identical(want, colnames(mat))) {
    mat <- mat[, want, drop = FALSE]
  }
  rm(blocks)
  invisible(gc)
  h <- write_sparse_hdf5(mat, filepath = dest_path, group = group, overwrite = TRUE)
  attr(h, "n_runs") <- n_runs
  h
}

#' Cells failing gene-expression OR filter (chunked)
#' @keywords internal
.ic2_gene_filter_fail_ids <- function(
    handle,
    genes,
    min_exp = 1,
    chunk.size = 200000L
) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  gene_idx <- match(unique(as.character(genes)), handle$feature_ids)
  if (any(is.na(gene_idx))) {
    missing <- unique(as.character(genes))[is.na(match(unique(as.character(genes)), handle$feature_ids))]
    stop(
      "Assay lacks filter.by.gene feature(s): ",
      paste(missing, collapse = ", ")
    )
  }
  n_cells <- handle$n_cells
  fail <- logical(n_cells)
  starts <- seq.int(1L, n_cells, by = as.integer(chunk.size))
  for (start in starts) {
    end <- min(start + as.integer(chunk.size) - 1L, n_cells)
    cols <- start:end
    block <- .ic2_read_cols(handle, cols)
    sub <- block[gene_idx, , drop = FALSE]
    ok <- Matrix::colSums(sub >= min_exp) > 0
    fail[cols] <- !as.logical(ok)
  }
  handle$cell_ids[fail]
}

#' Streaming column sums for a handle
#' @keywords internal
.ic2_handle_col_sums <- function(handle, row_idx = NULL, chunk.size = 200000L, use.rcpp = TRUE) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = is.null(row_idx))
  n_cells <- handle$n_cells
  out <- numeric(n_cells)
  starts <- seq.int(1L, n_cells, by = as.integer(chunk.size))
  for (start in starts) {
    end <- min(start + as.integer(chunk.size) - 1L, n_cells)
    cols <- start:end
    block <- .ic2_read_cols(handle, cols)
    if (!is.null(row_idx)) {
      block <- block[row_idx, , drop = FALSE]
    }
    if (inherits(block, "dgCMatrix") && isTRUE(use.rcpp) &&
        exists("cpp_sparse_col_stats", mode = "function")) {
      out[cols] <- as.numeric(cpp_sparse_col_stats(block)$col_sums)
    } else {
      out[cols] <- as.numeric(Matrix::colSums(block))
    }
  }
  out
}

#' Per-cell library sizes: reuse QC \code{UMIs} when available, else column sums
#'
#' @param cell_ids Character cell ids to return (order preserved).
#' @param meta \code{cell_data} (or similar) with optional \code{UMIs}.
#' @param handle Assay handle used if UMIs are missing / incomplete.
#' @param warn.fallback If \code{TRUE}, warn when falling back to streaming
#'   column sums (expensive on large 10x matrices).
#' @return Named numeric library sizes for \code{cell_ids}.
#' @keywords internal
.ic2_cell_library_sizes <- function(
    cell_ids,
    meta = NULL,
    handle = NULL,
    chunk.size = 200000L,
    use.rcpp = TRUE,
    warn.fallback = TRUE
) {
  cell_ids <- as.character(cell_ids)
  if (length(cell_ids) < 1L) {
    return(numeric(0))
  }

  umi <- .ic2_match_umis_from_meta(cell_ids, meta)
  if (!is.null(umi) && all(is.finite(umi) & umi > 0)) {
    attr(umi, "source") <- "UMIs"
    return(umi)
  }

  if (is.null(handle)) {
    stop("Library sizes unavailable: no usable UMIs in cell_data and no assay handle.")
  }
  if (isTRUE(warn.fallback)) {
    n_miss <- if (is.null(umi)) {
      length(cell_ids)
    } else {
      sum(!is.finite(umi) | umi <= 0)
    }
    warning(
      "QC UMIs not reusable for ", format(n_miss, big.mark = ","),
      " / ", format(length(cell_ids), big.mark = ","),
      " cells  -  streaming column sums (slow on large matrices). ",
      "Run qc.stats() before norm.data() and keep cell ids aligned.",
      call. = FALSE
    )
  }
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  lib_all <- .ic2_handle_col_sums(
    handle, chunk.size = chunk.size, use.rcpp = use.rcpp
  )
  names(lib_all) <- handle$cell_ids
  out <- as.numeric(lib_all[cell_ids])
  names(out) <- cell_ids
  if (any(!is.finite(out) | out <= 0)) {
    stop(
      "Non-positive library sizes for ",
      sum(!is.finite(out) | out <= 0),
      " cells after column sums."
    )
  }
  attr(out, "source") <- "colSums"
  out
}

#' Match QC UMIs to cell ids (rownames, cell_id col, sample_prefix strip/add)
#' @keywords internal
.ic2_match_umis_from_meta <- function(cell_ids, meta) {
  cell_ids <- as.character(cell_ids)
  if (!is.data.frame(meta) || !"UMIs" %in% colnames(meta) || nrow(meta) < 1L) {
    return(NULL)
  }

  umi_by_id <- NULL
  if (!is.null(rownames(meta)) && any(nzchar(rownames(meta)))) {
    umi_by_id <- stats::setNames(as.numeric(meta$UMIs), as.character(rownames(meta)))
  }
  if ("cell_id" %in% colnames(meta)) {
    cid <- as.character(meta$cell_id)
    u <- as.numeric(meta$UMIs)
    if (is.null(umi_by_id)) {
      umi_by_id <- stats::setNames(u, cid)
    } else {
      miss <- which(!cid %in% names(umi_by_id))
      if (length(miss)) {
        umi_by_id <- c(umi_by_id, stats::setNames(u[miss], cid[miss]))
      }
    }
  }
  if (is.null(umi_by_id) || !length(umi_by_id)) return(NULL)

  ## Direct match
  out <- unname(umi_by_id[cell_ids])
  names(out) <- cell_ids
  if (all(is.finite(out) & out > 0)) return(out)

  ## Strip first sample_ prefix from query ids (aggregation / make.obj style)
  sp <- regexpr("_", cell_ids, fixed = TRUE)
  has_pref <- sp > 0L
  if (any(has_pref)) {
    bare <- cell_ids
    bare[has_pref] <- substring(cell_ids[has_pref], sp[has_pref] + 1L)
    alt <- unname(umi_by_id[bare])
    fill <- which((!is.finite(out) | out <= 0) & is.finite(alt) & alt > 0)
    if (length(fill)) out[fill] <- alt[fill]
    if (all(is.finite(out) & out > 0)) {
      names(out) <- cell_ids
      return(out)
    }
  }

  ## Add common sample_ prefixes from meta onto bare query ids
  meta_ids <- names(umi_by_id)
  msp <- regexpr("_", meta_ids, fixed = TRUE)
  prefs <- unique(substr(meta_ids[msp > 0L], 1L, msp[msp > 0L] - 1L))
  prefs <- prefs[nzchar(prefs)]
  if (length(prefs)) {
    for (pref in prefs) {
      cand <- paste0(pref, "_", cell_ids)
      alt <- unname(umi_by_id[cand])
      fill <- which((!is.finite(out) | out <= 0) & is.finite(alt) & alt > 0)
      if (length(fill)) out[fill] <- alt[fill]
      if (all(is.finite(out) & out > 0)) {
        names(out) <- cell_ids
        return(out)
      }
    }
  }

  names(out) <- cell_ids
  if (all(is.finite(out) & out > 0)) {
    return(out)
  }
  NULL
}

#' Streaming row sums across cells (optionally a column subset)
#' @keywords internal
.ic2_handle_row_sums <- function(handle, chunk.size = 200000L, col_idx = NULL) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  n_feat <- handle$n_features
  out <- numeric(n_feat)
  if (is.null(col_idx)) {
    cols_all <- seq_len(handle$n_cells)
  } else {
    cols_all <- as.integer(col_idx)
    cols_all <- cols_all[!is.na(cols_all) & cols_all >= 1L & cols_all <= handle$n_cells]
  }
  if (length(cols_all) < 1L) {
    names(out) <- handle$feature_ids
    return(out)
  }
  runs <- .ic2_contig_runs(sort(unique(cols_all)))
  for (run in runs) {
    starts <- seq.int(1L, length(run), by = as.integer(chunk.size))
    for (s in starts) {
      e <- min(s + as.integer(chunk.size) - 1L, length(run))
      block <- .ic2_read_cols(handle, run[s:e])
      out <- out + as.numeric(Matrix::rowSums(block))
    }
  }
  names(out) <- handle$feature_ids
  out
}

#' Expand a dgCMatrix onto a larger feature universe (missing genes = 0)
#' @keywords internal
.ic2_expand_sparse_to_universe <- function(mat, feature_universe) {
  mat <- .ic2_to_dgCMatrix(mat)
  fu <- as.character(feature_universe)
  m <- match(rownames(mat), fu)
  if (anyNA(m)) {
    stop("Assay has features not present in the feature universe.")
  }
  mat@i <- as.integer(m[mat@i + 1L] - 1L)
  mat@Dim <- c(length(fu), ncol(mat))
  mat@Dimnames <- list(fu, colnames(mat))
  mat
}

#' Apply a simple expression transform to an assay block
#' @keywords internal
.ic2_assay_transform_block <- function(block, transform = c("none", "log2p1", "log1p")) {
  transform <- match.arg(transform)
  if (identical(transform, "none")) return(block)
  if (inherits(block, "Matrix")) {
    block <- .ic2_to_dgCMatrix(block)
    if (identical(transform, "log2p1")) {
      block@x <- log2(block@x + 1)
    } else {
      block@x <- log1p(block@x)
    }
    return(block)
  }
  block <- as.matrix(block)
  storage.mode(block) <- "double"
  if (identical(transform, "log2p1")) {
    block <- log2(block + 1)
  } else {
    block <- log1p(block)
  }
  block
}

#' Normalize one or more assay handles into a single CSC HDF5 (streamed)
#'
#' Columns are appended sample-by-sample. Features are aligned to
#' \code{feature_universe}. Does not densify the full matrix.
#'
#' @keywords internal
.ic2_normalize_merged_hdf5 <- function(
    handles,
    size_factors_list,
    dest_path,
    feature_universe,
    group = "counts",
    chunk.size = 200000L,
    round.num = FALSE,
    rounding.digits = 3,
    overwrite = TRUE,
    progress = NULL,
    transform = c("none", "log2p1", "log1p")
) {
  transform <- match.arg(transform)
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  if (length(handles) < 1L) stop("No assay handles to normalize.")
  if (length(size_factors_list) != length(handles)) {
    stop("size_factors_list length must match handles.")
  }
  feature_universe <- as.character(feature_universe)
  if (length(feature_universe) < 1L) stop("'feature_universe' is empty.")
  chunk.size <- as.integer(chunk.size)

  ## Enrich + plan columns to write (preserve handle / assay order)
  plans <- vector("list", length(handles))
  n_out <- 0L
  for (i in seq_along(handles)) {
    h <- .ic2_enrich_handle(handles[[i]], load_dimnames = TRUE)
    sf <- size_factors_list[[i]]
    keep <- names(sf)
    if (is.null(keep) || length(keep) < 1L) {
      stop("Size factors for handle ", i, " must be a named vector of cells to write.")
    }
    col_idx <- match(keep, h$cell_ids)
    if (anyNA(col_idx)) {
      stop("Size-factor cell ids missing from handle ", i, ".")
    }
    ord <- order(col_idx)
    col_idx <- col_idx[ord]
    keep <- keep[ord]
    sf <- as.numeric(sf[keep])
    if (any(!is.finite(sf) | sf <= 0)) {
      stop("Non-positive size factors in handle ", i, ".")
    }
    plans[[i]] <- list(handle = h, col_idx = col_idx, cell_ids = keep, sf = sf)
    n_out <- n_out + length(col_idx)
  }
  if (n_out < 1L) stop("No cells to write into merged normalized assay.")

  ## Read kept columns efficiently: when drops are sparse, read the covering
  ## contiguous span once and subset (avoids thousands of tiny H5 runs).
  .read_keep_cols <- function(handle, col_idx) {
    col_idx <- as.integer(col_idx)
    n <- length(col_idx)
    if (n < 1L) stop("Empty column index in normalize.")
    lo <- col_idx[[1L]]
    hi <- col_idx[[n]]
    span <- hi - lo + 1L
    n_drop <- span - n
    if (n_drop > 0L && (n / span >= 0.5 || n_drop <= max(64L, as.integer(0.02 * n)))) {
      raw <- .ic2_read_cols(handle, seq.int(lo, hi))
      return(raw[, col_idx - lo + 1L, drop = FALSE])
    }
    .ic2_read_cols(handle, col_idx)
  }

  ## Pass 1: count nnz (progress = cells over both passes)
  nnz_total <- 0
  cells_done <- 0L
  for (i in seq_along(plans)) {
    pl <- plans[[i]]
    starts <- seq.int(1L, length(pl$col_idx), by = chunk.size)
    for (s in starts) {
      e <- min(s + chunk.size - 1L, length(pl$col_idx))
      block <- .read_keep_cols(pl$handle, pl$col_idx[s:e])
      block <- .ic2_to_dgCMatrix(block)
      nnz_total <- nnz_total + length(block@x)
      cells_done <- cells_done + (e - s + 1L)
      if (is.function(progress)) progress(cells_done, n_out * 2L)
    }
  }
  nnz_store <- max(as.numeric(nnz_total), 1)

  if (file.exists(dest_path)) {
    if (!isTRUE(overwrite)) {
      stop("File already exists at '", dest_path, "'. Use overwrite = TRUE.")
    }
    unlink(dest_path)
  }
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)

  h5 <- hdf5r::H5File$new(dest_path, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5$create_group(group)
  g[["format"]] <- "iCellR2_csc_v1"
  g[["dims"]] <- as.integer(c(length(feature_universe), n_out))
  g[["rownames"]] <- feature_universe

  chunk_nx <- as.numeric(min(nnz_store, 5e6))
  if (chunk_nx < 1) chunk_nx <- 1
  ds_i <- g$create_dataset(
    "i", dtype = hdf5r::h5types$H5T_NATIVE_INT,
    dims = nnz_store, chunk_dims = chunk_nx
  )
  ds_x <- g$create_dataset(
    "x", dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = nnz_store, chunk_dims = chunk_nx
  )
  ds_p <- g$create_dataset(
    "p", dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = n_out + 1L,
    chunk_dims = as.integer(min(n_out + 1L, 1e6L))
  )
  .ic2_h5_write_1d(ds_p, start1 = 1, values = 0)

  out_nnz <- 0
  out_col <- 0L
  colnames_out <- character(n_out)

  for (i in seq_along(plans)) {
    pl <- plans[[i]]
    starts <- seq.int(1L, length(pl$col_idx), by = chunk.size)
    for (s in starts) {
      e <- min(s + chunk.size - 1L, length(pl$col_idx))
      cols <- pl$col_idx[s:e]
      block <- .read_keep_cols(pl$handle, cols)
      rn <- rownames(block)
      cn <- colnames(block)
      block <- block %*% Matrix::Diagonal(x = 1 / pl$sf[s:e])
      if (!identical(transform, "none")) {
        block <- .ic2_assay_transform_block(block, transform = transform)
      }
      if (isTRUE(round.num)) {
        block@x <- round(block@x, digits = rounding.digits)
      }
      block <- .ic2_to_dgCMatrix(block)
      rownames(block) <- rn
      colnames(block) <- cn
      block <- .ic2_expand_sparse_to_universe(block, feature_universe)

      bp <- as.numeric(block@p)
      bi <- as.integer(block@i)
      bx <- as.numeric(block@x)
      nc_b <- ncol(block)
      nadd <- length(bx)
      if (nadd > 0L) {
        .ic2_h5_write_1d(ds_i, start1 = out_nnz + 1, values = bi)
        .ic2_h5_write_1d(ds_x, start1 = out_nnz + 1, values = bx)
      }
      ## Column pointers for this chunk
      .ic2_h5_write_1d(
        ds_p,
        start1 = out_col + 2L,
        values = out_nnz + bp[-1L]
      )
      colnames_out[(out_col + 1L):(out_col + nc_b)] <- cn
      out_nnz <- out_nnz + nadd
      out_col <- out_col + nc_b
      if (is.function(progress)) progress(n_out + out_col, n_out * 2L)
    }
  }

  if (out_col != n_out) {
    stop("Internal error: wrote ", out_col, " columns, expected ", n_out, ".")
  }
  g[["colnames"]] <- colnames_out
  g[["nnz"]] <- as.numeric(out_nnz)
  h5$flush()
  h5$close_all()
  on.exit(NULL)

  .ic2_assay_handle(
    path = dest_path,
    backend = "hdf5_sparse",
    group = group,
    n_features = length(feature_universe),
    n_cells = n_out,
    feature_ids = NULL,
    cell_ids = NULL
  )
}

#' Materialize a feature x cell CSC HDF5 subset for fast PCA streaming
#'
#' Reads full assay columns once, keeps only \code{features} rows, writes a
#' skinny CSC suitable for repeated cov / project passes. Optionally applies
#' size factors + transform so callers can stream with transform = "none".
#'
#' @keywords internal
.ic2_materialize_feature_subset_hdf5 <- function(
    handles,
    features,
    cell_ids,
    dest_path,
    group = "counts",
    chunk.size = 200000L,
    size_factors = NULL,
    transform = c("none", "log2p1", "log1p"),
    overwrite = TRUE,
    progress = NULL
) {
  transform <- match.arg(transform)
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required.")
  }
  features <- unique(as.character(features))
  cell_ids <- as.character(cell_ids)
  if (length(features) < 1L || length(cell_ids) < 1L) {
    stop("feature subset materialize needs >= 1 feature and >= 1 cell.")
  }
  chunk.size <- as.integer(chunk.size)
  if (!is.finite(chunk.size) || chunk.size < 1L) chunk.size <- 200000L

  ## Plan columns across handles (preserve cell_ids order)
  remaining <- cell_ids
  plans <- list()
  for (h in handles) {
    if (length(remaining) < 1L) break
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    local <- intersect(remaining, h$cell_ids)
    if (length(local) < 1L) next
    row_idx <- match(features, h$feature_ids)
    if (all(is.na(row_idx))) next
    col_idx <- match(local, h$cell_ids)
    ord <- order(col_idx)
    local <- local[ord]
    col_idx <- col_idx[ord]
    plans[[length(plans) + 1L]] <- list(
      handle = h,
      col_idx = col_idx,
      cell_ids = local,
      row_idx = row_idx
    )
    remaining <- setdiff(remaining, local)
  }
  if (length(remaining) > 0L) {
    stop(
      "Could not locate ", length(remaining),
      " cells while building feature subset cache."
    )
  }
  n_out <- length(cell_ids)
  p <- length(features)

  ## Pass 1: nnz of feature-subsetted blocks
  nnz_total <- 0
  seen <- 0L
  for (pl in plans) {
    starts <- seq.int(1L, length(pl$col_idx), by = chunk.size)
    for (s in starts) {
      e <- min(s + chunk.size - 1L, length(pl$col_idx))
      raw <- .ic2_read_cols(pl$handle, pl$col_idx[s:e])
      present <- which(!is.na(pl$row_idx))
      if (length(present) < 1L) {
        seen <- seen + (e - s + 1L)
        next
      }
      sub <- raw[pl$row_idx[present], , drop = FALSE]
      sub <- .ic2_to_dgCMatrix(sub)
      ## nnz unchanged by SF / log transforms on existing nonzeros
      nnz_total <- nnz_total + length(sub@x)
      seen <- seen + ncol(sub)
      if (is.function(progress)) progress(seen, n_out * 2L)
    }
  }
  nnz_store <- max(as.numeric(nnz_total), 1)

  if (file.exists(dest_path)) {
    if (!isTRUE(overwrite)) stop("Destination exists: ", dest_path)
    unlink(dest_path)
  }
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)

  h5 <- hdf5r::H5File$new(dest_path, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- h5$create_group(group)
  g[["format"]] <- "iCellR2_csc_v1"
  g[["dims"]] <- as.integer(c(p, n_out))
  g[["rownames"]] <- features

  chunk_nx <- as.integer(min(max(nnz_store, 1), 5000000L))
  ds_i <- g$create_dataset(
    "i", dtype = hdf5r::h5types$H5T_NATIVE_INT,
    dims = nnz_store, chunk_dims = chunk_nx
  )
  ds_x <- g$create_dataset(
    "x", dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = nnz_store, chunk_dims = chunk_nx
  )
  ds_p <- g$create_dataset(
    "p", dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = n_out + 1L,
    chunk_dims = as.integer(min(n_out + 1L, 1e6L))
  )
  .ic2_h5_write_1d(ds_p, start1 = 1, values = 0)

  out_nnz <- 0
  out_col <- 0L
  colnames_out <- character(n_out)
  ## Re-order plan columns to match cell_ids (samples may interleave).
  all_col_handle <- integer(0)
  all_col_idx <- integer(0)
  all_row_idx <- list()
  all_cn <- character(0)
  for (pi in seq_along(plans)) {
    pl <- plans[[pi]]
    all_col_handle <- c(all_col_handle, rep(pi, length(pl$col_idx)))
    all_col_idx <- c(all_col_idx, pl$col_idx)
    all_cn <- c(all_cn, pl$cell_ids)
  }
  ord <- match(cell_ids, all_cn)
  if (anyNA(ord)) stop("Internal error: cell order mismatch in feature cache.")
  all_col_handle <- all_col_handle[ord]
  all_col_idx <- all_col_idx[ord]
  all_cn <- all_cn[ord]

  starts <- seq.int(1L, n_out, by = chunk.size)
  for (s in starts) {
    e <- min(s + chunk.size - 1L, n_out)
    ## Split chunk by handle to avoid mixing
    pos <- s
    while (pos <= e) {
      pi <- all_col_handle[pos]
      end_run <- pos
      while (end_run < e && all_col_handle[end_run + 1L] == pi) {
        end_run <- end_run + 1L
      }
      pl <- plans[[pi]]
      cols <- all_col_idx[pos:end_run]
      cn <- all_cn[pos:end_run]
      raw <- .ic2_read_cols(pl$handle, cols)
      present <- which(!is.na(pl$row_idx))
      block <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0), x = numeric(0),
        dims = c(p, length(cols)),
        dimnames = list(features, cn)
      )
      if (length(present) > 0L) {
        sub <- raw[pl$row_idx[present], , drop = FALSE]
        sub <- .ic2_to_dgCMatrix(sub)
        ## Place into full feature rows
        i_full <- as.integer(present[sub@i + 1L] - 1L)
        block <- Matrix::sparseMatrix(
          i = i_full + 1L,
          p = sub@p,
          x = as.numeric(sub@x),
          dims = c(p, length(cols)),
          dimnames = list(features, cn),
          index1 = TRUE
        )
      }
      block <- .ic2_scale_cols_by_sf(block, size_factors)
      block <- .ic2_assay_transform_block(block, transform = transform)
      block <- .ic2_to_dgCMatrix(block)
      bp <- as.numeric(block@p)
      bi <- as.integer(block@i)
      bx <- as.numeric(block@x)
      nc_b <- ncol(block)
      nadd <- length(bx)
      if (nadd > 0L) {
        .ic2_h5_write_1d(ds_i, start1 = out_nnz + 1, values = bi)
        .ic2_h5_write_1d(ds_x, start1 = out_nnz + 1, values = bx)
      }
      .ic2_h5_write_1d(
        ds_p,
        start1 = out_col + 2L,
        values = out_nnz + bp[-1L]
      )
      colnames_out[(out_col + 1L):(out_col + nc_b)] <- cn
      out_nnz <- out_nnz + nadd
      out_col <- out_col + nc_b
      pos <- end_run + 1L
    }
    if (is.function(progress)) progress(n_out + e, n_out * 2L)
  }

  if (out_col != n_out) {
    stop("Internal error: wrote ", out_col, " columns, expected ", n_out, ".")
  }
  ## Shrink datasets if nnz overestimated (zeros after transform rare)
  g[["colnames"]] <- colnames_out
  g[["nnz"]] <- as.numeric(out_nnz)
  h5$flush()
  h5$close_all()
  on.exit(NULL)

  .ic2_assay_handle(
    path = dest_path,
    backend = "hdf5_sparse",
    group = group,
    n_features = p,
    n_cells = n_out,
    feature_ids = NULL,
    cell_ids = NULL
  )
}

#' Normalize a handle to a new on-disk HDF5 assay (chunked; sparse throughout)
#'
#' @param keep_cell_ids Optional cell ids to write. If set, only those columns
#'   are normalized and written (raw source handle is not rewritten).
#' @keywords internal
.ic2_normalize_handle <- function(
    handle,
    dest_path,
    size_factors,
    group = "counts",
    chunk.size = 200000L,
    round.num = FALSE,
    rounding.digits = 3,
    keep_cell_ids = NULL,
    progress = NULL,
    transform = c("none", "log2p1", "log1p")
) {
  transform <- match.arg(transform)
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  if (is.null(names(size_factors))) {
    if (length(size_factors) != handle$n_cells) {
      stop("Unnamed size_factors length must equal n_cells in the assay handle.")
    }
    names(size_factors) <- handle$cell_ids
  }

  if (is.null(keep_cell_ids)) {
    keep_cell_ids <- handle$cell_ids
  } else {
    keep_cell_ids <- intersect(as.character(keep_cell_ids), handle$cell_ids)
  }
  if (length(keep_cell_ids) < 1L) {
    stop("No cells to write for normalize handle '", handle$path, "'.")
  }
  sf <- as.numeric(size_factors[keep_cell_ids])
  names(sf) <- keep_cell_ids
  if (any(!is.finite(sf) | sf <= 0)) {
    stop("All size factors for kept cells must be finite and > 0.")
  }

  ## Stream write via merged writer (single handle; no gene-universe expand)
  .ic2_normalize_merged_hdf5(
    handles = list(handle),
    size_factors_list = list(sf),
    dest_path = dest_path,
    feature_universe = as.character(handle$feature_ids),
    group = group,
    chunk.size = chunk.size,
    round.num = round.num,
    rounding.digits = rounding.digits,
    overwrite = TRUE,
    progress = progress,
    transform = transform
  )
}

#' Replace/update an assay slot with handle(s); update storage paths
#' @keywords internal
.ic2_assign_assay <- function(
    x,
    modality,
    assay.name,
    handles,
    virtual = FALSE,
    feature_universe = NULL,
    set.active = TRUE
) {
  if (isTRUE(virtual)) {
    sample_ids <- vapply(handles, function(h) {
      if (!is.null(h$sample_id) && nzchar(h$sample_id)) h$sample_id else NA_character_
    }, character(1))
    assay_ref <- list(
      aggregation_type = "virtual",
      sample_handles = handles,
      sample_assay_paths = lapply(handles, function(h) h$path),
      feature_universe = feature_universe
    )
    x@storage$assay_paths <- list(
      modality = modality,
      assay = assay.name,
      aggregation = "virtual",
      sample_paths = lapply(handles, function(h) h$path)
    )
  } else {
    h <- handles[[1L]]
    assay_ref <- list(
      aggregation_type = "single",
      handle = h
    )
    x@storage$assay_paths <- list(
      modality = modality,
      assay = assay.name,
      filepath = h$path,
      backend = h$backend,
      group = h$group
    )
  }
  x@modalities[[modality]]$assays[[assay.name]] <- assay_ref
  if (isTRUE(set.active)) {
    x@modalities[[modality]]$active_assay <- assay.name
  }
  x
}

#' Persist cell_data to SODA if a path is registered
#' @keywords internal
.ic2_persist_cell_data <- function(x) {
  path <- x@storage$metadata_paths$cell_data
  if (is.null(path) || !nzchar(as.character(path)[1L])) {
    root <- x@storage$root_dir
    if (!is.null(root) && nzchar(root)) {
      path <- file.path(root, "metadata", "cell_data")
    } else {
      return(invisible(x))
    }
  }
  stem <- sub("\\.(rds|parquet|tsv)$", "", path, ignore.case = TRUE)
  new_path <- .ic2_soda_write_table(x@cell_data, stem)
  x@storage$metadata_paths$cell_data <- new_path
  invisible(x)
}

#' Append a command log entry
#' @keywords internal
.ic2_log_command <- function(x, name, arguments, outputs, input_dir = NULL, output_dir = NULL) {
  if (is.null(x@commands$command_log)) {
    x@commands$command_log <- list()
  }
  ts <- as.character(Sys.time())
  soda_root <- tryCatch(x@storage$root_dir, error = function(e) NULL)

  ## Resolve I/O dirs for the human-readable log
  if (is.null(input_dir) && !is.null(arguments)) {
    for (key in c("input_dir", "input", "in.dir")) {
      if (!is.null(arguments[[key]]) && is.character(arguments[[key]]) &&
          length(arguments[[key]]) == 1L && nzchar(arguments[[key]])) {
        input_dir <- arguments[[key]]
        break
      }
    }
  }
  if (is.null(output_dir) && !is.null(arguments)) {
    for (key in c("output_dir", "out.dir", "h5.file", "dest")) {
      if (!is.null(arguments[[key]]) && is.character(arguments[[key]]) &&
          length(arguments[[key]]) == 1L && nzchar(arguments[[key]])) {
        output_dir <- arguments[[key]]
        break
      }
    }
  }
  if (is.null(output_dir) || !nzchar(as.character(output_dir)[1L])) {
    output_dir <- soda_root
  }
  if (is.null(input_dir) || !nzchar(as.character(input_dir)[1L])) {
    input_dir <- soda_root
  }
  ## Prefer directories when a file path was passed
  input_dir_disp <- as.character(input_dir)[1L]
  output_dir_disp <- as.character(output_dir)[1L]
  if (nzchar(input_dir_disp) && file.exists(input_dir_disp) && !dir.exists(input_dir_disp)) {
    input_dir_disp <- dirname(input_dir_disp)
  }
  if (nzchar(output_dir_disp) && file.exists(output_dir_disp) && !dir.exists(output_dir_disp)) {
    output_dir_disp <- dirname(output_dir_disp)
  }
  if (nzchar(input_dir_disp)) {
    input_dir_disp <- normalizePath(input_dir_disp, winslash = "/", mustWork = FALSE)
  }
  if (nzchar(output_dir_disp)) {
    output_dir_disp <- normalizePath(output_dir_disp, winslash = "/", mustWork = FALSE)
  }

  key <- paste0(name, "_", length(x@commands$command_log) + 1L)
  x@commands$command_log[[key]] <- list(
    function_name = name,
    timestamp = ts,
    input_dir = input_dir_disp,
    output_dir = output_dir_disp,
    arguments = arguments,
    outputs = outputs
  )
  x@validity$checked <- ts

  ## Mirror into human-readable SODA log (commands.rds unchanged)
  if (!is.null(soda_root) && nzchar(as.character(soda_root)[1L])) {
    out_txt <- if (is.null(outputs) || length(outputs) == 0L) {
      "  (none)"
    } else {
      .ic2_fmt_log_args(as.list(outputs))
    }
    .ic2_append_soda_log(soda_root, c(
      paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] COMMAND  ", name),
      paste0("  input_dir:  ", input_dir_disp %||% " - "),
      paste0("  output_dir: ", output_dir_disp %||% " - "),
      "  arguments:",
      .ic2_fmt_log_args(arguments),
      "  outputs:",
      out_txt,
      ""
    ))
  }
  x
}

#' Accumulate poscounts log-sum and n_pos per gene (DESeq-style)
#' @keywords internal
.ic2_handle_deseq_poscounts_accum <- function(handle, chunk.size = 200000L, col_idx = NULL) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  n_feat <- handle$n_features
  sum_log <- numeric(n_feat)
  n_pos <- integer(n_feat)
  if (is.null(col_idx)) {
    cols_all <- seq_len(handle$n_cells)
  } else {
    cols_all <- as.integer(col_idx)
    cols_all <- cols_all[!is.na(cols_all) & cols_all >= 1L & cols_all <= handle$n_cells]
  }
  if (length(cols_all) < 1L) {
    names(sum_log) <- handle$feature_ids
    names(n_pos) <- handle$feature_ids
    return(list(sum_log = sum_log, n_pos = n_pos))
  }
  runs <- .ic2_contig_runs(sort(unique(cols_all)))
  for (run in runs) {
    starts <- seq.int(1L, length(run), by = as.integer(chunk.size))
    for (s in starts) {
      e <- min(s + as.integer(chunk.size) - 1L, length(run))
      block <- .ic2_to_dgCMatrix(.ic2_read_cols(handle, run[s:e]))
      i <- block@i + 1L
      x <- as.numeric(block@x)
      keep <- is.finite(x) & x > 0
      if (!any(keep)) next
      i <- i[keep]
      lx <- log(x[keep])
      ## Accumulate without densifying
      for (k in seq_along(i)) {
        gi <- i[k]
        sum_log[gi] <- sum_log[gi] + lx[k]
      }
      n_pos <- n_pos + tabulate(i, nbins = n_feat)
    }
  }
  names(sum_log) <- handle$feature_ids
  names(n_pos) <- handle$feature_ids
  list(sum_log = sum_log, n_pos = n_pos)
}

#' DESeq median-of-ratios size factors given gene log-geometric means
#' @keywords internal
.ic2_handle_deseq_size_factors <- function(
    handle,
    loggeomeans,
    chunk.size = 200000L
) {
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  n_cells <- handle$n_cells
  sf <- rep(NA_real_, n_cells)
  loggeomeans <- as.numeric(loggeomeans)
  if (length(loggeomeans) != handle$n_features) {
    stop("loggeomeans length must equal n_features.")
  }
  finite_gene <- is.finite(loggeomeans)
  if (!any(finite_gene)) {
    return(sf)
  }
  starts <- seq.int(1L, n_cells, by = as.integer(chunk.size))
  for (start in starts) {
    end <- min(start + as.integer(chunk.size) - 1L, n_cells)
    cols <- start:end
    block <- .ic2_to_dgCMatrix(.ic2_read_cols(handle, cols))
    ## Per column: median(count / geoMean) over genes with count > 0 & finite geoMean
    p <- nrow(block)
    for (j in seq_along(cols)) {
      ## Extract sparse column j (1-based in block)
      p0 <- block@p[j] + 1L
      p1 <- block@p[j + 1L]
      if (p1 < p0) {
        sf[cols[j]] <- NA_real_
        next
      }
      rows <- block@i[p0:p1] + 1L
      vals <- as.numeric(block@x[p0:p1])
      ok <- vals > 0 & is.finite(vals) & finite_gene[rows]
      if (!any(ok)) {
        sf[cols[j]] <- NA_real_
        next
      }
      ratios <- vals[ok] / exp(loggeomeans[rows[ok]])
      sf[cols[j]] <- stats::median(ratios)
    }
  }
  sf
}
