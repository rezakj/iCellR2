#' Add a scATAC peak matrix on its own modality
#'
#' Reads a Cell Ranger ATAC or ARC peak matrix (\code{filtered_peak_bc_matrix.h5},
#' \code{filtered_feature_bc_matrix.h5}, or a 10x MTX directory) and writes
#' \code{modalities/atac/assays/counts}. The RNA assay and \code{reductions/pca}
#' are not modified. \code{fragments.tsv.gz} is recorded and not loaded.
#'
#' With \code{x = NULL} this creates an ATAC-only SODA archive. No RNA assay
#' is required. Counts stay on disk, the same way as \code{\link{load10x}}.
#' Two or more samples whose peak names differ are streamed onto one shared
#' set from \code{\link{merge.peaks}}. Pass \code{peaks} to use a BED merge
#' you already built. \code{fragments.tsv.gz} is recorded and not loaded.
#'
#' @param x An iCellR2 object, or \code{NULL} to start from ATAC alone (no RNA).
#' @param input One or more ATAC output directories, 10x HDF5 files, or MTX directories.
#'   The same formats as \code{\link{load10x}}. The matrix stays on disk.
#' @param out.dir SODA stem used only when \code{x} is \code{NULL}. The archive
#'   is \code{paste0(out.dir, ".soda")}.
#' @param overwrite Replace an existing SODA archive when \code{x} is \code{NULL}.
#' @param sample.ids Sample ids, same length as \code{input}. Default: folder name.
#'   Must match \code{cell_data$sample_id} when the object has more than one sample.
#' @param peaks Optional consensus peaks from \code{\link{merge.peaks}} (the
#'   data.frame, or the TSV written by the CLI). \code{NULL} merges the peak
#'   names in \code{input} when there are two or more samples.
#' @param feature.type Feature-type filter for a multiome H5. Default
#'   \code{"Peaks"}. Peak-only files that have no feature-type column are read whole.
#' @param assay.name Assay to write. Default \code{"counts"}.
#' @param verbose Progress style.
#' @return The iCellR2 object.
#'
#' @examples
#' \dontrun{
#' obj <- add.atac(input = "atac_outs", sample.ids = "donor1", out.dir = "atac_only")
#' obj <- add.atac(obj, input = "atac_outs", sample.ids = "donor1")
#'
#' shared <- merge.peaks(list(
#'   donor1 = "donor1_peaks.bed",
#'   donor2 = "donor2_peaks.bed"
#' ))
#' obj <- add.atac(
#'   obj,
#'   input = c("donor1_outs", "donor2_outs"),
#'   sample.ids = c("donor1", "donor2"),
#'   peaks = shared
#' )
#' }
#'
#' @export
add.atac <- function(
    x = NULL,
    input = NULL,
    sample.ids = NULL,
    peaks = NULL,
    feature.type = "Peaks",
    assay.name = "counts",
    out.dir = "iCellR2",
    overwrite = FALSE,
    verbose = TRUE
) {
  if (!is.null(x) && !inherits(x, "iCellR2")) {
    stop("x should be an iCellR2 object, or NULL to start from ATAC alone.")
  }
  if (is.null(input) || !length(input)) stop("Provide input (an ATAC directory, H5, or MTX folder).")
  input <- as.character(input)
  rep <- .ic2_new_reporter("add.atac", verbose)
  if (is.null(sample.ids)) {
    sample.ids <- vapply(input, .ic2_default_sample_id, character(1))
  }
  sample.ids <- as.character(sample.ids)
  if (length(sample.ids) != length(input)) {
    stop("'sample.ids' must have one id per input.")
  }

  consensus <- .ic2_atac_read_peaks_arg(peaks)
  samps <- vector("list", length(input))
  frag <- vector("list", length(input))
  for (i in seq_along(input)) {
    .ic2_step(rep, paste0("Attaching ", sample.ids[[i]], " (on disk)"))
    samps[[i]] <- .ic2_atac_virtual_sample(
      input[[i]], sample.ids[[i]], feature.type, verbose
    )
    frag[[i]] <- list(
      sample_id = sample.ids[[i]],
      fragments = samps[[i]]$fragments
    )
  }
  if (is.null(consensus) && length(samps) > 1L && !.ic2_atac_same_features(samps)) {
    consensus <- merge.peaks(.ic2_atac_feature_table(samps))
    .ic2_step(rep, paste0("Merged peaks: ", nrow(consensus), " consensus intervals"))
  }
  if (!is.null(consensus)) {
    object_cells <- if (is.null(x)) NULL else .ic2_adt_cell_ids(x)
    samps <- lapply(samps, function(s) {
      .ic2_atac_stream_consensus(s, consensus, object_cells, verbose)
    })
  }

  if (is.null(x)) {
    x <- .ic2_atac_make_object(samps, out.dir, overwrite, verbose)
  } else {
    x <- .ic2_atac_register_virtual(x, samps, assay.name, overwrite)
  }

  feat <- .ic2_enrich_handle(samps[[1L]]$handle, load_dimnames = TRUE)$feature_ids
  parsed <- .ic2_parse_peak_id(feat)
  fd <- data.frame(
    feature_id = feat,
    feature_name = feat,
    feature_type = "Peaks",
    chromosome = parsed$chromosome,
    start = parsed$start,
    end = parsed$end,
    stringsAsFactors = FALSE
  )
  root <- x@storage$root_dir
  fd_path <- .ic2_soda_write_table(fd, file.path(root, "metadata", "atac_feature_data"))
  frag_path <- file.path(root, "metadata", "atac_fragments.json")
  .ic2_soda_write_json(frag, frag_path)
  x <- .ic2_ensure_modality(x, "atac")
  x@modalities$atac$feature_data <- fd
  x@storage$metadata_paths$atac_feature_data <- fd_path
  x@storage$modality_paths$atac <- file.path(root, "modalities", "atac")
  x@params$atac <- list(
    assay.name = assay.name,
    n_features = length(feat),
    samples = sample.ids,
    fragments = frag,
    on_disk = TRUE
  )
  x <- .ic2_log_command(
    x, "add.atac",
    list(samples = sample.ids, n_features = length(feat), assay.name = assay.name),
    list(feature_data = fd_path, fragments = frag_path)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  .ic2_finish(rep, extra_lines = paste0(
    "Peaks: ", length(feat),
    "   |   samples: ", paste(sample.ids, collapse = ", "),
    if (is.null(consensus)) "   |   on-disk 10x handle" else "   |   streamed onto the shared peak set"
  ))
  x
}

#' Load one ATAC matrix the way \code{load10x} does: a handle, not a dense matrix.
#' @keywords internal
.ic2_atac_virtual_sample <- function(input, sample.id, feature.type, verbose) {
  found <- .ic2_atac_find(input)
  meta <- file.path(
    "iCellR2_import",
    paste0(sample.id, "_atac_tenx_meta.json")
  )
  samp <- tryCatch(
    load10x(
      found$matrix,
      sample.id = sample.id,
      feature.type = feature.type,
      output = "virtual",
      modality = "atac",
      h5.file = meta,
      overwrite = TRUE,
      verbose = verbose
    ),
    error = function(e) e
  )
  if (inherits(samp, "error") && !is.null(feature.type)) {
    samp <- load10x(
      found$matrix,
      sample.id = sample.id,
      feature.type = NULL,
      output = "virtual",
      modality = "atac",
      h5.file = meta,
      overwrite = TRUE,
      verbose = verbose
    )
  } else if (inherits(samp, "error")) {
    stop(samp$message)
  }
  samp$fragments <- found$fragments
  samp
}

#' @keywords internal
.ic2_atac_sample_features <- function(samp) {
  h <- if (!is.null(samp$handle)) samp$handle else samp
  h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
  as.character(h$feature_ids)
}

#' @keywords internal
.ic2_atac_same_features <- function(samps) {
  ids <- lapply(samps, .ic2_atac_sample_features)
  all(vapply(ids[-1L], function(z) identical(z, ids[[1L]]), logical(1)))
}

#' @keywords internal
.ic2_atac_feature_table <- function(samps) {
  tabs <- lapply(samps, function(s) {
    tab <- .ic2_parse_peak_id(.ic2_atac_sample_features(s))
    tab$sample <- s$sample_id %||% s$handle$sample_id
    tab[!is.na(tab$start), , drop = FALSE]
  })
  do.call(rbind, tabs)
}

#' @keywords internal
.ic2_atac_make_object <- function(samps, out.dir, overwrite, verbose) {
  if (length(samps) == 1L) {
    return(make.obj(
      samps[[1L]], out.dir = out.dir, modality = "atac",
      overwrite = overwrite, verbose = verbose
    ))
  }
  agg <- data.aggregation(
    samples = samps, modality = "atac", output = "virtual", verbose = verbose
  )
  make.obj(
    agg, out.dir = out.dir, modality = "atac",
    overwrite = overwrite, verbose = verbose
  )
}

#' Point an existing object at ATAC handles. Does not read the counts.
#' @keywords internal
.ic2_atac_register_virtual <- function(x, samps, assay.name, overwrite) {
  x <- .ic2_ensure_modality(x, "atac")
  object_cells <- .ic2_adt_cell_ids(x)
  handles <- lapply(samps, function(s) {
    h <- if (!is.null(s$handle)) s$handle else s
    sid <- s$sample_id %||% h$sample_id
    .ic2_atac_place_handle(x, h, sid, assay.name, object_cells, overwrite)
  })
  x <- .ic2_assign_assay(
    x, modality = "atac", assay.name = assay.name,
    handles = handles, virtual = length(handles) > 1L, set.active = TRUE
  )
  x
}

#' Copy a 10x or CSC handle into the SODA tree and align cell ids.
#' @keywords internal
.ic2_atac_place_handle <- function(x, handle, sample_id, assay.name, object_cells, overwrite) {
  dest <- .ic2_soda_assay_path(x, "atac", assay.name, sample_id = sample_id, ext = ".h5")
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  cells <- .ic2_atac_align_barcodes(handle$cell_ids, sample_id, object_cells)
  if (identical(handle$backend, "tenx_h5")) {
    meta <- .ic2_soda_read_tenx_meta(handle$meta_path)
    if (is.null(meta)) stop("ATAC handle is missing tenx meta.")
    n_src <- as.integer(meta$n_features_src %||% 0L)
    n_keep <- length(meta$row_keep %||% integer(0))
    mixed <- !is.null(meta$feature.type) && n_src > 0L && n_keep > 0L && n_keep < n_src
    if (isTRUE(mixed)) {
      ## Peaks get their own HDF5. The RNA file is not reused.
      if (length(cells) != length(meta$cell_ids %||% cells)) cells <- meta$cell_ids
      handle <- .ic2_split_10x_feature_type(
        src_path = handle$path,
        dest_path = dest,
        feature.type = meta$feature.type %||% "Peaks",
        cell.ids = cells,
        sample_id = sample_id,
        modality = "atac"
      )
      handle$cell_ids <- as.character(cells)
    } else {
      if (!identical(normalizePath(handle$path, mustWork = FALSE),
                     normalizePath(dest, mustWork = FALSE))) {
        if (file.exists(dest) && !isTRUE(overwrite)) {
          stop("ATAC file already exists at ", dest, ". Use overwrite = TRUE.")
        }
        if (file.exists(dest)) unlink(dest)
        if (!file.copy(handle$path, dest, overwrite = TRUE)) {
          stop("Failed to copy ATAC matrix to ", dest)
        }
      }
      meta_dest <- file.path(dirname(dest), "tenx_meta.json")
      meta$cell_ids <- cells
      meta$source <- normalizePath(dest, winslash = "/", mustWork = TRUE)
      meta_dest <- .ic2_soda_write_tenx_meta(meta, meta_dest)
      handle$path <- normalizePath(dest, winslash = "/", mustWork = TRUE)
      handle$meta_path <- meta_dest
      handle$cell_ids <- cells
    }
  } else {
    if (!identical(normalizePath(handle$path, mustWork = FALSE),
                   normalizePath(dest, mustWork = FALSE))) {
      if (file.exists(dest)) unlink(dest)
      if (!file.copy(handle$path, dest, overwrite = TRUE)) {
        stop("Failed to copy ATAC matrix to ", dest)
      }
      handle$path <- normalizePath(dest, winslash = "/", mustWork = TRUE)
    }
  }
  handle$cell_ids <- cells
  handle$sample_id <- sample_id
  handle$modality <- "atac"
  handle$n_cells <- length(cells)
  handle
}

#' Prefer the barcode form already used in cell_data.
#' @keywords internal
.ic2_atac_align_barcodes <- function(barcodes, sample_id, object_cells) {
  barcodes <- .ic2_clean_names(as.character(barcodes))
  pref <- paste(sample_id, barcodes, sep = "_")
  if (is.null(object_cells) || !length(object_cells)) return(barcodes)
  n_pref <- sum(pref %in% object_cells)
  n_raw <- sum(barcodes %in% object_cells)
  if (n_pref > n_raw) pref else barcodes
}

#' Stream one sample onto the consensus peaks. Counts that share an interval are summed.
#' @keywords internal
.ic2_atac_stream_consensus <- function(samp, consensus, object_cells, verbose) {
  h <- if (!is.null(samp$handle)) samp$handle else samp
  sid <- samp$sample_id %||% h$sample_id
  meta <- .ic2_soda_read_tenx_meta(h$meta_path)
  if (is.null(meta)) stop("Cannot remap peaks without a 10x HDF5 handle: ", sid)
  src_n <- as.integer(meta$n_features_src %||% meta$n_features)
  features <- as.character(meta$feature_ids)
  mapped <- .ic2_atac_map_rows(features, consensus)
  row.map <- rep(NA_integer_, src_n)
  if (!is.null(meta$row_keep) && length(meta$row_keep)) {
    row.map[as.integer(meta$row_keep)] <- ifelse(is.na(mapped), NA_integer_, as.integer(mapped) - 1L)
  } else {
    if (length(mapped) != src_n) {
      stop("Peak map length does not match the ATAC matrix for ", sid, ".")
    }
    row.map <- ifelse(is.na(mapped), NA_integer_, as.integer(mapped) - 1L)
  }
  cells <- .ic2_atac_align_barcodes(meta$cell_ids, sid, object_cells)
  dest <- tempfile(pattern = paste0("atac_", sid, "_"), fileext = ".h5")
  handle <- .ic2_convert_10x_h5_to_csc_hdf5(
    src_path = h$path,
    dest_path = dest,
    feature.type = NULL,
    overwrite = TRUE,
    row.map = row.map,
    dest.features = as.character(consensus$peak_id),
    cell.ids = cells
  )
  handle$sample_id <- sid
  handle$modality <- "atac"
  samp$handle <- handle
  samp
}

#' @keywords internal
.ic2_ensure_modality <- function(x, modality) {
  cur <- x@modalities[[modality]]
  if (is.null(cur) || !is.list(cur$assays)) {
    x@modalities[[modality]] <- list(
      assays = list(), active_assay = NULL, feature_data = NULL,
      reductions = list(), active_reduction = NULL, graphs = list(),
      sketches = list(), markers = list(), models = list(),
      params = list(), misc = list()
    )
  }
  x
}

#' @keywords internal
.ic2_atac_find <- function(input) {
  input <- as.character(input)[1L]
  if (!file.exists(input) && !dir.exists(input)) {
    stop("ATAC input not found: ", input)
  }
  matrix <- NULL
  root <- input
  if (!dir.exists(input)) {
    matrix <- input
    root <- dirname(input)
  } else {
    rel <- c(
      "filtered_peak_bc_matrix.h5",
      "filtered_feature_bc_matrix.h5",
      "raw_peak_bc_matrix.h5",
      file.path("outs", "filtered_peak_bc_matrix.h5"),
      file.path("outs", "filtered_feature_bc_matrix.h5"),
      "filtered_peak_bc_matrix",
      "filtered_feature_bc_matrix"
    )
    for (r in rel) {
      p <- file.path(input, r)
      if (file.exists(p) || dir.exists(p)) {
        matrix <- p
        break
      }
    }
    if (is.null(matrix)) {
      stop(
        "No peak matrix under ", input,
        ". Expected filtered_peak_bc_matrix.h5 or filtered_feature_bc_matrix.h5."
      )
    }
  }
  frag <- NULL
  peaks <- NULL
  for (dir in unique(c(root, dirname(matrix), input))) {
    if (!dir.exists(dir)) next
    f <- file.path(dir, "fragments.tsv.gz")
    if (file.exists(f)) frag <- normalizePath(f, winslash = "/", mustWork = FALSE)
    b <- file.path(dir, "peaks.bed")
    if (file.exists(b)) peaks <- normalizePath(b, winslash = "/", mustWork = FALSE)
  }
  list(matrix = matrix, fragments = frag, peaks = peaks)
}

#' @keywords internal
.ic2_atac_read_matrix <- function(path, feature.type) {
  read_once <- function(ftype) {
    load10x(
      path,
      sample.id = "atac",
      feature.type = ftype,
      output = "sparse",
      modality = "atac",
      verbose = FALSE
    )
  }
  mat <- tryCatch(read_once(feature.type), error = function(e) e)
  if (inherits(mat, "error") && !is.null(feature.type)) {
    mat <- tryCatch(read_once(NULL), error = function(e) e)
  }
  if (inherits(mat, "error")) stop(mat$message)
  if (!inherits(mat, "Matrix") && !is.matrix(mat)) {
    stop("ATAC matrix from ", path, " could not be read.")
  }
  mat
}

#' @keywords internal
.ic2_atac_read_peaks_arg <- function(peaks) {
  if (is.null(peaks)) return(NULL)
  if (is.character(peaks) && length(peaks) == 1L && file.exists(peaks)) {
    tab <- utils::read.delim(peaks, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    map_path <- sub("(\\.[^.]+)$", "_map\\1", peaks)
    if (file.exists(map_path)) {
      attr(tab, "map") <- utils::read.delim(map_path, sep = "\t", stringsAsFactors = FALSE)
    }
    peaks <- tab
  }
  if (!is.data.frame(peaks)) stop("'peaks' must be a merge.peaks() table or a TSV path.")
  need <- c("chromosome", "start", "end", "peak_id")
  if (!all(need %in% names(peaks))) {
    stop("'peaks' needs columns chromosome, start, end, and peak_id.")
  }
  peaks
}

#' @keywords internal
.ic2_atac_peaks_from_names <- function(ids) {
  parsed <- .ic2_parse_peak_id(ids)
  ok <- !is.na(parsed$start)
  if (!any(ok)) {
    return(data.frame(
      chromosome = rep(NA_character_, length(ids)),
      start = NA_integer_, end = NA_integer_,
      peak_id = as.character(ids),
      stringsAsFactors = FALSE
    ))
  }
  parsed$peak_id <- ifelse(ok, .ic2_peak_ids(parsed$chromosome, parsed$start, parsed$end), ids)
  parsed <- parsed[ok, , drop = FALSE]
  key <- paste(parsed$chromosome, parsed$start, parsed$end, sep = "|")
  parsed <- parsed[!duplicated(key), , drop = FALSE]
  parsed
}

#' Integer row in the consensus table for each sample peak. NA if none.
#' @keywords internal
.ic2_atac_map_rows <- function(peak_ids, consensus) {
  parsed <- .ic2_parse_peak_id(peak_ids)
  n <- length(peak_ids)
  out <- rep(NA_integer_, n)
  if (!nrow(consensus)) return(out)
  cons <- consensus
  cons$idx <- seq_len(nrow(cons))
  exact <- match(
    paste(parsed$chromosome, parsed$start, parsed$end, sep = "|"),
    paste(cons$chromosome, cons$start, cons$end, sep = "|")
  )
  out[!is.na(exact)] <- exact[!is.na(exact)]
  need <- which(is.na(out) & !is.na(parsed$start))
  if (!length(need)) {
    ## Unparsed names: match peak_id text
    miss <- which(is.na(out))
    if (length(miss)) out[miss] <- match(peak_ids[miss], cons$peak_id)
    return(out)
  }
  bin <- 100000L
  for (chr in unique(parsed$chromosome[need])) {
    pr <- need[parsed$chromosome[need] == chr]
    cr <- cons[cons$chromosome == chr, , drop = FALSE]
    if (!length(pr) || !nrow(cr)) next
    width <- max(1L, max(cr$end - cr$start))
    buckets <- split(seq_len(nrow(cr)), floor(cr$start / bin))
    for (ii in pr) {
      gs <- parsed$start[[ii]]
      ge <- parsed$end[[ii]]
      b1 <- floor((gs - width) / bin)
      b2 <- floor(ge / bin)
      cand <- integer(0)
      for (bb in seq.int(b1, b2)) {
        hit <- buckets[[as.character(bb)]]
        if (!is.null(hit)) cand <- c(cand, hit)
      }
      if (!length(cand)) next
      ov <- pmin(ge, cr$end[cand]) - pmax(gs, cr$start[cand]) + 1L
      ok <- which(ov > 0L)
      if (!length(ok)) next
      out[[ii]] <- cr$idx[[cand[ok][which.max(ov[ok])]]]
    }
  }
  miss <- which(is.na(out))
  if (length(miss)) {
    hit <- match(peak_ids[miss], cons$peak_id)
    out[miss][!is.na(hit)] <- hit[!is.na(hit)]
  }
  out
}

#' Column index in cell_data for each ATAC barcode. NA if unmatched.
#' @keywords internal
.ic2_atac_match_cells <- function(barcodes, cd, sample_id, cell_ids) {
  dot <- function(z) gsub("-", ".", as.character(z), fixed = TRUE)
  n <- length(cell_ids)
  sids <- if (!is.null(cd) && "sample_id" %in% colnames(cd)) {
    as.character(cd$sample_id)
  } else {
    rep(NA_character_, n)
  }
  bars <- if (!is.null(cd) && "barcode" %in% colnames(cd)) as.character(cd$barcode) else cell_ids
  use <- if (any(sids == sample_id, na.rm = TRUE)) which(sids == sample_id) else seq_len(n)
  keys <- c(
    dot(bars[use]),
    dot(cell_ids[use]),
    dot(paste(sample_id, bars[use], sep = "_"))
  )
  cols <- rep(use, 3L)
  ok <- !is.na(keys) & nzchar(keys)
  keys <- keys[ok]
  cols <- cols[ok]
  first <- !duplicated(keys)
  lut <- setNames(cols[first], keys[first])
  hit <- match(dot(barcodes), names(lut))
  out <- rep(NA_integer_, length(barcodes))
  out[!is.na(hit)] <- lut[hit[!is.na(hit)]]
  out
}
