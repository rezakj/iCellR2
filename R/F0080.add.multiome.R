#' Add RNA and ATAC from the same cells
#'
#' This is RNA + ATAC only. Gene expression and peaks stay in separate
#' modalities. Neither matrix is loaded into memory: each input is attached
#' the way \code{\link{load10x}} attaches a 10x directory, 10x HDF5 file, or
#' MTX folder, including multi-million-cell runs.
#'
#' \code{rna.input} and \code{atac.input} are paired. One sample is one RNA
#' file and one ATAC file. A Cell Ranger ARC folder can be passed as both
#' when the gene and peak matrices live in that same folder. With two or
#' more samples, pass \code{peaks} from \code{\link{merge.peaks}}, or leave
#' it empty and differing peak sets are merged from their names, then
#' streamed onto that shared set.
#'
#' @param x An iCellR2 object, or \code{NULL} to create one from \code{rna.input}.
#' @param rna.input RNA directories or count files. Same formats as \code{\link{load10x}}.
#' @param atac.input ATAC directories or peak matrices, one per RNA input.
#' @param input Optional. Used as both \code{rna.input} and \code{atac.input}
#'   when a folder contains both matrices.
#' @param sample.ids Sample ids, one per pair.
#' @param sample.id Alias of a single \code{sample.ids} value.
#' @param peaks Optional consensus peaks from \code{\link{merge.peaks}}.
#' @param out.dir SODA stem used only when \code{x} is \code{NULL}.
#' @param overwrite Replace an existing SODA archive when creating \code{x}.
#' @param verbose Progress style.
#' @return The iCellR2 object.
#'
#' @examples
#' \dontrun{
#' obj <- add.multiome(
#'   rna.input = "donor1_rna",
#'   atac.input = "donor1_atac",
#'   sample.id = "donor1",
#'   out.dir = "multiome"
#' )
#'
#' ## One ARC folder holds both matrices.
#' obj <- add.multiome(
#'   input = "donor1_arc",
#'   sample.id = "donor1",
#'   out.dir = "multiome"
#' )
#'
#' shared <- merge.peaks(list(
#'   donor1 = "donor1_peaks.bed",
#'   donor2 = "donor2_peaks.bed"
#' ))
#' obj <- add.multiome(
#'   rna.input = c("donor1_rna", "donor2_rna"),
#'   atac.input = c("donor1_atac", "donor2_atac"),
#'   sample.ids = c("donor1", "donor2"),
#'   peaks = shared,
#'   out.dir = "multiome"
#' )
#' }
#'
#' @export
add.multiome <- function(
    x = NULL,
    rna.input = NULL,
    atac.input = NULL,
    input = NULL,
    sample.ids = NULL,
    sample.id = NULL,
    peaks = NULL,
    out.dir = "iCellR2",
    overwrite = FALSE,
    verbose = TRUE
) {
  if (is.null(rna.input) && is.null(atac.input) && !is.null(input)) {
    rna.input <- input
    atac.input <- input
  }
  if (is.null(rna.input) || !length(rna.input) || is.null(atac.input) || !length(atac.input)) {
    stop("Provide rna.input and atac.input (or input when both matrices are in one folder).")
  }
  rna.input <- as.character(rna.input)
  atac.input <- as.character(atac.input)
  if (length(rna.input) != length(atac.input)) {
    stop("'rna.input' and 'atac.input' must be paired and the same length.")
  }
  if (is.null(sample.ids) && !is.null(sample.id)) sample.ids <- sample.id
  if (is.null(sample.ids)) {
    sample.ids <- vapply(rna.input, .ic2_default_sample_id, character(1))
  }
  sample.ids <- as.character(sample.ids)
  if (length(sample.ids) != length(rna.input)) {
    stop("'sample.ids' must have one id per RNA/ATAC pair.")
  }
  if (!is.null(x) && !inherits(x, "iCellR2")) {
    stop("x should be an iCellR2 object or NULL.")
  }

  if (is.null(x)) {
    samps <- lapply(seq_along(rna.input), function(i) {
      path <- .ic2_multiome_matrix(rna.input[[i]], which = "rna")
      .ic2_load10x_typed(
        path, sample.ids[[i]], "Gene Expression", "rna", verbose
      )
    })
    if (length(samps) == 1L) {
      x <- make.obj(
        samps[[1L]], out.dir = out.dir, modality = "rna",
        overwrite = overwrite, verbose = verbose
      )
    } else {
      agg <- data.aggregation(
        samples = samps, modality = "rna", output = "virtual", verbose = verbose
      )
      x <- make.obj(
        agg, out.dir = out.dir, modality = "rna",
        overwrite = overwrite, verbose = verbose
      )
    }
  }

  add.atac(
    x,
    input = vapply(atac.input, function(p) .ic2_multiome_matrix(p, which = "atac"), character(1)),
    sample.ids = sample.ids,
    peaks = peaks,
    feature.type = "Peaks",
    verbose = verbose
  )
}

#' Pick the RNA or ATAC matrix inside a 10x / ARC folder. A bare file is used as-is.
#' @keywords internal
.ic2_multiome_matrix <- function(input, which = c("rna", "atac")) {
  which <- match.arg(which)
  input <- as.character(input)[1L]
  if (!dir.exists(input)) {
    if (!file.exists(input)) stop("Input not found: ", input)
    return(normalizePath(input, winslash = "/", mustWork = TRUE))
  }
  if (identical(which, "atac")) {
    return(.ic2_atac_find(input)$matrix)
  }
  rel <- c(
    "filtered_feature_bc_matrix.h5",
    "raw_feature_bc_matrix.h5",
    file.path("outs", "filtered_feature_bc_matrix.h5"),
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix"
  )
  for (r in rel) {
    p <- file.path(input, r)
    if (file.exists(p) || dir.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "No RNA matrix under ", input,
    ". Expected filtered_feature_bc_matrix.h5. Pass that file as rna.input."
  )
}

#' \code{load10x} virtual handle. Falls back when the file has no feature-type column.
#' @keywords internal
.ic2_load10x_typed <- function(path, sample.id, feature.type, modality, verbose) {
  meta <- file.path(
    "iCellR2_import",
    paste0(sample.id, "_", modality, "_tenx_meta.json")
  )
  samp <- tryCatch(
    load10x(
      path,
      sample.id = sample.id,
      feature.type = feature.type,
      output = "virtual",
      modality = modality,
      h5.file = meta,
      overwrite = TRUE,
      verbose = verbose
    ),
    error = function(e) e
  )
  if (!inherits(samp, "error")) return(samp)
  if (is.null(feature.type)) stop(samp$message)
  load10x(
    path,
    sample.id = sample.id,
    feature.type = NULL,
    output = "virtual",
    modality = modality,
    h5.file = meta,
    overwrite = TRUE,
    verbose = verbose
  )
}
