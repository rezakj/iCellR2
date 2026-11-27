#' Add CITE-seq antibody counts as their own modality
#'
#' Writes \code{modalities/adt/assays/counts} aligned to \code{cell_data}.
#' The RNA assay is not modified. Antibody names are prefixed with \code{ADT_}.
#'
#' @param x An iCellR2 object.
#' @param adt.data Features-by-cells data.frame or matrix, or a TSV path
#'   (rownames = antibodies, colnames = cell ids).
#' @param assay.name Assay to write. Default \code{"counts"}.
#' @param overwrite Replace an existing ADT assay file.
#' @return The iCellR2 object.
#' @export
add.adt <- function(x = NULL, adt.data = NULL, assay.name = "counts", overwrite = TRUE) {
  if (!inherits(x, "iCellR2")) {
    stop("x should be an object of class iCellR2")
  }
  if (is.null(adt.data)) stop("Provide adt.data (matrix or TSV path).")
  if (is.character(adt.data) && length(adt.data) == 1L) {
    if (!file.exists(adt.data)) stop("ADT file not found: ", adt.data)
    adt.data <- utils::read.delim(
      adt.data, sep = "\t", check.names = FALSE,
      row.names = 1L, stringsAsFactors = FALSE
    )
  }
  if (is.data.frame(adt.data)) adt.data <- as.matrix(adt.data)
  if (!is.matrix(adt.data) && !inherits(adt.data, "Matrix")) {
    stop("adt.data must be a matrix, data.frame, or TSV path.")
  }
  storage.mode(adt.data) <- "double"
  adt.data[is.na(adt.data)] <- 0

  cell_ids <- .ic2_adt_cell_ids(x)
  adt.data <- .ic2_adt_orient(adt.data, cell_ids)
  if (is.null(rownames(adt.data)) || any(!nzchar(rownames(adt.data)))) {
    stop("ADT matrix needs antibody rownames.")
  }
  rownames(adt.data) <- .ic2_adt_feature_ids(rownames(adt.data))
  matched <- .ic2_adt_match_cells(colnames(adt.data), cell_ids)
  if (!any(!is.na(matched))) {
    stop("No ADT columns matched cell_data cell ids.")
  }
  n_miss <- sum(is.na(matched))
  if (n_miss > 0L) {
    message("add.adt: ", n_miss, " object cells have no ADT counts (stored as 0).")
  }

  present <- which(!is.na(matched))
  sm <- .ic2_to_dgCMatrix(adt.data[, matched[present], drop = FALSE])
  mat <- Matrix::sparseMatrix(
    i = sm@i + 1L,
    j = rep.int(present, diff(sm@p)),
    x = sm@x,
    dims = c(nrow(sm), length(cell_ids)),
    dimnames = list(rownames(sm), cell_ids),
    repr = "C"
  )

  dest <- .ic2_soda_assay_path(x, "adt", assay.name, ext = ".h5")
  handle <- write_sparse_hdf5(mat, filepath = dest, group = "counts", overwrite = overwrite)
  handle$modality <- "adt"

  if (is.null(x@modalities$adt) || !is.list(x@modalities$adt$assays)) {
    x@modalities$adt <- list(
      assays = list(), active_assay = NULL, feature_data = NULL,
      reductions = list(), active_reduction = NULL, graphs = list(),
      sketches = list(), markers = list(), models = list(),
      params = list(), misc = list()
    )
  }
  x <- .ic2_assign_assay(
    x, modality = "adt", assay.name = assay.name,
    handles = list(handle), virtual = FALSE, set.active = TRUE
  )
  fd <- data.frame(
    feature_id = rownames(mat),
    feature_name = rownames(mat),
    feature_type = "Antibody Capture",
    stringsAsFactors = FALSE
  )
  root <- x@storage$root_dir
  fd_path <- .ic2_soda_write_table(fd, file.path(root, "metadata", "adt_feature_data"))
  x@modalities$adt$feature_data <- fd
  x@storage$metadata_paths$adt_feature_data <- fd_path
  x@storage$modality_paths$adt <- file.path(root, "modalities", "adt")
  x@params$adt <- list(assay.name = assay.name, n_features = nrow(mat), n_cells = ncol(mat))
  x <- .ic2_log_command(
    x, "add.adt",
    list(assay.name = assay.name, n_features = nrow(mat), n_cells = ncol(mat)),
    list(counts = dest, feature_data = fd_path)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}

#' @keywords internal
.ic2_adt_cell_ids <- function(x) {
  cd <- x@cell_data
  if (!is.null(cd) && "cell_id" %in% colnames(cd)) {
    return(as.character(cd$cell_id))
  }
  if (!is.null(rownames(cd)) && length(rownames(cd)) > 0L) {
    return(as.character(rownames(cd)))
  }
  stop("cell_data has no cell_id column.")
}

#' @keywords internal
.ic2_adt_orient <- function(mat, cell_ids) {
  cn <- colnames(mat)
  rn <- rownames(mat)
  hit_col <- if (length(cn)) mean(cn %in% cell_ids) else 0
  hit_row <- if (length(rn)) mean(rn %in% cell_ids) else 0
  if (hit_row > hit_col) mat <- t(mat)
  mat
}

#' @keywords internal
.ic2_adt_feature_ids <- function(ids) {
  ids <- as.character(ids)
  ids <- gsub("/", "_or_", ids, fixed = TRUE)
  ids <- gsub(" ", "_", ids, fixed = TRUE)
  need <- !startsWith(ids, "ADT_")
  ids[need] <- paste0("ADT_", ids[need])
  make.unique(ids, sep = "_")
}

#' @keywords internal
.ic2_adt_match_cells <- function(adt_cells, cell_ids) {
  m <- match(cell_ids, adt_cells)
  if (anyNA(m)) {
    m2 <- match(
      gsub("-", ".", cell_ids, fixed = TRUE),
      gsub("-", ".", adt_cells, fixed = TRUE)
    )
    m[is.na(m)] <- m2[is.na(m)]
  }
  m
}

#' Per-cell CLR on a features-by-cells sparse matrix. Zeros stay implicit.
#' Seurat-style: geo-mean of log1p over all features, then log1p(x / geo).
#' @keywords internal
.ic2_adt_clr <- function(mat) {
  mat <- .ic2_to_dgCMatrix(mat)
  n_feat <- nrow(mat)
  if (n_feat < 1L) return(mat)
  x <- mat@x
  p <- mat@p
  for (col in seq_len(ncol(mat))) {
    a <- p[col] + 1L
    b <- p[col + 1L]
    if (a > b) next
    vals <- x[a:b]
    pos <- vals > 0
    if (!any(pos)) next
    geo <- exp(sum(log1p(vals[pos])) / n_feat)
    x[a:b] <- log1p(vals / geo)
  }
  mat@x <- x
  mat
}
