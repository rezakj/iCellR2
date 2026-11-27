#' Attach Xenium cell centroids
#'
#' Reads \code{cells.parquet} or \code{cells.csv.gz} and matches \code{cell_id}
#' to the RNA barcodes. Centroids are stored in microns
#' (\code{spatial_x}, \code{spatial_y}, and \code{spatial_x_um}). The morphology
#' image and \code{transcripts.parquet} are recorded and not loaded.
#'
#' @param x An iCellR2 object whose cells are the Xenium cell-feature matrix.
#' @param input Xenium output directory.
#' @param condition.names Sample id. Default: folder name.
#' @return The object.
#' @examples
#' \dontrun{
#' obj <- add.xenium(obj, input = "xenium_outs", condition.names = "xen1")
#' }
#' @export
add.xenium <- function(x = NULL, input = NULL, condition.names = NULL) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  if (is.null(input) || !dir.exists(input)) stop("Provide a Xenium output directory.")
  input <- normalizePath(input, winslash = "/", mustWork = TRUE)
  cells_path <- .ic2_xenium_find(input, c("cells.parquet", "cells.csv.gz", "cells.csv"))
  if (is.na(cells_path)) stop("No cells.parquet or cells.csv.gz in ", input)
  cells <- .ic2_xenium_read_cells(cells_path)
  sample_id <- if (!is.null(condition.names)) as.character(condition.names)[1L] else basename(input)
  um <- .ic2_xenium_pixel_size(input)
  img <- .ic2_xenium_find(input, c(
    "morphology_mip.ome.tif", "morphology_focus.ome.tif", "morphology.ome.tif"
  ))
  tx <- .ic2_xenium_find(input, c("transcripts.parquet", "transcripts.csv.gz"))
  root <- x@storage$root_dir
  dest <- file.path(root, "images", sample_id)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  .ic2_vdj_write_tsv(cells, file.path(dest, "cells.tsv"))
  .ic2_soda_write_json(list(
    platform = "xenium",
    um_per_px = um,
    image = if (!is.na(img)) img else NA_character_,
    transcripts = if (!is.na(tx)) tx else NA_character_,
    cells = cells_path
  ), file.path(dest, "xenium.json"))
  cell_ids <- .ic2_adt_cell_ids(x)
  cd <- x@cell_data
  cd$spatial_x <- NA_real_
  cd$spatial_y <- NA_real_
  cd$spatial_x_um <- NA_real_
  cd$spatial_y_um <- NA_real_
  cd$spatial_um_per_px <- NA_real_
  cd$spatial_sample <- NA_character_
  cd$spatial_platform <- NA_character_
  cd$in_tissue <- NA_integer_
  hit <- .ic2_vdj_match_cells(cells$cell_id, cell_ids, sample_id)
  ok <- !is.na(hit)
  if (!any(ok)) stop("add.xenium: no cell_id values matched the RNA barcodes.")
  cd$spatial_x[hit[ok]] <- cells$x_centroid[ok]
  cd$spatial_y[hit[ok]] <- cells$y_centroid[ok]
  cd$spatial_x_um[hit[ok]] <- cells$x_centroid[ok]
  cd$spatial_y_um[hit[ok]] <- cells$y_centroid[ok]
  if (is.finite(um)) cd$spatial_um_per_px[hit[ok]] <- um
  cd$spatial_sample[hit[ok]] <- sample_id
  cd$spatial_platform[hit[ok]] <- "xenium"
  cd$in_tissue[hit[ok]] <- 1L
  x@cell_data <- cd
  x <- .ic2_persist_cell_data(x)
  info <- list(
    dir = dest, n_matched = sum(ok), platform = "xenium",
    um_per_px = um, image = img, transcripts = tx
  )
  x@misc$spatial <- info
  x@storage$image_paths <- info
  x <- .ic2_log_command(
    x, "add.xenium",
    list(sample = sample_id, n_matched = sum(ok)),
    list(xenium = dest)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}

#' @keywords internal
.ic2_xenium_find <- function(dir, names) {
  hits <- list.files(dir, recursive = TRUE, full.names = TRUE)
  base <- basename(hits)
  for (nm in names) {
    i <- match(nm, base)
    if (!is.na(i)) return(normalizePath(hits[[i]], winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

#' @keywords internal
.ic2_xenium_read_cells <- function(path) {
  df <- if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (requireNamespace("arrow", quietly = TRUE)) {
      as.data.frame(arrow::read_parquet(path))
    } else if (requireNamespace("nanoparquet", quietly = TRUE)) {
      as.data.frame(nanoparquet::read_parquet(path))
    } else {
      stop("Reading ", path, " needs package 'arrow' or 'nanoparquet'.")
    }
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  names(df) <- tolower(names(df))
  if (!"cell_id" %in% names(df)) stop("Xenium cells table has no cell_id: ", path)
  if (!all(c("x_centroid", "y_centroid") %in% names(df))) {
    stop("Xenium cells table has no x_centroid / y_centroid: ", path)
  }
  df$cell_id <- as.character(df$cell_id)
  df$x_centroid <- as.numeric(df$x_centroid)
  df$y_centroid <- as.numeric(df$y_centroid)
  df[, c("cell_id", "x_centroid", "y_centroid"), drop = FALSE]
}

#' @keywords internal
.ic2_xenium_pixel_size <- function(dir) {
  path <- .ic2_xenium_find(dir, "experiment.xenium")
  if (is.na(path)) return(NA_real_)
  obj <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  .ic2_xenium_find_number(obj, c("pixel_size", "microns_per_pixel"))
}

#' @keywords internal
.ic2_xenium_find_number <- function(obj, keys) {
  if (is.null(obj)) return(NA_real_)
  if (is.list(obj)) {
    nm <- names(obj)
    for (k in keys) {
      if (!is.null(nm) && k %in% nm && is.numeric(obj[[k]]) && length(obj[[k]]) == 1L) {
        return(as.numeric(obj[[k]]))
      }
    }
    if (!is.null(nm) && "pixels_per_um" %in% nm && is.numeric(obj$pixels_per_um) && obj$pixels_per_um > 0) {
      return(1 / as.numeric(obj$pixels_per_um))
    }
    for (el in obj) {
      hit <- .ic2_xenium_find_number(el, keys)
      if (is.finite(hit)) return(hit)
    }
  }
  NA_real_
}
