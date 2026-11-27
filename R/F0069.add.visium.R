#' Attach Visium images and spot positions to an iCellR2 object
#'
#' Copies each spatial folder under \code{images/<sample>/} and adds
#' \code{spatial_x}, \code{spatial_y}, and \code{in_tissue} to \code{cell_data}.
#' Barcodes are matched to RNA cells and are not rewritten. Pass one directory,
#' several directories, or the list returned by \code{\link{capture.image.10x}}.
#'
#' @param x An iCellR2 object.
#' @param input Directory, character vector of directories, or a capture list.
#'   \code{image.data.list} is the same argument.
#' @param image.data.list Alias of \code{input}.
#' @param condition.names Sample ids, one per directory. Default: folder name.
#' @param bin.size Visium HD bin in microns when \code{input} is a Space Ranger
#'   \code{outs} folder. Default \code{8}. Use \code{2} or \code{16} to override.
#' @return The object. Coordinates are full-resolution pixels in
#'   \code{spatial_x}/\code{spatial_y} and microns in \code{spatial_x_um}/
#'   \code{spatial_y_um}. \code{spatial_platform} is \code{visium} or \code{visium_hd}.
#'   \code{add.10x.image} is the same function.
#' @examples
#' \dontrun{
#' obj <- add.visium(obj, input = "outs/spatial", condition.names = "sample1")
#' obj <- add.visium(obj, input = "spaceranger/outs", bin.size = 8, condition.names = "hd1")
#' }
#' @rdname add.visium
#' @export
add.visium <- function(
    x = NULL,
    input = NULL,
    image.data.list = NULL,
    condition.names = NULL,
    bin.size = 8
) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  src <- if (!is.null(input)) input else image.data.list
  if (is.null(src)) stop("Provide input: a spatial directory or a list of them.")
  items <- .ic2_spatial_as_items(src, condition.names, bin.size)
  root <- x@storage$root_dir
  cell_ids <- .ic2_adt_cell_ids(x)
  cd <- x@cell_data
  cd$spatial_x <- NA_real_
  cd$spatial_y <- NA_real_
  cd$spatial_x_um <- NA_real_
  cd$spatial_y_um <- NA_real_
  cd$spatial_um_per_px <- NA_real_
  cd$array_row <- NA_integer_
  cd$array_col <- NA_integer_
  cd$in_tissue <- NA_integer_
  cd$spatial_sample <- NA_character_
  cd$spatial_platform <- NA_character_
  cd$spatial_bin_um <- NA_real_
  info <- list()
  for (i in seq_along(items)) {
    it <- items[[i]]
    sample_id <- it$sample_id
    dest <- file.path(root, "images", sample_id)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    pos_out <- .ic2_vdj_write_tsv(it$position, file.path(dest, "tissue_positions.tsv"))
    pics <- it$image_paths
    if (is.null(pics) && !is.null(it$image_path) && !is.na(it$image_path)) pics <- it$image_path
    for (pic in pics) {
      if (!is.na(pic) && file.exists(pic)) {
        file.copy(pic, file.path(dest, basename(pic)), overwrite = TRUE)
      }
    }
    if (!is.null(it$scale)) {
      .ic2_soda_write_json(it$scale, file.path(dest, "scalefactors.json"))
    }
    hit <- .ic2_vdj_match_cells(it$position$barcode, cell_ids, sample_id)
    ok <- !is.na(hit)
    if (!any(ok)) {
      message("add.visium: no barcodes in ", sample_id, " matched cell_data.")
    } else {
      cd$spatial_x[hit[ok]] <- as.numeric(it$position$pxl_col_in_fullres[ok])
      cd$spatial_y[hit[ok]] <- as.numeric(it$position$pxl_row_in_fullres[ok])
      if ("array_row" %in% colnames(it$position)) {
        cd$array_row[hit[ok]] <- as.integer(it$position$array_row[ok])
      }
      if ("array_col" %in% colnames(it$position)) {
        cd$array_col[hit[ok]] <- as.integer(it$position$array_col[ok])
      }
      if ("in_tissue" %in% colnames(it$position)) {
        cd$in_tissue[hit[ok]] <- as.integer(it$position$in_tissue[ok])
      }
      cd$spatial_sample[hit[ok]] <- sample_id
      plat <- it$platform %||% "visium"
      cd$spatial_platform[hit[ok]] <- plat
      if (is.finite(it$bin_um %||% NA_real_)) cd$spatial_bin_um[hit[ok]] <- it$bin_um
      um <- it$um_per_px
      if (is.finite(um %||% NA_real_)) {
        cd$spatial_um_per_px[hit[ok]] <- um
        cd$spatial_x_um[hit[ok]] <- cd$spatial_x[hit[ok]] * um
        cd$spatial_y_um[hit[ok]] <- cd$spatial_y[hit[ok]] * um
      }
    }
    info[[sample_id]] <- list(
      dir = dest, positions = pos_out, n_matched = sum(ok),
      platform = it$platform %||% "visium", bin_um = it$bin_um, um_per_px = it$um_per_px
    )
  }
  x@cell_data <- cd
  x <- .ic2_persist_cell_data(x)
  x@misc$spatial <- info
  x@storage$image_paths <- info
  x <- .ic2_log_command(
    x, "add.visium",
    list(samples = names(info)),
    list(images = file.path(root, "images"))
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}

#' @rdname add.visium
#' @export
add.10x.image <- add.visium

#' @keywords internal
.ic2_spatial_as_items <- function(src, condition.names, bin.size = 8) {
  if (is.character(src)) {
    src <- as.list(src)
  }
  if (is.list(src) && !is.null(src$position)) src <- list(src)
  if (!is.list(src) || length(src) < 1L) stop("input must be a directory or a capture list.")
  if (!is.null(condition.names) && length(condition.names) != length(src)) {
    stop("condition.names must have one id per spatial input.")
  }
  lapply(seq_along(src), function(i) {
    el <- src[[i]]
    if (is.character(el) && length(el) == 1L) el <- capture.image.10x(el, bin.size = bin.size)
    if (is.null(el$position)) stop("Spatial input ", i, " has no position table.")
    sid <- if (!is.null(condition.names)) as.character(condition.names[[i]]) else {
      p <- el$positions_path %||% el$image_path %||% paste0("sample", i)
      basename(dirname(p))
    }
    if (!nzchar(sid) || identical(sid, ".")) sid <- paste0("sample", i)
    el$sample_id <- sid
    el
  })
}
