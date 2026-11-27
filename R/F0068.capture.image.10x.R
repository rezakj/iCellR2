#' Read a 10x Visium spatial folder
#'
#' Looks for \code{tissue_positions.csv} or \code{tissue_positions_list.csv},
#' \code{scalefactors_json.json}, and a tissue image. Barcodes are not rewritten.
#' The image file is not decoded; the returned list holds its path. Spot
#' tables stay a data.frame (one row per barcode), which is the right size
#' even for Visium HD.
#'
#' @param dir.10x Spatial directory, a Visium HD \code{square_*um} folder, or a
#'   Space Ranger \code{outs} folder.
#' @param bin.size Visium HD bin in microns when \code{dir.10x} contains
#'   \code{binned_outputs}. One of 8 (default), 16, or 2. Ignored for classic Visium.
#' @return A list with \code{position}, \code{scale}, \code{image_path},
#'   \code{platform} (\code{visium} or \code{visium_hd}), \code{bin_um}, and
#'   \code{um_per_px}.
#' @examples
#' \dontrun{
#' img <- capture.image.10x("outs/spatial")
#' hd <- capture.image.10x("spaceranger/outs", bin.size = 8)
#' }
#' @export
capture.image.10x <- function(dir.10x = NULL, bin.size = 8) {
  if (is.null(dir.10x) || !dir.exists(dir.10x)) {
    stop("Provide a 10x spatial directory (tissue positions, scalefactors, image).")
  }
  dir.10x <- .ic2_spatial_hd_dir(dir.10x, bin.size)
  pos_path <- .ic2_spatial_find(dir.10x, c(
    "tissue_positions.parquet", "tissue_positions.csv", "tissue_positions_list.csv"
  ))
  if (is.na(pos_path)) stop("No tissue_positions.parquet or tissue_positions.csv in ", dir.10x)
  scale_path <- .ic2_spatial_find(dir.10x, "scalefactors_json.json")
  imgs <- .ic2_spatial_find_all(dir.10x, c(
    "tissue_lowres_image.png", "tissue_hires_image.png"
  ))
  pos <- .ic2_spatial_read_positions(pos_path)
  scale <- if (!is.na(scale_path)) {
    jsonlite::fromJSON(scale_path, simplifyVector = TRUE)
  } else {
    NULL
  }
  platform <- .ic2_spatial_platform(dir.10x, scale)
  list(
    position = pos,
    scale = scale,
    image_path = if (length(imgs)) imgs[[1L]] else NA_character_,
    image_paths = imgs,
    positions_path = pos_path,
    scale_path = if (!is.na(scale_path)) scale_path else NA_character_,
    platform = platform,
    bin_um = .ic2_spatial_bin_um(dir.10x, scale),
    um_per_px = .ic2_spatial_um_per_px(scale, platform),
    spatial_dir = dir.10x
  )
}

#' @keywords internal
.ic2_spatial_find <- function(dir.10x, names) {
  names <- as.character(names)
  roots <- c(dir.10x, file.path(dir.10x, "spatial"))
  for (nm in names) {
    for (root in roots) {
      p <- file.path(root, nm)
      if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  NA_character_
}

#' Folder that holds one Visium HD bin size, or \code{path} when it is not HD.
#' @keywords internal
.ic2_spatial_hd_dir <- function(path, bin.size = 8) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (grepl("^square_[0-9]+um$", basename(path))) return(path)
  binned <- file.path(path, "binned_outputs")
  if (!dir.exists(binned)) return(path)
  um <- as.integer(bin.size %||% 8L)
  if (!um %in% c(2L, 8L, 16L)) stop("bin.size must be 2, 8, or 16.")
  dest <- file.path(binned, sprintf("square_%03dum", um))
  if (!dir.exists(dest)) stop("Visium HD bin folder not found: ", dest)
  normalizePath(dest, winslash = "/", mustWork = TRUE)
}

#' Count matrix for \code{load10x}: the HD \code{.h5} when \code{binned_outputs} exists.
#' @keywords internal
.ic2_spatial_resolve_matrix <- function(path, bin.size = 8) {
  if (is.null(path) || !dir.exists(path)) return(path)
  base <- .ic2_spatial_hd_dir(path, bin.size)
  if (identical(base, normalizePath(path, winslash = "/", mustWork = TRUE)) &&
      !grepl("^square_[0-9]+um$", basename(base))) {
    return(path)
  }
  h5 <- file.path(base, "filtered_feature_bc_matrix.h5")
  if (file.exists(h5)) return(normalizePath(h5, winslash = "/", mustWork = TRUE))
  base
}

#' @keywords internal
.ic2_spatial_find_all <- function(dir.10x, names) {
  names <- as.character(names)
  roots <- c(dir.10x, file.path(dir.10x, "spatial"))
  found <- character(0)
  for (nm in names) {
    for (root in roots) {
      p <- file.path(root, nm)
      if (file.exists(p)) found <- c(found, normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  unique(found)
}

#' @keywords internal
.ic2_spatial_platform <- function(path, scale) {
  if (grepl("square_[0-9]+um", path)) return("visium_hd")
  if (!is.null(scale) && !is.null(scale$bin_size_um)) return("visium_hd")
  "visium"
}

#' Full-resolution microns per pixel. Classic Visium uses the 55 um spot diameter.
#' @keywords internal
.ic2_spatial_um_per_px <- function(scale, platform) {
  if (is.null(scale)) return(NA_real_)
  if (!is.null(scale$microns_per_pixel) && is.finite(scale$microns_per_pixel)) {
    return(as.numeric(scale$microns_per_pixel))
  }
  diam <- scale$spot_diameter_fullres
  if (identical(platform, "visium") && !is.null(diam) && is.finite(diam) && diam > 0) {
    return(55 / as.numeric(diam))
  }
  NA_real_
}

#' @keywords internal
.ic2_spatial_bin_um <- function(path, scale) {
  m <- regexpr("square_([0-9]+)um", path, perl = TRUE)
  if (m > 0L) {
    hit <- regmatches(path, m)
    return(as.numeric(sub("^square_0*", "", sub("um$", "", hit))))
  }
  if (!is.null(scale$bin_size_um) && is.finite(scale$bin_size_um)) {
    return(as.numeric(scale$bin_size_um))
  }
  if (identical(.ic2_spatial_platform(path, scale), "visium")) return(55)
  NA_real_
}

#' @keywords internal
.ic2_spatial_read_positions <- function(path) {
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    df <- if (requireNamespace("arrow", quietly = TRUE)) {
      as.data.frame(arrow::read_parquet(path))
    } else if (requireNamespace("nanoparquet", quietly = TRUE)) {
      as.data.frame(nanoparquet::read_parquet(path))
    } else {
      stop("Reading ", path, " needs package 'arrow' or 'nanoparquet'.")
    }
    names(df) <- tolower(names(df))
    if (!"barcode" %in% names(df)) stop("Positions file has no barcode column: ", path)
    df$barcode <- as.character(df$barcode)
    if (!"pxl_col_in_fullres" %in% names(df) && "pxl_col" %in% names(df)) {
      df$pxl_col_in_fullres <- df$pxl_col
    }
    if (!"pxl_row_in_fullres" %in% names(df) && "pxl_row" %in% names(df)) {
      df$pxl_row_in_fullres <- df$pxl_row
    }
    return(df)
  }
  line <- readLines(path, n = 1L, warn = FALSE)
  header <- grepl("barcode", line, ignore.case = TRUE)
  if (header) {
    df <- utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- utils::read.csv(path, header = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- c(
      "barcode", "in_tissue", "array_row", "array_col",
      "pxl_row_in_fullres", "pxl_col_in_fullres"
    )[seq_len(ncol(df))]
  }
  names(df) <- tolower(names(df))
  if (!"barcode" %in% names(df)) stop("Positions file has no barcode column: ", path)
  df$barcode <- as.character(df$barcode)
  if (!"pxl_col_in_fullres" %in% names(df) && "pxl_col" %in% names(df)) {
    df$pxl_col_in_fullres <- df$pxl_col
  }
  if (!"pxl_row_in_fullres" %in% names(df) && "pxl_row" %in% names(df)) {
    df$pxl_row_in_fullres <- df$pxl_row
  }
  df
}
