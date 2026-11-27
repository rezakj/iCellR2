#' Spatial neighbor graph
#'
#' Neighbor rule follows the assay. \code{method = "auto"} uses hex for classic
#' Visium, a square lattice for Visium HD, and a micron radius for Xenium.
#' Bins with \code{in_tissue == 0} are left out. Edges are one undirected list
#' at \code{graphs/spatial_knn/edges.tsv}. No spot-by-spot distance matrix.
#'
#' @param x An iCellR2 object after \code{\link{add.visium}} or \code{\link{add.xenium}}.
#' @param method \code{"auto"}, \code{"hex"}, \code{"square"}, \code{"knn"}, or \code{"radius"}.
#' @param k Neighbors for \code{knn}. For \code{square}, 4 is orthogonal and 8
#'   includes diagonals. For \code{radius}, this is the maximum neighbors kept
#'   inside the radius. Default \code{6}.
#' @param radius.um Radius in microns for \code{method = "radius"}. Default \code{30}.
#'   Distances use \code{spatial_x_um} / \code{spatial_y_um}.
#' @return The object.
#' @examples
#' \dontrun{
#' obj <- spatial.graph(obj)                 # hex, square, or 30 um radius
#' obj <- spatial.graph(obj, method = "square", k = 8)
#' }
#' @export
spatial.graph <- function(x = NULL, method = c("auto", "hex", "square", "knn", "radius"), k = 6L, radius.um = 30) {
  method <- match.arg(method)
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  cd <- x@cell_data
  if (!all(c("spatial_x", "spatial_y") %in% colnames(cd))) {
    stop("No spatial coordinates. Run add.visium() or add.xenium() first.")
  }
  ids <- if ("cell_id" %in% colnames(cd)) as.character(cd$cell_id) else rownames(cd)
  ok <- is.finite(cd$spatial_x) & is.finite(cd$spatial_y)
  if ("in_tissue" %in% colnames(cd)) ok <- ok & (is.na(cd$in_tissue) | cd$in_tissue == 1L)
  cd <- cd[ok, , drop = FALSE]
  ids <- ids[ok]
  if (nrow(cd) < 2L) stop("Need at least two spots.")
  samples <- if ("spatial_sample" %in% colnames(cd)) as.character(cd$spatial_sample) else rep("sample", nrow(cd))
  samples[is.na(samples)] <- "sample"
  from <- character(0)
  to <- character(0)
  used <- character(0)
  for (s in unique(samples)) {
    sel <- which(samples == s)
    if (length(sel) < 2L) next
    how <- .ic2_spatial_method(method, cd[sel, , drop = FALSE])
    used <- c(used, how)
    pairs <- if (identical(how, "hex")) {
      .ic2_spatial_hex(ids[sel], cd$array_row[sel], cd$array_col[sel])
    } else if (identical(how, "square")) {
      .ic2_spatial_square(ids[sel], cd$array_row[sel], cd$array_col[sel], as.integer(k))
    } else if (identical(how, "radius")) {
      .ic2_spatial_radius(ids[sel], cd[sel, , drop = FALSE], radius.um, as.integer(k))
    } else {
      .ic2_spatial_knn(ids[sel], cd$spatial_x[sel], cd$spatial_y[sel], as.integer(k))
    }
    from <- c(from, pairs$from)
    to <- c(to, pairs$to)
  }
  used <- unique(used)
  if (length(from) < 1L) stop("Spatial graph has no edges.")
  root <- x@storage$root_dir
  gdir <- file.path(root, "graphs", "spatial_knn")
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  path <- .ic2_vdj_write_tsv(
    data.frame(from = from, to = to, stringsAsFactors = FALSE),
    file.path(gdir, "edges.tsv")
  )
  .ic2_soda_write_json(list(
    name = "spatial_knn", method = paste(used, collapse = "+"),
    k = as.integer(k), radius_um = as.numeric(radius.um),
    n_edges = length(from), path = path
  ), file.path(gdir, "graph.json"))
  x@graphs$spatial_knn <- path
  x <- .ic2_log_command(
    x, "spatial.graph",
    list(method = paste(used, collapse = "+"), n_edges = length(from)),
    list(edges = path)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}
#' @keywords internal
.ic2_spatial_method <- function(method, cd) {
  if (!identical(method, "auto")) return(method)
  plat <- if ("spatial_platform" %in% colnames(cd)) {
    unique(stats::na.omit(as.character(cd$spatial_platform)))
  } else {
    character(0)
  }
  if (any(plat == "xenium")) return("radius")
  if (any(plat == "visium_hd")) return("square")
  if (all(c("array_row", "array_col") %in% colnames(cd)) &&
      any(is.finite(cd$array_row) & is.finite(cd$array_col))) {
    return("hex")
  }
  "knn"
}

#' @keywords internal
.ic2_spatial_square <- function(ids, row, col, k = 4L) {
  r <- as.integer(row)
  c <- as.integer(col)
  offs <- if (as.integer(k) >= 8L) {
    rbind(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L),
          c(-1L, -1L), c(-1L, 1L), c(1L, -1L), c(1L, 1L))
  } else {
    rbind(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L))
  }
  nr <- sapply(seq_len(nrow(offs)), function(i) r + offs[i, 1L])
  nc <- sapply(seq_len(nrow(offs)), function(i) c + offs[i, 2L])
  .ic2_spatial_lookup(ids, r, c, nr, nc)
}

#' @keywords internal
.ic2_spatial_lookup <- function(ids, r, c, nr, nc) {
  key <- paste(r, c, sep = ",")
  lookup <- stats::setNames(ids, key)
  nkey <- paste(as.vector(nr), as.vector(nc), sep = ",")
  hit <- lookup[nkey]
  self <- rep(ids, ncol(as.matrix(nr)))
  keep <- !is.na(hit) & !is.na(self) & self < hit
  list(from = self[keep], to = unname(hit[keep]))
}

#' @keywords internal
.ic2_spatial_radius <- function(ids, cd, radius.um, k) {
  if (!all(c("spatial_x_um", "spatial_y_um") %in% colnames(cd)) ||
      !any(is.finite(cd$spatial_x_um))) {
    stop("method = 'radius' needs spatial_x_um / spatial_y_um. Run add.visium() or add.xenium().")
  }
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Package 'FNN' is required for a radius graph.")
  xy <- cbind(as.numeric(cd$spatial_x_um), as.numeric(cd$spatial_y_um))
  kk <- min(max(as.integer(k), 15L), nrow(xy) - 1L)
  kn <- FNN::get.knn(xy, k = kk)
  self <- rep(ids, each = kk)
  nbr <- ids[as.vector(t(kn$nn.index))]
  dist <- as.vector(t(kn$nn.dist))
  keep <- self < nbr & is.finite(dist) & dist <= as.numeric(radius.um)
  list(from = self[keep], to = nbr[keep])
}

#' @keywords internal
.ic2_spatial_hex <- function(ids, row, col) {
  r <- as.integer(row)
  c <- as.integer(col)
  even <- !is.na(r) & (r %% 2L == 0L)
  nr <- cbind(r, r, r - 1L, r - 1L, r + 1L, r + 1L)
  shift <- ifelse(even, -1L, 1L)
  nc <- cbind(c - 1L, c + 1L, c, c + shift, c, c + shift)
  .ic2_spatial_lookup(ids, r, c, nr, nc)
}

#' @keywords internal
.ic2_spatial_knn <- function(ids, x, y, k) {
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Package 'FNN' is required for coordinate kNN.")
  xy <- cbind(as.numeric(x), as.numeric(y))
  kk <- min(as.integer(k), nrow(xy) - 1L)
  kn <- FNN::get.knn(xy, k = kk)
  self <- rep(ids, each = kk)
  nbr <- ids[as.vector(t(kn$nn.index))]
  keep <- self < nbr
  pair <- unique(paste(self[keep], nbr[keep], sep = "\r"))
  sp <- strsplit(pair, "\r", fixed = TRUE)
  list(from = vapply(sp, `[[`, character(1), 1L), to = vapply(sp, `[[`, character(1), 2L))
}

