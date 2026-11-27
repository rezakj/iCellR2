#' Plot spots on Visium coordinates
#'
#' Uses \code{spatial_x} and \code{spatial_y} from \code{cell_data}. Classic
#' Visium and Xenium are drawn as points, down-sampled to \code{max.spots}.
#' Visium HD is drawn as a binned raster so millions of bins are not points.
#' The whole slide uses the low-resolution image. \code{xlim}/\code{ylim} switch
#' to the hires image and crop. A gene color reads that one gene in cell chunks.
#'
#' @param x An iCellR2 object after \code{\link{add.visium}} or \code{\link{add.xenium}}.
#' @param col.by \code{"clusters"}, \code{"gene"}, \code{"cc"}, \code{"nGenes"},
#'   \code{"UMIs"}, or \code{"mt"}.
#' @param gene Gene id when \code{col.by = "gene"}.
#' @param conds.to.plot \code{spatial_sample} or \code{condition_id} values to keep.
#' @param max.spots Maximum points, or maximum raster bins. Default \code{200000}.
#' @param mark \code{"auto"} (points, or a raster for Visium HD), \code{"point"},
#'   or \code{"raster"}.
#' @param xlim,ylim Optional crop in the plotted coordinate system (image pixels
#'   for Visium, microns for Xenium). A crop uses the hires image when it exists.
#' @param cell.size,cell.transparency,cell.colors,back.col Plot settings.
#' @param interactive Write an HTML widget of the sample.
#' @param out.name PNG, or HTML when \code{interactive = TRUE}.
#' @param assay.name Assay for \code{col.by = "gene"}. Default \code{"data"}.
#' @param scaleValue,min.scale,max.scale,anno.clust,anno.size,anno.col,data.type
#'   Classic arguments. Cluster labels are drawn when \code{anno.clust = TRUE}.
#' @return A ggplot, or the HTML path when interactive.
#' @examples
#' \dontrun{
#' spatial.plot(obj, col.by = "clusters", out.name = "spatial.png")
#' spatial.plot(obj, col.by = "gene", gene = "CD3E",
#'              xlim = c(2000, 4000), ylim = c(2000, 4000))
#' }
#' @export
spatial.plot <- function(
    x = NULL,
    col.by = "clusters",
    gene = NULL,
    conds.to.plot = NULL,
    max.spots = 200000L,
    mark = c("auto", "point", "raster"),
    xlim = NULL,
    ylim = NULL,
    cell.size = 0.6,
    cell.transparency = 0.9,
    cell.colors = c("grey80", "#c23b22"),
    back.col = "white",
    interactive = FALSE,
    out.name = NULL,
    assay.name = "data",
    scaleValue = FALSE,
    min.scale = NA_real_,
    max.scale = NA_real_,
    anno.clust = TRUE,
    anno.size = 3,
    anno.col = "black",
    data.type = NULL
) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  cd <- x@cell_data
  if (!all(c("spatial_x", "spatial_y") %in% colnames(cd))) {
    stop("No spatial coordinates. Run add.visium() or add.xenium() first.")
  }
  ok <- is.finite(cd$spatial_x) & is.finite(cd$spatial_y)
  if ("in_tissue" %in% colnames(cd)) ok <- ok & (is.na(cd$in_tissue) | cd$in_tissue == 1L)
  cd <- cd[ok, , drop = FALSE]
  if (!is.null(conds.to.plot)) {
    key <- if ("spatial_sample" %in% colnames(cd)) cd$spatial_sample else cd$condition_id
    cd <- cd[key %in% conds.to.plot, , drop = FALSE]
  }
  if (nrow(cd) < 1L) stop("No spots left to plot.")
  mark <- match.arg(mark)
  if (identical(mark, "auto")) {
    plat <- if ("spatial_platform" %in% colnames(cd)) cd$spatial_platform else NA_character_
    mark <- if (any(plat == "visium_hd", na.rm = TRUE)) "raster" else "point"
  }
  laid <- .ic2_spatial_lay_on_image(x, cd, hires = !is.null(xlim) || !is.null(ylim))
  cd <- laid$cd
  if (!is.null(xlim)) cd <- cd[cd$spatial_x >= xlim[[1L]] & cd$spatial_x <= xlim[[2L]], , drop = FALSE]
  if (!is.null(ylim)) cd <- cd[cd$spatial_y >= ylim[[1L]] & cd$spatial_y <= ylim[[2L]], , drop = FALSE]
  if (nrow(cd) < 1L) stop("No spots left in the crop.")
  color <- .ic2_spatial_color(x, cd, col.by, gene, assay.name)
  if (isTRUE(scaleValue) && is.numeric(color)) {
    if (is.finite(min.scale)) color[color < min.scale] <- min.scale
    if (is.finite(max.scale)) color[color > max.scale] <- max.scale
  }
  max.spots <- as.integer(max.spots)
  if (identical(mark, "raster")) {
    df <- .ic2_spatial_bin_raster(cd$spatial_x, cd$spatial_y, color, max.spots)
  } else {
    idx <- seq_len(nrow(cd))
    if (nrow(cd) > max.spots) idx <- sort(sample.int(nrow(cd), max.spots))
    df <- data.frame(
      x = cd$spatial_x[idx],
      y = cd$spatial_y[idx],
      color = if (is.numeric(color)) color[idx] else as.character(color[idx]),
      stringsAsFactors = FALSE
    )
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = back.col, colour = NA))
  if (!is.null(laid$image)) {
    im <- laid$image
    p <- p + ggplot2::annotation_raster(
      im, xmin = 0, xmax = ncol(im), ymin = nrow(im), ymax = 0
    )
  } else {
    p <- p + ggplot2::scale_y_reverse()
  }
  if (identical(mark, "raster")) {
    aes_fill <- ggplot2::aes(fill = color)
    p <- p + ggplot2::geom_tile(aes_fill, width = df$width[[1L]], height = df$height[[1L]], alpha = cell.transparency)
    if (is.numeric(df$color)) {
      p <- p + ggplot2::scale_fill_gradient(low = cell.colors[[1L]], high = cell.colors[[min(2L, length(cell.colors))]])
    }
  } else if (is.numeric(df$color)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(colour = color), size = cell.size, alpha = cell.transparency
    ) + ggplot2::scale_colour_gradient(low = cell.colors[[1L]], high = cell.colors[[min(2L, length(cell.colors))]])
  } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(colour = color), size = cell.size, alpha = cell.transparency
    )
  }
  if (isTRUE(anno.clust) && identical(col.by, "clusters") && "cluster_id" %in% colnames(cd)) {
    lab <- cd
    if (nrow(lab) > max.spots) lab <- lab[sample.int(nrow(lab), max.spots), , drop = FALSE]
    cents <- stats::aggregate(cbind(spatial_x, spatial_y) ~ cluster_id, data = lab, FUN = mean)
    p <- p + ggplot2::geom_text(
      data = cents,
      ggplot2::aes(x = spatial_x, y = spatial_y, label = cluster_id),
      inherit.aes = FALSE, size = anno.size, colour = anno.col
    )
  }
  if (!is.null(out.name) && !isTRUE(interactive)) {
    ggplot2::ggsave(out.name, p, width = 7, height = 7, dpi = 150)
  }
  if (isTRUE(interactive)) {
    if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
      stop("Packages 'plotly' and 'htmlwidgets' are required for interactive = TRUE.")
    }
    html <- if (is.null(out.name)) "spatial.plot.html" else {
      if (grepl("\\.html$", out.name, ignore.case = TRUE)) out.name else paste0(out.name, ".html")
    }
    htmlwidgets::saveWidget(plotly::ggplotly(p), html)
    return(invisible(html))
  }
  p
}

#' Scale full-resolution pixels onto the low-resolution tissue image.
#' @keywords internal
.ic2_spatial_lay_on_image <- function(x, cd, hires = FALSE) {
  samples <- if ("spatial_sample" %in% colnames(cd)) unique(stats::na.omit(cd$spatial_sample)) else character()
  plat <- if ("spatial_platform" %in% colnames(cd)) unique(stats::na.omit(as.character(cd$spatial_platform))) else character()
  if (length(samples) != 1L) return(list(cd = cd, image = NULL))
  if (any(plat == "xenium")) return(list(cd = cd, image = NULL))
  root <- x@storage$root_dir
  dir <- file.path(root, "images", samples[[1L]])
  js <- file.path(dir, "scalefactors.json")
  sf <- 1
  if (file.exists(js)) {
    sc <- jsonlite::fromJSON(js, simplifyVector = TRUE)
    img_names <- list.files(dir, pattern = "tissue_.*image\\.png$", full.names = TRUE)
    use_hi <- isTRUE(hires) || (length(img_names) == 1L && grepl("hires", basename(img_names[[1L]])))
    sf <- if (use_hi) sc$tissue_hires_scalef %||% sc$tissue_lowres_scalef else {
      sc$tissue_lowres_scalef %||% sc$tissue_hires_scalef
    }
    if (is.null(sf) || !is.finite(sf)) sf <- 1
  }
  cd$spatial_x <- cd$spatial_x * sf
  cd$spatial_y <- cd$spatial_y * sf
  img <- NULL
  pics <- if (isTRUE(hires)) {
    list.files(dir, pattern = "tissue_hires_image\\.png$", full.names = TRUE)
  } else {
    character(0)
  }
  if (length(pics) < 1L) pics <- list.files(dir, pattern = "tissue_lowres_image\\.png$", full.names = TRUE)
  if (length(pics) < 1L) pics <- list.files(dir, pattern = "tissue_.*image\\.png$", full.names = TRUE)
  if (length(pics) >= 1L && requireNamespace("png", quietly = TRUE)) {
    img <- png::readPNG(pics[[1L]])
  }
  list(cd = cd, image = img)
}

#' Bin spots into at most \code{max.bins} tiles. Numeric colors are means.
#' @keywords internal
.ic2_spatial_bin_raster <- function(x, y, color, max.bins) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  color <- color[ok]
  xr <- range(x)
  yr <- range(y)
  if (diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
  if (diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
  nx <- max(1L, min(1200L, as.integer(sqrt(max.bins))))
  ny <- nx
  ix <- pmax(1L, pmin(nx, as.integer((x - xr[[1L]]) / diff(xr) * (nx - 1L)) + 1L))
  iy <- pmax(1L, pmin(ny, as.integer((y - yr[[1L]]) / diff(yr) * (ny - 1L)) + 1L))
  key <- paste(ix, iy, sep = ",")
  dx <- diff(xr) / nx
  dy <- diff(yr) / ny
  if (is.numeric(color)) {
    mean_c <- tapply(color, key, mean)
    data.frame(
      x = xr[[1L]] + (as.integer(sub(",.*", "", names(mean_c))) - 0.5) * dx,
      y = yr[[1L]] + (as.integer(sub(".*,", "", names(mean_c))) - 0.5) * dy,
      color = as.numeric(mean_c),
      width = dx, height = dy,
      stringsAsFactors = FALSE
    )
  } else {
    mode_c <- tapply(as.character(color), key, function(z) {
      tt <- tabulate(match(z, unique(z)))
      unique(z)[which.max(tt)]
    })
    data.frame(
      x = xr[[1L]] + (as.integer(sub(",.*", "", names(mode_c))) - 0.5) * dx,
      y = yr[[1L]] + (as.integer(sub(".*,", "", names(mode_c))) - 0.5) * dy,
      color = as.character(mode_c),
      width = dx, height = dy,
      stringsAsFactors = FALSE
    )
  }
}

#' @keywords internal
.ic2_spatial_color <- function(x, cd, col.by, gene, assay.name) {
  if (identical(col.by, "gene")) {
    if (is.null(gene) || !nzchar(as.character(gene)[1L])) stop("Provide gene when col.by = 'gene'.")
    vals <- .ic2_spatial_gene(x, as.character(gene)[1L], assay.name)
    ids <- if ("cell_id" %in% colnames(cd)) as.character(cd$cell_id) else rownames(cd)
    out <- vals[ids]
    out[is.na(out)] <- 0
    return(as.numeric(out))
  }
  if (identical(col.by, "cc")) {
    if (!"Phase" %in% colnames(cd)) stop("No Phase column. Run cc() first.")
    return(as.character(cd$Phase))
  }
  if (identical(col.by, "nGenes")) {
    col <- intersect(c("nFeature_RNA", "n_genes", "nGenes"), colnames(cd))
    if (length(col) < 1L) stop("No gene-count column on cell_data.")
    return(as.numeric(cd[[col[[1L]]]]))
  }
  if (identical(col.by, "UMIs")) {
    col <- intersect(c("nCount_RNA", "n_umis", "nUMI", "UMIs"), colnames(cd))
    if (length(col) < 1L) stop("No UMI column on cell_data.")
    return(as.numeric(cd[[col[[1L]]]]))
  }
  if (identical(col.by, "mt")) {
    col <- intersect(c("percent_mito", "percent.mt", "mito"), colnames(cd))
    if (length(col) < 1L) stop("No mito column on cell_data.")
    return(as.numeric(cd[[col[[1L]]]]))
  }
  if (!"cluster_id" %in% colnames(cd)) stop("No cluster_id. Run clustering or set col.by.")
  as.character(cd$cluster_id)
}

#' One gene, chunked. Named by cell id.
#' @keywords internal
.ic2_spatial_gene <- function(x, gene, assay.name) {
  handles <- .ic2_get_assay_handles(x, "rna", assay.name)
  vals <- numeric(0)
  ids <- character(0)
  found <- FALSE
  chunk <- 50000L
  for (h in handles) {
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    ri <- match(tolower(gene), tolower(h$feature_ids))
    if (is.na(ri)) next
    found <- TRUE
    starts <- seq.int(1L, h$n_cells, by = chunk)
    for (s in starts) {
      e <- min(s + chunk - 1L, h$n_cells)
      block <- .ic2_read_cols(h, s:e)
      vals <- c(vals, as.numeric(block[ri, ]))
      ids <- c(ids, h$cell_ids[s:e])
    }
  }
  if (!found) stop("Gene not in assay '", assay.name, "': ", gene)
  names(vals) <- ids
  vals
}
