#' Pseudotime stream / arrows on a 2D embedding
#'
#' Overlay a velocity-style \strong{stream}, straight \strong{arrows}, or
#' \strong{both} side by side on a 2D embedding. Follows the local increase
#' in Pseudotime (UMAP or another reduction). No RNA velocity required.
#'
#' @param x An \code{iCellR2} object after \code{\link{run.pseudotime}}.
#' @param modality Modality. Default \code{"rna"}.
#' @param plot.type Embedding for the canvas. Default \code{"umap"}.
#' @param dims Two dimensions. Default \code{1:2}.
#' @param col.by Point colors (same as \code{\link{cluster.plot}}).
#'   Default \code{"clusters"}.
#' @param pt.col Metadata column for the gradient. Default \code{"Pseudotime"}.
#' @param style \code{"both"} (stream | grid side by side, default),
#'   \code{"stream"} (curvy paths), or \code{"arrows"} (grid quiver).
#' @param grid.n Bins per axis for the vector field. Default \code{50}
#'   for streams, \code{45} for arrows when left \code{NULL}.
#' @param min.cells Minimum cells in a neighborhood to estimate a vector.
#'   Default \code{12}.
#' @param smooth Neighborhood radius in bins for the local fit. Default \code{1}.
#' @param field.smooth Extra passes of neighbor averaging on the vector
#'   field (makes streams curvier / less noisy). Default \code{3}.
#' @param n.streams Approximate number of streamlines (style = stream).
#'   Default \code{220}. Seeds prefer denser bins.
#' @param stream.step Integration step as a fraction of axis span.
#'   Default \code{0.006}.
#' @param stream.steps Max forward integration steps per stream. Default \code{55}.
#' @param stream.backward Backward steps from each seed (adds origin -> flow look).
#'   Default \code{20}.
#' @param stream.width Base line width for streams. Default \code{0.25}.
#' @param stream.taper If \code{TRUE}, fade / thin streams toward the tip.
#'   Default \code{TRUE}.
#' @param arrow.scale Multiplier on straight-arrow length (\code{style = "arrows"}).
#'   Default \code{1}.
#' @param arrow.color Stream / arrow colour. Default \code{"grey10"}.
#' @param arrow.alpha Stream / arrow alpha. Default \code{0.85}.
#' @param arrow.size Line width for \code{style = "arrows"}. Default \code{0.35}.
#' @param arrow.every Place triangular heads every N path steps along each
#'   streamline (scVelo streamplot style). Default \code{8}. Set \code{0}
#'   to disable heads. Legacy alias: \code{tip.arrows = FALSE} turns heads off.
#' @param tip.arrows Deprecated; use \code{arrow.every}. If \code{FALSE},
#'   disables heads.
#' @param arrow.head Size of mid-stream arrowheads in cm. Default \code{0.08}.
#' @param cell.size,cell.transparency,cell.color,back.col,max.cells,anno.clust,anno.size,seed,verbose
#'   Passed through like \code{\link{cluster.plot}}. Cluster labels skip
#'   groups with fewer than 10 plotted cells.
#'
#' @return A \code{ggplot} (single style) or a printable arranged grob
#'   when \code{style = "both"}.
#'
#' @examples
#' \dontrun{
#' obj <- run.pseudotime(obj, root.cluster = 1)
#' pseudotime.arrows(obj)                       # stream | grid side by side
#' pseudotime.arrows(obj, style = "stream")
#' pseudotime.arrows(obj, style = "arrows")
#' }
#'
#' @importFrom rlang .data
#' @export
pseudotime.arrows <- function(
    x,
    modality = "rna",
    plot.type = c("umap", "pca", "cgmap", "cgpca", "cgumap", "ptumap"),
    dims = 1:2,
    col.by = "clusters",
    pt.col = "Pseudotime",
    style = c("both", "stream", "arrows"),
    grid.n = NULL,
    min.cells = 12L,
    smooth = 1L,
    field.smooth = 3L,
    n.streams = 220L,
    stream.step = 0.006,
    stream.steps = 55L,
    stream.backward = 20L,
    stream.width = 0.25,
    stream.taper = TRUE,
    arrow.scale = 1,
    arrow.color = "grey10",
    arrow.alpha = 0.85,
    arrow.size = 0.35,
    arrow.every = 8L,
    arrow.head = 0.08,
    tip.arrows = TRUE,
    cell.size = 1,
    cell.transparency = 0.45,
    cell.color = "grey20",
    back.col = "white",
    max.cells = 500000L,
    anno.clust = TRUE,
    anno.size = 4,
    seed = 1,
    verbose = TRUE
) {
  .ic2_require_object(x)
  plot.type <- match.arg(plot.type)
  style <- match.arg(style)
  modality <- tolower(modality)
  dims <- as.integer(dims)
  if (length(dims) != 2L) stop("'dims' must have length 2 for arrow plots.")

  ## Side-by-side: stream + short grid arrows
  if (identical(style, "both")) {
    rep <- .ic2_new_reporter("pseudotime.arrows", verbose, soda_root = x@storage$root_dir)
    .ic2_step(rep, "[1/2] Stream panel")
    p_stream <- pseudotime.arrows(
      x,
      modality = modality,
      plot.type = plot.type,
      dims = dims,
      col.by = col.by,
      pt.col = pt.col,
      style = "stream",
      grid.n = grid.n,
      min.cells = min.cells,
      smooth = smooth,
      field.smooth = field.smooth,
      n.streams = n.streams,
      stream.step = stream.step,
      stream.steps = stream.steps,
      stream.backward = stream.backward,
      stream.width = stream.width,
      stream.taper = stream.taper,
      arrow.scale = arrow.scale,
      arrow.color = arrow.color,
      arrow.alpha = arrow.alpha,
      arrow.size = arrow.size,
      arrow.every = arrow.every,
      arrow.head = arrow.head,
      tip.arrows = tip.arrows,
      cell.size = cell.size,
      cell.transparency = cell.transparency,
      cell.color = cell.color,
      back.col = back.col,
      max.cells = max.cells,
      anno.clust = anno.clust,
      anno.size = anno.size,
      seed = seed,
      verbose = FALSE
    )
    .ic2_step(rep, "[2/2] Grid-arrow panel")
    p_arrows <- pseudotime.arrows(
      x,
      modality = modality,
      plot.type = plot.type,
      dims = dims,
      col.by = col.by,
      pt.col = pt.col,
      style = "arrows",
      grid.n = grid.n,
      min.cells = min.cells,
      smooth = smooth,
      field.smooth = field.smooth,
      n.streams = n.streams,
      stream.step = stream.step,
      stream.steps = stream.steps,
      stream.backward = stream.backward,
      stream.width = stream.width,
      stream.taper = stream.taper,
      arrow.scale = arrow.scale,
      arrow.color = arrow.color,
      arrow.alpha = arrow.alpha,
      arrow.size = arrow.size,
      arrow.every = arrow.every,
      arrow.head = arrow.head,
      tip.arrows = tip.arrows,
      cell.size = cell.size,
      cell.transparency = cell.transparency,
      cell.color = cell.color,
      back.col = back.col,
      max.cells = max.cells,
      anno.clust = anno.clust,
      anno.size = anno.size,
      seed = seed,
      verbose = FALSE
    )
    p_stream <- p_stream + ggplot2::labs(title = paste0(toupper(plot.type), " + stream"))
    p_arrows <- p_arrows + ggplot2::labs(title = paste0(toupper(plot.type), " + grid arrows"))
    g <- gridExtra::arrangeGrob(p_stream, p_arrows, ncol = 2L)
    class(g) <- c("pseudotime_arrows_both", class(g))
    rep <- .ic2_banner_start(rep, c(
      paste0("Reduction: ", plot.type, "   |   style=both   |   ", pt.col),
      "Panels: stream | grid arrows"
    ))
    .ic2_finish(rep, extra_lines = "arranged grob   |   style=both")
    return(invisible(g))
  }

  if (is.null(grid.n)) grid.n <- if (identical(style, "stream")) 50L else 45L
  grid.n <- as.integer(grid.n)[1L]
  min.cells <- as.integer(min.cells)[1L]
  smooth <- as.integer(smooth)[1L]
  field.smooth <- as.integer(field.smooth)[1L]
  n.streams <- as.integer(n.streams)[1L]
  stream.steps <- as.integer(stream.steps)[1L]
  stream.backward <- as.integer(stream.backward)[1L]
  arrow.every <- as.integer(arrow.every)[1L]
  arrow.head <- as.numeric(arrow.head)[1L]
  if (isFALSE(tip.arrows)) arrow.every <- 0L
  if (grid.n < 5L) stop("'grid.n' must be >= 5.")
  if (min.cells < 3L) stop("'min.cells' must be >= 3.")
  if (smooth < 0L) stop("'smooth' must be >= 0.")
  if (field.smooth < 0L) stop("'field.smooth' must be >= 0.")
  pt.col <- as.character(pt.col)[1L]
  arrow.scale <- as.numeric(arrow.scale)[1L]
  stream.step <- as.numeric(stream.step)[1L]
  stream.width <- as.numeric(stream.width)[1L]
  if (!is.finite(arrow.scale) || arrow.scale <= 0) {
    stop("'arrow.scale' must be a positive number.")
  }
  if (!is.finite(stream.step) || stream.step <= 0) {
    stop("'stream.step' must be a positive number.")
  }

  rep <- .ic2_new_reporter("pseudotime.arrows", verbose, soda_root = x@storage$root_dir)

  if (!pt.col %in% colnames(x@cell_data)) {
    stop("No '", pt.col, "' in cell_data. Run run.pseudotime() first.")
  }

  .ic2_step(rep, paste0("[1/4] Loading ", plot.type, " + ", pt.col))
  emb_full <- .ic2_get_plot_embedding(x, modality, plot.type)
  if (max(dims) > ncol(emb_full)) {
    stop(
      "Requested dims up to ", max(dims),
      " but ", plot.type, " only has ", ncol(emb_full), " columns."
    )
  }
  meta <- x@cell_data
  common <- intersect(rownames(emb_full), rownames(meta))
  if (length(common) < 10L) stop("Embedding rownames do not match cell_data.")
  emb_full <- emb_full[common, dims, drop = FALSE]
  meta <- meta[common, , drop = FALSE]
  pt_all <- as.numeric(meta[[pt.col]])
  ok <- is.finite(pt_all) & is.finite(emb_full[, 1L]) & is.finite(emb_full[, 2L])
  if (sum(ok) < 10L) stop("Too few finite Pseudotime / embedding values.")
  emb_full <- emb_full[ok, , drop = FALSE]
  meta <- meta[ok, , drop = FALSE]
  pt_all <- pt_all[ok]
  storage.mode(emb_full) <- "double"
  colnames(emb_full) <- c("Dim1", "Dim2")

  .ic2_step(rep, "[2/4] Estimating smoothed Pseudotime vector field")
  keep <- seq_len(nrow(emb_full))
  if (!is.null(max.cells) && is.finite(max.cells) && nrow(emb_full) > as.integer(max.cells)) {
    set.seed(seed)
    keep <- sort(sample.int(nrow(emb_full), as.integer(max.cells)))
    .ic2_step(rep, paste0(
      "  Downsampled ", .ic2_fmt_num(length(keep)),
      " / ", .ic2_fmt_num(nrow(emb_full)), " cells"
    ))
  }
  emb <- emb_full[keep, , drop = FALSE]
  meta_s <- meta[keep, , drop = FALSE]
  pt <- pt_all[keep]

  xr <- range(emb[, 1L], na.rm = TRUE)
  yr <- range(emb[, 2L], na.rm = TRUE)
  xspan <- diff(xr); yspan <- diff(yr)
  if (!is.finite(xspan) || xspan <= 0) xspan <- 1
  if (!is.finite(yspan) || yspan <= 0) yspan <- 1
  pad <- 1e-8
  gx <- pmin(grid.n, pmax(1L, as.integer(floor((emb[, 1L] - xr[1L]) / (xspan + pad) * grid.n)) + 1L))
  gy <- pmin(grid.n, pmax(1L, as.integer(floor((emb[, 2L] - yr[1L]) / (yspan + pad) * grid.n)) + 1L))
  bin_id <- gx + (gy - 1L) * grid.n
  bin_list <- split(seq_along(bin_id), bin_id)

  ## Unit-direction field on grid + density mask
  vx <- matrix(NA_real_, grid.n, grid.n)
  vy <- matrix(NA_real_, grid.n, grid.n)
  dens <- matrix(0L, grid.n, grid.n)
  for (ix in seq_len(grid.n)) {
    for (iy in seq_len(grid.n)) {
      ids <- integer()
      for (dx in (-smooth):smooth) {
        for (dy in (-smooth):smooth) {
          jx <- ix + dx; jy <- iy + dy
          if (jx < 1L || jx > grid.n || jy < 1L || jy > grid.n) next
          key <- as.character(jx + (jy - 1L) * grid.n)
          if (!is.null(bin_list[[key]])) ids <- c(ids, bin_list[[key]])
        }
      }
      center_ids0 <- bin_list[[as.character(ix + (iy - 1L) * grid.n)]]
      dens[ix, iy] <- if (is.null(center_ids0)) 0L else length(center_ids0)
      if (length(ids) < min.cells) next
      X <- cbind(1, emb[ids, 1L], emb[ids, 2L])
      fit <- tryCatch(stats::lm.fit(X, pt[ids]), error = function(e) NULL)
      if (is.null(fit) || anyNA(fit$coefficients)) next
      ax <- as.numeric(fit$coefficients[2L])
      ay <- as.numeric(fit$coefficients[3L])
      if (!is.finite(ax) || !is.finite(ay)) next
      nrm <- sqrt(ax * ax + ay * ay)
      if (!is.finite(nrm) || nrm < .Machine$double.eps) next
      vx[ix, iy] <- ax / nrm
      vy[ix, iy] <- ay / nrm
    }
  }

  ## Neighbor averaging (keeps directions unit-ish)
  if (field.smooth > 0L) {
    for (pass in seq_len(field.smooth)) {
      vx2 <- vx; vy2 <- vy
      for (ix in seq_len(grid.n)) {
        for (iy in seq_len(grid.n)) {
          xs <- max(1L, ix - 1L):min(grid.n, ix + 1L)
          ys <- max(1L, iy - 1L):min(grid.n, iy + 1L)
          block_x <- vx[xs, ys, drop = FALSE]
          block_y <- vy[xs, ys, drop = FALSE]
          okb <- is.finite(block_x) & is.finite(block_y)
          if (!any(okb)) next
          sx <- mean(block_x[okb]); sy <- mean(block_y[okb])
          nrm <- sqrt(sx * sx + sy * sy)
          if (!is.finite(nrm) || nrm < .Machine$double.eps) next
          vx2[ix, iy] <- sx / nrm
          vy2[ix, iy] <- sy / nrm
        }
      }
      vx <- vx2; vy <- vy2
    }
  }

  dens_cut <- max(1L, as.integer(stats::quantile(dens[dens > 0], probs = 0.15, names = FALSE)))
  mask <- dens >= dens_cut & is.finite(vx) & is.finite(vy)

  ## Grid center coordinates
  xc <- xr[1L] + (seq_len(grid.n) - 0.5) / grid.n * xspan
  yc <- yr[1L] + (seq_len(grid.n) - 0.5) / grid.n * yspan

  .ic2_interp_field <- function(px, py) {
    ## bilinear sample of unit field; NA outside mask
    if (px < xr[1L] || px > xr[2L] || py < yr[1L] || py > yr[2L]) {
      return(c(NA_real_, NA_real_))
    }
    fx <- (px - xr[1L]) / (xspan + pad) * grid.n
    fy <- (py - yr[1L]) / (yspan + pad) * grid.n
    i0 <- max(1L, min(grid.n - 1L, as.integer(floor(fx))))
    j0 <- max(1L, min(grid.n - 1L, as.integer(floor(fy))))
    i1 <- i0 + 1L; j1 <- j0 + 1L
    if (!mask[i0, j0] && !mask[i1, j0] && !mask[i0, j1] && !mask[i1, j1]) {
      return(c(NA_real_, NA_real_))
    }
    tx <- fx - i0; ty <- fy - j0
    ## Pull NA vectors as 0 for interpolation, then renormalize
    grab <- function(i, j) {
      if (!mask[i, j]) return(c(0, 0, 0))
      c(vx[i, j], vy[i, j], 1)
    }
    a <- grab(i0, j0); b <- grab(i1, j0); c0 <- grab(i0, j1); d <- grab(i1, j1)
    w00 <- (1 - tx) * (1 - ty); w10 <- tx * (1 - ty); w01 <- (1 - tx) * ty; w11 <- tx * ty
    wx <- w00 * a[1] + w10 * b[1] + w01 * c0[1] + w11 * d[1]
    wy <- w00 * a[2] + w10 * b[2] + w01 * c0[2] + w11 * d[2]
    ww <- w00 * a[3] + w10 * b[3] + w01 * c0[3] + w11 * d[3]
    if (!is.finite(ww) || ww < 1e-6) return(c(NA_real_, NA_real_))
    nrm <- sqrt(wx * wx + wy * wy)
    if (!is.finite(nrm) || nrm < .Machine$double.eps) return(c(NA_real_, NA_real_))
    c(wx / nrm, wy / nrm)
  }

  stream_df <- NULL
  tip_df <- NULL
  arr <- NULL

  if (identical(style, "arrows")) {
    .ic2_step(rep, "[3/4] Building straight arrows")
    step <- 0.55 * min(xspan, yspan) / grid.n * arrow.scale
    rows <- list(); n_arr <- 0L
    for (ix in seq_len(grid.n)) {
      for (iy in seq_len(grid.n)) {
        if (!mask[ix, iy]) next
        n_arr <- n_arr + 1L
        rows[[n_arr]] <- c(
          x = xc[ix], y = yc[iy],
          xend = xc[ix] + vx[ix, iy] * step,
          yend = yc[iy] + vy[ix, iy] * step
        )
      }
    }
    if (n_arr < 1L) stop("No arrows could be estimated.")
    arr <- as.data.frame(do.call(rbind, rows))
    colnames(arr) <- c("x", "y", "xend", "yend")
  } else {
    .ic2_step(rep, "[3/4] Integrating curvy streamlines")
    set.seed(seed)
    ## Density-aware seeds: prefer occupied bins with more cells
    seed_ix <- which(mask, arr.ind = TRUE)
    if (nrow(seed_ix) < 1L) stop("No vector field  -  cannot draw streams.")
    seed_w <- dens[seed_ix]
    seed_w[!is.finite(seed_w) | seed_w < 1] <- 1
    seed_w <- as.numeric(seed_w)
    n_seed <- max(1L, min(as.integer(n.streams), max(1L, nrow(seed_ix))))
    pick <- sample.int(
      nrow(seed_ix), n_seed,
      replace = n_seed > nrow(seed_ix),
      prob = seed_w
    )
    step_len <- stream.step * min(xspan, yspan)

    paths <- vector("list", n_seed)
    tips <- list()
    n_path <- 0L
    for (s in seq_len(n_seed)) {
      ix <- seed_ix[pick[s], 1L]; iy <- seed_ix[pick[s], 2L]
      x0 <- xc[ix] + stats::runif(1, -0.4, 0.4) * (xspan / grid.n)
      y0 <- yc[iy] + stats::runif(1, -0.4, 0.4) * (yspan / grid.n)

      ## integrate backward then forward so stream covers both sides of seed
      trace <- function(x, y, steps, sign = 1) {
        xs <- x; ys <- y
        for (k in seq_len(steps)) {
          vv <- .ic2_interp_field(xs[length(xs)], ys[length(ys)])
          if (anyNA(vv)) break
          ## RK2
          mid <- .ic2_interp_field(
            xs[length(xs)] + sign * 0.5 * step_len * vv[1],
            ys[length(ys)] + sign * 0.5 * step_len * vv[2]
          )
          if (anyNA(mid)) {
            xs <- c(xs, xs[length(xs)] + sign * step_len * vv[1])
            ys <- c(ys, ys[length(ys)] + sign * step_len * vv[2])
            break
          }
          xs <- c(xs, xs[length(xs)] + sign * step_len * mid[1])
          ys <- c(ys, ys[length(ys)] + sign * step_len * mid[2])
        }
        list(x = xs, y = ys)
      }
      back <- if (stream.backward > 0L) trace(x0, y0, stream.backward, sign = -1) else list(x = x0, y = y0)
      fwd <- trace(x0, y0, stream.steps, sign = 1)
      ## stitch: reverse back (excluding seed) + forward
      bx <- rev(back$x); by <- rev(back$y)
      if (length(bx) > 1L) {
        bx <- bx[-length(bx)]; by <- by[-length(by)]
      } else {
        bx <- numeric(); by <- numeric()
      }
      xs <- c(bx, fwd$x); ys <- c(by, fwd$y)
      if (length(xs) < 3L) next
      n_path <- n_path + 1L
      t01 <- seq(0, 1, length.out = length(xs))
      paths[[n_path]] <- data.frame(
        x = xs, y = ys, id = n_path, t = t01,
        stringsAsFactors = FALSE
      )
      ## scVelo-like: arrowheads at regular intervals along the whole stream
      if (arrow.every > 0L && length(xs) >= (arrow.every + 2L)) {
        idx <- seq.int(arrow.every, length(xs) - 1L, by = arrow.every)
        for (ii in idx) {
          i0 <- max(1L, ii - 2L)
          tips[[length(tips) + 1L]] <- c(
            x = xs[i0], y = ys[i0],
            xend = xs[ii], yend = ys[ii]
          )
        }
      }
    }
    if (n_path < 1L) stop("Failed to integrate streamlines.")
    stream_df <- do.call(rbind, paths[seq_len(n_path)])
    if (length(tips)) {
      tip_df <- as.data.frame(do.call(rbind, tips))
      colnames(tip_df) <- c("x", "y", "xend", "yend")
    } else {
      tip_df <- NULL
    }
    .ic2_step(rep, paste0(
      "  Streams: ", .ic2_fmt_num(n_path),
      if (is.null(tip_df)) "" else paste0("   |   heads: ", .ic2_fmt_num(nrow(tip_df)))
    ))
  }

  .ic2_step(rep, "[4/4] Plotting")
  if (identical(col.by, "clusters") && !"cluster_id" %in% colnames(meta_s)) {
    if ("condition_id" %in% colnames(meta_s)) col.by <- "conditions" else col.by <- "monochrome"
  }
  color_info <- .ic2_plot_color_vector(meta_s, col.by = col.by)
  df <- data.frame(
    Dim1 = emb[, 1L], Dim2 = emb[, 2L], .color = color_info$values,
    stringsAsFactors = FALSE
  )
  ## Grid arrows default to solid arrow.color (black). Set
  ## options(iCellR2.pt.arrows.cluster = TRUE) to colour by clusters.
  color_by_cluster <- isTRUE(getOption("iCellR2.pt.arrows.cluster", FALSE)) &&
    isTRUE(color_info$discrete) && !identical(col.by, "monochrome")

  ## Map each stream / arrowhead to the nearest cell's cluster colour
  .ic2_nearest_color <- function(px, py) {
    ## bin lookup first (fast), else nearest of a small sample
    ix <- pmin(grid.n, pmax(1L, as.integer(floor((px - xr[1L]) / (xspan + pad) * grid.n)) + 1L))
    iy <- pmin(grid.n, pmax(1L, as.integer(floor((py - yr[1L]) / (yspan + pad) * grid.n)) + 1L))
    key <- as.character(ix + (iy - 1L) * grid.n)
    ids <- bin_list[[key]]
    if (is.null(ids) || length(ids) < 1L) {
      ## search 3x3 neighborhood of bins
      ids <- integer()
      for (dx in -1:1) for (dy in -1:1) {
        jx <- ix + dx; jy <- iy + dy
        if (jx < 1L || jx > grid.n || jy < 1L || jy > grid.n) next
        kk <- as.character(jx + (jy - 1L) * grid.n)
        if (!is.null(bin_list[[kk]])) ids <- c(ids, bin_list[[kk]])
      }
    }
    if (length(ids) < 1L) {
      ## fallback: nearest among up to 2k cells
      n_s <- min(2000L, nrow(emb))
      set.seed(seed)
      samp <- if (nrow(emb) > n_s) sort(sample.int(nrow(emb), n_s)) else seq_len(nrow(emb))
      d2 <- (emb[samp, 1L] - px)^2 + (emb[samp, 2L] - py)^2
      return(df$.color[samp[which.min(d2)]])
    }
    ## majority colour in the bin
    vv <- df$.color[ids]
    tb <- sort(table(vv), decreasing = TRUE)
    vv_levels <- names(tb)[1L]
    if (is.factor(df$.color)) {
      factor(vv_levels, levels = levels(df$.color))
    } else {
      vv_levels
    }
  }

  if (color_by_cluster && identical(style, "stream") && !is.null(stream_df)) {
    ## one colour per stream id (from stream midpoint)
    mids <- stats::aggregate(
      cbind(x, y) ~ id, data = stream_df, FUN = mean, na.rm = TRUE
    )
    mids$.color <- vector("list", nrow(mids))
    for (i in seq_len(nrow(mids))) {
      mids$.color[[i]] <- .ic2_nearest_color(mids$x[i], mids$y[i])
    }
    mids$.color <- do.call(c, mids$.color)
    if (is.factor(df$.color)) {
      mids$.color <- factor(as.character(mids$.color), levels = levels(df$.color))
    }
    stream_df$.color <- mids$.color[match(stream_df$id, mids$id)]
    if (!is.null(tip_df) && nrow(tip_df) > 0L) {
      tip_df$.color <- vector("list", nrow(tip_df))
      for (i in seq_len(nrow(tip_df))) {
        tip_df$.color[[i]] <- .ic2_nearest_color(
          (tip_df$x[i] + tip_df$xend[i]) / 2,
          (tip_df$y[i] + tip_df$yend[i]) / 2
        )
      }
      tip_df$.color <- do.call(c, tip_df$.color)
      if (is.factor(df$.color)) {
        tip_df$.color <- factor(as.character(tip_df$.color), levels = levels(df$.color))
      }
    }
  }
  if (color_by_cluster && identical(style, "arrows") && !is.null(arr)) {
    arr$.color <- vector("list", nrow(arr))
    for (i in seq_len(nrow(arr))) {
      arr$.color[[i]] <- .ic2_nearest_color(arr$x[i], arr$y[i])
    }
    arr$.color <- do.call(c, arr$.color)
    if (is.factor(df$.color)) {
      arr$.color <- factor(as.character(arr$.color), levels = levels(df$.color))
    }
  }

  title <- paste0(toupper(plot.type), " + ", pt.col, " ", style)
  xlab <- paste0(plot.type, dims[1L]); ylab <- paste0(plot.type, dims[2L])

  if (identical(col.by, "monochrome")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2)) +
      ggplot2::geom_point(
        size = cell.size, alpha = cell.transparency, colour = cell.color, stroke = 0
      )
  } else if (isTRUE(color_info$discrete)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2, colour = .data$.color)) +
      ggplot2::geom_point(
        size = cell.size, alpha = cell.transparency, stroke = 0
      ) +
      ggplot2::scale_colour_discrete(name = color_info$title) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 5.5, alpha = 1, stroke = 0),
        title.position = "top"
      ))
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2, colour = .data$.color)) +
      ggplot2::geom_point(
        size = cell.size, alpha = cell.transparency, stroke = 0
      ) +
      ggplot2::scale_colour_viridis_c(name = color_info$title, option = "C") +
      ggplot2::guides(colour = ggplot2::guide_colorbar(
        barwidth = 0.55, barheight = 6, ticks = FALSE, frame.colour = NA, title.position = "top"
      ))
  }

  if (identical(style, "arrows")) {
    if (isTRUE(color_by_cluster) && ".color" %in% names(arr)) {
      p <- p + ggplot2::geom_segment(
        data = arr,
        mapping = ggplot2::aes(
          x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
          colour = .data$.color
        ),
        inherit.aes = FALSE,
        alpha = arrow.alpha, linewidth = arrow.size,
        arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed", angle = 22),
        lineend = "round", linejoin = "mitre",
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_segment(
        data = arr,
        mapping = ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        inherit.aes = FALSE,
        colour = arrow.color, alpha = arrow.alpha, linewidth = arrow.size,
        arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed", angle = 22),
        lineend = "round", linejoin = "mitre"
      )
    }
  } else {
    ## Segment form so we can taper alpha / width along each stream
    seg <- stream_df
    seg$xend <- c(seg$x[-1], NA)
    seg$yend <- c(seg$y[-1], NA)
    seg$t2 <- c(seg$t[-1], NA)
    seg$id2 <- c(seg$id[-1], NA)
    if (".color" %in% names(seg)) {
      seg$.color2 <- c(as.character(seg$.color[-1]), NA)
    }
    seg <- seg[!is.na(seg$xend) & seg$id == seg$id2, , drop = FALSE]
    if (isTRUE(stream.taper)) {
      seg$alpha <- arrow.alpha * (0.25 + 0.75 * seg$t2)
      seg$lw <- stream.width * (0.35 + 0.65 * seg$t2)
    } else {
      seg$alpha <- arrow.alpha
      seg$lw <- stream.width
    }
    if (isTRUE(color_by_cluster) && ".color" %in% names(seg)) {
      p <- p + ggplot2::geom_segment(
        data = seg,
        mapping = ggplot2::aes(
          x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
          colour = .data$.color, alpha = .data$alpha, linewidth = .data$lw
        ),
        inherit.aes = FALSE,
        lineend = "round",
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_segment(
        data = seg,
        mapping = ggplot2::aes(
          x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
          alpha = .data$alpha, linewidth = .data$lw
        ),
        inherit.aes = FALSE,
        colour = arrow.color,
        lineend = "round",
        show.legend = FALSE
      )
    }
    p <- p +
      ggplot2::scale_alpha_identity() +
      ggplot2::scale_linewidth_identity()
    if (!is.null(tip_df) && nrow(tip_df) > 0L && arrow.every > 0L) {
      if (isTRUE(color_by_cluster) && ".color" %in% names(tip_df)) {
        p <- p + ggplot2::geom_segment(
          data = tip_df,
          mapping = ggplot2::aes(
            x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
            colour = .data$.color
          ),
          inherit.aes = FALSE,
          alpha = arrow.alpha,
          linewidth = max(0.15, stream.width * 0.55),
          arrow = grid::arrow(
            length = grid::unit(arrow.head, "cm"),
            type = "closed",
            angle = 25
          ),
          lineend = "butt", linejoin = "mitre",
          show.legend = FALSE
        )
      } else {
        p <- p + ggplot2::geom_segment(
          data = tip_df,
          mapping = ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
          inherit.aes = FALSE,
          colour = arrow.color, alpha = arrow.alpha,
          linewidth = max(0.15, stream.width * 0.55),
          arrow = grid::arrow(
            length = grid::unit(arrow.head, "cm"),
            type = "closed",
            angle = 25
          ),
          lineend = "butt", linejoin = "mitre"
        )
      }
    }
  }


  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    .ic2_theme_cluster_plot(back.col = back.col) +
    ggplot2::coord_equal()

  if (isTRUE(anno.clust) && identical(col.by, "clusters") && isTRUE(color_info$discrete)) {
    counts <- table(df$.color)
    keep_lab <- names(counts)[as.integer(counts) >= 10L]
    if (length(keep_lab) > 0L) {
      df_lab <- df[as.character(df$.color) %in% keep_lab, , drop = FALSE]
      cents <- stats::aggregate(
        cbind(Dim1, Dim2) ~ .color, data = df_lab, FUN = mean, na.rm = TRUE
      )
      p <- p + ggplot2::annotate(
        "text", x = cents$Dim1, y = cents$Dim2, label = cents$.color,
        size = anno.size, fontface = "bold", colour = "black"
      )
    }
  }

  n_show <- if (identical(style, "stream")) length(unique(stream_df$id)) else nrow(arr)
  rep <- .ic2_banner_start(rep, c(
    paste0("Reduction: ", plot.type, "   |   style=", style, "   |   ", pt.col),
    paste0("Paths/arrows: ", .ic2_fmt_num(n_show), "   |   grid=", grid.n)
  ))
  .ic2_finish(rep, extra_lines = paste0(
    "ggplot   |   cells: ", .ic2_fmt_num(nrow(df)),
    "   |   style=", style, "   |   n=", .ic2_fmt_num(n_show)
  ))
  invisible(p)
}

#' Print a side-by-side pseudotime arrow grob
#'
#' @param x A \code{pseudotime_arrows_both} grob.
#' @param newpage Start a new grid page. Default \code{TRUE}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.pseudotime_arrows_both <- function(x, newpage = TRUE, ...) {
  if (isTRUE(newpage)) grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}
