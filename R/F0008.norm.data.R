#' Normalize an iCellR2 assay (chunked, on-disk)
#'
#' Normalizes a counts assay using streaming column chunks. By default
#' (\code{merge = NULL} / auto) the normalized assay is always
#' \strong{one sparse HDF5} when materialized:
#' \itemize{
#'   \item \strong{one sample:} write that sample's normalized file (no merge).
#'   \item \strong{multiple samples:} merge into a single combined file.
#' }
#' Set \code{merge = FALSE} to keep one normalized HDF5 per sample (virtual
#' assay). Set \code{merge = TRUE} to force a single merged file.
#'
#' The full matrix is never densified with \code{as.matrix()}; size factors
#' and normalized values stay sparse. Raw counts under SODA are left untouched.
#'
#' By default only cells present in \code{x@cell_data} are written to the
#' normalized assay. After \code{\link{cell.filter}}, that means raw counts
#' stay intact in SODA while \code{"data"} contains filtered cells only.
#' After \code{\link{cell.filter}}, the fast 10x copy path is skipped so the
#' written \code{"data"} assay contains only kept cells (CSC rebuild / merge).
#' Prefer \code{\link{run.pca}} on \code{"scale.data"} when available
#' (gene.model-sized), else \code{"data"} with on-the-fly \code{log2p1}.
#'
#' Default \code{materialize = TRUE}: write a portable on-disk \code{"data"}
#' assay. Unfiltered native 10x H5 sources may use a fast copy + in-place
#' \code{matrix/data} scale. Set \code{materialize = FALSE} for lazy norm
#' (size factors only; applied on read) when you want to skip the assay write.
#'
#' For virtual multi-sample projects, feature ranks / size-factor scaling are
#' still computed globally across samples.
#'
#' Supported methods:
#' \describe{
#'   \item{\code{global.glsf}}{Library-size factors from all features,
#'     scaled to mean 1. Default method. Reuses QC \code{UMIs} when
#'     present; otherwise streams column sums.}
#'   \item{\code{rpm}}{Relative counts:
#'     \code{counts / library_size * rpm.factor}.
#'     Default \code{rpm.factor = 1e4}. Reuses QC \code{UMIs} when present.}
#'   \item{\code{LogNormalize}}{Library-size relative counts with
#'     \code{scale.factor} (default \code{10000}), intended for use with
#'     natural \code{log1p} downstream (or on materialize). Reuses QC
#'     \code{UMIs} when present.}
#'   \item{\code{ranked.glsf}}{Global scaling factors from the top
#'     \code{top.rank} features by global mean (classic iCellR).}
#'   \item{\code{deseq}}{DESeq2-style median-of-ratios size factors
#'     (poscounts geometric means). Heavier than library-size methods;
#'     prefer for smaller objects or bulk-like matrices.}
#'   \item{\code{spike.in}}{User-supplied per-cell size factors.}
#'   \item{\code{no.norm}}{Copy assay to a new handle without scaling.}
#' }
#'
#' For multimodal projects, set \code{modality} (e.g. \code{"atac"}) instead
#' of the old \code{ATAC.data} flag. Optional \code{match.cells.to.modality}
#' further intersects cell ids with another modality.
#'
#' @param x An object of class \code{iCellR2}.
#' @param modality Modality to normalize. Default \code{"rna"}.
#' @param assay.name Source assay; default active assay.
#' @param new.assay.name Name of the normalized assay to create.
#'   Default \code{"data"}.
#' @param norm.method Normalization method (see details).
#'   Default \code{"global.glsf"}.
#' @param top.rank Number of top features for \code{ranked.glsf}.
#'   Default \code{500}.
#' @param spike.in.factors Numeric size factors (same order as cells in
#'   the assay) when \code{norm.method = "spike.in"}.
#' @param rpm.factor Library-size multiplier for \code{rpm}
#'   (\code{counts / lib * rpm.factor}). Default \code{1e4}.
#' @param scale.factor Library-size multiplier for \code{LogNormalize}
#'   (\code{counts / lib * scale.factor}). Default \code{10000}.
#' @param rounding.digits Digits for optional rounding of nonzero values.
#' @param round.num If \code{TRUE}, round nonzero normalized values.
#' @param match.cells.to.modality Optional modality whose cell ids should
#'   be intersected with this assay (replacement for old ATAC.filter).
#' @param use.cell.data If \code{TRUE} (default), write only cells in
#'   \code{rownames(x@cell_data)} (post-\code{cell.filter} barcodes).
#' @param merge \code{NULL} (default / auto): one sample  ->  single per-sample
#'   HDF5; multiple samples  ->  one merged HDF5. \code{TRUE} forces merge;
#'   \code{FALSE} keeps one file per sample (virtual assay).
#' @param chunk.size Cells per disk chunk. Default \code{500000}.
#' @param use.rcpp Use Rcpp column sums when available.
#' @param set.active Set \code{new.assay.name} as the active assay.
#' @param materialize If \code{TRUE} (default), write a full normalized
#'   assay. Native 10x H5 sources use a fast copy + in-place
#'   \code{matrix/data} scale (script-speed); other backends rebuild CSC
#'   HDF5. Set \code{FALSE} to write size factors only and register a lazy
#'   \code{"data"} assay over raw counts.
#' @param overwrite Overwrite an existing normalized assay file.
#' @param verbose Progress style: \code{TRUE}/\code{"full"}, \code{"compact"},
#'   or \code{FALSE}/\code{"silent"}.
#'
#' @return An updated \code{iCellR2} object with a new on-disk assay handle.
#'
#' @examples
#' \dontrun{
#' obj <- norm.data(obj)  # global.glsf + materialize (default; reuses QC UMIs)
#' obj <- norm.data(obj, materialize = FALSE)  # lazy size factors only
#' obj <- norm.data(obj, norm.method = "rpm")
#' obj <- norm.data(obj, norm.method = "LogNormalize")  # scale.factor = 1e4
#' obj <- norm.data(obj, norm.method = "ranked.glsf", top.rank = 500)
#' obj <- norm.data(obj, norm.method = "deseq")
#' obj <- run.pca(obj, assay.name = "data", transform = "log2p1")
#' }
#'
#' @export
norm.data <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    new.assay.name = "data",
    norm.method = c(
      "global.glsf", "rpm", "LogNormalize", "ranked.glsf",
      "deseq", "spike.in", "no.norm"
    ),
    top.rank = 500,
    spike.in.factors = NULL,
    rpm.factor = 1e4,
    scale.factor = 10000,
    rounding.digits = 3,
    round.num = FALSE,
    match.cells.to.modality = NULL,
    use.cell.data = TRUE,
    merge = NULL,
    chunk.size = 500000L,
    use.rcpp = TRUE,
    set.active = TRUE,
    materialize = TRUE,
    overwrite = TRUE,
    verbose = TRUE
) {
  .ic2_require_object(x)
  norm.method <- match.arg(norm.method)
  modality <- tolower(modality)
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  chunk.size <- as.integer(chunk.size)
  top.rank <- as.integer(top.rank)
  rep <- .ic2_new_reporter("norm.data", verbose, soda_root = x@storage$root_dir)

  if (!is.character(new.assay.name) || length(new.assay.name) != 1L || !nzchar(new.assay.name)) {
    stop("'new.assay.name' must be a non-empty character string.")
  }

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  ## Auto: one sample  ->  no merge; multi-sample  ->  single merged HDF5
  if (is.null(merge)) {
    merge <- length(handles) > 1L
  } else {
    merge <- isTRUE(merge)
  }
  rep <- .ic2_register_assay_inputs(rep, x, modality, assay.name, label = paste0("Assay[", assay.name, "]"))
  is_virtual <- .ic2_assay_is_virtual(x, modality, assay.name)
  n_cells_guess <- sum(vapply(handles, function(h) {
    as.integer(.ic2_enrich_handle(h, FALSE)$n_cells %||% 0L)
  }, integer(1)))

  ## Default: write only barcodes still in cell_data (post cell.filter)
  keep_ids <- NULL
  if (isTRUE(use.cell.data)) {
    keep_ids <- rownames(x@cell_data)
    if (is.null(keep_ids) || length(keep_ids) < 1L) {
      stop("'use.cell.data = TRUE' but x@cell_data has no cell ids.")
    }
  }
  if (!is.null(match.cells.to.modality)) {
    match.cells.to.modality <- tolower(match.cells.to.modality)
    ref_ids <- rownames(x@cell_data)
    keep_ids <- if (is.null(keep_ids)) ref_ids else intersect(keep_ids, ref_ids)
  }

  n_write <- if (!is.null(keep_ids)) length(keep_ids) else n_cells_guess
  rep <- .ic2_banner_start(rep, c(
    paste0(
      "Method: ", norm.method,
      if (identical(norm.method, "ranked.glsf")) {
        paste0(" (top.rank = ", top.rank, ")")
      } else if (identical(norm.method, "rpm")) {
        paste0(" (rpm.factor = ", rpm.factor, ")")
      } else if (identical(norm.method, "LogNormalize")) {
        paste0(" (scale.factor = ", scale.factor, ")")
      } else {
        ""
      },
      if (!isTRUE(materialize)) "   |   lazy" else "   |   materialize"
    ),
    paste0("Assay: ", assay.name, "  ->  ", new.assay.name),
    paste0(
      "Source cells: ", .ic2_fmt_num(n_cells_guess),
      " | Writing: ", .ic2_fmt_num(n_write),
      " | Handles: ", length(handles),
      if (is_virtual) " (virtual)" else "",
      if (!is.null(keep_ids)) " | filtered via cell_data" else ""
    )
  ))
  if (!is.null(keep_ids) && n_write < n_cells_guess) {
    .ic2_step(rep, paste0(
      "Restricting to ", .ic2_fmt_num(n_write), " / ", .ic2_fmt_num(n_cells_guess),
      " cells from cell_data (raw counts unchanged)."
    ))
  }
  if (!is.null(match.cells.to.modality)) {
    .ic2_step(rep, paste0(
      "match.cells.to.modality = '", match.cells.to.modality, "' applied."
    ))
  }

  ## Compute size factors (ranked.glsf / deseq need a global feature pass)
  pooled_row_sums <- NULL
  deseq_loggeomeans <- NULL
  if (identical(norm.method, "ranked.glsf")) {
    .ic2_step(rep, "[1/3] Estimating feature ranks (streaming row sums)")
    for (i in seq_along(handles)) {
      h <- .ic2_enrich_handle(handles[[i]], load_dimnames = TRUE)
      col_idx <- NULL
      if (!is.null(keep_ids)) {
        col_idx <- which(h$cell_ids %in% keep_ids)
        if (length(col_idx) < 1L) {
          stop("No cells left after cell_data filter for handle ", i, ".")
        }
      }
      rs <- .ic2_handle_row_sums(h, chunk.size = chunk.size, col_idx = col_idx)
      if (is.null(pooled_row_sums)) {
        pooled_row_sums <- rs
      } else {
        common <- intersect(names(pooled_row_sums), names(rs))
        only1 <- setdiff(names(pooled_row_sums), names(rs))
        only2 <- setdiff(names(rs), names(pooled_row_sums))
        pooled_row_sums <- c(
          pooled_row_sums[common] + rs[common],
          pooled_row_sums[only1],
          rs[only2]
        )
      }
      rep <- .ic2_progress_bar(rep, i, length(handles), prefix = "  rank pass")
    }
    ord <- order(pooled_row_sums, decreasing = TRUE)
    top_n <- min(top.rank, length(ord))
    top_features <- names(pooled_row_sums)[ord[seq_len(top_n)]]
  } else if (identical(norm.method, "deseq")) {
    .ic2_step(rep, "[1/3] DESeq poscounts geometric means (streaming)")
    sum_log <- NULL
    n_pos <- NULL
    for (i in seq_along(handles)) {
      h <- .ic2_enrich_handle(handles[[i]], load_dimnames = TRUE)
      col_idx <- NULL
      if (!is.null(keep_ids)) {
        col_idx <- which(h$cell_ids %in% keep_ids)
        if (length(col_idx) < 1L) {
          stop("No cells left after cell_data filter for handle ", i, ".")
        }
      }
      acc <- .ic2_handle_deseq_poscounts_accum(
        h, chunk.size = chunk.size, col_idx = col_idx
      )
      if (is.null(sum_log)) {
        sum_log <- acc$sum_log
        n_pos <- acc$n_pos
      } else {
        common <- intersect(names(sum_log), names(acc$sum_log))
        only1 <- setdiff(names(sum_log), names(acc$sum_log))
        only2 <- setdiff(names(acc$sum_log), names(sum_log))
        sum_log <- c(
          sum_log[common] + acc$sum_log[common],
          sum_log[only1],
          acc$sum_log[only2]
        )
        n_pos <- c(
          n_pos[common] + acc$n_pos[common],
          n_pos[only1],
          acc$n_pos[only2]
        )
      }
      rep <- .ic2_progress_bar(rep, i, length(handles), prefix = "  deseq geo")
    }
    deseq_loggeomeans <- rep(-Inf, length(sum_log))
    names(deseq_loggeomeans) <- names(sum_log)
    ok <- n_pos > 0L & is.finite(sum_log)
    deseq_loggeomeans[ok] <- sum_log[ok] / n_pos[ok]
    if (!any(is.finite(deseq_loggeomeans))) {
      stop("DESeq size factors: no gene with positive counts across cells.")
    }
    top_features <- NULL
  } else {
    top_features <- NULL
    .ic2_step(rep, "[1/3] Preparing size-factor method")
  }

  .ic2_step(rep, "[2/3] Computing size factors")
  all_sf <- list()
  lib_sources <- character(0)

  for (i in seq_along(handles)) {
    h <- .ic2_enrich_handle(handles[[i]], load_dimnames = TRUE)
    local_keep <- h$cell_ids
    if (!is.null(keep_ids)) {
      local_keep <- intersect(h$cell_ids, keep_ids)
      if (length(local_keep) < 1L) {
        stop("No cells left after cell_data filter for handle ", i, ".")
      }
    }

    if (identical(norm.method, "no.norm")) {
      sf <- rep(1, length(local_keep))
      names(sf) <- local_keep
    } else if (identical(norm.method, "spike.in")) {
      if (is.null(spike.in.factors)) {
        stop("'spike.in.factors' must be provided for norm.method = 'spike.in'.")
      }
      if (is.null(names(spike.in.factors))) {
        if (length(spike.in.factors) != h$n_cells) {
          stop("Unnamed spike.in.factors length must equal n_cells for each handle.")
        }
        sf <- as.numeric(spike.in.factors)
        names(sf) <- h$cell_ids
        sf <- sf[local_keep]
      } else {
        sf <- as.numeric(spike.in.factors[local_keep])
        names(sf) <- local_keep
        if (any(is.na(sf))) stop("Missing spike.in.factors for some cell ids.")
      }
    } else if (norm.method %in% c("global.glsf", "rpm", "LogNormalize")) {
      sf <- .ic2_cell_library_sizes(
        local_keep,
        meta = x@cell_data,
        handle = h,
        chunk.size = chunk.size,
        use.rcpp = use.rcpp
      )
      lib_sources <- c(lib_sources, attr(sf, "source") %||% "colSums")
    } else if (identical(norm.method, "ranked.glsf")) {
      top_idx <- match(top_features, h$feature_ids)
      top_idx <- top_idx[!is.na(top_idx)]
      if (length(top_idx) < 1L) {
        stop("None of the ranked top features were found in assay handle ", i, ".")
      }
      sf_all <- .ic2_handle_col_sums(
        h, row_idx = top_idx, chunk.size = chunk.size, use.rcpp = use.rcpp
      )
      names(sf_all) <- h$cell_ids
      sf <- sf_all[local_keep]
      ## cells with zero top-gene counts fall back to full library size
      if (any(sf <= 0 | !is.finite(sf))) {
        lib <- .ic2_cell_library_sizes(
          local_keep,
          meta = x@cell_data,
          handle = h,
          chunk.size = chunk.size,
          use.rcpp = use.rcpp
        )
        bad <- which(sf <= 0 | !is.finite(sf))
        sf[bad] <- lib[bad]
        lib_sources <- c(lib_sources, paste0("ranked+", attr(lib, "source") %||% "colSums"))
      } else {
        lib_sources <- c(lib_sources, "ranked.top")
      }
    } else if (identical(norm.method, "deseq")) {
      lg <- deseq_loggeomeans[match(h$feature_ids, names(deseq_loggeomeans))]
      names(lg) <- h$feature_ids
      lg[is.na(lg)] <- -Inf
      sf_all <- .ic2_handle_deseq_size_factors(
        h, loggeomeans = lg, chunk.size = chunk.size
      )
      names(sf_all) <- h$cell_ids
      sf <- sf_all[local_keep]
      ## Empty / failed cells: fall back to library size
      if (any(sf <= 0 | !is.finite(sf))) {
        lib <- .ic2_cell_library_sizes(
          local_keep,
          meta = x@cell_data,
          handle = h,
          chunk.size = chunk.size,
          use.rcpp = use.rcpp
        )
        bad <- which(sf <= 0 | !is.finite(sf))
        sf[bad] <- lib[bad]
        lib_sources <- c(lib_sources, paste0("deseq+", attr(lib, "source") %||% "colSums"))
      } else {
        lib_sources <- c(lib_sources, "deseq")
      }
    } else {
      stop("Unsupported norm.method.")
    }

    ## Keep size factors only for cells we will write
    if (any(sf <= 0 | !is.finite(sf))) {
      stop(
        "Non-positive size factors for ",
        sum(sf <= 0 | !is.finite(sf)),
        " cells in handle ", i, ". Check empty cells / filtering."
      )
    }

    all_sf[[i]] <- sf
    rep <- .ic2_progress_bar(rep, i, length(handles), prefix = "  size factors")
  }

  if (length(lib_sources) > 0L) {
    src_tab <- table(lib_sources)
    .ic2_step(rep, paste0(
      "Library sizes: ",
      paste(sprintf("%s x %s", names(src_tab), as.integer(src_tab)), collapse = ", ")
    ))
    if (any(grepl("colSums", names(src_tab), fixed = TRUE))) {
      .ic2_step(rep, paste0(
        "WARNING: streamed column sums for library sizes  -  slow on large 10x H5. ",
        "Prefer qc.stats() UMIs (Library sizes: UMIs x N)."
      ))
    }
  }

  ## Scale factors relative to global mean (glsf / deseq)  -  filtered cells only
  if (norm.method %in% c("global.glsf", "ranked.glsf", "deseq")) {
    pooled <- unlist(all_sf, use.names = TRUE)
    mu <- mean(pooled)
    if (!is.finite(mu) || mu <= 0) stop("Invalid mean size factor.")
    all_sf <- lapply(all_sf, function(sf) sf / mu)
  } else if (identical(norm.method, "rpm")) {
    all_sf <- lapply(all_sf, function(sf) sf / rpm.factor)
  } else if (identical(norm.method, "LogNormalize")) {
    all_sf <- lapply(all_sf, function(sf) sf / scale.factor)
  }

  .ic2_step(rep, if (!isTRUE(materialize)) {
    "[3/3] Writing size factors only (lazy norm; raw counts unchanged)"
  } else if (isTRUE(merge)) {
    "[3/3] Writing single merged normalized HDF5 (raw untouched)"
  } else {
    "[3/3] Writing per-sample normalized HDF5 (virtual assay; raw untouched)"
  })

  ## Feature universe from actual assay dimnames (guarantees row match)
  feature_universe <- character(0)
  for (h in handles) {
    h2 <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    feature_universe <- union(feature_universe, as.character(h2$feature_ids))
  }
  feature_universe <- as.character(feature_universe)
  if (length(feature_universe) < 1L) stop("No features found across assay handles.")

  ## Size factors first (needed for both lazy and materialized)
  sf_all <- unlist(all_sf, use.names = TRUE)
  sf_dir <- file.path(
    x@storage$root_dir, "modalities", modality, "assays", new.assay.name
  )
  dir.create(sf_dir, recursive = TRUE, showWarnings = FALSE)
  sf_ids <- names(sf_all)
  if (is.null(sf_ids) || !length(sf_ids)) {
    sf_ids <- as.character(seq_along(sf_all))
  }
  utils::write.table(
    data.frame(
      cell_id = sf_ids,
      size_factor = as.numeric(sf_all),
      stringsAsFactors = FALSE
    ),
    file = file.path(sf_dir, "size_factors.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  sf_path <- file.path(sf_dir, "size_factors.tsv")
  .ic2_soda_unlink_rds(file.path(sf_dir, "size_factors.rds"))

  new_handles <- list()
  dest_paths <- character(0)
  mat_transform <- if (identical(norm.method, "LogNormalize")) "log1p" else "none"

  if (!isTRUE(materialize)) {
    x <- .ic2_assign_lazy_norm_assay(
      x,
      modality = modality,
      assay.name = new.assay.name,
      source_assay = assay.name,
      size_factors_path = sf_path,
      feature_universe = feature_universe,
      n_cells = length(sf_all),
      set.active = set.active
    )
  } else {
    ## Fast path: native 10x H5  ->  copy + patch matrix/data (script-speed).
    ## Disabled when cell.filter dropped cells  -  full copy would leave dropped
    ## barcodes in "data" and force gappy PCA column reads.
    cells_filtered <- !is.null(keep_ids) && as.integer(n_write) < as.integer(n_cells_guess)
    if (isTRUE(cells_filtered)) {
      .ic2_step(rep, paste0(
        "Fast 10x copy skipped (filtered cells: writing ",
        .ic2_fmt_num(n_write), " / ", .ic2_fmt_num(n_cells_guess),
        " via CSC rebuild)."
      ))
    }
    tenx_ok <- !isTRUE(cells_filtered) &&
      norm.method %in% c("global.glsf", "rpm", "LogNormalize", "no.norm") &&
      length(handles) >= 1L &&
      (!isTRUE(merge) || length(handles) == 1L) &&
      all(vapply(seq_along(handles), function(i) {
        .ic2_tenx_inplace_norm_eligible(
          handles[[i]], norm.method, feature_universe = feature_universe
        )
      }, logical(1)))

    if (isTRUE(tenx_ok)) {
      .ic2_step(rep, paste0(
        "[3/3] Fast 10x materialize: copy H5 + in-place scale matrix/data (",
        norm.method, ")"
      ))
      new_handles <- vector("list", length(handles))
      dest_paths <- character(length(handles))
      for (i in seq_along(handles)) {
        h <- .ic2_enrich_handle(handles[[i]], load_dimnames = TRUE)
        sid <- h$sample_id
        if (is.null(sid) || !nzchar(as.character(sid)[1L])) {
          sid <- paste0("sample", i)
        }
        dest <- .ic2_soda_assay_path(
          x, modality, new.assay.name, sample_id = sid, ext = ".h5"
        )
        if (file.exists(dest) && !overwrite) {
          stop("Normalized assay exists at '", dest, "'. Use overwrite = TRUE.")
        }
        .ic2_step(rep, paste0("  Sample ", sid, " (", i, "/", length(handles), ")"))
        nh <- .ic2_materialize_tenx_inplace(
          handle = h,
          dest_path = dest,
          size_factors = all_sf[[i]],
          chunk.size = chunk.size,
          transform = mat_transform,
          round.num = round.num,
          rounding.digits = rounding.digits,
          overwrite = overwrite,
          progress = function(cur, tot) {
            rep <<- .ic2_progress_bar(rep, cur, tot, prefix = paste0("  ", sid))
          }
        )
        nh$sample_id <- sid
        nh$modality <- modality
        new_handles[[i]] <- nh
        dest_paths[i] <- dest
      }
      x <- .ic2_assign_assay(
        x,
        modality = modality,
        assay.name = new.assay.name,
        handles = new_handles,
        virtual = length(new_handles) > 1L,
        feature_universe = if (length(new_handles) > 1L) feature_universe else NULL,
        set.active = set.active
      )
    } else if (isTRUE(merge)) {
      dest <- .ic2_soda_assay_path(x, modality, new.assay.name, sample_id = NULL, ext = ".h5")
      if (file.exists(dest) && !overwrite) {
        stop("Normalized assay exists at '", dest, "'. Use overwrite = TRUE.")
      }
      .ic2_step(rep, "[3/3] Writing merged normalized HDF5 (CSC rebuild)")
      merged <- .ic2_normalize_merged_hdf5(
        handles = handles,
        size_factors_list = all_sf,
        dest_path = dest,
        feature_universe = feature_universe,
        chunk.size = chunk.size,
        round.num = round.num,
        rounding.digits = rounding.digits,
        overwrite = overwrite,
        transform = mat_transform,
        progress = function(cur, tot) {
          rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  write assay")
        }
      )
      merged$modality <- modality
      new_handles <- list(merged)
      dest_paths <- dest
      x <- .ic2_assign_assay(
        x,
        modality = modality,
        assay.name = new.assay.name,
        handles = new_handles,
        virtual = FALSE,
        feature_universe = NULL,
        set.active = set.active
      )
    } else {
      .ic2_step(rep, "[3/3] Writing per-sample normalized HDF5 (CSC rebuild)")
      new_handles <- vector("list", length(handles))
      dest_paths <- character(length(handles))
      for (i in seq_along(handles)) {
        h <- .ic2_enrich_handle(handles[[i]], load_dimnames = TRUE)
        sid <- h$sample_id
        if (is.null(sid) || !nzchar(as.character(sid)[1L])) {
          sid <- paste0("sample", i)
        }
        dest <- .ic2_soda_assay_path(
          x, modality, new.assay.name, sample_id = sid, ext = ".h5"
        )
        if (file.exists(dest) && !overwrite) {
          stop("Normalized assay exists at '", dest, "'. Use overwrite = TRUE.")
        }
        .ic2_step(rep, paste0("  Sample ", sid, " (", i, "/", length(handles), ")"))
        nh <- .ic2_normalize_handle(
          handle = h,
          dest_path = dest,
          size_factors = all_sf[[i]],
          chunk.size = chunk.size,
          round.num = round.num,
          rounding.digits = rounding.digits,
          keep_cell_ids = names(all_sf[[i]]),
          transform = mat_transform,
          progress = function(cur, tot) {
            rep <<- .ic2_progress_bar(rep, cur, tot, prefix = paste0("  ", sid))
          }
        )
        nh$sample_id <- sid
        nh$modality <- modality
        new_handles[[i]] <- nh
        dest_paths[i] <- dest
      }
      x <- .ic2_assign_assay(
        x,
        modality = modality,
        assay.name = new.assay.name,
        handles = new_handles,
        virtual = length(new_handles) > 1L,
        feature_universe = if (length(new_handles) > 1L) feature_universe else NULL,
        set.active = set.active
      )
    }
  }

  ## Store compact size-factor summary (not necessarily all factors if huge)
  x@params$normalization[[modality]] <- list(
    assay.name = assay.name,
    new.assay.name = new.assay.name,
    norm.method = norm.method,
    merge = isTRUE(merge),
    materialize = isTRUE(materialize),
    lazy_norm = !isTRUE(materialize),
    top.rank = if (identical(norm.method, "ranked.glsf")) top.rank else NA_integer_,
    rpm.factor = if (identical(norm.method, "rpm")) rpm.factor else NA_real_,
    scale.factor = if (identical(norm.method, "LogNormalize")) scale.factor else NA_real_,
    preferred_transform = if (identical(norm.method, "LogNormalize")) "log1p" else NA_character_,
    n_cells = length(sf_all),
    size_factor_mean = mean(sf_all),
    size_factor_median = stats::median(sf_all),
    chunk.size = chunk.size,
    size_factors_path = sf_path
  )

  x <- .ic2_log_command(
    x,
    "norm.data",
    arguments = list(
      modality = modality,
      assay.name = assay.name,
      new.assay.name = new.assay.name,
      norm.method = norm.method,
      top.rank = top.rank,
      chunk.size = chunk.size,
      materialize = isTRUE(materialize)
    ),
    outputs = list(
      n_cells = length(sf_all),
      size_factors_path = sf_path,
      lazy_norm = !isTRUE(materialize)
    )
  )

  try(.ic2_soda_write_manifest(x), silent = TRUE)

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  for (p in dest_paths) rep <- .ic2_add_output(rep, "Norm assay", p)
  rep <- .ic2_add_output(rep, "Size factors", sf_path)
  rep <- .ic2_add_output(rep, "Active assay", new.assay.name)

  rep <- .ic2_add_scale(rep, c(
    paste0(
      "Mode:     ",
      if (!isTRUE(materialize)) {
        "lazy size factors only (no full data HDF5)"
      } else if (isTRUE(merge)) {
        "sparse column scaling  ->  one merged HDF5"
      } else {
        "sparse column scaling  ->  per-sample HDF5 (virtual)"
      }
    ),
    "Passes:   stream  ->  size factors",
    if (isTRUE(materialize)) "            ->  write assay" else "            ->  skip assay write",
    "Memory:   chunk-bound, not dataset-bound",
    "Next:     gene.stats / run.pca(..., transform = \"log2p1\") on assay \"data\""
  ))
  wall <- as.numeric(difftime(Sys.time(), rep$wall0, units = "secs"))
  if (is.finite(wall) && wall > 0) {
    rep$stats$cells_per_sec <- length(sf_all) / wall
  }
  n_feat_out <- length(feature_universe)
  n_cell_out <- length(sf_all)
  if (isTRUE(materialize) && length(new_handles) > 0L) {
    n_feat_out <- if (isTRUE(merge)) {
      new_handles[[1L]]$n_features
    } else {
      length(feature_universe)
    }
    n_cell_out <- sum(vapply(new_handles, function(h) as.integer(h$n_cells), integer(1)))
  }
  .ic2_finish(rep, extra_lines = paste0(
    if (!isTRUE(materialize)) {
      "Lazy norm: "
    } else if (isTRUE(merge)) {
      "Merged assay: "
    } else {
      "Norm assay: "
    },
    .ic2_fmt_num(n_feat_out), "  x  ", .ic2_fmt_num(n_cell_out),
    if (isTRUE(materialize) && !isTRUE(merge)) paste0("   |   samples: ", length(new_handles)) else "",
    "   |   size factors mean ", sprintf("%.2f", mean(sf_all)),
    "  median ", sprintf("%.2f", stats::median(sf_all))
  ))
  x
}
