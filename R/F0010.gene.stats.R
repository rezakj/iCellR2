#' Compute per-gene statistics (chunked, on-disk)
#'
#' Streams an assay handle and computes per-feature statistics without
#' densifying the matrix. Results are stored in the modality and written
#' under the SODA archive for use by \code{\link{make.gene.model}}.
#'
#' For a single native 10x HDF5 assay (\code{tenx_h5}), uses a bulk CSC
#' path (full \code{data}/\code{indices}/\code{indptr} + \code{tabulate}/
#' \code{rowsum}), matching the speed of the standalone
#' \code{gene_stats_10x_h5.R} script. Multi-sample / non-tenx backends
#' still stream in chunks.
#'
#' By default summarizes normalized \code{"data"} (\code{use.raw = FALSE}).
#' Run \code{\link{norm.data}} first. Set \code{use.raw = TRUE} for raw
#' \code{"counts"}. When \code{use.raw = FALSE}, size factors are applied
#' \emph{only} if that assay is still lazy (not materialized). A
#' materialized \code{"data"} matrix is used as-is with no extra scaling.
#'
#' @param x An \code{iCellR2} object.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Assay to summarize. Default \code{NULL}: uses
#'   \code{"data"} when \code{use.raw = FALSE}, or \code{"counts"} when
#'   \code{use.raw = TRUE}.
#' @param use.raw If \code{FALSE} (default), use normalized \code{"data"}
#'   (materialized or lazy). If \code{TRUE}, prefer raw \code{"counts"}.
#' @param group.by Optional metadata column (e.g. \code{"condition_id"})
#'   to also compute per-group stats. Default \code{NULL} (all cells only).
#' @param chunk.size Cells per disk chunk. Default \code{200000}.
#' @param nnz.budget Max nonzeros per HDF5 hyperslab in the tenx bulk path.
#'   Default \code{1e9}. Raise on big-RAM / local-SSD servers; lower if RAM-bound.
#' @param use.rcpp Use Rcpp row accumulators when available.
#' @param use.cell.data If \code{TRUE} (default), only cells in
#'   \code{rownames(x@cell_data)} contribute (post-\code{cell.filter}).
#' @param update.object If \code{TRUE}, store results on \code{x}.
#' @param verbose Progress style: \code{TRUE}/\code{"full"}, \code{"compact"},
#'   or \code{FALSE}/\code{"silent"}.
#'
#' @return Updated \code{iCellR2} object (or a data frame if
#'   \code{update.object = FALSE}).
#'
#' @examples
#' \dontrun{
#' obj <- gene.stats(obj)                         # normalized "data" (default)
#' obj <- gene.stats(obj, use.raw = TRUE)         # raw counts
#' obj <- gene.stats(obj, group.by = "condition_id")
#' }
#'
#' @export
gene.stats <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    use.raw = FALSE,
    group.by = NULL,
    chunk.size = 200000L,
    nnz.budget = 1e9,
    use.rcpp = TRUE,
    use.cell.data = TRUE,
    update.object = TRUE,
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  assay.name <- .ic2_resolve_gene_stats_assay(x, modality, assay.name, use.raw)
  chunk.size <- as.integer(chunk.size)
  rep <- .ic2_new_reporter("gene.stats", verbose, soda_root = x@storage$root_dir)

  if (!is.null(group.by)) {
    if (!group.by %in% colnames(x@cell_data)) {
      stop("'group.by' column '", group.by, "' not found in x@cell_data.")
    }
  }

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  ## Lazy SF only when use.raw = FALSE and norm.data was not materialized
  size_factors <- NULL
  if (!isTRUE(use.raw) && .ic2_assay_is_lazy_norm(x, modality, assay.name)) {
    size_factors <- .ic2_assay_size_factors(x, modality, assay.name)
  }
  rep <- .ic2_register_assay_inputs(rep, x, modality, assay.name, label = paste0("Assay[", assay.name, "]"))
  cell_mask <- if (isTRUE(use.cell.data)) rownames(x@cell_data) else NULL
  n_cells_guess <- if (!is.null(cell_mask)) {
    length(cell_mask)
  } else {
    sum(vapply(handles, function(h) {
      as.integer(.ic2_enrich_handle(h, FALSE)$n_cells %||% 0L)
    }, integer(1)))
  }

  rep <- .ic2_banner_start(rep, c(
    paste0(
      "Modality: ", modality, " | Assay: ", assay.name,
      if (isTRUE(use.raw)) {
        " (raw)"
      } else if (!is.null(size_factors)) {
        " (lazy norm)"
      } else {
        " (normalized, materialized)"
      }
    ),
    paste0(
      "Cells: ", .ic2_fmt_num(n_cells_guess),
      " | Handles: ", length(handles),
      if (!is.null(group.by)) paste0(" | group.by: ", group.by) else " | group.by: all",
      if (!is.null(cell_mask)) " | masked to cell_data" else ""
    )
  ))
  if (!is.null(size_factors)) {
    .ic2_step(rep, paste0(
      "Lazy norm: scaling ", .ic2_fmt_num(length(size_factors)),
      " cells with size factors (norm.data not materialized)"
    ))
  } else if (!isTRUE(use.raw) && !.ic2_assay_is_lazy_norm(x, modality, assay.name)) {
    .ic2_step(rep, "Materialized 'data' assay  -  using values as-is (no SF rescale)")
  }

  ## Slabbed CSC for tenx_h5: read indptr once, then hyperslab data/indices
  ## in nnz-bounded chunks (avoids loading ~2.6B nnz into RAM / mem.maxVSize).
  .ic2_accum_tenx_bulk <- function(
      handle,
      col_idx,
      size_factors = NULL,
      nnz.budget = 1e9
  ) {
    handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
    if (!identical(handle$backend, "tenx_h5")) return(NULL)
    if (is.null(handle$meta_path) || !file.exists(handle$meta_path)) return(NULL)
    meta <- .ic2_soda_read_tenx_meta(handle$meta_path)
    if (is.null(meta)) return(NULL)
    if (!requireNamespace("hdf5r", quietly = TRUE)) return(NULL)

    path <- handle$path
    group <- as.character(meta$group %||% handle$group %||% "matrix")[1L]
    n_cells <- as.integer(meta$n_cells %||% handle$n_cells)
    n_feat <- as.integer(handle$n_features)
    genes <- handle$feature_ids
    col_idx <- sort(unique(as.integer(col_idx)))
    if (!length(col_idx) || anyNA(col_idx)) return(NULL)
    nnz.budget <- as.numeric(nnz.budget)[1L]
    if (!is.finite(nnz.budget) || nnz.budget < 1e6) nnz.budget <- 1e9

    h5 <- hdf5r::H5File$new(path, mode = "r")
    on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
    g <- h5[[group]]
    indptr <- as.numeric(g[["indptr"]]$read())
    if (length(indptr) != n_cells + 1L) {
      stop("Invalid tenx CSC: length(indptr) != n_cells + 1")
    }
    ds_data <- g[["data"]]
    ds_idx <- g[["indices"]]

    ## Optional feature filter from tenx meta (row_keep)
    remap <- NULL
    row_keep <- meta$row_keep
    if (!is.null(row_keep) && length(row_keep)) {
      n_src <- as.integer(meta$n_features_src %||% n_feat)
      remap <- rep(NA_integer_, n_src)
      remap[as.integer(row_keep)] <- seq_along(row_keep) - 1L
    }

    sf_by_col <- NULL
    if (!is.null(size_factors)) {
      sf_by_col <- as.numeric(size_factors[handle$cell_ids])
      sf_by_col[!is.finite(sf_by_col) | is.na(sf_by_col)] <- 1
    }

    sumv <- numeric(n_feat)
    sumsq <- numeric(n_feat)
    nnzv <- integer(n_feat)
    n_use <- length(col_idx)
    done <- 0L

    ## Keep mask: after cell.filter, dropped barcodes punch holes in col_idx.
    ## Contig-only runs  ->  thousands of tiny H5 reads (slow on NFS/server disks).
    ## Prefer covering spans of the original matrix; skip dropped columns when
    ## accumulating (same totals, far fewer hyperslabs).
    keep_l <- logical(n_cells)
    keep_l[col_idx] <- TRUE
    span_lo <- col_idx[[1L]]
    span_hi <- col_idx[[n_use]]
    dens <- n_use / max(1, span_hi - span_lo + 1L)
    use_cover <- dens >= 0.5 || (span_hi - span_lo + 1L - n_use) <= max(64L, as.integer(0.05 * n_use))

    .accum_slab <- function(c_start, c_end, cols_keep) {
      if (length(cols_keep) < 1L) return(invisible(NULL))
      off0 <- indptr[c_start]
      n_x <- as.numeric(indptr[c_end + 1L] - off0)
      if (n_x < 1) return(invisible(NULL))
      vals <- as.numeric(.ic2_h5_read_1d(ds_data, start1 = off0 + 1, count = n_x))
      gene_idx0 <- as.integer(.ic2_h5_read_1d(ds_idx, start1 = off0 + 1, count = n_x))
      cols_here <- c_start:c_end
      ## Drop non-kept columns inside the covering span
      if (length(cols_keep) < length(cols_here)) {
        keep_col <- keep_l[cols_here]
        if (!all(keep_col)) {
          reps <- as.integer(indptr[cols_here + 1L] - indptr[cols_here])
          take <- rep.int(keep_col, times = reps)
          vals <- vals[take]
          gene_idx0 <- gene_idx0[take]
          cols_here <- cols_here[keep_col]
          rm(take, keep_col, reps)
        }
      }
      if (!length(vals)) return(invisible(NULL))
      if (!is.null(sf_by_col)) {
        reps <- as.integer(indptr[cols_here + 1L] - indptr[cols_here])
        vals <- vals * rep.int(sf_by_col[cols_here], times = reps)
      }
      if (!is.null(remap)) {
        mapped <- remap[gene_idx0 + 1L]
        keep <- !is.na(mapped)
        gene_idx0 <- mapped[keep]
        vals <- vals[keep]
        rm(mapped, keep)
      }
      gene_idx <- gene_idx0 + 1L
      nnzv <<- nnzv + as.integer(tabulate(gene_idx, nbins = n_feat))
      if (length(vals)) {
        sums <- rowsum(vals, group = gene_idx, reorder = FALSE)
        ix <- as.integer(rownames(sums))
        sumv[ix] <<- sumv[ix] + sums[, 1L]
        sums_sq <- rowsum(vals * vals, group = gene_idx, reorder = FALSE)
        ix2 <- as.integer(rownames(sums_sq))
        sumsq[ix2] <<- sumsq[ix2] + sums_sq[, 1L]
      }
      invisible(NULL)
    }

    if (isTRUE(use_cover)) {
      ## Walk original column index; grow slabs by nnz budget
      pos <- span_lo
      while (pos <= span_hi) {
        c_start <- pos
        c_end <- c_start
        ## Grow contiguous original span until nnz budget
        while (c_end < span_hi) {
          nnz_try <- indptr[c_end + 2L] - indptr[c_start]
          if (nnz_try > nnz.budget) break
          c_end <- c_end + 1L
        }
        ## Ensure at least one column
        if (c_end < c_start) c_end <- c_start
        cols_keep <- which(keep_l[c_start:c_end]) + (c_start - 1L)
        .accum_slab(c_start, c_end, cols_keep)
        done <- done + length(cols_keep)
        rep <<- .ic2_progress_bar(rep, done, n_use, prefix = "  bulk")
        pos <- c_end + 1L
        ## Skip leading dropped columns quickly
        while (pos <= span_hi && !keep_l[pos]) pos <- pos + 1L
      }
    } else {
      runs <- .ic2_contig_runs(col_idx)
      for (run in runs) {
        i <- 1L
        nr <- length(run)
        while (i <= nr) {
          c_start <- run[i]
          j <- i
          while (j < nr) {
            nnz_try <- indptr[run[j + 1L] + 1L] - indptr[c_start]
            if (nnz_try > nnz.budget) break
            j <- j + 1L
          }
          c_end <- run[j]
          .accum_slab(c_start, c_end, c_start:c_end)
          done <- done + (j - i + 1L)
          rep <<- .ic2_progress_bar(rep, done, n_use, prefix = "  bulk")
          i <- j + 1L
        }
      }
    }

    list(
      genes = genes,
      sum = sumv,
      sumsq = sumsq,
      nnz = as.integer(nnzv),
      n_cells = n_use,
      method = "bulk_csc"
    )
  }

  .ic2_accum_handle <- function(handle, cell_mask = NULL) {
    handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
    n_feat <- handle$n_features
    n_cells <- handle$n_cells
    idx <- seq_len(n_cells)
    if (!is.null(cell_mask)) {
      idx <- which(handle$cell_ids %in% cell_mask)
      if (length(idx) < 1L) {
        return(list(
          genes = handle$feature_ids,
          sum = numeric(n_feat),
          sumsq = numeric(n_feat),
          nnz = integer(n_feat),
          n_cells = 0L
        ))
      }
    }

    ## Fast path: native 10x H5 CSC (single-file / per-sample handle)
    if (identical(handle$backend, "tenx_h5")) {
      bulk <- tryCatch(
        .ic2_accum_tenx_bulk(handle, idx, size_factors, nnz.budget = nnz.budget),
        error = function(e) {
          .ic2_step(rep, paste0(
            "Bulk CSC unavailable (", conditionMessage(e), "); using chunked path"
          ))
          NULL
        }
      )
      if (!is.null(bulk)) return(bulk)
    }

    sumv <- numeric(n_feat)
    sumsq <- numeric(n_feat)
    nnzv <- integer(n_feat)
    n_use <- length(idx)

    ## contiguous runs of idx for fast HDF5 reads
    runs <- .ic2_contig_runs(sort(idx))
    done_cells <- 0L
    for (run in runs) {
      starts <- seq.int(1L, length(run), by = chunk.size)
      for (s in starts) {
        e <- min(s + chunk.size - 1L, length(run))
        cols <- run[s:e]
        block <- .ic2_read_cols(handle, cols)
        if (!is.null(size_factors)) {
          block <- .ic2_scale_cols_by_sf(block, size_factors)
        }
        if (inherits(block, "dgCMatrix") && isTRUE(use.rcpp) &&
            exists("cpp_sparse_row_accum", mode = "function")) {
          cpp_sparse_row_accum(block, sumv, sumsq, nnzv)
        } else {
          ## R fallback
          sumv <- sumv + as.numeric(Matrix::rowSums(block))
          nnzv <- nnzv + as.integer(Matrix::rowSums(block != 0))
          ## sumsq
          if (inherits(block, "dgCMatrix")) {
            for (j in seq_len(ncol(block))) {
              if (block@p[j] < block@p[j + 1L]) {
                ii <- (block@p[j] + 1L):block@p[j + 1L]
                rows <- block@i[ii] + 1L
                vals <- block@x[ii]
                sumsq[rows] <- sumsq[rows] + vals * vals
              }
            }
          } else {
            sumsq <- sumsq + as.numeric(Matrix::rowSums(block * block))
          }
        }
        done_cells <- done_cells + length(cols)
        rep <<- .ic2_progress_bar(rep, done_cells, n_use, prefix = "  genes")
      }
    }
    list(
      genes = handle$feature_ids,
      sum = sumv,
      sumsq = sumsq,
      nnz = nnzv,
      n_cells = n_use,
      method = "chunked"
    )
  }

  .ic2_stats_table <- function(acc, condition = "all") {
    n <- acc$n_cells
    meanExp <- if (n > 0) acc$sum / n else rep(0, length(acc$sum))
    ## population SD across cells (including zeros)
    varv <- if (n > 1) {
      pmax(0, (acc$sumsq / n) - (meanExp * meanExp)) * (n / (n - 1))
    } else {
      rep(0, length(acc$sum))
    }
    SDs <- sqrt(varv)
    data.frame(
      genes = acc$genes,
      numberOfCells = as.integer(acc$nnz),
      totalNumberOfCells = as.integer(n),
      percentOfCells = if (n > 0) (as.numeric(acc$nnz) / n) * 100 else 0,
      meanExp = as.numeric(meanExp),
      SDs = as.numeric(SDs),
      condition = condition,
      stringsAsFactors = FALSE
    )
  }

  ## Global stats across all handles
  bulk_eligible <- length(handles) == 1L &&
    identical(handles[[1L]]$backend %||% NA_character_, "tenx_h5")
  .ic2_step(rep, paste0(
    "[1/2] ",
    if (bulk_eligible) {
      "Slabbed CSC gene stats (tenx_h5; covering spans, low RAM)"
    } else {
      "Streaming assay for global gene stats"
    }
  ))
  global_acc <- NULL
  for (i in seq_along(handles)) {
    acc <- .ic2_accum_handle(handles[[i]], cell_mask = cell_mask)
    if (identical(acc$method, "bulk_csc")) {
      .ic2_step(rep, paste0(
        "Handle ", i, "/", length(handles), ": bulk CSC  |  cells ",
        .ic2_fmt_num(acc$n_cells)
      ))
    }
    if (is.null(global_acc)) {
      global_acc <- acc
    } else {
      ## align by gene name
      if (!identical(global_acc$genes, acc$genes)) {
        common <- intersect(global_acc$genes, acc$genes)
        m1 <- match(common, global_acc$genes)
        m2 <- match(common, acc$genes)
        global_acc$sum[m1] <- global_acc$sum[m1] + acc$sum[m2]
        global_acc$sumsq[m1] <- global_acc$sumsq[m1] + acc$sumsq[m2]
        global_acc$nnz[m1] <- global_acc$nnz[m1] + acc$nnz[m2]
        global_acc$n_cells <- global_acc$n_cells + acc$n_cells
      } else {
        global_acc$sum <- global_acc$sum + acc$sum
        global_acc$sumsq <- global_acc$sumsq + acc$sumsq
        global_acc$nnz <- global_acc$nnz + acc$nnz
        global_acc$n_cells <- global_acc$n_cells + acc$n_cells
      }
    }
  }
  tables <- list(.ic2_stats_table(global_acc, "all"))

  if (!is.null(group.by)) {
    .ic2_step(rep, paste0("[2/2] Per-group stats by ", group.by))
    levels <- unique(as.character(x@cell_data[[group.by]]))
    levels <- levels[!is.na(levels) & nzchar(levels)]
    for (gi in seq_along(levels)) {
      g <- levels[gi]
      cells_g <- rownames(x@cell_data)[as.character(x@cell_data[[group.by]]) == g]
      acc_g <- NULL
      for (h in handles) {
        acc <- .ic2_accum_handle(h, cell_mask = cells_g)
        if (is.null(acc_g)) {
          acc_g <- acc
        } else if (identical(acc_g$genes, acc$genes)) {
          acc_g$sum <- acc_g$sum + acc$sum
          acc_g$sumsq <- acc_g$sumsq + acc$sumsq
          acc_g$nnz <- acc_g$nnz + acc$nnz
          acc_g$n_cells <- acc_g$n_cells + acc$n_cells
        }
      }
      tables[[length(tables) + 1L]] <- .ic2_stats_table(acc_g, g)
      rep <- .ic2_progress_bar(rep, gi, length(levels), prefix = "  groups")
    }
  } else {
    .ic2_step(rep, "[2/2] No group.by  -  global stats only")
  }

  Table <- do.call(rbind, tables)
  rownames(Table) <- NULL

  if (!update.object) {
    .ic2_finish(rep)
    return(Table)
  }

  ## Store on object
  x@modalities[[modality]]$gene.data <- Table
  ## also compact columns onto feature_data when condition == all
  all_tab <- Table[Table$condition == "all", , drop = FALSE]
  if (!is.null(x@modalities[[modality]]$feature_data) && nrow(all_tab) > 0L) {
    fd <- x@modalities[[modality]]$feature_data
    mm <- match(rownames(fd), all_tab$genes)
    fd$numberOfCells <- NA_integer_
    fd$percentOfCells <- NA_real_
    fd$meanExp <- NA_real_
    fd$SDs <- NA_real_
    ok <- !is.na(mm)
    fd$numberOfCells[ok] <- all_tab$numberOfCells[mm[ok]]
    fd$percentOfCells[ok] <- all_tab$percentOfCells[mm[ok]]
    fd$meanExp[ok] <- all_tab$meanExp[mm[ok]]
    fd$SDs[ok] <- all_tab$SDs[mm[ok]]
    x@modalities[[modality]]$feature_data <- fd
  }

  stats_stem <- file.path(
    x@storage$root_dir, "metadata", paste0(modality, "_gene_stats")
  )
  dir.create(dirname(stats_stem), recursive = TRUE, showWarnings = FALSE)
  stats_tsv <- .ic2_soda_write_table(Table, stats_stem)
  stats_path <- stats_tsv
  if (!is.null(x@storage$metadata_paths)) {
    x@storage$metadata_paths$gene_stats <- stats_path
  }

  x@params$feature_selection[[modality]] <- list(
    gene.stats = list(
      assay.name = assay.name,
      group.by = group.by,
      n_genes = length(unique(Table$genes[Table$condition == "all"])),
      chunk.size = chunk.size,
      path = stats_path
    )
  )

  x <- .ic2_log_command(
    x, "gene.stats",
    arguments = list(
      modality = modality, assay.name = assay.name,
      use.raw = use.raw, group.by = group.by, chunk.size = chunk.size, nnz.budget = nnz.budget
    ),
    outputs = list(path = stats_path, n_rows = nrow(Table))
  )

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  rep <- .ic2_add_output(rep, "Gene stats", stats_path)
  wall <- as.numeric(difftime(Sys.time(), rep$wall0, units = "secs"))
  if (is.finite(wall) && wall > 0) {
    rep$stats$cells_per_sec <- n_cells_guess / wall
  }
  .ic2_finish(rep, extra_lines = paste0(
    "Genes: ", .ic2_fmt_num(length(unique(Table$genes))),
    "   |   rows: ", .ic2_fmt_num(nrow(Table))
  ))
  x
}
