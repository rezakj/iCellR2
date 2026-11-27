#' Filter HTO assignments and export per-sample RNA
#'
#' Companion to \code{\link{hto.anno}}. Applies classic tutorial filters
#' (\code{percent.match > min.percent}, optional coverage / good-assignment
#' rules), then writes one RNA matrix per hashtag  -  same layout as the
#' classic \code{RNAs_<tag>_<n>/} loop, with optional 10x H5 output.
#'
#' \strong{Two call styles}
#' \itemize{
#'   \item \strong{Classic:} \code{hto.demux(anno = data, rna.data = MyRNAs)}.
#'     \code{rna.data} may be a matrix or anything \code{\link{load10x}}
#'     accepts: a virtual sample, a 10x HDF5 file, or a 10x directory.
#'   \item \strong{SODA:} \code{hto.demux(x = obj)} streams RNA columns per
#'     assignment (never densifies the full matrix).
#' }
#' \code{out.format = "10x_hdf5"} (aliases \code{"h5"}, \code{"tenx_h5"})
#' writes one HDF5 file per hashtag by streaming cell chunks.
#'
#' @param x Optional \code{iCellR2} object already annotated by
#'   \code{hto.anno} (needs \code{hto_assignment} / \code{assignment.annotation}
#'   in \code{cell_data}).
#' @param anno Classic annotation \code{data.frame} from \code{hto.anno}
#'   (cells  x  columns including \code{assignment.annotation}). Required when
#'   \code{x} is missing.
#' @param rna.data RNA counts used when \code{x} is missing. A genes  x  cells
#'   matrix, a \code{load10x} virtual sample, a 10x HDF5 path, or a 10x
#'   directory. Paths and virtual samples are streamed.
#' @param modality Modality for object mode. Default \code{"rna"}.
#' @param assay.name Source assay for RNA export (default active / counts).
#' @param min.percent Keep cells with \code{percent.match > min.percent}
#'   (classic uses strict \code{>}). Default \code{70}. Set \code{NULL} to skip.
#' @param min.coverage Keep cells with \code{coverage > min.coverage}.
#'   Default \code{NULL} (off).
#' @param keep.good.only If \code{TRUE}, also require
#'   \code{assignment.threshold == "good.assignment"}.
#' @param drop.low.cov If \code{TRUE}, drop cells flagged \code{low.cov}.
#' @param assignments Optional character vector of HTO tags to export.
#'   Default: all unique assignments after filtering.
#' @param subset.cells If \code{TRUE} (object mode), keep only passing cells
#'   in \code{cell_data}. Default \code{FALSE}.
#' @param write If \code{TRUE}, write per-HTO RNA matrices. Default \code{TRUE}.
#' @param out.dir Output directory. Classic default \code{"."}; SODA default
#'   \code{file.path(soda, "hto_split")}.
#' @param out.format \code{"10x_dir"} (MTX + barcodes + genes) or
#'   \code{"10x_hdf5"}. Aliases: \code{"h5"}, \code{"tenx_h5"}, \code{"mtx"}.
#' @param feature.type Feature type kept when \code{rna.data} is a path.
#'   Default \code{"Gene Expression"}. Ignored for a matrix or a virtual
#'   sample that is already loaded.
#' @param prefix Folder / file prefix. Default \code{"RNAs"}
#'    ->  \code{RNAs_<tag>_<n>}.
#' @param adt.pattern Features matching this regex are dropped from RNA when
#'   splitting a SODA object (HTO tags themselves are also dropped). Default
#'   \code{"TotalFF|TotalSeq|adt_|ADT"}. Set \code{NULL} to skip ADT drop.
#' @param hto.pattern Extra regex of features to drop from RNA (object mode).
#'   Default \code{NULL}; HTO features already annotated are always dropped.
#' @param write.tsv Write filtered annotation TSV. Default \code{TRUE}.
#' @param out.name Stem for filtered TSV. Default \code{"hto_demux"}.
#' @param verbose Progress messages.
#'
#' @return Classic: filtered annotation \code{data.frame} (invisible list
#'   attribute \code{split_paths} when written). Object: updated \code{x}
#'   with \code{hto_pass} and optional subset.
#'
#' @seealso \code{\link{hto.anno}}
#'
#' @examples
#' \dontrun{
#' anno <- hto.anno(hto.data = MyHTOs, plot = TRUE, plot.format = "pdf")
#' anno2 <- hto.demux(
#'   anno = anno,
#'   rna.data = MyRNAs,
#'   min.percent = 70,
#'   out.format = "10x_dir"
#' )
#' ## or H5:
#' hto.demux(
#'   anno = anno,
#'   rna.data = "filtered_feature_bc_matrix.h5",
#'   min.percent = 80,
#'   min.coverage = 100,
#'   out.format = "h5"
#' )
#'
#' obj <- hto.anno(obj, hto.pattern = "^CMO", feature.type = NULL)
#' obj <- hto.demux(obj, min.percent = 70, subset.cells = TRUE, out.format = "tenx_h5")
#' }
#'
#' @export
hto.demux <- function(
    x = NULL,
    anno = NULL,
    rna.data = NULL,
    modality = "rna",
    assay.name = NULL,
    min.percent = 70,
    min.coverage = NULL,
    keep.good.only = FALSE,
    drop.low.cov = FALSE,
    assignments = NULL,
    subset.cells = FALSE,
    write = TRUE,
    out.dir = NULL,
    out.format = "10x_dir",
    feature.type = "Gene Expression",
    prefix = "RNAs",
    adt.pattern = "TotalFF|TotalSeq|adt_|ADT",
    hto.pattern = NULL,
    write.tsv = TRUE,
    out.name = "hto_demux",
    verbose = TRUE
) {
  out.format <- .ic2_hto_out_format(out.format)

  ## ---- Classic path ----
  if (is.null(x)) {
    if (is.null(anno)) {
      stop("Provide annotation via 'anno=' (from hto.anno) or an iCellR2 object 'x'.")
    }
    anno <- .ic2_hto_normalize_anno(anno)
    n0 <- nrow(anno)
    pass <- .ic2_hto_filter_pass(
      anno,
      min.percent = min.percent,
      min.coverage = min.coverage,
      keep.good.only = keep.good.only,
      drop.low.cov = drop.low.cov
    )
    anno_f <- anno[pass, , drop = FALSE]
    if (!is.null(assignments) && length(assignments)) {
      keep_a <- as.character(anno_f$assignment.annotation) %in% as.character(assignments)
      anno_f <- anno_f[keep_a, , drop = FALSE]
    }
    if (isTRUE(verbose)) {
      message(
        "hto.demux: ", .ic2_fmt_num(n0), "  ->  ", .ic2_fmt_num(nrow(anno_f)),
        " cells after filters"
      )
    }
    if (nrow(anno_f) < 1L) {
      stop("No cells passed HTO demux filters.")
    }

    od <- out.dir %||% "."
    if (isTRUE(write.tsv)) {
      tsv <- if (grepl("[/\\\\]|\\.tsv$", out.name)) out.name else paste0(out.name, ".tsv")
      if (!grepl("\\.tsv$", tsv, ignore.case = TRUE)) tsv <- paste0(tsv, ".tsv")
      tsv <- if (identical(dirname(tsv), ".")) file.path(od, basename(tsv)) else tsv
      .ic2_hto_write_anno_tsv(anno_f, tsv)
      if (isTRUE(verbose)) message("Wrote ", normalizePath(tsv, winslash = "/", mustWork = FALSE))
    }

    split_paths <- character(0)
    if (isTRUE(write)) {
      if (is.null(rna.data)) {
        stop("Classic demux with write=TRUE needs rna.data= (matrix, virtual sample, H5, or 10x directory).")
      }
      if (.ic2_hto_rna_is_ondisk(rna.data)) {
        handles <- .ic2_hto_rna_handles(rna.data, feature.type = feature.type)
        split_paths <- .ic2_hto_split_ondisk(
          handles = handles,
          anno = anno_f,
          out.dir = od,
          format = out.format,
          prefix = prefix,
          adt.pattern = adt.pattern,
          hto.pattern = hto.pattern,
          verbose = verbose
        )
      } else {
        split_paths <- .ic2_hto_split_rna_matrix(
          rna.data = rna.data,
          anno = anno_f,
          out.dir = od,
          format = if (identical(out.format, "10x_hdf5")) "tenx_h5" else "10x_dir",
          prefix = prefix,
          verbose = verbose
        )
      }
      if (isTRUE(verbose)) {
        message("Wrote ", length(split_paths), " sample(s) under ", normalizePath(od, winslash = "/", mustWork = FALSE))
      }
    }
    attr(anno_f, "split_paths") <- split_paths
    return(anno_f)
  }

  ## ---- SODA object path ----
  .ic2_require_object(x)
  modality <- tolower(modality)
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  rep <- .ic2_new_reporter("hto.demux", verbose, soda_root = x@storage$root_dir)

  cd <- x@cell_data
  if (!("hto_assignment" %in% names(cd) || "assignment.annotation" %in% names(cd))) {
    stop("Run hto.anno(x) first so cell_data has HTO assignments.")
  }
  ## Build anno frame aligned to cell_data
  anno <- data.frame(
    assignment.annotation = as.character(
      if ("hto_assignment" %in% names(cd)) cd$hto_assignment else cd$assignment.annotation
    ),
    percent.match = as.numeric(
      if ("hto_percent_match" %in% names(cd)) cd$hto_percent_match else cd$percent.match
    ),
    coverage = as.numeric(
      if ("hto_coverage" %in% names(cd)) cd$hto_coverage else cd$coverage
    ),
    low.cov = as.logical(
      if ("hto_low_cov" %in% names(cd)) cd$hto_low_cov else cd$low.cov
    ),
    assignment.threshold = as.character(
      if ("hto_threshold" %in% names(cd)) cd$hto_threshold else cd$assignment.threshold
    ),
    row.names = rownames(cd),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  n0 <- nrow(anno)
  pass <- .ic2_hto_filter_pass(
    anno,
    min.percent = min.percent,
    min.coverage = min.coverage,
    keep.good.only = keep.good.only,
    drop.low.cov = drop.low.cov
  )
  anno$hto_pass <- pass
  if (!is.null(assignments) && length(assignments)) {
    anno$hto_pass <- anno$hto_pass &
      (as.character(anno$assignment.annotation) %in% as.character(assignments))
  }

  x@cell_data$hto_pass <- as.logical(anno$hto_pass)
  n_pass <- sum(x@cell_data$hto_pass, na.rm = TRUE)
  if (n_pass < 1L) stop("No cells passed HTO demux filters.")

  .ic2_step(rep, paste0(
    "[1/3] Filters: ", .ic2_fmt_num(n0), "  ->  ", .ic2_fmt_num(n_pass), " cells"
  ))
  rep <- .ic2_banner_start(rep, c(
    paste0("Pass: ", .ic2_fmt_num(n_pass), " / ", .ic2_fmt_num(n0)),
    paste0(
      "min.percent=", min.percent %||% "NULL",
      "  min.coverage=", min.coverage %||% "NULL",
      "  good.only=", keep.good.only,
      "  drop.low.cov=", drop.low.cov
    ),
    paste0("Export: ", if (isTRUE(write)) out.format else "none",
           "  subset.cells=", isTRUE(subset.cells))
  ))

  od <- out.dir %||% file.path(x@storage$root_dir, "hto_split")
  meta_dir <- file.path(x@storage$root_dir, "metadata", "hto")
  dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

  anno_f <- anno[anno$hto_pass, , drop = FALSE]
  tsv_path <- NULL
  if (isTRUE(write.tsv)) {
    tsv_path <- file.path(meta_dir, paste0(out.name, ".filtered.tsv"))
    .ic2_hto_write_anno_tsv(anno_f, tsv_path)
  }

  if (isTRUE(subset.cells)) {
    .ic2_step(rep, "[2/3] Subsetting cell_data to passing cells")
    keep <- which(x@cell_data$hto_pass)
    x@cell_data <- x@cell_data[keep, , drop = FALSE]
    anno_f <- anno_f[rownames(x@cell_data), , drop = FALSE]
  } else {
    .ic2_step(rep, "[2/3] Keeping full cell_data (hto_pass flagged)")
  }
  .ic2_persist_cell_data(x)

  split_paths <- character(0)
  if (isTRUE(write)) {
    .ic2_step(rep, paste0("[3/3] Writing per-HTO RNA (", out.format, ")"))
    handles <- .ic2_get_assay_handles(x, modality, assay.name)
    all_feats <- character(0)
    for (h in handles) {
      h2 <- .ic2_enrich_handle(h, load_dimnames = TRUE)
      all_feats <- union(all_feats, as.character(h2$feature_ids))
    }
    drop <- character(0)
    ## Prefer stored HTO feature list from hto.anno params
    hto_feats <- tryCatch(
      as.character(x@params$hto[[modality]]$hto.features),
      error = function(e) character(0)
    )
    if (length(hto_feats)) drop <- union(drop, hto_feats)
    if (!is.null(adt.pattern) && nzchar(adt.pattern)) {
      drop <- union(drop, grep(adt.pattern, all_feats, value = TRUE, ignore.case = TRUE))
    }
    if (!is.null(hto.pattern) && nzchar(hto.pattern)) {
      drop <- union(drop, grep(hto.pattern, all_feats, value = TRUE, ignore.case = TRUE))
    }
    ## Always drop assignment tag names if they appear as features
    drop <- union(drop, unique(as.character(anno_f$assignment.annotation)))
    rna_feats <- setdiff(all_feats, drop)
    if (length(rna_feats) < 1L) {
      stop("No RNA features left after excluding HTO/ADT; check patterns.")
    }
    split_paths <- .ic2_hto_split_ondisk(
      handles = handles,
      anno = anno_f,
      out.dir = od,
      format = out.format,
      prefix = prefix,
      adt.pattern = NULL,
      hto.pattern = NULL,
      features = rna_feats,
      verbose = verbose,
      progress = function(cur, tot) {
        rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  split")
      }
    )
    for (p in split_paths) rep <- .ic2_add_output(rep, "HTO demux", p)
  } else {
    .ic2_step(rep, "[3/3] Skip RNA export (write = FALSE)")
  }

  x@params$hto_demux <- x@params$hto_demux %||% list()
  x@params$hto_demux[[modality]] <- list(
    assay.name = assay.name,
    min.percent = min.percent,
    min.coverage = min.coverage,
    keep.good.only = isTRUE(keep.good.only),
    drop.low.cov = isTRUE(drop.low.cov),
    subset.cells = isTRUE(subset.cells),
    out.format = if (isTRUE(write)) out.format else NULL,
    out.dir = if (isTRUE(write)) od else NULL,
    n_cells_in = n0,
    n_cells_pass = n_pass,
    n_cells_out = nrow(x@cell_data),
    tsv = tsv_path,
    split_paths = split_paths
  )

  x <- .ic2_log_command(
    x,
    "hto.demux",
    arguments = list(
      modality = modality,
      min.percent = min.percent,
      min.coverage = min.coverage,
      subset.cells = isTRUE(subset.cells),
      write = isTRUE(write),
      out.format = out.format
    ),
    outputs = list(
      n_pass = n_pass,
      n_cells = nrow(x@cell_data),
      tsv = tsv_path,
      n_splits = length(split_paths)
    )
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)

  if (!is.null(tsv_path)) rep <- .ic2_add_output(rep, "Filtered TSV", tsv_path)
  rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  rep <- .ic2_add_scale(rep, c(
    paste0("Pass:     ", .ic2_fmt_num(n_pass), " / ", .ic2_fmt_num(n0)),
    paste0("Samples:  ", length(split_paths)),
    paste0("Cells:    ", .ic2_fmt_num(nrow(x@cell_data)))
  ))
  .ic2_finish(rep, extra_lines = paste0(
    "HTO demux: ", .ic2_fmt_num(n_pass), " cells  ->  ",
    length(split_paths), " sample(s)"
  ))
  x
}


#' Normalize classic / object-style HTO annotation tables
#' @keywords internal
.ic2_hto_normalize_anno <- function(anno) {
  if (is.null(anno) || !nrow(anno)) stop("'anno' is empty.")
  anno <- as.data.frame(anno, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(rownames(anno)) || !length(rownames(anno))) {
    if ("ID" %in% names(anno)) {
      rownames(anno) <- as.character(anno$ID)
    } else {
      stop("Annotation needs rownames (cell ids) or an ID column.")
    }
  }
  if (!("assignment.annotation" %in% names(anno)) && "hto_assignment" %in% names(anno)) {
    anno$assignment.annotation <- as.character(anno$hto_assignment)
  }
  if (!("percent.match" %in% names(anno)) && "hto_percent_match" %in% names(anno)) {
    anno$percent.match <- as.numeric(anno$hto_percent_match)
  }
  if (!("coverage" %in% names(anno)) && "hto_coverage" %in% names(anno)) {
    anno$coverage <- as.numeric(anno$hto_coverage)
  }
  if (!("low.cov" %in% names(anno)) && "hto_low_cov" %in% names(anno)) {
    anno$low.cov <- as.logical(anno$hto_low_cov)
  }
  if (!("assignment.threshold" %in% names(anno)) && "hto_threshold" %in% names(anno)) {
    anno$assignment.threshold <- as.character(anno$hto_threshold)
  }
  if (!("assignment.annotation" %in% names(anno))) {
    stop("Annotation needs assignment.annotation (run hto.anno first).")
  }
  anno
}


#' Classic HTO filter mask (strict > for percent / coverage)
#' @keywords internal
.ic2_hto_filter_pass <- function(
    anno,
    min.percent = NULL,
    min.coverage = NULL,
    keep.good.only = FALSE,
    drop.low.cov = FALSE
) {
  n <- nrow(anno)
  pass <- rep(TRUE, n)
  if (!is.null(min.percent)) {
    pm <- as.numeric(anno$percent.match)
    pass <- pass & (pm > as.numeric(min.percent))
  }
  if (!is.null(min.coverage)) {
    cv <- as.numeric(anno$coverage)
    pass <- pass & (cv > as.numeric(min.coverage))
  }
  if (isTRUE(keep.good.only)) {
    thr <- as.character(anno$assignment.threshold)
    pass <- pass & (thr == "good.assignment")
  }
  if (isTRUE(drop.low.cov)) {
    pass <- pass & !as.logical(anno$low.cov)
  }
  pass[is.na(pass)] <- FALSE
  pass
}


#' Normalize hto.demux output format names.
#' @keywords internal
.ic2_hto_out_format <- function(out.format) {
  key <- tolower(as.character(out.format)[1L])
  if (key %in% c("10x_hdf5", "tenx_h5", "h5", "hdf5", "10x_h5")) return("10x_hdf5")
  if (key %in% c("10x_dir", "mtx", "dir")) return("10x_dir")
  stop("'out.format' must be \"10x_dir\" or \"10x_hdf5\" (aliases: h5, tenx_h5, mtx).")
}


#' TRUE when rna.data is a path or a load10x handle rather than a matrix.
#' @keywords internal
.ic2_hto_rna_is_ondisk <- function(rna.data) {
  inherits(rna.data, "iCellR2_sample") ||
    inherits(rna.data, "iCellR2_assay_handle") ||
    (is.character(rna.data) && length(rna.data) == 1L && !is.na(rna.data) &&
       file.exists(rna.data))
}


#' Resolve rna.data to assay handles. Paths go through load10x virtual.
#' @keywords internal
.ic2_hto_rna_handles <- function(rna.data, feature.type = "Gene Expression") {
  if (inherits(rna.data, "iCellR2_assay_handle")) return(list(rna.data))
  if (inherits(rna.data, "iCellR2_sample") || (is.list(rna.data) && !is.null(rna.data$handle))) {
    return(list(rna.data$handle))
  }
  if (is.character(rna.data) && length(rna.data) == 1L && file.exists(rna.data)) {
    samp <- load10x(
      input = rna.data,
      sample.id = "hto_demux_rna",
      output = "virtual",
      feature.type = feature.type,
      overwrite = TRUE,
      verbose = FALSE
    )
    return(list(samp$handle))
  }
  stop("'rna.data' must be a matrix, a load10x virtual sample, an HDF5 file, or a 10x directory.")
}


#' Stream one RNA file per HTO assignment.
#' @keywords internal
.ic2_hto_split_ondisk <- function(
    handles,
    anno,
    out.dir = ".",
    format = "10x_hdf5",
    prefix = "RNAs",
    adt.pattern = NULL,
    hto.pattern = NULL,
    features = NULL,
    verbose = TRUE,
    progress = NULL
) {
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  handles <- lapply(handles, function(h) .ic2_enrich_handle(h, load_dimnames = TRUE))
  if (is.null(features)) {
    features <- character(0)
    for (h in handles) features <- union(features, as.character(h$feature_ids))
    drop <- character(0)
    if (!is.null(adt.pattern) && nzchar(adt.pattern)) {
      drop <- union(drop, grep(adt.pattern, features, value = TRUE, ignore.case = TRUE))
    }
    if (!is.null(hto.pattern) && nzchar(hto.pattern)) {
      drop <- union(drop, grep(hto.pattern, features, value = TRUE, ignore.case = TRUE))
    }
    drop <- union(drop, unique(as.character(anno$assignment.annotation)))
    features <- setdiff(features, drop)
  }
  features <- as.character(features)
  if (length(features) < 1L) stop("No RNA features left to export.")

  assign_chr <- as.character(anno$assignment.annotation)
  tags <- unique(assign_chr)
  tags <- tags[!is.na(tags) & nzchar(tags)]
  if (length(tags) < 1L) stop("No HTO assignments to demux.")
  cells_by_tag <- split(rownames(anno), factor(assign_chr, levels = tags))
  out_paths <- character(0)
  for (i in seq_along(tags)) {
    tag <- tags[[i]]
    if (is.function(progress)) progress(i - 1L, length(tags))
    cells <- cells_by_tag[[tag]]
    if (isTRUE(verbose)) message(" getting sample ", tag, " ...")
    path <- .ic2_hto_stream_tag(
      handles = handles,
      cells = cells,
      features = features,
      out.dir = out.dir,
      prefix = prefix,
      tag = tag,
      format = format,
      verbose = verbose
    )
    if (is.null(path)) {
      if (isTRUE(verbose)) message("  (no overlapping cells; skip)")
      next
    }
    out_paths <- c(out_paths, path)
  }
  if (is.function(progress)) progress(length(tags), length(tags))
  out_paths
}


#' Remap a sparse block onto a fixed feature order (drop non-matching rows).
#' @keywords internal
.ic2_hto_remap_features <- function(raw, features, feature_ids = NULL) {
  raw <- .ic2_to_dgCMatrix(raw)
  rn <- rownames(raw)
  if (is.null(rn)) rn <- as.character(feature_ids)
  features <- as.character(features)
  nc <- ncol(raw)
  if (identical(as.character(rn), features)) return(raw)
  dest_row <- match(as.character(rn), features)
  ii <- raw@i + 1L
  keep <- !is.na(dest_row[ii])
  if (!any(keep)) {
    return(Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(length(features), nc),
      dimnames = list(features, colnames(raw))
    ))
  }
  jj <- rep.int(seq_len(nc), diff(raw@p))
  Matrix::sparseMatrix(
    i = dest_row[ii[keep]],
    j = jj[keep],
    x = raw@x[keep],
    dims = c(length(features), nc),
    dimnames = list(features, colnames(raw)),
    index1 = TRUE
  )
}


#' Stream selected cells from assay handles into one 10x directory or HDF5.
#'
#' Single-pass read: accumulate remapped CSC pieces, then one bulk write.
#' Avoids the old nnz double-read and per-column HDF5 hyperslab writes.
#' @keywords internal
.ic2_hto_stream_tag <- function(
    handles, cells, features, out.dir, prefix, tag, format, chunk.cells = 4000L, verbose = TRUE
) {
  cells <- unique(as.character(cells))
  cells <- cells[!is.na(cells) & nzchar(cells)]
  features <- as.character(features)
  plan <- list()
  left <- cells
  for (h in handles) {
    if (length(left) < 1L) break
    ids <- as.character(h$cell_ids)
    hit <- match(left, ids)
    if (all(is.na(hit))) hit <- match(.ic2_clean_names(left), .ic2_clean_names(ids))
    take <- which(!is.na(hit))
    if (length(take) < 1L) next
    cols <- as.integer(hit[take])
    ord <- order(cols)
    plan[[length(plan) + 1L]] <- list(handle = h, cols = cols[ord], ids = ids[cols[ord]])
    left <- left[-take]
  }
  if (length(plan) < 1L) return(NULL)
  chunk.cells <- as.integer(max(1L, chunk.cells))

  barcodes <- unlist(lapply(plan, function(p) p$ids), use.names = FALSE)
  nc <- length(barcodes)
  tag_clean <- .ic2_clean_names(tag)[1L]
  dest <- file.path(out.dir, paste(prefix, tag_clean, nc, sep = "_"))
  if (isTRUE(verbose)) {
    message(" number of cells ", nc)
    message(" writing sample ", tag, " ...")
  }

  ## Single pass: read each selected column once, remap features, stash CSC pieces.
  i_parts <- list()
  x_parts <- list()
  p_out <- numeric(nc + 1L)
  p_out[1L] <- 0
  col0 <- 0L
  nnz_acc <- 0
  part_n <- 0L
  for (p in plan) {
    starts <- seq.int(1L, length(p$cols), by = chunk.cells)
    for (s in starts) {
      e <- min(s + chunk.cells - 1L, length(p$cols))
      raw <- .ic2_read_cols(p$handle, p$cols[s:e])
      sub <- .ic2_hto_remap_features(raw, features, feature_ids = p$handle$feature_ids)
      n_sub <- ncol(sub)
      if (length(sub@x)) {
        part_n <- part_n + 1L
        i_parts[[part_n]] <- as.integer(sub@i)
        x_parts[[part_n]] <- as.numeric(sub@x)
        dp <- as.numeric(diff(sub@p))
        p_out[(col0 + 2L):(col0 + n_sub + 1L)] <- nnz_acc + cumsum(dp)
        nnz_acc <- nnz_acc + sum(dp)
      } else if (n_sub > 0L) {
        p_out[(col0 + 2L):(col0 + n_sub + 1L)] <- nnz_acc
      }
      col0 <- col0 + n_sub
    }
  }
  i_all <- if (part_n) unlist(i_parts, use.names = FALSE) else integer(0)
  x_all <- if (part_n) unlist(x_parts, use.names = FALSE) else numeric(0)

  if (identical(format, "10x_dir")) {
    mat <- Matrix::sparseMatrix(
      i = i_all,
      p = as.integer(p_out),
      x = x_all,
      dims = c(length(features), nc),
      dimnames = list(features, barcodes),
      index1 = FALSE
    )
    .ic2_write_10x_mtx_dir(mat, dest)
    return(dest)
  }

  dest_h5 <- paste0(dest, ".h5")
  .ic2_write_tenx_hdf5(
    filepath = dest_h5,
    i = i_all,
    p = as.integer(p_out),
    x = x_all,
    barcodes = as.character(barcodes),
    feature_ids = features,
    feature_names = features,
    feature_types = "Gene Expression",
    genome = "unknown",
    overwrite = TRUE
  )
  dest_h5
}
