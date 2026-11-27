#' Annotate hashtag oligos (HTOs)
#'
#' Classic iCellR assignment rule, SODA-aware for \code{iCellR2}:
#' for each cell, take the dominant HTO count as a percent of total HTO
#' coverage; mark \code{good.assignment} when that percent is at least
#' \code{assignment.thr}, and \code{low.cov} when total coverage is below
#' \code{mean(coverage) / cov.thr}.
#'
#' This function \strong{only annotates} (and optionally writes TSV / QC plot).
#' Filtering and per-sample RNA export live in \code{\link{hto.demux}}.
#'
#' \strong{Two call styles}
#' \itemize{
#'   \item \strong{Classic matrix:} \code{hto.anno(hto.data = MyHTOs)} with
#'     HTOs  x  cells, or a \code{load10x} virtual sample
#'     (\code{feature.type = "Antibody Capture"}). Returns a data frame.
#'   \item \strong{iCellR2 object:} \code{hto.anno(x = obj)} streams only the
#'     HTO features from the on-disk assay, writes assignment columns into
#'     \code{cell_data}, and saves a TSV under SODA \code{metadata/}.
#' }
#'
#' HTO features are resolved in order: explicit \code{hto.features}, else
#' 10x \code{feature.type} (default Antibody / Multiplexing Capture), else
#' \code{hto.pattern} (CMO / TotalSeq / Hashtag / HTO).
#'
#' @param x Optional \code{iCellR2} object.
#' @param hto.data Optional HTO counts (features  x  cells), or a
#'   \code{load10x(..., output = "virtual")} sample. Required when \code{x}
#'   is missing.
#' @param modality Modality for object mode. Default \code{"rna"}.
#' @param assay.name Source assay (default active / \code{"counts"}).
#' @param hto.features Character vector of HTO feature names. Optional.
#' @param hto.pattern Regex used when \code{hto.features} is \code{NULL}
#'   and feature types are unavailable.
#' @param feature.type 10x feature type(s) to treat as HTOs. Set \code{NULL}
#'   to skip and use \code{hto.pattern} only.
#' @param cov.thr Divide mean HTO coverage by this to set the low-coverage
#'   threshold. Default \code{10}.
#' @param assignment.thr Percent match (max/sum x 100) at or above which a
#'   cell is \code{good.assignment}. Default \code{80}.
#' @param write.tsv Write annotation table (classic: \code{ID} + columns;
#'   uses \code{data.table::fwrite} when available).
#' @param out.name Stem for the TSV / plot. Classic example:
#'   \code{"HTOs_annotated_HSThigh"}.
#' @param plot If \code{TRUE}, write classic two-panel QC plot.
#' @param plot.format \code{"png"} (default) or \code{"pdf"}.
#' @param max.plot.points Max jitter points in QC plots. Default \code{25000}.
#' @param verbose Progress style.
#'
#' @return If \code{x} is missing: a data.frame (cells  x  HTO counts +
#'   annotation columns). If \code{x} is an \code{iCellR2} object: the
#'   updated object with HTO columns in \code{cell_data}.
#'
#' @seealso \code{\link{hto.demux}} for filter + per-HTO RNA export.
#'
#' @examples
#' \dontrun{
#' ## Classic annotate
#' anno <- hto.anno(
#'   hto.data = MyHTOs,
#'   cov.thr = 10,
#'   assignment.thr = 80,
#'   out.name = "HTOs_annotated_HSThigh",
#'   plot = TRUE,
#'   plot.format = "pdf"
#' )
#'
#' ## Then filter + split
#' hto.demux(anno = anno, rna.data = MyRNAs, min.percent = 70, out.format = "10x_dir")
#'
#' ## SODA
#' obj <- hto.anno(obj, hto.pattern = "^CMO", feature.type = NULL, plot = TRUE)
#' obj <- hto.demux(obj, min.percent = 70, subset.cells = TRUE, out.format = "tenx_h5")
#' }
#'
#' @export
hto.anno <- function(
    x = NULL,
    hto.data = NULL,
    modality = "rna",
    assay.name = NULL,
    hto.features = NULL,
    hto.pattern = "^CMO|CMO|TotalSeq|Hashtag|hashtag|HTO",
    feature.type = c("Antibody Capture", "Multiplexing Capture"),
    cov.thr = 10,
    assignment.thr = 80,
    write.tsv = TRUE,
    out.name = "hto_anno",
    plot = FALSE,
    plot.format = c("png", "pdf"),
    max.plot.points = 25000L,
    verbose = TRUE
) {
  plot.format <- match.arg(plot.format)

  ## ---- Classic matrix path (annotate only) ----
  if (is.null(x)) {
    if (is.null(hto.data)) {
      stop("Provide either an iCellR2 object 'x' or classic 'hto.data' (HTOs  x  cells).")
    }
    hto.data <- .ic2_hto_as_matrix(
      hto.data,
      hto.features = hto.features,
      hto.pattern = hto.pattern
    )
    anno <- .ic2_hto_anno_core(
      hto.data = hto.data,
      cov.thr = cov.thr,
      assignment.thr = assignment.thr
    )
    if (isTRUE(write.tsv)) {
      out <- if (!is.null(out.name) && grepl("[/\\\\]|\\.tsv$", out.name)) {
        out.name
      } else {
        paste0(out.name %||% "HTOs_annotated", ".tsv")
      }
      if (!grepl("\\.tsv$", out, ignore.case = TRUE)) out <- paste0(out, ".tsv")
      .ic2_hto_write_anno_tsv(anno, out)
      if (isTRUE(verbose)) message("Wrote ", normalizePath(out, winslash = "/", mustWork = FALSE))
    }
    if (isTRUE(plot)) {
      pfile <- paste0(
        sub("\\.tsv$", "", out.name %||% "HTO_stats", ignore.case = TRUE),
        "_stats.", plot.format
      )
      if (identical(out.name, "hto_anno") || is.null(out.name)) {
        pfile <- paste0("HTO_stats.", plot.format)
      }
      .ic2_hto_anno_plot(anno, pfile, format = plot.format, max.points = max.plot.points)
      if (isTRUE(verbose)) message("Wrote ", normalizePath(pfile, winslash = "/", mustWork = FALSE))
    }
    return(anno)
  }

  .ic2_require_object(x)
  modality <- tolower(modality)
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  rep <- .ic2_new_reporter("hto.anno", verbose, soda_root = x@storage$root_dir)

  handles <- .ic2_get_assay_handles(x, modality, assay.name)
  rep <- .ic2_register_assay_inputs(
    rep, x, modality, assay.name, label = paste0("Assay[", assay.name, "]")
  )

  cell_ids <- rownames(x@cell_data)
  if (is.null(cell_ids) || length(cell_ids) < 1L) {
    stop("x@cell_data has no cell ids.")
  }

  .ic2_step(rep, "[1/3] Resolving HTO features")
  hto_feats <- .ic2_resolve_hto_features(
    handles = handles,
    hto.features = hto.features,
    hto.pattern = hto.pattern,
    feature.type = feature.type
  )
  rep <- .ic2_banner_start(rep, c(
    paste0("HTOs: ", length(hto_feats), "   |   cells: ", .ic2_fmt_num(length(cell_ids))),
    paste0("Assay: ", assay.name, "   |   cov.thr = ", cov.thr, "   |   assignment.thr = ", assignment.thr),
    paste0("Features: ", paste(utils::head(hto_feats, 5L), collapse = ", "),
           if (length(hto_feats) > 5L) ", ..." else "")
  ))

  .ic2_step(rep, paste0(
    "[2/3] Streaming ", length(hto_feats), " HTO feature(s)  x  ",
    .ic2_fmt_num(length(cell_ids)), " cells"
  ))
  hto_cx <- .ic2_load_features_for_cells(
    handles = handles,
    features = hto_feats,
    cell_ids = cell_ids,
    chunk.size = 200000L,
    transform = "none",
    size_factors = NULL
  )
  hto_mat <- t(hto_cx)
  rownames(hto_mat) <- colnames(hto_cx)
  colnames(hto_mat) <- rownames(hto_cx)

  anno <- .ic2_hto_anno_core(
    hto.data = hto_mat,
    cov.thr = cov.thr,
    assignment.thr = assignment.thr
  )
  anno <- anno[cell_ids, , drop = FALSE]

  .ic2_step(rep, "[3/3] Writing TSV + QC plot + cell_data")
  out_dir <- file.path(x@storage$root_dir, "metadata", "hto")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tsv_path <- NULL
  if (isTRUE(write.tsv)) {
    tsv_path <- file.path(out_dir, paste0(out.name, ".tsv"))
    .ic2_hto_write_anno_tsv(anno, tsv_path)
  }

  plot_path <- NULL
  if (isTRUE(plot)) {
    plot_path <- .ic2_hto_anno_plot(
      anno,
      file.path(out_dir, paste0(out.name, "_stats.", plot.format)),
      format = plot.format,
      max.points = max.plot.points
    )
  }

  x@cell_data$hto_assignment <- as.character(anno$assignment.annotation)
  x@cell_data$hto_percent_match <- as.numeric(anno$percent.match)
  x@cell_data$hto_coverage <- as.numeric(anno$coverage)
  x@cell_data$hto_low_cov <- as.logical(anno$low.cov)
  x@cell_data$hto_threshold <- as.character(anno$assignment.threshold)
  x@cell_data$assignment.annotation <- x@cell_data$hto_assignment
  x@cell_data$percent.match <- x@cell_data$hto_percent_match
  x@cell_data$assignment.threshold <- x@cell_data$hto_threshold
  .ic2_persist_cell_data(x)

  x@params$hto <- x@params$hto %||% list()
  x@params$hto[[modality]] <- list(
    assay.name = assay.name,
    hto.features = hto_feats,
    cov.thr = cov.thr,
    assignment.thr = assignment.thr,
    n_cells = nrow(x@cell_data),
    n_good = sum(as.character(anno$assignment.threshold) == "good.assignment", na.rm = TRUE),
    tsv = tsv_path,
    plot = plot_path
  )

  x <- .ic2_log_command(
    x,
    "hto.anno",
    arguments = list(
      modality = modality,
      assay.name = assay.name,
      cov.thr = cov.thr,
      assignment.thr = assignment.thr,
      n_hto = length(hto_feats)
    ),
    outputs = list(
      n_cells = nrow(x@cell_data),
      tsv = tsv_path,
      plot = plot_path
    )
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)

  if (!is.null(tsv_path)) rep <- .ic2_add_output(rep, "HTO TSV", tsv_path)
  if (!is.null(plot_path)) rep <- .ic2_add_output(rep, "HTO plot", plot_path)
  rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  rep <- .ic2_add_scale(rep, c(
    paste0("HTOs:     ", length(hto_feats)),
    paste0(
      "Assign:   good=", sum(anno$assignment.threshold == "good.assignment", na.rm = TRUE),
      "  unsure=", sum(anno$assignment.threshold == "unsure", na.rm = TRUE),
      "  low.cov=", sum(anno$low.cov, na.rm = TRUE)
    ),
    paste0("Cells:    ", .ic2_fmt_num(nrow(x@cell_data)))
  ))
  .ic2_finish(rep, extra_lines = paste0(
    "HTO annotate: ", length(hto_feats), " tags  x  ",
    .ic2_fmt_num(nrow(x@cell_data)), " cells"
  ))
  x
}


#' HTOs  x  cells matrix from a matrix or a load10x virtual sample.
#' @keywords internal
.ic2_hto_as_matrix <- function(hto.data, hto.features = NULL, hto.pattern = NULL) {
  if (inherits(hto.data, "Matrix") || is.matrix(hto.data) || is.data.frame(hto.data)) {
    return(hto.data)
  }
  handle <- hto.data
  if (inherits(hto.data, "iCellR2_sample") || is.list(hto.data) && !is.null(hto.data$handle)) {
    handle <- hto.data$handle
  }
  if (!inherits(handle, "iCellR2_assay_handle") &&
      (is.null(handle$path) || is.null(handle$backend))) {
    stop("'hto.data' must be an HTOs  x  cells matrix or a load10x virtual sample.")
  }
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  feats <- as.character(handle$feature_ids)
  cells <- as.character(handle$cell_ids)
  if (!is.null(hto.features) && length(hto.features)) {
    feats <- intersect(feats, as.character(hto.features))
  } else if (!is.null(hto.pattern) && nzchar(hto.pattern)) {
    hit <- grep(hto.pattern, feats, value = TRUE, ignore.case = TRUE)
    if (length(hit) > 0L) feats <- hit
  }
  if (length(feats) < 1L || length(cells) < 1L) {
    stop("No HTO features found in 'hto.data'.")
  }
  cx <- .ic2_load_features_for_cells(
    handles = list(handle),
    features = feats,
    cell_ids = cells,
    transform = "none"
  )
  out <- t(cx)
  rownames(out) <- colnames(cx)
  colnames(out) <- rownames(cx)
  out
}

#' Classic HTO assignment on an HTOs  x  cells matrix
#' @keywords internal
.ic2_hto_anno_core <- function(hto.data, cov.thr = 10, assignment.thr = 80) {
  data <- hto.data
  if (inherits(data, "Matrix")) {
    data <- as.matrix(data)
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(data) < 1L || ncol(data) < 1L) {
    stop("'hto.data' must be a non-empty HTOs  x  cells matrix.")
  }
  ## Ensure numeric
  for (j in seq_len(ncol(data))) {
    data[[j]] <- as.numeric(data[[j]])
  }
  data[is.na(data)] <- 0

  ## cells  x  HTOs
  data_t <- as.data.frame(t(data), stringsAsFactors = FALSE, check.names = FALSE)
  mySum <- as.numeric(rowSums(data_t))
  neg <- mean(mySum) / as.numeric(cov.thr)
  myCov <- mySum < neg
  myMax <- as.numeric(apply(data_t, 1L, max))
  myRatio <- ifelse(mySum > 0, (myMax / mySum) * 100, 0)
  thr <- ifelse(myRatio >= as.numeric(assignment.thr), "good.assignment", "unsure")
  myWhich <- colnames(data_t)[max.col(as.matrix(data_t), ties.method = "first")]

  out <- cbind(
    data_t,
    assignment.annotation = myWhich,
    percent.match = myRatio,
    coverage = mySum,
    low.cov = myCov,
    assignment.threshold = thr,
    stringsAsFactors = FALSE
  )
  out
}


#' Resolve HTO feature names from handles / pattern / feature_type
#' @keywords internal
.ic2_resolve_hto_features <- function(
    handles,
    hto.features = NULL,
    hto.pattern = "^CMO|CMO|TotalSeq|Hashtag|hashtag|HTO",
    feature.type = c("Antibody Capture", "Multiplexing Capture")
) {
  if (!is.null(hto.features) && length(hto.features)) {
    return(unique(as.character(hto.features)))
  }

  all_feats <- character(0)
  typed <- character(0)
  ftypes <- feature.type
  if (!is.null(ftypes) && length(ftypes)) {
    ftypes <- as.character(ftypes)
    ftypes <- ftypes[nzchar(ftypes) & !is.na(ftypes)]
  } else {
    ftypes <- character(0)
  }

  for (h in handles) {
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    all_feats <- union(all_feats, as.character(h$feature_ids))
    if (length(ftypes) && identical(h$backend, "tenx_h5") &&
        !is.null(h$path) && file.exists(h$path)) {
      for (ft in ftypes) {
        info <- tryCatch(
          .ic2_tenx_h5_inspect(h$path, feature.field = "name", feature.type = ft),
          error = function(e) NULL
        )
        if (!is.null(info) && length(info$features)) {
          typed <- union(typed, as.character(info$features))
        }
      }
    }
  }
  if (length(typed) > 0L) {
    hit <- intersect(typed, all_feats)
    if (length(hit) > 0L) return(hit)
    return(typed)
  }

  if (is.null(hto.pattern) || !nzchar(hto.pattern)) {
    stop(
      "Could not resolve HTO features. Pass hto.features=, or ensure the ",
      "assay has 10x feature_type (Antibody Capture / Multiplexing Capture), ",
      "or set hto.pattern= (e.g. \"^CMO\")."
    )
  }
  hit <- grep(hto.pattern, all_feats, value = TRUE, ignore.case = TRUE)
  if (length(hit) < 1L) {
    stop(
      "No features matched hto.pattern = '", hto.pattern, "'. ",
      "For CellPlex use hto.pattern = \"^CMO\"; for TotalSeq use \"TotalSeq|Hashtag\"."
    )
  }
  hit
}


#' Write one RNA matrix per HTO assignment (classic split loop)
#' @keywords internal
.ic2_hto_split_rna_samples <- function(
    handles,
    anno,
    rna.features,
    out.dir = "hto_split",
    format = c("10x_dir", "tenx_h5"),
    prefix = "RNAs",
    progress = NULL
) {
  format <- match.arg(format)
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  tags <- unique(as.character(anno$assignment.annotation))
  tags <- tags[!is.na(tags) & nzchar(tags)]
  if (length(tags) < 1L) stop("No HTO assignments to split.")

  out_paths <- character(0)
  for (i in seq_along(tags)) {
    tag <- tags[[i]]
    cells <- rownames(anno)[as.character(anno$assignment.annotation) == tag]
    cells <- as.character(cells)
    if (length(cells) < 1L) next
    if (is.function(progress)) progress(i - 1L, length(tags))

    mat <- .ic2_gather_assay_cols(handles, cells, rna.features)
    mat <- .ic2_to_dgCMatrix(mat)
    tag_clean <- .ic2_clean_names(tag)[1L]
    name <- paste(prefix, tag_clean, ncol(mat), sep = "_")
    dest <- file.path(out.dir, name)

    if (identical(format, "10x_dir")) {
      .ic2_write_10x_mtx_dir(mat, dest)
      out_paths <- c(out_paths, dest)
    } else {
      dest_h5 <- paste0(dest, ".h5")
      .ic2_write_matrix_tenx_h5(
        mat,
        filepath = dest_h5,
        sample_id = tag_clean,
        modality = "rna",
        feature_types = "Gene Expression",
        overwrite = TRUE
      )
      out_paths <- c(out_paths, dest_h5)
    }
  }
  if (is.function(progress)) progress(length(tags), length(tags))
  out_paths
}


#' Classic split from an in-memory RNA matrix (features  x  cells)
#' @keywords internal
.ic2_hto_split_rna_matrix <- function(
    rna.data,
    anno,
    out.dir = ".",
    format = c("10x_dir", "tenx_h5"),
    prefix = "RNAs",
    progress = NULL,
    verbose = TRUE
) {
  format <- match.arg(format)
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

  ## Keep sparse: never densify a data.frame / Matrix via as.matrix().
  rna <- .ic2_to_dgCMatrix(rna.data)
  if (is.null(colnames(rna)) || is.null(rownames(rna))) {
    stop("'rna.data' must have rownames (genes) and colnames (cells).")
  }

  assign_col <- if ("assignment.annotation" %in% names(anno)) {
    "assignment.annotation"
  } else if ("hto_assignment" %in% names(anno)) {
    "hto_assignment"
  } else {
    stop("Annotation needs assignment.annotation (or hto_assignment).")
  }
  cell_ids <- rownames(anno)
  if (is.null(cell_ids) || !length(cell_ids)) {
    if ("ID" %in% names(anno)) cell_ids <- as.character(anno$ID)
    else stop("Annotation needs rownames or an ID column.")
  }
  tags <- unique(as.character(anno[[assign_col]]))
  tags <- tags[!is.na(tags) & nzchar(tags)]
  if (length(tags) < 1L) stop("No HTO assignments to demux.")

  out_paths <- character(0)
  cn <- colnames(rna)
  for (i in seq_along(tags)) {
    tag <- tags[[i]]
    if (is.function(progress)) progress(i - 1L, length(tags))
    sample_cells <- cell_ids[as.character(anno[[assign_col]]) == tag]
    sample_cells <- intersect(as.character(sample_cells), cn)
    if (isTRUE(verbose)) {
      message(" getting sample ", tag, " ...")
    }
    if (length(sample_cells) < 1L) {
      if (isTRUE(verbose)) message("  (no overlapping cells; skip)")
      next
    }
    sample <- rna[, sample_cells, drop = FALSE]
    if (isTRUE(verbose)) {
      message(" number of cells ", ncol(sample))
      message(" writing sample ", tag, " ...")
    }
    mat <- .ic2_to_dgCMatrix(sample)
    tag_clean <- .ic2_clean_names(tag)[1L]
    name <- paste(prefix, tag_clean, ncol(mat), sep = "_")
    dest <- file.path(out.dir, name)
    if (identical(format, "10x_dir")) {
      .ic2_write_10x_mtx_dir(mat, dest)
      out_paths <- c(out_paths, dest)
    } else {
      dest_h5 <- paste0(dest, ".h5")
      .ic2_write_matrix_tenx_h5(
        mat,
        filepath = dest_h5,
        sample_id = tag_clean,
        modality = "rna",
        feature_types = "Gene Expression",
        overwrite = TRUE
      )
      out_paths <- c(out_paths, dest_h5)
    }
  }
  if (is.function(progress)) progress(length(tags), length(tags))
  out_paths
}


#' Write a classic 10x MTX directory (matrix.mtx + barcodes + genes)
#' @keywords internal
.ic2_write_10x_mtx_dir <- function(mat, dir_path) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required to write 10x_dir splits.")
  }
  mat <- .ic2_to_dgCMatrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Matrix must have rownames and colnames for 10x_dir write.")
  }
  if (dir.exists(dir_path)) unlink(dir_path, recursive = TRUE)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  Matrix::writeMM(mat, file.path(dir_path, "matrix.mtx"))
  utils::write.table(
    colnames(mat),
    file = gzfile(file.path(dir_path, "barcodes.tsv.gz")),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  ## genes.tsv.gz: id + name (classic two-column)
  genes <- cbind(rownames(mat), rownames(mat))
  utils::write.table(
    genes,
    file = gzfile(file.path(dir_path, "genes.tsv.gz")),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  invisible(dir_path)
}


#' Fast TSV write: ID + annotation (classic cbind + write.table)
#' @keywords internal
.ic2_hto_write_anno_tsv <- function(anno, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tab <- data.frame(
    ID = rownames(anno),
    anno,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(tab, file = path, sep = "\t", quote = FALSE)
  } else {
    utils::write.table(tab, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  invisible(path)
}


#' Classic-style HTO QC plot (optimized for large n)
#'
#' Same two panels as the iCellR tutorial. For large cell counts, jitter uses
#' a random subsample (\code{max.points}) while boxplots use a larger sample
#' so the PDF/PNG stays fast.
#'
#' @keywords internal
.ic2_hto_anno_plot <- function(
    anno,
    path,
    format = c("png", "pdf"),
    max.points = 25000L
) {
  format <- match.arg(format)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping HTO plot.")
    return(NULL)
  }
  n <- nrow(anno)
  max.points <- as.integer(max.points)
  if (!is.finite(max.points) || max.points < 500L) max.points <- 25000L

  df <- data.frame(
    assignment.annotation = as.character(anno$assignment.annotation),
    percent.match = as.numeric(anno$percent.match),
    low.cov = factor(as.logical(anno$low.cov), levels = c(FALSE, TRUE)),
    stringsAsFactors = FALSE
  )

  ## Boxplot sample (up to 100k) + jitter subsample
  set.seed(1L)
  box_n <- min(n, max(max.points, 100000L))
  jit_n <- min(n, max.points)
  idx_box <- if (box_n < n) sample.int(n, box_n) else seq_len(n)
  idx_jit <- if (jit_n < n) sample.int(n, jit_n) else seq_len(n)
  df_box <- df[idx_box, , drop = FALSE]
  df_jit <- df[idx_jit, , drop = FALSE]

  subtitle <- if (jit_n < n) {
    paste0("jitter n=", format(jit_n, big.mark = ","), " / ", format(n, big.mark = ","))
  } else {
    NULL
  }

  A <- ggplot2::ggplot(df_box, ggplot2::aes(assignment.annotation, percent.match)) +
    ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA, coef = 1.5) +
    ggplot2::geom_jitter(
      data = df_jit,
      alpha = 0.2, color = "blue", width = 0.2, height = 0, size = 0.6
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(x = "assignment.annotation", y = "percent.match", subtitle = subtitle)

  B <- ggplot2::ggplot(df_box, ggplot2::aes(low.cov, percent.match)) +
    ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA, coef = 1.5) +
    ggplot2::geom_jitter(
      data = df_jit,
      alpha = 0.2, color = "blue", width = 0.2, height = 0, size = 0.6
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(x = "low.cov", y = "percent.match", subtitle = subtitle)

  dir.create(dirname(normalizePath(path, winslash = "/", mustWork = FALSE)),
             recursive = TRUE, showWarnings = FALSE)
  ## dirname of relative path may be "." 
  dn <- dirname(path)
  if (!identical(dn, ".") && !identical(dn, "")) {
    dir.create(dn, recursive = TRUE, showWarnings = FALSE)
  }

  if (identical(format, "pdf")) {
    grDevices::pdf(path, width = 8, height = 8, useDingbats = FALSE)
  } else {
    grDevices::png(path, width = 8, height = 8, units = "in", res = 300)
  }
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(A, B, ncol = 2)
  } else {
    print(A)
  }
  grDevices::dev.off()
  on.exit(NULL)
  path
}
