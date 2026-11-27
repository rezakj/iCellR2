#' Filter cells in an iCellR2 project
#'
#' Filters cells using QC columns in \code{x@cell_data} and optional
#' gene-expression / cell-id rules. By default this is **metadata-only**:
#' SODA keeps the raw counts assay unchanged; good cell ids are retained in
#' \code{x@cell_data}. Downstream \code{\link{norm.data}} writes the
#' normalized assay with filtered cells only; \code{\link{data.scale}} and
#' \code{\link{run.pca}} then use that filtered normalized / scaled data.
#'
#' Set \code{update.assay = TRUE} only if you explicitly want a physical
#' filtered counts copy on disk (slow for atlas-scale data).
#'
#' Filtering criteria (combined with AND across categories):
#' \itemize{
#'   \item mitochondrial fraction (\code{mito.percent})
#'   \item genes per cell (\code{nGenes})
#'   \item UMIs per cell (\code{UMIs})
#'   \item explicit cell ids to drop or keep
#'   \item optional gene expression gate (OR across genes): keep cells with
#'         expression \eqn{\ge} \code{filter.by.gene.exp.min} for at least
#'         one listed gene (evaluated in chunks on disk)
#' }
#'
#' Run \code{\link{qc.stats}} before mito / nGenes / UMI filters.
#'
#' @param x An object of class \code{iCellR2}.
#' @param modality Modality to filter. Default \code{"rna"}.
#' @param assay.name Source assay used only for optional gene-expression
#'   gates (and for optional physical rewrite). Default: active assay.
#' @param min.mito,max.mito Allowed mitochondrial fraction range.
#' @param min.genes,max.genes Allowed genes-per-cell range.
#' @param min.umis,max.umis Allowed UMI-per-cell range.
#' @param filter.by.cell.id Optional character vector of cell ids to remove.
#' @param keep.cell.id Optional character vector of cell ids to retain
#'   (mutually exclusive with \code{filter.by.cell.id}).
#' @param filter.by.gene Optional character vector of feature ids used for
#'   an expression gate (OR logic across genes).
#' @param filter.by.gene.exp.min Minimum expression for the gene gate.
#'   Default \code{1}.
#' @param update.assay If \code{FALSE} (default), only subset
#'   \code{x@cell_data} and keep the raw assay as-is. If \code{TRUE}, also
#'   write a filtered on-disk counts assay (slow at scale).
#' @param new.assay.name Name for the filtered assay when
#'   \code{update.assay = TRUE}. Default \code{paste0(assay.name, "_filtered")}.
#' @param set.active If \code{TRUE} and \code{update.assay = TRUE}, set
#'   \code{new.assay.name} as active.
#' @param chunk.size Cells per disk chunk for gene-gate / optional rewrite.
#' @param overwrite Overwrite existing filtered assay files.
#' @param verbose Progress style: \code{TRUE}/\code{"full"}, \code{"compact"},
#'   or \code{FALSE}/\code{"silent"}.
#'
#' @return An updated \code{iCellR2} object.
#'
#' @examples
#' \dontrun{
#' obj <- qc.stats(obj)
#' obj <- cell.filter(
#'   obj,
#'   max.mito = 0.1,
#'   min.genes = 200,
#'   max.genes = 6000,
#'   min.umis = 500
#' )
#' obj <- norm.data(obj)   # writes "data" with filtered cells only
#' }
#'
#' @export
cell.filter <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    min.mito = 0,
    max.mito = 1,
    min.genes = 0,
    max.genes = Inf,
    min.umis = 0,
    max.umis = Inf,
    filter.by.cell.id = NULL,
    keep.cell.id = NULL,
    filter.by.gene = NULL,
    filter.by.gene.exp.min = 1,
    update.assay = FALSE,
    new.assay.name = NULL,
    set.active = TRUE,
    chunk.size = 200000L,
    overwrite = TRUE,
    verbose = TRUE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  assay.name <- .ic2_resolve_assay_name(x, modality, assay.name)
  chunk.size <- as.integer(chunk.size)
  rep <- .ic2_new_reporter("cell.filter", verbose, soda_root = x@storage$root_dir)

  if (!is.null(filter.by.cell.id) && !is.null(keep.cell.id)) {
    stop("Set only one of 'filter.by.cell.id' or 'keep.cell.id'.")
  }

  cd <- x@cell_data
  if (is.null(cd) || !is.data.frame(cd) || nrow(cd) < 1L) {
    stop("'x@cell_data' is empty. Run make.obj() / qc.stats() first.")
  }
  cell_ids <- rownames(cd)
  n0 <- length(cell_ids)
  n_feat <- if (!is.null(x@misc$n_features)) x@misc$n_features else NA_integer_
  keep <- rep(TRUE, n0)
  names(keep) <- cell_ids

  if (is.null(new.assay.name)) {
    new.assay.name <- paste0(assay.name, "_filtered")
  }

  rep <- .ic2_banner_start(rep, c(
    paste0(
      "Modality: ", modality, " | Assay: ", assay.name,
      if (isTRUE(update.assay)) paste0("  ->  ", new.assay.name) else " (IDs only; raw unchanged)"
    ),
    paste0(
      "Cells in: ", .ic2_fmt_num(n0),
      if (is.finite(n_feat)) paste0(" | Features: ", .ic2_fmt_num(n_feat)) else ""
    )
  ))

  need_qc <- any(c(min.mito > 0, max.mito < 1, min.genes > 0, is.finite(max.genes),
                   min.umis > 0, is.finite(max.umis)))
  qc_cols <- c("mito.percent", "nGenes", "UMIs")
  if (need_qc && !all(qc_cols %in% colnames(cd))) {
    stop(
      "QC columns missing in 'x@cell_data' (",
      paste(setdiff(qc_cols, colnames(cd)), collapse = ", "),
      "). Run qc.stats() first."
    )
  }

  if (need_qc) {
    keep <- keep &
      cd$mito.percent >= min.mito & cd$mito.percent <= max.mito &
      cd$nGenes >= min.genes & cd$nGenes <= max.genes &
      cd$UMIs >= min.umis & cd$UMIs <= max.umis
    .ic2_step(rep, paste0(
      "[QC gates] mito [", min.mito, ", ", max.mito, "]  genes [",
      min.genes, ", ", max.genes, "]  UMIs [", min.umis, ", ", max.umis, "]"
    ))
  }

  if (!is.null(filter.by.cell.id)) {
    filter.by.cell.id <- unique(as.character(filter.by.cell.id))
    missing <- setdiff(filter.by.cell.id, cell_ids)
    if (length(missing) > 0L) {
      stop("filter.by.cell.id not found: ", paste(utils::head(missing, 10), collapse = ", "))
    }
    keep[filter.by.cell.id] <- FALSE
    .ic2_step(rep, paste0("Removed ", length(filter.by.cell.id), " cell id(s) via filter.by.cell.id."))
  }

  if (!is.null(keep.cell.id)) {
    keep.cell.id <- unique(as.character(keep.cell.id))
    missing <- setdiff(keep.cell.id, cell_ids)
    if (length(missing) > 0L) {
      stop("keep.cell.id not found: ", paste(utils::head(missing, 10), collapse = ", "))
    }
    keep[] <- FALSE
    keep[keep.cell.id] <- TRUE
    .ic2_step(rep, paste0("Retained ", length(keep.cell.id), " cell id(s) via keep.cell.id."))
  }

  if (!is.null(filter.by.gene) && length(filter.by.gene) > 0L) {
    .ic2_step(rep, paste0("[Gene gate] streaming expression for: ",
                          paste(filter.by.gene, collapse = ", ")))
    handles <- .ic2_get_assay_handles(x, modality, assay.name)
    rep <- .ic2_register_assay_inputs(rep, x, modality, assay.name, label = paste0("Assay[", assay.name, "]"))
    gene_fail <- character(0)
    for (h in handles) {
      gene_fail <- c(
        gene_fail,
        .ic2_gene_filter_fail_ids(
          h,
          genes = filter.by.gene,
          min_exp = filter.by.gene.exp.min,
          chunk.size = chunk.size
        )
      )
    }
    gene_fail <- unique(gene_fail)
    gene_fail <- intersect(gene_fail, names(keep)[keep])
    keep[gene_fail] <- FALSE
    .ic2_step(rep, paste0(
      "  Gene gate (>= ", filter.by.gene.exp.min, "): removed ",
      .ic2_fmt_num(length(gene_fail)), " cells."
    ))
  }

  keep_ids <- names(keep)[keep]
  n_keep <- length(keep_ids)
  n_drop <- n0 - n_keep
  if (n_keep < 1L) {
    stop("All cells were filtered out. Relax thresholds and try again.")
  }
  .ic2_step(rep, paste0(
    "  Keep: ", .ic2_fmt_num(n_keep), " / ", .ic2_fmt_num(n0),
    "  (", .ic2_fmt_num(n_drop), " removed)"
  ))
  .ic2_msg(rep, "")

  x@cell_data <- cd[keep_ids, , drop = FALSE]
  x@misc$n_cells <- n_keep
  .ic2_persist_cell_data(x)

  ## Persist kept barcodes for downstream / audit (raw assay stays intact)
  if (!is.null(x@storage$root_dir)) {
    filt_dir <- file.path(x@storage$root_dir, "modalities", modality, "filters")
    dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
    keep_tsv <- file.path(filt_dir, "kept_cell_ids.tsv")
    utils::write.table(
      data.frame(cell_id = as.character(keep_ids), stringsAsFactors = FALSE),
      file = keep_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )
    .ic2_soda_unlink_rds(file.path(filt_dir, "kept_cell_ids.rds"))
    x@storage$metadata_paths$kept_cell_ids <- keep_tsv
    .ic2_step(rep, paste0("Saved kept cell ids  ->  ", keep_tsv))
  }

  n_runs_total <- 0L
  dest_paths <- character(0)

  if (isTRUE(update.assay)) {
    handles <- .ic2_get_assay_handles(x, modality, assay.name)
    rep <- .ic2_register_assay_inputs(rep, x, modality, assay.name, label = paste0("Assay[", assay.name, "]"))
    is_virtual <- .ic2_assay_is_virtual(x, modality, assay.name)
    new_handles <- vector("list", length(handles))
    .ic2_step(rep, paste0(
      "[optional] Subsetting on-disk counts assay",
      if (length(handles) > 1L) paste0("  -  ", length(handles), " samples") else ""
    ))

    for (i in seq_along(handles)) {
      h <- .ic2_enrich_handle(handles[[i]], load_dimnames = TRUE)
      local_keep <- intersect(keep_ids, h$cell_ids)
      if (length(local_keep) < 1L) {
        stop("Sample/assay handle '", h$path, "' has no cells left after filtering.")
      }
      sid <- if (!is.null(h$sample_id) && nzchar(as.character(h$sample_id))) {
        h$sample_id
      } else if (is_virtual) {
        paste0("sample", i)
      } else {
        NULL
      }
      dest <- .ic2_soda_assay_path(x, modality, new.assay.name, sample_id = sid, ext = ".h5")
      if (file.exists(dest) && !overwrite) {
        stop("Filtered assay exists at '", dest, "'. Use overwrite = TRUE.")
      }
      prefix <- if (length(handles) > 1L) paste0("  sample ", i, "/", length(handles)) else " "
      new_handles[[i]] <- .ic2_subset_handle_cells(
        handle = h,
        keep_cell_ids = local_keep,
        dest_path = dest,
        chunk.size = chunk.size,
        progress_cb = function(cur, tot) {
          rep <<- .ic2_progress_bar(rep, cur, tot, prefix = prefix)
        }
      )
      n_runs_total <- n_runs_total + as.integer(attr(new_handles[[i]], "n_runs") %||% 0L)
      dest_paths <- c(dest_paths, dest)
      new_handles[[i]]$sample_id <- sid
      new_handles[[i]]$modality <- modality
      new_handles[[i]]$feature_ids <- NULL
      new_handles[[i]]$cell_ids <- NULL
    }

    .ic2_step(rep, "Registering filtered assay handles")
    feature_universe <- NULL
    if (is_virtual) {
      old <- x@modalities[[modality]]$assays[[assay.name]]
      feature_universe <- old$feature_universe
    }

    x <- .ic2_assign_assay(
      x,
      modality = modality,
      assay.name = new.assay.name,
      handles = new_handles,
      virtual = is_virtual,
      feature_universe = feature_universe,
      set.active = set.active
    )
  } else {
    .ic2_step(rep, "Raw assay unchanged  -  norm.data() will write filtered cells only.")
  }

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  if (!is.null(x@storage$metadata_paths$kept_cell_ids)) {
    rep <- .ic2_add_output(rep, "Kept IDs", x@storage$metadata_paths$kept_cell_ids)
  }
  for (p in dest_paths) rep <- .ic2_add_output(rep, "Assay", p)
  if (!is.null(x@storage$metadata_paths$cell_data)) {
    rep <- .ic2_add_output(rep, "Metadata", x@storage$metadata_paths$cell_data)
  }

  rep <- .ic2_add_scale(rep, c(
    if (isTRUE(update.assay)) {
      paste0("Mode:     physical counts rewrite (", n_runs_total, " HDF5 run(s))")
    } else {
      "Mode:     cell-ID filter only (raw counts left intact in SODA)"
    },
    paste0("Kept:     ", .ic2_fmt_num(n_keep), " / ", .ic2_fmt_num(n0), " barcodes in cell_data"),
    "Next:     norm.data()  ->  filtered 'data'; data.scale / PCA use that"
  ))

  wall <- as.numeric(difftime(Sys.time(), rep$wall0, units = "secs"))
  if (is.finite(wall) && wall > 0) {
    rep$stats$cells_per_sec <- n_keep / wall
  }

  filter_summary <- list(
    n_before = n0,
    n_after = n_keep,
    n_removed = n_drop,
    min.mito = min.mito,
    max.mito = max.mito,
    min.genes = min.genes,
    max.genes = max.genes,
    min.umis = min.umis,
    max.umis = max.umis,
    filter.by.gene = filter.by.gene,
    filter.by.gene.exp.min = filter.by.gene.exp.min,
    update.assay = update.assay,
    new.assay.name = if (isTRUE(update.assay)) new.assay.name else NA_character_
  )
  x@params$qc[[paste0(modality, "_filter")]] <- filter_summary
  x@misc$my.filters <- filter_summary

  x <- .ic2_log_command(
    x,
    "cell.filter",
    arguments = list(
      modality = modality,
      assay.name = assay.name,
      new.assay.name = new.assay.name,
      update.assay = update.assay,
      chunk.size = chunk.size
    ),
    outputs = list(n_before = n0, n_after = n_keep, n_removed = n_drop)
  )

  .ic2_finish(rep)
  x
}
