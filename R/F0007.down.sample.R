#' Down-sample cells to balance groups
#'
#' Randomly down-samples cells so each level of a metadata column
#' (default \code{condition_id}) has the same number of cells. Uses
#' \code{x@cell_data} only for sampling decisions; optionally rewrites the
#' active assay on disk for retained cells (handle-only / HDF5 path).
#'
#' @param x An object of class \code{iCellR2}.
#' @param modality Modality. Default \code{"rna"}.
#' @param assay.name Source assay; default active assay.
#' @param group.by Metadata column in \code{x@cell_data} used to balance
#'   groups. Default \code{"condition_id"}.
#' @param n.cells Target cells per group. If \code{NULL} (default), uses
#'   the size of the smallest group.
#' @param seed Integer random seed. Default \code{1}.
#' @param update.assay If \code{FALSE} (default), only subset
#'   \code{x@cell_data}. Normalized / scaled assays written later will use
#'   the retained barcodes. Set \code{TRUE} only for a physical assay rewrite.
#' @param new.assay.name Name for the down-sampled assay. Default
#'   \code{paste0(assay.name, "_downsampled")}.
#' @param set.active If \code{TRUE}, set the new assay as active.
#' @param chunk.size Cells per disk chunk when rewriting the assay.
#' @param overwrite Overwrite existing down-sampled assay files.
#' @param verbose Progress style: \code{TRUE}/\code{"full"}, \code{"compact"},
#'   or \code{FALSE}/\code{"silent"}.
#'
#' @return An updated \code{iCellR2} object.
#'
#' @examples
#' \dontrun{
#' obj <- down.sample(obj, group.by = "condition_id")
#' obj <- down.sample(obj, group.by = "sample_id", n.cells = 1000, seed = 42)
#' }
#'
#' @export
down.sample <- function(
    x,
    modality = "rna",
    assay.name = NULL,
    group.by = "condition_id",
    n.cells = NULL,
    seed = 1,
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
  rep <- .ic2_new_reporter("down.sample", verbose, soda_root = x@storage$root_dir)

  cd <- x@cell_data
  if (is.null(cd) || nrow(cd) < 1L) {
    stop("'x@cell_data' is empty.")
  }
  if (!group.by %in% colnames(cd)) {
    stop("'group.by' column '", group.by, "' not found in x@cell_data.")
  }

  groups <- as.character(cd[[group.by]])
  if (any(is.na(groups) | !nzchar(groups))) {
    stop("'group.by' contains missing or empty values.")
  }
  tab <- table(groups)
  if (length(tab) < 2L) {
    stop("Need at least two levels in '", group.by, "' to down-sample.")
  }

  target <- if (is.null(n.cells)) {
    as.integer(min(tab))
  } else {
    as.integer(n.cells)
  }
  if (is.na(target) || target < 1L) {
    stop("'n.cells' must be a positive integer.")
  }
  if (target > min(tab)) {
    stop(
      "'n.cells' = ", target, " exceeds smallest group size (", min(tab),
      "). Reduce n.cells or omit it."
    )
  }

  if (is.null(new.assay.name)) {
    new.assay.name <- paste0(assay.name, "_downsampled")
  }

  before <- paste0(names(tab), " (", .ic2_fmt_num(as.integer(tab)), ")", collapse = ", ")
  rep <- .ic2_banner_start(rep, c(
    paste0("Modality: ", modality, " | Assay: ", assay.name,
           if (isTRUE(update.assay)) paste0("  ->  ", new.assay.name) else ""),
    paste0("Balance by: ", group.by, " | target/group: ", .ic2_fmt_num(target)),
    paste0("Before: ", before)
  ))

  set.seed(seed)
  keep_ids <- unlist(lapply(names(tab), function(g) {
    ids <- rownames(cd)[groups == g]
    if (length(ids) == target) return(ids)
    ids[sample.int(length(ids), size = target, replace = FALSE)]
  }), use.names = FALSE)
  keep_ids <- keep_ids[order(match(keep_ids, rownames(cd)))]

  x@cell_data <- cd[keep_ids, , drop = FALSE]
  x@misc$n_cells <- length(keep_ids)
  .ic2_persist_cell_data(x)

  after_tab <- table(as.character(x@cell_data[[group.by]]))
  after <- paste0(names(after_tab), " (", .ic2_fmt_num(as.integer(after_tab)), ")", collapse = ", ")
  .ic2_step(rep, paste0("After:  ", after))
  .ic2_msg(rep, "")

  n_runs_total <- 0L
  dest_paths <- character(0)

  if (isTRUE(update.assay)) {
    handles <- .ic2_get_assay_handles(x, modality, assay.name)
    is_virtual <- .ic2_assay_is_virtual(x, modality, assay.name)
    new_handles <- vector("list", length(handles))
    .ic2_step(rep, "[1/2] Writing down-sampled on-disk assay (contiguous HDF5 runs)")

    for (i in seq_along(handles)) {
      h <- .ic2_enrich_handle(handles[[i]], load_dimnames = TRUE)
      local_keep <- intersect(keep_ids, h$cell_ids)
      if (length(local_keep) < 1L) {
        stop("Handle '", h$path, "' has no retained cells after down-sampling.")
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
        stop("Down-sampled assay exists at '", dest, "'. Use overwrite = TRUE.")
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

    .ic2_step(rep, "[2/2] Registering down-sampled assay handles")
    feature_universe <- NULL
    if (is_virtual) {
      feature_universe <- x@modalities[[modality]]$assays[[assay.name]]$feature_universe
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
    .ic2_step(rep, "update.assay = FALSE: cell_data updated; on-disk assays unchanged.")
  }

  if (!is.null(x@storage$root_dir)) {
    rep <- .ic2_add_output(rep, "Archive", x@storage$root_dir)
  }
  for (p in dest_paths) rep <- .ic2_add_output(rep, "Assay", p)
  if (!is.null(x@storage$metadata_paths$cell_data)) {
    rep <- .ic2_add_output(rep, "Metadata", x@storage$metadata_paths$cell_data)
  }

  rep <- .ic2_add_scale(rep, c(
    "Mode:     metadata sampling first, then sparse on-disk rewrite",
    paste0("I/O:      ", n_runs_total, " contiguous HDF5 run(s)"),
    "Memory:   no dense merge  -  S4 stays handle-only",
    paste0("Balanced: ", length(tab), " groups  x  ", .ic2_fmt_num(target), " cells")
  ))

  wall <- as.numeric(difftime(Sys.time(), rep$wall0, units = "secs"))
  if (is.finite(wall) && wall > 0) {
    rep$stats$cells_per_sec <- length(keep_ids) / wall
  }

  x@params$misc$down.sample <- list(
    group.by = group.by,
    n.cells = target,
    seed = seed,
    n_before = nrow(cd),
    n_after = length(keep_ids),
    counts_before = as.list(tab),
    counts_after = as.list(after_tab),
    new.assay.name = if (isTRUE(update.assay)) new.assay.name else NA_character_
  )

  x <- .ic2_log_command(
    x,
    "down.sample",
    arguments = list(
      modality = modality,
      assay.name = assay.name,
      group.by = group.by,
      n.cells = target,
      seed = seed,
      update.assay = update.assay
    ),
    outputs = list(n_before = nrow(cd), n_after = length(keep_ids))
  )

  .ic2_finish(rep)
  x
}
