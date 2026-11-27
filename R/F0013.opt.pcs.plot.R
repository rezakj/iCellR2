#' Plot PCA standard deviations and suggest optimal PCs
#'
#' Reads PCA stdev from the SODA archive / object (written by
#' \code{\link{run.pca}}) and plots an elbow. Suggested \code{opt.pcs}
#' uses a knee/elbow detector on PC variance (typically near the visual
#' bend, often ~10 for scRNA-seq).
#'
#' @param x An \code{iCellR2} object (after \code{run.pca}).
#' @param modality Modality. Default \code{"rna"}.
#' @param pcs.in.plot Number of PCs to show. Default \code{50}.
#' @param method Heuristic passed to the internal suggester:
#'   \code{"knee"} (default), \code{"legacy"}, or \code{"variance"}.
#' @param verbose Progress style (usually quiet for a plot helper).
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline
#'   scale_x_continuous scale_y_continuous scale_colour_gradient
#'   theme_bw xlab ylab ggtitle
#' @export
opt.pcs.plot <- function(
    x = NULL,
    modality = "rna",
    pcs.in.plot = 50,
    method = c("knee", "legacy", "variance"),
    verbose = FALSE
) {
  .ic2_require_object(x)
  modality <- tolower(modality)
  method <- match.arg(method)
  rep <- .ic2_new_reporter("opt.pcs.plot", verbose, soda_root = x@storage$root_dir)

  pca <- .ic2_get_pca_ref(x, modality)
  sdev <- .ic2_get_pca_stdev(pca)
  if (length(sdev) < 2L) {
    stop("PCA stdev not found or too short. Run run.pca() first.")
  }

  OPTpcs <- .ic2_suggest_opt_pcs(sdev, method = method)
  stored <- pca$opt.pcs %||% x@misc$opt.pcs
  titel <- paste0("Optimal number of PCs (1:", OPTpcs, "; ", method, ")")

  pcs.in.plot <- as.integer(min(as.integer(pcs.in.plot), length(sdev)))
  DATA <- data.frame(
    PCs = seq_len(pcs.in.plot),
    SDs = as.numeric(sdev[seq_len(pcs.in.plot)]),
    stringsAsFactors = FALSE
  )
  DATA$LogSDs <- log2(pmax(DATA$SDs, .Machine$double.eps))

  rep <- .ic2_banner_start(rep, c(
    paste0("Modality: ", modality),
    paste0(
      "PCs available: ", length(sdev),
      "   |   suggested opt.pcs: ", OPTpcs,
      if (!is.null(stored)) paste0("   |   stored: ", stored) else ""
    )
  ))
  .ic2_finish(rep, extra_lines = paste0("Plotting first ", pcs.in.plot, " PCs"))

  p <- ggplot2::ggplot(DATA, ggplot2::aes(y = SDs, x = PCs, col = LogSDs)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = pmax(1, DATA$LogSDs)) +
    ggplot2::scale_x_continuous(trans = "log1p") +
    ggplot2::scale_y_continuous(trans = "log1p") +
    ggplot2::geom_vline(xintercept = OPTpcs, linetype = "dotted") +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::xlab("Principal Components") +
    ggplot2::ylab("PC Standard Deviations") +
    ggplot2::scale_colour_gradient(low = "gray", high = "red", name = "log2 (SDs)") +
    ggplot2::ggtitle(titel)

  if (!is.null(stored) && is.finite(stored) && as.integer(stored) != OPTpcs) {
    p <- p + ggplot2::geom_vline(xintercept = as.integer(stored), linetype = "dashed", colour = "steelblue")
  }
  p
}

#' Resolve PCA reduction reference from an iCellR2 object
#' @keywords internal
.ic2_get_pca_ref <- function(x, modality = "rna") {
  modality <- tolower(modality)
  pca <- x@modalities[[modality]]$reductions$pca
  if (is.null(pca)) {
    stop("No PCA found for modality '", modality, "'. Run run.pca() first.")
  }
  pca
}

#' Load PCA stdev vector
#' @keywords internal
.ic2_get_pca_stdev <- function(pca) {
  if (!is.null(pca$stdev) && length(pca$stdev) > 0L) {
    return(as.numeric(pca$stdev))
  }
  if (!is.null(pca$paths$stdev) && file.exists(pca$paths$stdev)) {
    return(.ic2_read_stdev_file(pca$paths$stdev))
  }
  ## path may be stored on the ref as stdev file beside embeddings
  if (is.character(pca$cell_embeddings) && length(pca$cell_embeddings) == 1L) {
    for (cand in c(
      file.path(dirname(pca$cell_embeddings), "stdev.tsv"),
      file.path(dirname(pca$cell_embeddings), "stdev.rds")
    )) {
      if (file.exists(cand)) return(.ic2_read_stdev_file(cand))
    }
  }
  if (!is.null(pca$meta_path) && file.exists(pca$meta_path)) {
    meta <- .ic2_soda_read_meta(pca$meta_path)
    if (!is.null(meta$paths$stdev) && file.exists(meta$paths$stdev)) {
      return(.ic2_read_stdev_file(meta$paths$stdev))
    }
  }
  numeric(0)
}

.ic2_read_stdev_file <- function(path) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    df <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE)
    if ("stdev" %in% names(df)) return(as.numeric(df$stdev))
    return(as.numeric(df[[ncol(df)]]))
  }
  as.numeric(readRDS(path))
}

#' Load PCA cell embeddings matrix (cells  x  PCs)
#' @keywords internal
.ic2_get_pca_embeddings <- function(pca) {
  path <- pca$cell_embeddings
  if (is.null(path) && !is.null(pca$paths$embeddings)) {
    path <- pca$paths$embeddings
  }
  if (is.null(path) && !is.null(pca$meta_path) && file.exists(pca$meta_path)) {
    meta <- .ic2_soda_read_meta(pca$meta_path)
    path <- meta$paths$embeddings
  }
  if (is.null(path) || !file.exists(path)) {
    stop("PCA cell embeddings not found on disk. Run run.pca() first.")
  }
  emb <- if (grepl("\\.(h5|hdf5)$", path, ignore.case = TRUE)) {
    .ic2_h5_read_dense(path)
  } else {
    readRDS(path)
  }
  if (!is.matrix(emb)) emb <- as.matrix(emb)
  emb
}
