#' Score cells by gene sets (streamed)
#'
#' Classic \code{i.score}: each named gene set is scored per cell, then the
#' winning set is assigned with the same rule as \code{\link{hto.anno}}
#' (dominant score as a percent of the total).
#'
#' Only the signature genes (and, for \code{tirosh}, expression-matched
#' control genes) are read from the on-disk assay, in cell chunks. The full
#' matrix is never densified. \code{gsva}, \code{ssgsea}, and \code{plage}
#' need a dense matrix and are not used here. Use \code{tirosh} (default),
#' \code{mean}, \code{sum}, or \code{zscore}.
#'
#' Gene sets are a named list. A character vector of names still works when
#' those objects exist in the calling environment (classic style).
#'
#' Scores and the assignment are written to \code{cell_data} and
#' \code{metadata/scoring/}.
#'
#' @param object An \code{iCellR2} object. Alias: \code{x}.
#' @param x Alias for \code{object}.
#' @param data.type Classic assay alias: \code{"main.data"} (normalized
#'   \code{data}), \code{"raw.data"} (\code{counts}), \code{"imputed.data"}.
#' @param assay.name Assay to read. Overrides \code{data.type}.
#' @param scoring.List Named list of gene ids, or a character vector of
#'   object names that hold those vectors. At least two sets.
#' @param scoring.method \code{"tirosh"} (default), \code{"mean"},
#'   \code{"sum"}, or \code{"zscore"}.
#' @param return.stats If \code{TRUE}, return the annotation data.frame.
#'   If \code{FALSE} (default), return the updated object. The object and
#'   TSV are updated either way.
#' @param ctrl.size Control genes per signature gene for \code{tirosh}.
#'   Default: size of the smallest set.
#' @param n.bin Expression bins for \code{tirosh} controls. Default \code{25}.
#' @param cov.thr,assignment.thr Passed to the HTO-style assignment.
#'   Defaults \code{10} and \code{50} (classic \code{i.score}).
#' @param chunk.size Cells per HDF5 read. Default \code{200000}.
#' @param modality Modality. Default \code{"rna"}.
#' @param seed RNG seed for control-gene sampling.
#' @param verbose Progress style.
#'
#' @return The updated object, or the annotation data.frame when
#'   \code{return.stats = TRUE}.
#'
#' @examples
#' \dontrun{
#' sets <- list(G1S = g1s_genes, G2M = g2m_genes, S = s_genes)
#' obj <- i.score(obj, scoring.List = sets)
#' tab <- i.score(obj, scoring.List = sets, return.stats = TRUE)
#' }
#'
#' @export
i.score <- function(
    object = NULL,
    data.type = "main.data",
    scoring.List = NULL,
    return.stats = FALSE,
    scoring.method = c("tirosh", "mean", "sum", "zscore"),
    x = NULL,
    assay.name = NULL,
    ctrl.size = NULL,
    n.bin = 25L,
    cov.thr = 10,
    assignment.thr = 50,
    chunk.size = 200000L,
    modality = "rna",
    seed = 1,
    verbose = TRUE
) {
  if (is.null(object) && !is.null(x)) object <- x
  .ic2_require_object(object)
  modality <- tolower(modality)
  scoring.method <- match.arg(scoring.method)
  chunk.size <- as.integer(chunk.size)
  if (!is.finite(chunk.size) || chunk.size < 1L) chunk.size <- 200000L

  sets <- .ic2_score_as_list(scoring.List)
  if (length(sets) < 2L) {
    stop("scoring.List needs at least 2 gene sets (named list, or names of objects).")
  }

  if (is.null(assay.name) || !nzchar(as.character(assay.name)[1L])) {
    assay.name <- switch(
      data.type,
      main.data = "data",
      raw.data = "counts",
      imputed.data = "imputed",
      as.character(data.type)[1L]
    )
  }
  assay.name <- .ic2_resolve_assay_name(object, modality, assay.name)
  handles <- .ic2_get_assay_handles(object, modality, assay.name)
  feature_ids <- character(0)
  cell_ids <- character(0)
  for (h in handles) {
    h2 <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    feature_ids <- union(feature_ids, as.character(h2$feature_ids))
    cell_ids <- union(cell_ids, as.character(h2$cell_ids))
  }
  cd_ids <- rownames(object@cell_data)
  if (!is.null(cd_ids) && length(cd_ids)) {
    cell_ids <- intersect(cell_ids, cd_ids)
  }
  if (length(cell_ids) < 2L) stop("Need at least 2 cells to score.")

  sets <- lapply(sets, function(g) .ic2_score_match_genes(g, feature_ids))
  empty <- names(sets)[vapply(sets, length, 1L) < 1L]
  if (length(empty)) {
    stop("No genes matched for: ", paste(empty, collapse = ", "), ".")
  }

  rep <- .ic2_new_reporter("i.score", verbose, soda_root = object@storage$root_dir)
  rep <- .ic2_banner_start(rep, c(
    paste0("Method: ", scoring.method, "   |   assay: ", assay.name),
    paste0("Sets: ", paste(sprintf("%s (%d)", names(sets), lengths(sets)), collapse = ", "))
  ))

  ctrl <- NULL
  if (identical(scoring.method, "tirosh")) {
    .ic2_step(rep, "[1/3] Expression bins for control genes")
    means <- .ic2_score_gene_means(
      handles, feature_ids, chunk.size, object = object, modality = modality
    )
    ctrl <- .ic2_score_controls(sets, means, n.bin = as.integer(n.bin), ctrl.size = ctrl.size, seed = seed)
  } else {
    .ic2_step(rep, "[1/3] Gene sets matched")
  }

  .ic2_step(rep, "[2/3] Streaming scores")
  mu <- NULL
  sdv <- NULL
  if (identical(scoring.method, "zscore")) {
    mom <- .ic2_score_moments(handles, cell_ids, unique(unlist(sets)), chunk.size)
    mu <- mom$mu
    sdv <- mom$sd
  }
  scores <- .ic2_score_stream(
    handles, cell_ids, sets, ctrl = ctrl,
    method = scoring.method, mu = mu, sdv = sdv, chunk.size = chunk.size,
    progress = function(cur, tot) {
      rep <<- .ic2_progress_bar(rep, cur, tot, prefix = "  score")
    }
  )

  .ic2_step(rep, "[3/3] Assigning the winning set")
  anno <- .ic2_hto_anno_core(
    hto.data = t(scores),
    cov.thr = cov.thr,
    assignment.thr = assignment.thr
  )
  anno <- anno[cell_ids, , drop = FALSE]

  for (nm in colnames(scores)) {
    object@cell_data[[paste0("score_", nm)]] <- NA_real_
    object@cell_data[rownames(anno), paste0("score_", nm)] <- anno[[nm]]
  }
  object@cell_data$score_assignment <- NA_character_
  object@cell_data[rownames(anno), "score_assignment"] <- as.character(anno$assignment.annotation)
  object@cell_data$score_percent <- NA_real_
  object@cell_data[rownames(anno), "score_percent"] <- anno$percent.match
  object <- .ic2_persist_cell_data(object)

  score_dir <- file.path(object@storage$root_dir, "metadata", "scoring")
  dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
  out_tsv <- file.path(score_dir, "scores.tsv")
  utils::write.table(
    data.frame(cell_id = rownames(anno), anno, stringsAsFactors = FALSE, check.names = FALSE),
    file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE
  )

  object <- .ic2_log_command(
    object, "i.score",
    arguments = list(
      scoring.method = scoring.method, assay.name = assay.name,
      sets = names(sets), assignment.thr = assignment.thr
    ),
    outputs = list(scores = out_tsv, assignment = "cell_data$score_assignment")
  )
  .ic2_finish(rep, extra_lines = paste0(
    "Assignment  ->  cell_data$score_assignment   |   ", normalizePath(out_tsv, winslash = "/", mustWork = FALSE)
  ))
  if (isTRUE(return.stats)) anno else object
}


#' @keywords internal
.ic2_score_as_list <- function(scoring.List) {
  if (is.null(scoring.List)) {
    stop("Provide scoring.List: a named list of gene vectors, or names of those objects.")
  }
  if (is.list(scoring.List)) {
    if (is.null(names(scoring.List)) || any(!nzchar(names(scoring.List)))) {
      stop("scoring.List must be a named list (names are the score columns).")
    }
    return(lapply(scoring.List, function(g) unique(as.character(g))))
  }
  nms <- as.character(scoring.List)
  env <- parent.frame(2L)
  out <- lapply(nms, function(nm) {
    if (!exists(nm, envir = env, inherits = TRUE)) {
      stop("Gene set '", nm, "' was not found. Pass a named list instead of names.")
    }
    unique(as.character(get(nm, envir = env, inherits = TRUE)))
  })
  names(out) <- nms
  out
}


#' @keywords internal
.ic2_score_match_genes <- function(genes, feature_ids) {
  key <- tolower(feature_ids)
  hit <- match(tolower(genes), key)
  unique(feature_ids[hit[!is.na(hit)]])
}


#' Gene means for Tirosh bins. Uses gene.stats meanExp when it covers the assay.
#' @keywords internal
.ic2_score_gene_means <- function(handles, feature_ids, chunk.size, object = NULL, modality = "rna") {
  fd <- if (!is.null(object)) object@modalities[[modality]]$feature_data else NULL
  if (!is.null(fd) && "meanExp" %in% colnames(fd)) {
    rn <- rownames(fd)
    if (is.null(rn) && "gene" %in% colnames(fd)) rn <- fd$gene
    if (!is.null(rn)) {
      m <- match(feature_ids, rn)
      if (mean(!is.na(m)) > 0.9 && any(is.finite(fd$meanExp[m]))) {
        means <- rep(NA_real_, length(feature_ids))
        names(means) <- feature_ids
        ok <- !is.na(m)
        means[ok] <- as.numeric(fd$meanExp[m[ok]])
        means[!is.finite(means)] <- 0
        return(means)
      }
    }
  }
  sums <- numeric(length(feature_ids))
  names(sums) <- feature_ids
  n <- 0L
  for (h in handles) {
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    idx <- match(feature_ids, h$feature_ids)
    present <- which(!is.na(idx))
    n <- n + length(h$cell_ids)
    runs <- .ic2_contig_runs(seq_along(h$cell_ids))
    done <- 0L
    for (run in runs) {
      for (s in seq.int(1L, length(run), by = chunk.size)) {
        e <- min(s + chunk.size - 1L, length(run))
        raw <- .ic2_to_dgCMatrix(.ic2_read_cols(h, run[s:e]))
        part <- numeric(length(feature_ids))
        if (length(present)) {
          part[present] <- Matrix::rowSums(raw[idx[present], , drop = FALSE])
        }
        sums <- sums + part
        done <- done + (e - s + 1L)
      }
    }
  }
  if (n < 1L) stop("No cells while computing gene means.")
  sums / n
}


#' @keywords internal
.ic2_score_controls <- function(sets, means, n.bin = 25L, ctrl.size = NULL, seed = 1) {
  means <- means[is.finite(means)]
  means <- sort(means)
  n.bin <- max(2L, as.integer(n.bin))
  bins <- as.integer(cut(seq_along(means), breaks = n.bin, labels = FALSE))
  names(bins) <- names(means)
  if (is.null(ctrl.size)) ctrl.size <- min(lengths(sets))
  ctrl.size <- max(1L, as.integer(ctrl.size))
  set.seed(seed)
  lapply(sets, function(genes) {
    ctrl <- character(0)
    for (g in genes) {
      b <- bins[match(tolower(g), tolower(names(bins)))]
      if (is.na(b)) next
      pool <- names(bins)[bins == b]
      pool <- setdiff(pool, genes)
      if (!length(pool)) next
      take <- min(ctrl.size, length(pool))
      ctrl <- c(ctrl, sample(pool, take, replace = FALSE))
    }
    unique(ctrl)
  })
}


#' One streamed pass: cells  x  sets score matrix.
#' @keywords internal
.ic2_score_stream <- function(handles, cell_ids, sets, ctrl = NULL, method = "mean", mu = NULL, sdv = NULL, chunk.size = 200000L, progress = NULL) {
  set_names <- names(sets)
  need <- unique(unlist(c(sets, ctrl), use.names = FALSE))
  scores <- matrix(0, nrow = length(cell_ids), ncol = length(set_names),
                   dimnames = list(cell_ids, set_names))
  seen <- 0L
  for (h in handles) {
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    local <- intersect(cell_ids, h$cell_ids)
    if (!length(local)) next
    row_idx <- match(need, h$feature_ids)
    names(row_idx) <- need
    col_idx <- match(local, h$cell_ids)
    ord <- order(col_idx)
    local <- local[ord]
    col_idx <- col_idx[ord]
    runs <- .ic2_contig_runs(col_idx)
    for (run in runs) {
      for (s in seq.int(1L, length(run), by = chunk.size)) {
        e <- min(s + chunk.size - 1L, length(run))
        cols <- run[s:e]
        raw <- .ic2_to_dgCMatrix(.ic2_read_cols(h, cols))
        cn <- h$cell_ids[cols]
        present_rows <- row_idx[!is.na(row_idx)]
        block <- matrix(0, nrow = length(need), ncol = length(cols), dimnames = list(need, cn))
        if (length(present_rows)) {
          sub <- as.matrix(raw[present_rows, , drop = FALSE])
          block[names(present_rows), ] <- sub
        }
        for (nm in set_names) {
          g <- intersect(sets[[nm]], need)
          if (!length(g)) next
          if (identical(method, "sum")) {
            scores[cn, nm] <- colSums(block[g, , drop = FALSE])
          } else if (identical(method, "zscore")) {
            z <- block[g, , drop = FALSE]
            for (gi in seq_len(nrow(z))) {
              m <- mu[g[gi]]
              s <- sdv[g[gi]]
              if (!is.finite(m)) m <- 0
              if (!is.finite(s) || s == 0) s <- 1
              z[gi, ] <- (z[gi, ] - m) / s
            }
            scores[cn, nm] <- colMeans(z)
          } else {
            scores[cn, nm] <- colMeans(block[g, , drop = FALSE])
          }
          if (!is.null(ctrl) && length(ctrl[[nm]])) {
            cg <- intersect(ctrl[[nm]], need)
            if (length(cg)) {
              scores[cn, nm] <- scores[cn, nm] - colMeans(block[cg, , drop = FALSE])
            }
          }
        }
        seen <- seen + length(cols)
        if (is.function(progress)) progress(seen, length(cell_ids))
      }
    }
  }
  scores
}


#' Mean and SD of selected genes, one streamed pass.
#' @keywords internal
.ic2_score_moments <- function(handles, cell_ids, genes, chunk.size) {
  genes <- unique(as.character(genes))
  sum1 <- numeric(length(genes))
  sum2 <- numeric(length(genes))
  names(sum1) <- genes
  names(sum2) <- genes
  n <- 0L
  for (h in handles) {
    h <- .ic2_enrich_handle(h, load_dimnames = TRUE)
    local <- intersect(cell_ids, h$cell_ids)
    if (!length(local)) next
    row_idx <- match(genes, h$feature_ids)
    names(row_idx) <- genes
    col_idx <- match(local, h$cell_ids)
    ord <- order(col_idx)
    col_idx <- col_idx[ord]
    runs <- .ic2_contig_runs(col_idx)
    for (run in runs) {
      for (s in seq.int(1L, length(run), by = chunk.size)) {
        e <- min(s + chunk.size - 1L, length(run))
        raw <- .ic2_to_dgCMatrix(.ic2_read_cols(h, run[s:e]))
        present <- row_idx[!is.na(row_idx)]
        if (!length(present)) next
        sub <- as.matrix(raw[present, , drop = FALSE])
        sum1[names(present)] <- sum1[names(present)] + rowSums(sub)
        sum2[names(present)] <- sum2[names(present)] + rowSums(sub * sub)
        n <- n + ncol(sub)
      }
    }
  }
  if (n < 2L) stop("Need at least 2 cells for zscore.")
  mu <- sum1 / n
  v <- (sum2 - (sum1 * sum1) / n) / (n - 1L)
  v[v < 0] <- 0
  list(mu = mu, sd = sqrt(v))
}
