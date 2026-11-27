#' SODA portable I/O + open.soda()
#'
#' Professional layout (schema 1.0.0)  -  \strong{no \code{.rds} in the archive}:
#' \itemize{
#'   \item Assays / embeddings: HDF5 (\code{.h5})
#'   \item Tables (\code{cell_data}, maps, \ldots): TSV (\code{.tsv})
#'   \item Nested metadata / manifests: JSON (\code{.json}) + \code{soda.info}
#'   \item Logs / KPIs: \code{.txt} / \code{.log}
#'   \item UMAP models: \code{.uwot} via \code{uwot::save_uwot}
#'   \item Optional interchange: Parquet (\code{.parquet}) when \pkg{arrow} is available
#' }
#' Legacy \code{.rds} files are still \emph{readable} once, then should be removed.
#'
#' @name soda_portable
NULL

#' Unlink legacy .rds path(s) under a SODA archive
#' @keywords internal
.ic2_soda_unlink_rds <- function(...) {
  paths <- unlist(list(...), use.names = FALSE)
  paths <- as.character(paths)
  paths <- paths[nzchar(paths) & !is.na(paths)]
  if (!length(paths)) return(invisible(0L))
  ## Only delete paths that look like .rds files (never wipe dirs)
  paths <- paths[grepl("\\.rds$", paths, ignore.case = TRUE) & file.exists(paths)]
  if (!length(paths)) return(invisible(0L))
  unlink(paths, force = TRUE)
  invisible(length(paths))
}

#' Write tenx_h5 sidecar meta as JSON (never .rds)
#' @keywords internal
.ic2_soda_write_tenx_meta <- function(meta, path) {
  path <- as.character(path)[1L]
  path <- sub("\\.rds$", ".json", path, ignore.case = TRUE)
  if (!grepl("\\.json$", path, ignore.case = TRUE)) {
    path <- paste0(tools::file_path_sans_ext(path), ".json")
  }
  .ic2_soda_write_json(meta, path)
}

#' Read tenx_h5 sidecar meta (JSON preferred; legacy .rds fallback)
#' @keywords internal
.ic2_soda_read_tenx_meta <- function(path) {
  if (is.null(path) || !nzchar(as.character(path)[1L])) return(NULL)
  path <- as.character(path)[1L]
  if (!file.exists(path)) {
    alt <- c(
      sub("\\.rds$", ".json", path, ignore.case = TRUE),
      sub("\\.json$", ".rds", path, ignore.case = TRUE)
    )
    path <- alt[file.exists(alt)][1L]
    if (is.na(path)) return(NULL)
  }
  if (grepl("\\.json$", path, ignore.case = TRUE)) {
    return(jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE))
  }
  tryCatch(readRDS(path), error = function(e) NULL)
}

#' Write a portable JSON metadata file
#' @keywords internal
.ic2_soda_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!grepl("\\.json$", path, ignore.case = TRUE)) {
    path <- paste0(path, ".json")
  }
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  ## Drop legacy RDS sibling
  .ic2_soda_unlink_rds(sub("\\.json$", ".rds", path, ignore.case = TRUE))
  invisible(path)
}

#' Persist a uwot model with save_uwot (not RDS)
#' @keywords internal
.ic2_soda_write_uwot_model <- function(model, path) {
  if (is.null(model)) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ## Official uwot archive (Annoy index + embedding params); not R serialize
  if (file.exists(path)) unlink(path, force = TRUE)
  uwot::save_uwot(model, file = path, unload = TRUE, verbose = FALSE)
  ## Drop legacy RDS sibling if present
  unlink(sub("\\.uwot$", ".rds", path, ignore.case = TRUE), force = TRUE)
  unlink(file.path(dirname(path), "umap_model.rds"), force = TRUE)
  invisible(path)
}

#' Load a uwot model (prefers .uwot; falls back to legacy .rds)
#' @keywords internal
.ic2_soda_load_uwot_model <- function(path) {
  if (is.null(path) || !nzchar(as.character(path)[1L])) return(NULL)
  path <- as.character(path)[1L]
  if (!file.exists(path)) {
    alt <- c(
      sub("\\.rds$", ".uwot", path, ignore.case = TRUE),
      sub("\\.uwot$", ".rds", path, ignore.case = TRUE),
      file.path(dirname(path), "umap_model.uwot"),
      file.path(dirname(path), "umap_model.rds")
    )
    path <- alt[file.exists(alt)][1L]
    if (is.na(path)) return(NULL)
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }
  uwot::load_uwot(file = path, verbose = FALSE)
}

#' Read reduction / clustering meta (JSON preferred, legacy RDS fallback)
#' @keywords internal
.ic2_soda_read_meta <- function(path) {
  if (is.null(path) || !nzchar(as.character(path)[1L]) || !file.exists(path)) {
    return(NULL)
  }
  if (grepl("\\.json$", path, ignore.case = TRUE)) {
    return(jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE))
  }
  readRDS(path)
}

#' Resolve meta path: prefer reduction.json / *.json over *_meta.rds
#' @keywords internal
.ic2_soda_find_meta <- function(..., dir = NULL) {
  cands <- c(...)
  if (!is.null(dir)) {
    cands <- c(
      file.path(dir, "reduction.json"),
      file.path(dir, "clustering.json"),
      list.files(dir, pattern = "\\.json$", full.names = TRUE),
      list.files(dir, pattern = "_meta\\.rds$", full.names = TRUE),
      cands
    )
  }
  cands <- unique(as.character(cands))
  hit <- cands[nzchar(cands) & file.exists(cands)]
  if (length(hit)) hit[[1L]] else NA_character_
}

.ic2_soda_find_table <- function(stem) {
  cands <- c(
    paste0(stem, ".tsv"),
    paste0(stem, ".parquet"),
    paste0(stem, ".rds"), ## legacy read-only
    stem
  )
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[[1L]] else NA_character_
}

.ic2_soda_try_table <- function(stem) {
  path <- .ic2_soda_find_table(stem)
  if (is.na(path)) return(NULL)
  tryCatch(.ic2_soda_read_table(path), error = function(e) NULL)
}

.ic2_soda_read_table <- function(path_or_stem) {
  path <- if (file.exists(path_or_stem)) {
    path_or_stem
  } else {
    .ic2_soda_find_table(path_or_stem)
  }
  if (is.na(path) || !file.exists(path)) stop("Table not found: ", path_or_stem)
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Legacy Parquet table found but package 'arrow' is not installed: ", path,
           ". Re-save as TSV or install arrow.")
    }
    df <- as.data.frame(arrow::read_parquet(path))
  } else if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    df <- readRDS(path)
    if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!is.data.frame(df)) stop("Expected a table in: ", path)
  } else {
    df <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }
  if ("cell_id" %in% names(df)) {
    rn <- as.character(df$cell_id)
    if (!anyDuplicated(rn) && !anyNA(rn)) rownames(df) <- rn
  } else if ("feature_id" %in% names(df)) {
    rownames(df) <- as.character(df$feature_id)
  }
  attr(df, "soda_path") <- path
  df
}

.ic2_soda_write_table <- function(df, path_stem, id_col = NULL) {
  dir.create(dirname(path_stem), recursive = TRUE, showWarnings = FALSE)
  out <- df
  if (!is.null(rownames(df)) && !"cell_id" %in% names(out) &&
      !"feature_id" %in% names(out) && !"gene" %in% names(out)) {
    if ("feature_name" %in% names(out) || "feature_type" %in% names(out)) {
      out <- cbind(feature_id = rownames(df), out, stringsAsFactors = FALSE)
    } else {
      out <- cbind(cell_id = rownames(df), out, stringsAsFactors = FALSE)
    }
  }
  path <- paste0(path_stem, ".tsv")
  utils::write.table(
    out, file = path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
  )
  ## Never keep .rds twins in SODA; parquet optional later
  .ic2_soda_unlink_rds(paste0(path_stem, ".rds"))
  invisible(path)
}

#' Write a character / id vector as a one-column TSV (no .rds)
#' @keywords internal
.ic2_soda_write_ids <- function(ids, path_stem, col = "id") {
  dir.create(dirname(path_stem), recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(stats::setNames(list(as.character(ids)), col), stringsAsFactors = FALSE)
  .ic2_soda_write_table(df, path_stem)
}

#' Flatten a list of per-sample map data.frames and write as one table
#' @keywords internal
.ic2_soda_write_map_list <- function(map_list, path_stem) {
  if (is.null(map_list) || length(map_list) < 1L) {
    return(invisible(NA_character_))
  }
  ## already a single data.frame
  if (is.data.frame(map_list)) {
    return(.ic2_soda_write_table(map_list, path_stem))
  }
  parts <- lapply(seq_along(map_list), function(i) {
    df <- map_list[[i]]
    if (is.null(df) || !is.data.frame(df) || nrow(df) < 1L) return(NULL)
    if (!"map_index" %in% names(df)) {
      df$map_index <- as.integer(i)
    }
    ## preserve list name as sample hint when sample_id missing
    nm <- names(map_list)[i]
    if (!is.null(nm) && nzchar(nm) && !"sample_id" %in% names(df)) {
      df$sample_id <- nm
    }
    df
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) < 1L) return(invisible(NA_character_))
  ## align columns
  cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(df) {
    miss <- setdiff(cols, names(df))
    for (m in miss) df[[m]] <- NA
    df[cols]
  })
  big <- do.call(rbind, parts)
  rownames(big) <- NULL
  .ic2_soda_write_table(big, path_stem)
}

#' Read a flattened map table; optionally split back into a list by sample_id
#' @keywords internal
.ic2_soda_read_map_list <- function(path_stem, split = TRUE) {
  df <- .ic2_soda_try_table(path_stem)
  if (is.null(df)) return(NULL)
  if (!isTRUE(split)) return(df)
  if ("sample_id" %in% names(df)) {
    ids <- unique(as.character(df$sample_id))
    out <- lapply(ids, function(sid) {
      sub <- df[as.character(df$sample_id) == sid, , drop = FALSE]
      rownames(sub) <- NULL
      if ("map_index" %in% names(sub)) sub$map_index <- NULL
      sub
    })
    names(out) <- ids
    return(out)
  }
  if ("map_index" %in% names(df)) {
    ids <- sort(unique(as.integer(df$map_index)))
    out <- lapply(ids, function(i) {
      sub <- df[as.integer(df$map_index) == i, , drop = FALSE]
      rownames(sub) <- NULL
      sub$map_index <- NULL
      sub
    })
    return(out)
  }
  list(df)
}

.ic2_soda_try_json_or_rds <- function(stem) {
  cands <- c(paste0(stem, ".json"), paste0(stem, ".rds"), stem)
  path <- cands[file.exists(cands)][1L]
  if (is.na(path)) return(NULL)
  tryCatch({
    if (grepl("\\.json$", path, ignore.case = TRUE)) {
      jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
    } else {
      ## legacy .rds  -  read once
      readRDS(path)
    }
  }, error = function(e) NULL)
}

.ic2_write_kv_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  keys <- names(x)
  vals <- vapply(keys, function(k) {
    v <- x[[k]]
    if (is.null(v)) return("null")
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.list(v) || length(v) > 1L) {
      return(as.character(jsonlite::toJSON(v, auto_unbox = TRUE, null = "null")))
    }
    as.character(v)
  }, character(1))
  utils::write.table(
    data.frame(key = keys, value = vals, stringsAsFactors = FALSE),
    file = path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  invisible(path)
}

.ic2_read_kv_tsv <- function(path) {
  df <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("key", "value") %in% names(df))) stop("Expected key/value TSV: ", path)
  out <- as.list(stats::setNames(df$value, df$key))
  for (k in names(out)) {
    v <- out[[k]]
    if (identical(v, "null") || identical(v, "NULL")) out[[k]] <- NULL
    else if (identical(v, "true")) out[[k]] <- TRUE
    else if (identical(v, "false")) out[[k]] <- FALSE
    else if (grepl("^\\[|^\\{", v)) {
      out[[k]] <- tryCatch(jsonlite::fromJSON(v), error = function(e) v)
    } else {
      num <- suppressWarnings(as.numeric(v))
      if (!is.na(num) && identical(as.character(num), v)) out[[k]] <- num
    }
  }
  out
}

.ic2_h5_write_dense <- function(mat, path, group = "data", overwrite = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && isTRUE(overwrite)) file.remove(path)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- if (h5$exists(group)) h5[[group]] else h5$create_group(group)
  ## Explicit 2-D dataset (avoid compound / ambiguous assign).
  ## hdf5r stores R matrices with dims reversed vs h5py; Python write_dense_h5
  ## matches that on-disk layout. Readers auto-orient via row_ids/col_ids.
  g$create_dataset(
    "matrix",
    robj = mat,
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    chunk_dims = pmin(dim(mat), c(10000L, ncol(mat)))
  )
  if (!is.null(rownames(mat))) {
    g[["row_ids"]] <- as.character(rownames(mat))
  }
  if (!is.null(colnames(mat))) {
    g[["col_ids"]] <- as.character(colnames(mat))
  }
  h5$flush()
  invisible(path)
}

.ic2_h5_read_dense <- function(path, group = "data") {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(as.matrix(readRDS(path)))
  }
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  g <- if (nzchar(group) && h5$exists(group)) h5[[group]] else h5
  if (!g$exists("matrix")) {
    if (h5$exists("matrix")) g <- h5 else stop("No 'matrix' dataset in: ", path)
  }
  ds <- g[["matrix"]]
  ## hdf5r: 2-D datasets must be subset with [,] not [] ("1 vs 2" error)
  nd <- length(ds$dims)
  if (nd >= 2L) {
    mat <- ds[, ]
  } else if (nd == 1L) {
    mat <- ds[]
  } else {
    mat <- ds$read()
  }
  mat <- as.matrix(mat)
  rn <- NULL
  cn <- NULL
  if (g$exists("row_ids")) {
    rid <- g[["row_ids"]]
    rn <- as.character(if (length(rid$dims) >= 1L) rid[] else rid$read())
  }
  if (g$exists("col_ids")) {
    cid <- g[["col_ids"]]
    cn <- as.character(if (length(cid$dims) >= 1L) cid[] else cid$read())
  }
  ## Orient so nrow == length(row_ids), ncol == length(col_ids)
  ## (handles legacy archives + Python-written files).
  if (!is.null(rn) && !is.null(cn)) {
    if (nrow(mat) == length(cn) && ncol(mat) == length(rn) &&
        !(nrow(mat) == length(rn) && ncol(mat) == length(cn))) {
      mat <- t(mat)
    }
  } else if (!is.null(rn) && nrow(mat) != length(rn) && ncol(mat) == length(rn)) {
    mat <- t(mat)
  }
  if (!is.null(rn) && length(rn) == nrow(mat)) rownames(mat) <- rn
  if (!is.null(cn) && length(cn) == ncol(mat)) colnames(mat) <- cn
  mat
}


#' Make a path relative to the SODA root (portable soda.json)
#' @keywords internal
.ic2_soda_path_rel <- function(root, path) {
  if (is.null(path) || length(path) < 1L) return(NULL)
  path <- as.character(path)[[1L]]
  if (!nzchar(path) || identical(path, "null")) return(NULL)
  root <- as.character(root)[[1L]]
  if (!nzchar(root)) return(path)
  root_n <- normalizePath(root, winslash = "/", mustWork = FALSE)
  ## Already relative?
  is_abs <- grepl("^(~/|/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (!is_abs) {
    return(gsub("\\\\", "/", path, fixed = FALSE))
  }
  path_n <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
  root_slash <- paste0(root_n, "/")
  if (identical(path_n, root_n)) return(".")
  if (startsWith(path_n, root_slash)) {
    return(substring(path_n, nchar(root_slash) + 1L))
  }
  ## Outside archive  -  keep absolute (rare)
  path_n
}

#' Resolve a soda.json path against the SODA root
#' @keywords internal
.ic2_soda_path_abs <- function(root, path) {
  if (is.null(path) || length(path) < 1L) return(NULL)
  path <- as.character(path)[[1L]]
  if (!nzchar(path)) return(NULL)
  if (identical(path, ".")) {
    return(normalizePath(root, winslash = "/", mustWork = FALSE))
  }
  is_abs <- grepl("^(~/|/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (is_abs) {
    return(normalizePath(path.expand(path), winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(root, path), winslash = "/", mustWork = FALSE)
}

#' Recursively convert character paths in a list to SODA-relative
#' @keywords internal
.ic2_soda_paths_make_rel <- function(root, x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L) {
    return(.ic2_soda_path_rel(root, x))
  }
  if (is.list(x)) {
    return(lapply(x, function(z) .ic2_soda_paths_make_rel(root, z)))
  }
  x
}

#' Recursively resolve character paths in a list against SODA root
#' @keywords internal
.ic2_soda_paths_make_abs <- function(root, x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L) {
    ## Only treat as path if it looks like one
    if (grepl("[/\\\\]|\\.tsv$|\\.json$|\\.h5$|\\.rds$|\\.uwot$|metadata|modalities|reductions", x)) {
      return(.ic2_soda_path_abs(root, x))
    }
    return(x)
  }
  if (is.list(x)) {
    return(lapply(x, function(z) .ic2_soda_paths_make_abs(root, z)))
  }
  x
}




#' Package version string that never errors (dev / uninstalled OK)
#' @keywords internal
.ic2_package_version <- function() {
  ver <- tryCatch(
    as.character(utils::packageVersion("iCellR2")),
    error = function(e) NA_character_
  )
  if (is.character(ver) && length(ver) == 1L && !is.na(ver) && nzchar(ver)) {
    return(ver)
  }
  desc <- ""
  ns_path <- tryCatch({
    if (isNamespaceLoaded("iCellR2")) {
      as.character(getNamespaceInfo("iCellR2", "path"))[1L]
    } else {
      NA_character_
    }
  }, error = function(e) NA_character_)
  if (is.character(ns_path) && length(ns_path) == 1L && !is.na(ns_path) && nzchar(ns_path)) {
    cand <- file.path(ns_path, "DESCRIPTION")
    if (file.exists(cand)) desc <- cand
  }
  if (!nzchar(desc)) {
    desc <- tryCatch(
      system.file("DESCRIPTION", package = "iCellR2"),
      error = function(e) ""
    )
  }
  if (!nzchar(desc) || !file.exists(desc)) {
    for (p in c(
      file.path(getwd(), "DESCRIPTION"),
      file.path(getwd(), "iCellR2", "DESCRIPTION"),
      file.path(dirname(getwd()), "iCellR2", "DESCRIPTION")
    )) {
      if (file.exists(p)) {
        desc <- p
        break
      }
    }
  }
  if (nzchar(desc) && file.exists(desc)) {
    fields <- tryCatch(
      read.dcf(desc, fields = "Version")[1L, "Version"],
      error = function(e) NA_character_
    )
    if (is.character(fields) && length(fields) == 1L && !is.na(fields) && nzchar(fields)) {
      return(unname(as.character(fields)))
    }
  }
  "0.0.0-dev"
}

#' Empty multimodal slot (RNA+ADT, RNA+ATAC, ...). Not batch alignment.
#' @keywords internal
.ic2_empty_multimodal <- function() {
  list(
    rna_atac = NULL, rna_adt = NULL, rna_spatial = NULL,
    rna_visium = NULL, rna_xenium = NULL, multiome = NULL, misc = list()
  )
}

#' Compact multimodal entries for soda.json (character paths made relative).
#' @keywords internal
.ic2_multimodal_manifest <- function(x) {
  root <- x@storage$root_dir
  out <- list()
  mm <- x@multimodal %||% list()
  for (nm in names(mm)) {
    z <- mm[[nm]]
    if (is.null(z) || length(z) < 1L || identical(nm, "misc")) next
    out[[nm]] <- .ic2_multimodal_rel(z, root)
  }
  out
}

.ic2_multimodal_rel <- function(z, root) {
  if (is.character(z) && length(z) == 1L) {
    return(.ic2_soda_path_rel(root, z))
  }
  if (!is.list(z)) return(z)
  lapply(z, function(el) .ic2_multimodal_rel(el, root))
}

#' Read multimodal/*/multimodal.json. Falls back to a legacy integrations/ folder.
#' @keywords internal
.ic2_multimodal_from_disk <- function(root, slots = .ic2_empty_multimodal()) {
  base <- file.path(root, "multimodal")
  legacy <- file.path(root, "integrations")
  if (!dir.exists(base) && dir.exists(legacy)) base <- legacy
  if (!dir.exists(base)) return(slots)
  kids <- list.dirs(base, full.names = FALSE, recursive = FALSE)
  for (k in kids) {
    meta <- file.path(base, k, "multimodal.json")
    if (file.exists(meta)) {
      parsed <- tryCatch(
        jsonlite::fromJSON(meta, simplifyVector = FALSE),
        error = function(e) NULL
      )
      slots[[k]] <- parsed %||% list(dir = file.path(base, k))
    } else {
      slots[[k]] <- list(dir = file.path(base, k))
    }
  }
  slots
}

#' Write soda.json (+ soda.info) from an iCellR2 object
#'
#' Paths inside the manifest are stored \strong{relative to the SODA root}
#' so the archive can be moved or shared. Also refreshes README.md and a
#' light MANIFEST of core files.
#'
#' @keywords internal
.ic2_soda_write_manifest <- function(x) {
  root <- x@storage$root_dir
  if (is.null(root) || !nzchar(root)) return(invisible(NULL))
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  mods <- list()
  for (mod in names(x@modalities %||% list())) {
    m <- x@modalities[[mod]]
    assays <- list()
    for (an in names(m$assays %||% list())) {
      ar <- m$assays[[an]]
      apath <- ar$handle$path %||% ar$path %||% NULL
      assays[[an]] <- list(
        active = identical(m$active_assay, an),
        aggregation = ar$aggregation_type %||% "single",
        path = .ic2_soda_path_rel(root, apath)
      )
    }
    reductions <- list()
    for (rn in names(m$reductions %||% list())) {
      rr <- m$reductions[[rn]]
      reductions[[rn]] <- list(
        embeddings = .ic2_soda_path_rel(root, rr$cell_embeddings %||% NULL),
        loadings = .ic2_soda_path_rel(root, rr$feature_loadings %||% NULL)
      )
    }
    mods[[mod]] <- list(
      active_assay = m$active_assay %||% NULL,
      active_reduction = m$active_reduction %||% NULL,
      assays = assays,
      reductions = reductions
    )
  }

  meta_paths <- x@storage$metadata_paths %||% list()
  meta_rel <- lapply(meta_paths, function(p) .ic2_soda_path_rel(root, p))

  man <- list(
    format = "SODA",
    format_name = "Single Omics Data Archive",
    soda_version = "1.0.0",
    schema_version = "1.0.0",
    package_version = .ic2_package_version(),
    root_dir = ".",
    project_name = x@project_info$project_name %||% NA_character_,
    modalities = mods,
    metadata_paths = meta_rel,
    multimodal = .ic2_multimodal_manifest(x),
    created = as.character(Sys.time()),
    portable = TRUE
  )
  path <- file.path(root, "soda.json")
  jsonlite::write_json(man, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  .ic2_write_kv_tsv(list(
    format = "SODA",
    soda_version = man$soda_version,
    schema_version = man$schema_version,
    package_version = man$package_version,
    project_name = man$project_name,
    root_dir = ".",
    portable = "TRUE"
  ), file.path(root, "soda.info"))
  ## Drop legacy twin if present
  legacy <- file.path(root, "soda.tsv")
  if (file.exists(legacy)) unlink(legacy)
  .ic2_soda_write_architecture(root)
  try({
    inv <- .ic2_soda_write_inventory(root, man = man, x = x)
    .ic2_soda_write_readme(root, man = man, x = x)
    .ic2_soda_write_html_index(root, man = man, x = x, inventory = inv)
    .ic2_soda_write_file_manifest(root)
  }, silent = TRUE)
  invisible(path)
}


#' HTML-escape text for SODA index pages
#' @keywords internal
.ic2_html_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

#' Turn UTF-8 / C-locale dumps into safe HTML entities (arrows, dashes, ...)
#' @keywords internal
.ic2_html_fix_utf8 <- function(x) {
  x <- paste(as.character(x %||% ""), collapse = "\n")
  ## Literal C-locale dumps (often appear after paste/write under C locale)
  lit <- c(
    "<e2><86><92>" = "&rarr;",
    "<e2><86><90>" = "&larr;",
    "<e2><86><94>" = "&harr;",
    "<e2><80><94>" = "&mdash;",
    "<e2><80><93>" = "&ndash;",
    "<e2><80><a6>" = "&hellip;",
    "<e2><95><90>" = "&#9552;",
    "<e2><94><80>" = "&#9472;",
    "<e2><94><82>" = "&#9474;",
    "<e2><94><9c>" = "&#9500;",
    "<e2><94><94>" = "&#9492;",
    "<c2><b7>" = "&middot;",
    "<c3><97>" = "&times;",
    ## Escaped dumps (if html_escape ran first)
    "&lt;e2&gt;&lt;86&gt;&lt;92&gt;" = "&rarr;",
    "&lt;e2&gt;&lt;86&gt;&lt;90&gt;" = "&larr;",
    "&lt;e2&gt;&lt;86&gt;&lt;94&gt;" = "&harr;",
    "&lt;e2&gt;&lt;80&gt;&lt;94&gt;" = "&mdash;",
    "&lt;e2&gt;&lt;80&gt;&lt;93&gt;" = "&ndash;",
    "&lt;e2&gt;&lt;80&gt;&lt;a6&gt;" = "&hellip;",
    "&lt;e2&gt;&lt;95&gt;&lt;90&gt;" = "&#9552;",
    "&lt;c2&gt;&lt;b7&gt;" = "&middot;",
    "&lt;c3&gt;&lt;97&gt;" = "&times;",
    ## R <U+XXXX> forms
    "<U+2192>" = "&rarr;",
    "<U+2190>" = "&larr;",
    "<U+2194>" = "&harr;",
    "<U+2014>" = "&mdash;",
    "<U+2013>" = "&ndash;",
    "<U+2026>" = "&hellip;",
    "<U+00B7>" = "&middot;",
    "<U+00D7>" = "&times;",
    "<U+2550>" = "&#9552;",
    "&lt;U+2192&gt;" = "&rarr;",
    "&lt;U+2190&gt;" = "&larr;",
    "&lt;U+2194&gt;" = "&harr;",
    "&lt;U+2014&gt;" = "&mdash;",
    "&lt;U+2013&gt;" = "&ndash;",
    "&lt;U+2026&gt;" = "&hellip;",
    "&lt;U+00B7&gt;" = "&middot;",
    "&lt;U+00D7&gt;" = "&times;",
    "&lt;U+2550&gt;" = "&#9552;"
  )
  for (nm in names(lit)) x <- gsub(nm, lit[[nm]], x, fixed = TRUE)
  x <- tryCatch({
    x <- gsub(intToUtf8(0x2192L), "&rarr;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2190L), "&larr;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2194L), "&harr;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2014L), "&mdash;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2013L), "&ndash;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2026L), "&hellip;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x00B7L), "&middot;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x00D7L), "&times;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2550L), "&#9552;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x0394L), "&Delta;", x, fixed = TRUE, useBytes = TRUE)
    x
  }, error = function(e) x)
  x
}

#' Human labels for well-known SODA root documents
#' @keywords internal
.ic2_soda_doc_catalog <- function() {
  list(
    list(file = "soda.json", label = "Project metadata", blurb = "Machine-readable archive manifest (paths, versions, modalities)."),
    list(file = "soda.info", label = "Project metadata (text)", blurb = "Key=value twin of soda.json for quick skim."),
    list(file = "README.md", label = "README", blurb = "Human intro and how to open this archive."),
    list(file = "SODA_SPEC.md", label = "Format specification", blurb = "Stable SODA layout contract."),
    list(file = "MANIFEST.txt", label = "File inventory (light)", blurb = "Core path list with sizes."),
    list(file = "inventory.json", label = "Detailed inventory", blurb = "Per-file roles, dimensions, and HDF5 structure."),
    list(file = "index.html", label = "Archive home page", blurb = "This browser-friendly overview."),
    list(file = "benchmark.txt", label = "Timing / benchmark", blurb = "Wall / CPU / RAM KPIs from the run."),
    list(file = "iCellR2.log", label = "Run log", blurb = "Tool command log for this project."),
    list(file = "architecture.txt", label = "Folder layout", blurb = "On-disk tree and slot mapping notes.")
  )
}

#' Infer a short human role from a path inside a SODA root
#' @keywords internal
.ic2_soda_file_role <- function(rel) {
  rel <- as.character(rel)[1L]
  bn <- basename(rel)
  catalog <- .ic2_soda_doc_catalog()
  for (z in catalog) {
    if (identical(bn, z$file) || identical(rel, z$file)) {
      return(list(label = z$label, blurb = z$blurb))
    }
  }
  if (grepl("^metadata/cell_data\\.tsv$", rel)) {
    return(list(label = "Cell metadata table", blurb = "One row per cell (QC, clusters, annotations)."))
  }
  if (grepl("^metadata/sample_data\\.tsv$", rel)) {
    return(list(label = "Sample metadata table", blurb = "One row per sample / donor / batch."))
  }
  if (grepl("^metadata/.*feature_data\\.tsv$", rel)) {
    return(list(label = "Feature metadata table", blurb = "Gene / feature annotations for a modality."))
  }
  if (grepl("^metadata/.*gene_stats\\.tsv$", rel)) {
    return(list(label = "Gene statistics", blurb = "Per-gene means / sparsity used for modeling."))
  }
  if (grepl("^metadata/.*gene_model\\.tsv$", rel)) {
    return(list(label = "Gene model / HVG list", blurb = "Selected features for PCA and friends."))
  }
  if (grepl("^metadata/project_info\\.json$", rel)) {
    return(list(label = "Project info", blurb = "Project name, organism, genome, technologies."))
  }
  if (grepl("^metadata/params\\.json$", rel)) {
    return(list(label = "Analysis parameters", blurb = "Recorded tool parameters."))
  }
  if (grepl("^metadata/commands\\.json$", rel)) {
    return(list(label = "Command history", blurb = "Commands that built this archive."))
  }
  if (grepl("^metadata/version\\.json$", rel)) {
    return(list(label = "Version stamps", blurb = "Object / schema / package versions."))
  }
  if (grepl("^metadata/session_stats\\.json$", rel)) {
    return(list(label = "Session stats", blurb = "Progress / RAM session snapshot."))
  }
  if (grepl("^metadata/clustering/", rel)) {
    return(list(label = "Clustering output", blurb = "Cluster IDs and related tables."))
  }
  if (grepl("/assays/[^/]+/.+\\.h5$", rel, ignore.case = TRUE)) {
    assay <- sub("^.*/assays/([^/]+)/.*$", "\\1", rel)
    mod <- if (grepl("modalities/", rel)) sub("^modalities/([^/]+)/.*$", "\\1", rel) else "rna"
    return(list(
      label = sprintf("%s / %s matrix", mod, assay),
      blurb = sprintf("HDF5 expression matrix for assay '%s' (modality '%s').", assay, mod)
    ))
  }
  if (grepl("cell_embeddings\\.h5$", rel, ignore.case = TRUE)) {
    red <- basename(dirname(rel))
    return(list(label = sprintf("%s cell embeddings", red), blurb = "Dense cell  x  component embedding matrix."))
  }
  if (grepl("feature_loadings\\.h5$", rel, ignore.case = TRUE)) {
    red <- basename(dirname(rel))
    return(list(label = sprintf("%s feature loadings", red), blurb = "Dense feature  x  component loadings matrix."))
  }
  if (grepl("reduction\\.json$", rel)) {
    red <- basename(dirname(rel))
    return(list(label = sprintf("%s reduction metadata", red), blurb = "Parameters and provenance for this reduction."))
  }
  if (grepl("assay\\.json$", rel)) {
    return(list(label = "Assay sidecar", blurb = "Assay-level metadata (dims, source, lazy flags)."))
  }
  if (grepl("tenx_meta\\.json$", rel)) {
    return(list(label = "10x pointer metadata", blurb = "Native 10x HDF5 group / feature field notes."))
  }
  if (grepl("\\.uwot$", rel, ignore.case = TRUE)) {
    return(list(label = "UMAP model", blurb = "Saved uwot model for transform / reopen."))
  }
  if (grepl("^images/", rel)) {
    return(list(label = "Plot / image", blurb = "Exported figure for this project."))
  }
  ext <- tolower(tools::file_ext(bn))
  label <- switch(ext,
    "h5" = "HDF5 data file",
    "hdf5" = "HDF5 data file",
    "tsv" = "Table",
    "json" = "JSON metadata",
    "txt" = "Text document",
    "log" = "Log file",
    "md" = "Markdown document",
    "parquet" = "Parquet table",
    "html" = "HTML page",
    "png" = "PNG image",
    "Data file"
  )
  list(label = label, blurb = paste0("File inside the SODA archive (", bn, ")."))
}

#' Peek HDF5 structure without loading full matrices
#' @keywords internal
.ic2_soda_inspect_h5 <- function(path) {
  out <- list(
    kind = "unknown",
    format = NA_character_,
    group = NA_character_,
    dims = NULL,
    nnz = NA_real_,
    dtype = NA_character_,
    datasets = character(0),
    notes = character(0)
  )
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    out$notes <- "hdf5r not available"
    return(out)
  }
  if (!file.exists(path)) {
    out$notes <- "missing"
    return(out)
  }

  .safe_dims <- function(obj) {
    ## numeric: sparse data vectors can exceed .Machine$integer.max
    tryCatch(as.numeric(obj$dims), error = function(e) NULL)
  }
  .safe_dtype <- function(obj) {
    tryCatch(as.character(obj$get_type()$get_class()), error = function(e) NA_character_)
  }
  .ds_info <- function(obj, name) {
    tryCatch({
      if (inherits(obj, "H5D")) {
        d <- .safe_dims(obj)
        sprintf("%s [%s] dtype=%s", name, paste(d %||% "?", collapse = "x"), .safe_dtype(obj) %||% "?")
      } else if (inherits(obj, "H5Group")) {
        kids <- tryCatch(names(obj), error = function(e) character(0))
        sprintf("%s/ (group: %s)", name, paste(kids, collapse = ", "))
      } else name
    }, error = function(e) name)
  }

  h5 <- tryCatch(hdf5r::H5File$new(path, mode = "r"), error = function(e) NULL)
  if (is.null(h5)) {
    out$notes <- "unreadable"
    return(out)
  }
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  top <- tryCatch(names(h5), error = function(e) character(0))

  ## Dense reduction layout: /data/matrix (+ row_ids, col_ids)
  if ("data" %in% top) {
    g <- tryCatch(h5[["data"]], error = function(e) NULL)
    if (inherits(g, "H5Group") && "matrix" %in% names(g)) {
      ds <- g[["matrix"]]
      out$kind <- "dense"
      out$format <- "iCellR2_dense_matrix"
      out$group <- "data/matrix"
      out$dims <- .safe_dims(ds)
      out$dtype <- .safe_dtype(ds)
      out$datasets <- vapply(names(g), function(k) .ds_info(g[[k]], paste0("data/", k)), character(1))
      return(out)
    }
  }

  ## 10x native
  if ("matrix" %in% top) {
    g <- tryCatch(h5[["matrix"]], error = function(e) NULL)
    if (!is.null(g) && inherits(g, "H5Group")) {
      nms <- names(g)
      if (all(c("data", "indices", "indptr") %in% nms)) {
        out$kind <- "tenx_h5"
        out$format <- "tenx_h5"
        out$group <- "matrix"
        if ("shape" %in% nms) {
          sh <- tryCatch(as.numeric(g[["shape"]][]), error = function(e) NULL)
          if (length(sh)) out$dims <- sh
        }
        if ("data" %in% nms) {
          dd <- .safe_dims(g[["data"]])
          out$nnz <- if (length(dd)) as.numeric(dd[1L]) else NA_real_
          out$dtype <- .safe_dtype(g[["data"]])
        }
        out$datasets <- vapply(nms, function(nm) .ds_info(g[[nm]], paste0("matrix/", nm)), character(1))
        return(out)
      }
    }
  }

  ## iCellR2 CSC / dense-with-dims groups
  for (nm in top) {
    g <- tryCatch(h5[[nm]], error = function(e) NULL)
    if (is.null(g) || !inherits(g, "H5Group")) next
    nms <- names(g)
    fmt <- NA_character_
    if ("format" %in% nms) {
      fmt <- tryCatch(as.character(g[["format"]][][1L]), error = function(e) NA_character_)
    }
    if (all(c("i", "p", "x", "dims") %in% nms) || identical(fmt, "iCellR2_csc_v1")) {
      out$kind <- "icellr2_csc"
      out$format <- if (!is.na(fmt) && nzchar(fmt)) fmt else "iCellR2_csc_v1"
      out$group <- nm
      dims <- tryCatch(as.numeric(g[["dims"]][]), error = function(e) NULL)
      if (length(dims)) out$dims <- dims
      dd <- .safe_dims(g[["x"]])
      out$nnz <- if (length(dd)) as.numeric(dd[1L]) else NA_real_
      out$dtype <- .safe_dtype(g[["x"]])
      out$datasets <- vapply(nms, function(k) .ds_info(g[[k]], paste0(nm, "/", k)), character(1))
      return(out)
    }
    if ("matrix" %in% nms && inherits(tryCatch(g[["matrix"]], error = function(e) NULL), "H5D")) {
      ds <- g[["matrix"]]
      out$kind <- "dense"
      out$format <- if (!is.na(fmt) && nzchar(fmt)) fmt else "iCellR2_dense_matrix"
      out$group <- paste0(nm, "/matrix")
      out$dims <- .safe_dims(ds)
      out$dtype <- .safe_dtype(ds)
      out$datasets <- vapply(nms, function(k) .ds_info(g[[k]], paste0(nm, "/", k)), character(1))
      return(out)
    }
    if (("x" %in% nms && "dims" %in% nms) || identical(fmt, "iCellR2_dense_v1")) {
      out$kind <- "icellr2_dense"
      out$format <- if (!is.na(fmt) && nzchar(fmt)) fmt else "iCellR2_dense_v1"
      out$group <- nm
      if ("dims" %in% nms) {
        dims <- tryCatch(as.numeric(g[["dims"]][]), error = function(e) NULL)
        if (length(dims)) out$dims <- dims
      } else if ("x" %in% nms) {
        out$dims <- .safe_dims(g[["x"]])
      }
      if ("x" %in% nms) out$dtype <- .safe_dtype(g[["x"]])
      out$datasets <- vapply(nms, function(k) .ds_info(g[[k]], paste0(nm, "/", k)), character(1))
      return(out)
    }
  }

  out$datasets <- vapply(top, function(nm) .ds_info(h5[[nm]], nm), character(1))
  for (nm in top) {
    obj <- tryCatch(h5[[nm]], error = function(e) NULL)
    if (inherits(obj, "H5D")) {
      out$kind <- "dense"
      out$format <- "dense_dataset"
      out$group <- nm
      out$dims <- .safe_dims(obj)
      out$dtype <- .safe_dtype(obj)
      return(out)
    }
  }
  out
}

#' Lightweight non-H5 file summary (rows/cols/keys)
#' @keywords internal
.ic2_soda_inspect_tableish <- function(path, rel) {
  ext <- tolower(tools::file_ext(path))
  info <- list(rows = NA_integer_, cols = NA_integer_, columns = character(0), keys = character(0))
  skip_txt <- grepl("(benchmark|architecture|MANIFEST|iCellR2\\.log|\\.info\\.txt)$", basename(path), ignore.case = TRUE)
  if (ext %in% c("tsv", "csv") || (ext == "txt" && !skip_txt)) {
    tryCatch({
      hdr <- readLines(path, n = 1L, warn = FALSE)
      if (length(hdr) == 1L && nzchar(hdr)) {
        sep <- if (grepl("\t", hdr, fixed = TRUE)) "\t" else ","
        cols <- strsplit(hdr, sep, fixed = TRUE)[[1L]]
        info$columns <- cols
        info$cols <- length(cols)
      }
      sz <- file.info(path)$size %||% 0
      if (!is.na(sz) && sz <= 25e6) {
        nlines <- length(readLines(path, warn = FALSE))
        if (nlines >= 1L) info$rows <- as.integer(nlines - 1L)
      } else {
        ## Large tables: approximate via wc when available
        wc_line <- tryCatch(
          system2("wc", c("-l", path), stdout = TRUE, stderr = FALSE),
          error = function(e) NA_character_
        )
        if (length(wc_line) >= 1L && !is.na(wc_line[1L])) {
          n <- suppressWarnings(as.integer(strsplit(trimws(wc_line[1L]), "[[:space:]]+")[[1L]][1L]))
          if (!is.na(n) && n >= 1L) info$rows <- as.integer(n - 1L)
        }
      }
    }, error = function(e) NULL)
  }
  if (ext == "json") {
    js <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
    if (is.list(js)) info$keys <- names(js)
  }
  info
}

#' Format one inventory entry as human .info.txt text
#' @keywords internal
.ic2_soda_format_info_txt <- function(entry) {
  lines <- c(
    paste0("# ", entry$label %||% entry$path),
    paste0("path: ", entry$path),
    paste0("role: ", entry$label %||% ""),
    paste0("description: ", entry$blurb %||% ""),
    paste0("size_bytes: ", entry$size_bytes %||% ""),
    paste0("modified: ", entry$mtime %||% "")
  )
  if (!is.null(entry$h5)) {
    h <- entry$h5
    lines <- c(
      lines,
      "type: hdf5",
      paste0("h5_kind: ", h$kind %||% ""),
      paste0("h5_format: ", h$format %||% ""),
      paste0("h5_group: ", h$group %||% ""),
      paste0("dimensions: ", if (length(h$dims)) paste(h$dims, collapse = " x ") else ""),
      paste0("nnz: ", if (!is.null(h$nnz) && !is.na(h$nnz)) h$nnz else ""),
      paste0("dtype: ", h$dtype %||% ""),
      "datasets:"
    )
    if (length(h$datasets)) {
      lines <- c(lines, paste0("  - ", h$datasets))
    } else {
      lines <- c(lines, "  - (none listed)")
    }
  } else if (!is.null(entry$table)) {
    t <- entry$table
    lines <- c(
      lines,
      paste0("type: ", entry$ext %||% "table"),
      paste0("rows: ", t$rows %||% ""),
      paste0("cols: ", t$cols %||% "")
    )
    if (length(t$columns)) {
      show <- utils::head(t$columns, 40L)
      lines <- c(lines, paste0("columns: ", paste(show, collapse = ", ")),
                 if (length(t$columns) > 40L) paste0("  ... +", length(t$columns) - 40L, " more") else NULL)
    }
    if (length(t$keys)) {
      lines <- c(lines, paste0("json_keys: ", paste(t$keys, collapse = ", ")))
    }
  } else {
    lines <- c(lines, paste0("type: ", entry$ext %||% "file"))
  }
  lines <- c(lines, paste0("generated: ", entry$generated %||% as.character(Sys.time())), "")
  lines
}

#' Build inventory.json for a SODA root (no per-file *.info.txt sidecars)
#' @keywords internal
.ic2_soda_write_inventory <- function(root, man = NULL, x = NULL) {
  root <- as.character(root)[1L]
  if (!nzchar(root) || !dir.exists(root)) return(invisible(NULL))
  root_n <- normalizePath(root, winslash = "/", mustWork = TRUE)
  s <- .ic2_soda_summary(root_n, man = man, x = x)

  skip_re <- "(^|/)(\\.DS_Store|inventory\\.json|.*\\.info\\.txt|^view/)"
  allowed_re <- "\\.(json|h5|hdf5|parquet|tsv|txt|log|info|uwot|md|html|png|svg|pdf)$"
  all_files <- list.files(root_n, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  all_files <- all_files[!file.info(all_files)$isdir %in% TRUE]
  rels <- vapply(all_files, function(f) .ic2_soda_path_rel(root_n, f) %||% basename(f), character(1))
  keep <- !grepl(skip_re, rels) & grepl(allowed_re, rels, ignore.case = TRUE)
  rels <- unique(c(rels[keep], vapply(.ic2_soda_doc_catalog(), `[[`, character(1), "file")))
  ## Generated HTML chrome (keep source files + index.html only)
  viewer_html <- c(
    "docs.html",
    "soda.info.html", "README.html", "SODA_SPEC.html", "architecture.html",
    "iCellR2.log.html", "benchmark.html", "MANIFEST.html"
  )
  rels <- setdiff(rels, viewer_html)
  rels <- rels[nzchar(rels) & file.exists(file.path(root_n, rels))]
  rels <- sort(unique(rels[!grepl(skip_re, rels)]))

  generated <- as.character(Sys.time())
  entries <- vector("list", length(rels))
  for (i in seq_along(rels)) {
    rel <- rels[[i]]
    abs <- file.path(root_n, rel)
    fi <- file.info(abs)
    role <- .ic2_soda_file_role(rel)
    ext <- tolower(tools::file_ext(rel))
    entry <- list(
      path = rel,
      label = role$label,
      blurb = role$blurb,
      ext = ext,
      size_bytes = as.numeric(fi$size %||% NA_real_),
      mtime = as.character(fi$mtime %||% NA),
      generated = generated
    )
    if (ext %in% c("h5", "hdf5")) {
      entry$h5 <- tryCatch(.ic2_soda_inspect_h5(abs), error = function(e) {
        list(kind = "unknown", format = NA_character_, group = NA_character_,
             dims = NULL, nnz = NA_real_, dtype = NA_character_,
             datasets = character(0), notes = conditionMessage(e))
      })
    } else if (ext %in% c("tsv", "csv", "json", "txt", "md", "log", "info")) {
      if (!(ext %in% c("log") && isTRUE(fi$size > 5e6))) {
        entry$table <- .ic2_soda_inspect_tableish(abs, rel)
      }
    }
    entries[[i]] <- entry
  }

  ## Drop legacy per-file sidecars (inventory.json is the single source)
  old_infos <- list.files(root_n, pattern = "\\.info\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(old_infos)) try(unlink(old_infos), silent = TRUE)

  inv <- list(
    format = "SODA_inventory",
    inventory_version = "1.0.0",
    generated = generated,
    project_name = s$project,
    soda_basename = s$basename,
    n_cells = s$n_cells,
    soda_version = s$soda_version,
    schema_version = s$schema_version,
    package_version = s$package_version,
    n_files = length(entries),
    files = entries
  )
  jsonlite::write_json(
    inv,
    file.path(root_n, "inventory.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = NA
  )
  invisible(inv)
}

#' Collect lightweight SODA summary for README / HTML index
#' @keywords internal
.ic2_soda_summary <- function(root, man = NULL, x = NULL) {
  root <- as.character(root)[1L]
  if (is.null(man) && file.exists(file.path(root, "soda.json"))) {
    man <- tryCatch(
      jsonlite::fromJSON(file.path(root, "soda.json"), simplifyVector = TRUE, simplifyDataFrame = FALSE),
      error = function(e) NULL
    )
  }
  proj <- as.character(man$project_name %||% basename(root))[1L]
  ver <- as.character(man$soda_version %||% "1.0.0")[1L]
  schema <- as.character(man$schema_version %||% "1.0.0")[1L]
  pkg <- as.character(man$package_version %||% NA_character_)[1L]
  created <- as.character(man$created %||% NA_character_)[1L]

  n_cells <- NA_character_
  if (!is.null(x) && !is.null(x@cell_data)) {
    n_cells <- as.character(nrow(x@cell_data))
  } else {
    cd <- file.path(root, "metadata", "cell_data.tsv")
    if (file.exists(cd)) {
      n_cells <- tryCatch({
        as.character(max(0L, length(readLines(cd, warn = FALSE)) - 1L))
      }, error = function(e) NA_character_)
    }
  }

  mods <- names(man$modalities %||% list())
  if (length(mods) < 1L && dir.exists(file.path(root, "modalities"))) {
    mods <- list.dirs(file.path(root, "modalities"), full.names = FALSE, recursive = FALSE)
  }

  assay_rows <- list()
  for (mod in mods) {
    assay_root <- file.path(root, "modalities", mod, "assays")
    if (!dir.exists(assay_root)) next
    for (an in list.dirs(assay_root, full.names = FALSE, recursive = FALSE)) {
      adir <- file.path(assay_root, an)
      h5s <- list.files(adir, pattern = "\\.(h5|hdf5)$", recursive = TRUE, ignore.case = TRUE)
      assay_rows[[length(assay_rows) + 1L]] <- list(
        modality = mod,
        assay = an,
        n_h5 = length(h5s)
      )
    }
  }

  red_names <- character(0)
  for (rd in c(file.path(root, "reductions"), file.path(root, "modalities", "rna", "reductions"))) {
    if (dir.exists(rd)) {
      red_names <- c(red_names, list.dirs(rd, full.names = FALSE, recursive = FALSE))
    }
  }
  red_names <- sort(unique(red_names[nzchar(red_names)]))

  images <- character(0)
  img_dir <- file.path(root, "images")
  if (dir.exists(img_dir)) {
    images <- list.files(img_dir, pattern = "\\.(png|html|pdf|svg)$", ignore.case = TRUE, recursive = TRUE)
  }

  links <- lapply(.ic2_soda_doc_catalog(), function(z) {
    list(label = z$label, href = z$file, blurb = z$blurb, file = z$file)
  })
  links <- Filter(function(z) file.exists(file.path(root, z$href)), links)

  list(
    root = root,
    basename = basename(root),
    project = proj,
    soda_version = ver,
    schema_version = schema,
    package_version = pkg,
    created = created,
    n_cells = n_cells,
    modalities = mods,
    assays = assay_rows,
    reductions = red_names,
    images = images,
    links = links,
    generated = as.character(Sys.time())
  )
}

#' Minimal Markdown  ->  HTML for embedding README in index.html
#' @keywords internal
.ic2_md_inline <- function(x) {
  x <- as.character(x %||% "")
  x <- .ic2_html_escape(x)
  ## After escape so '&' in entities is not re-escaped
  x <- .ic2_html_fix_utf8(x)
  x <- gsub("`([^`]+)`", "<code>\\1</code>", x, perl = TRUE)
  x <- gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", x, perl = TRUE)
  x <- gsub("(?<!\\*)\\*([^*]+)\\*(?!\\*)", "<em>\\1</em>", x, perl = TRUE)
  x
}

#' @keywords internal
.ic2_md_to_html <- function(md) {
  if (is.null(md) || !nzchar(paste(md, collapse = ""))) {
    return("<p class=\"muted\">No README.md found.</p>")
  }
  lines <- if (length(md) == 1L) strsplit(md, "\n", fixed = TRUE)[[1L]] else as.character(md)
  out <- character(0)
  i <- 1L
  in_ul <- FALSE
  close_ul <- function() {
    if (in_ul) {
      out <<- c(out, "</ul>")
      in_ul <<- FALSE
    }
  }
  while (i <= length(lines)) {
    ln <- lines[[i]]
    if (startsWith(ln, "```")) {
      close_ul()
      i <- i + 1L
      buf <- character(0)
      while (i <= length(lines) && !startsWith(lines[[i]], "```")) {
        buf <- c(buf, lines[[i]])
        i <- i + 1L
      }
      out <- c(out, paste0("<pre class=\"md-code\"><code>", .ic2_html_escape(paste(buf, collapse = "\n")), "</code></pre>"))
      i <- i + 1L
      next
    }
    if (grepl("^###\\s+", ln)) {
      close_ul()
      out <- c(out, paste0("<h4>", .ic2_md_inline(sub("^###\\s+", "", ln)), "</h4>"))
      i <- i + 1L
      next
    }
    if (grepl("^##\\s+", ln)) {
      close_ul()
      out <- c(out, paste0("<h3>", .ic2_md_inline(sub("^##\\s+", "", ln)), "</h3>"))
      i <- i + 1L
      next
    }
    if (grepl("^#\\s+", ln)) {
      close_ul()
      out <- c(out, paste0("<h2 class=\"md-title\">", .ic2_md_inline(sub("^#\\s+", "", ln)), "</h2>"))
      i <- i + 1L
      next
    }
    if (grepl("^[-*]\\s+", ln)) {
      if (!in_ul) {
        out <- c(out, "<ul class=\"md-list\">")
        in_ul <- TRUE
      }
      out <- c(out, paste0("<li>", .ic2_md_inline(sub("^[-*]\\s+", "", ln)), "</li>"))
      i <- i + 1L
      next
    }
    if (!nzchar(trimws(ln))) {
      close_ul()
      i <- i + 1L
      next
    }
    if (grepl("^_.*_$", trimws(ln))) {
      close_ul()
      out <- c(out, paste0("<p class=\"md-meta\">", .ic2_md_inline(gsub("^_|_$", "", trimws(ln))), "</p>"))
      i <- i + 1L
      next
    }
    close_ul()
    out <- c(out, paste0("<p>", .ic2_md_inline(ln), "</p>"))
    i <- i + 1L
  }
  close_ul()
  paste(out, collapse = "\n")
}

#' Render soda.info key/value twin as an in-page HTML panel
#' @keywords internal
.ic2_soda_info_to_html <- function(path) {
  if (!file.exists(path)) {
    return("<p class=\"muted\">No soda.info found.</p>")
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  rows <- character(0)
  for (ln in lines) {
    if (!nzchar(trimws(ln)) || startsWith(trimws(ln), "#")) next
    if (grepl("^key\\tvalue$", ln, ignore.case = TRUE)) next
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) {
      if (grepl("=", ln, fixed = TRUE)) {
        parts <- strsplit(ln, "=", fixed = TRUE)[[1L]]
        parts <- c(parts[[1L]], paste(parts[-1L], collapse = "="))
      } else next
    }
    k <- trimws(parts[[1L]])
    v <- trimws(paste(parts[-1L], collapse = "\t"))
    rows <- c(rows, sprintf(
      "<div class=\"kv\"><span class=\"k\">%s</span><span class=\"v\">%s</span></div>",
      .ic2_html_escape(k), .ic2_html_escape(v)
    ))
  }
  if (!length(rows)) return("<p class=\"muted\">soda.info is empty.</p>")
  paste0("<div class=\"kv-grid\">", paste(rows, collapse = "\n"), "</div>")
}

#' Capture iCellR2 show() text for SODA index.html
#' @keywords internal
#' HTML-escape console show() text; keep box-drawing / checks as entities
#' @keywords internal
.ic2_soda_console_html <- function(txt) {
  x <- paste(as.character(txt %||% ""), collapse = "\n")
  x <- .ic2_html_escape(x)
  ## C-locale capture stores literal <U+XXXX>; after escape that is &lt;U+XXXX&gt;
  uni <- c(
    "&lt;U+2550&gt;" = "&#9552;",
    "&lt;U+2500&gt;" = "&#9472;",
    "&lt;U+2713&gt;" = "&#10003;",
    "&lt;U+2190&gt;" = "&larr;",
    "&lt;U+00B7&gt;" = "&middot;",
    "&lt;U+2014&gt;" = "&mdash;",
    "&lt;U+2013&gt;" = "&ndash;",
    "<U+2550>" = "&#9552;",
    "<U+2500>" = "&#9472;",
    "<U+2713>" = "&#10003;",
    "<U+2190>" = "&larr;",
    "<U+00B7>" = "&middot;"
  )
  for (nm in names(uni)) x <- gsub(nm, uni[[nm]], x, fixed = TRUE)
  x <- gsub("&lt;e2&gt;&lt;95&gt;&lt;90&gt;", "&#9552;", x, fixed = TRUE)
  x <- gsub("&lt;e2&gt;&lt;94&gt;&lt;80&gt;", "&#9472;", x, fixed = TRUE)
  x <- gsub("&lt;e2&gt;&lt;9c&gt;&lt;93&gt;", "&#10003;", x, fixed = TRUE)
  x <- gsub("<e2><95><90>", "&#9552;", x, fixed = TRUE)
  x <- gsub("<e2><94><80>", "&#9472;", x, fixed = TRUE)
  x <- gsub("<e2><9c><93>", "&#10003;", x, fixed = TRUE)
  x <- gsub("<c2><b7>", "&middot;", x, fixed = TRUE)
  x <- tryCatch({
    x <- gsub(intToUtf8(0x2550L), "&#9552;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2500L), "&#9472;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2713L), "&#10003;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2190L), "&larr;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2192L), "&rarr;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x00B7L), "&middot;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2014L), "&mdash;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2013L), "&ndash;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x2026L), "&hellip;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x00D7L), "&times;", x, fixed = TRUE, useBytes = TRUE)
    x <- gsub(intToUtf8(0x0394L), "&Delta;", x, fixed = TRUE, useBytes = TRUE)
    x
  }, error = function(e) x)
  heavy <- paste(rep("=", 48L), collapse = "")
  light <- paste(rep("-", 48L), collapse = "")
  if (grepl(heavy, x, fixed = TRUE)) {
    x <- gsub(heavy, paste(rep("&#9552;", 48L), collapse = ""), x, fixed = TRUE)
  }
  if (grepl(light, x, fixed = TRUE)) {
    x <- gsub(light, paste(rep("&#9472;", 48L), collapse = ""), x, fixed = TRUE)
  }
  x <- gsub("[x]", "[&#10003;]", x, fixed = TRUE)
  x <- gsub(" &lt;- active", " &larr; active", x, fixed = TRUE)
  x <- gsub(" <- active", " &larr; active", x, fixed = TRUE)
  x <- gsub("  |  ", " &middot; ", x, fixed = TRUE)
  ## Longer ASCII rules (architecture.txt uses width 80)
  while (grepl("====", x, fixed = TRUE)) {
    x <- gsub("====", "&#9552;&#9552;&#9552;&#9552;", x, fixed = TRUE)
  }
  while (grepl("----", x, fixed = TRUE)) {
    x <- gsub("----", "&#9472;&#9472;&#9472;&#9472;", x, fixed = TRUE)
  }
  .ic2_html_fix_utf8(x)
}

.ic2_soda_show_text <- function(root, man = NULL, x = NULL) {
  root <- as.character(root)[1L]
  obj <- x
  if (is.null(obj) && nzchar(root) && dir.exists(root)) {
    obj <- tryCatch(
      open.soda(root, verbose = FALSE),
      error = function(e) NULL
    )
  }
  if (!is.null(obj) && inherits(obj, "iCellR2")) {
    op <- options(iCellR2.show.utf8 = TRUE)
    on.exit(options(op), add = TRUE)
    txt <- tryCatch(
      utils::capture.output(methods::show(obj), type = "message"),
      error = function(e) character(0)
    )
    if (length(txt) >= 1L) {
      return(paste(txt, collapse = "\n"))
    }
  }
  ## Disk fallback (no object / show failed)
  s <- .ic2_soda_summary(root, man = man, x = x)
  lines <- c(
    paste(rep("=", 48L), collapse = ""),
    "  iCellR2",
    "  Massive-scale single-cell analysis",
    "  with a familiar R interface",
    paste(rep("=", 48L), collapse = ""),
    sprintf("  %-11s %s", "Project", s$project),
    sprintf(
      "  %-11s %s   |   schema %s   |   SODA %s",
      "Version", s$package_version %||% "NA", s$schema_version, s$soda_version
    ),
    sprintf("  %-11s %s", "Archive", s$basename),
    "",
    sprintf("  %-11s %s", "Cells", s$n_cells %||% "NA"),
    sprintf(
      "  %-11s %s",
      "Modalities",
      if (length(s$modalities)) paste(s$modalities, collapse = ", ") else "none"
    ),
    paste(rep("-", 48L), collapse = ""),
    "  Assays on disk",
    paste(rep("-", 48L), collapse = "")
  )
  if (length(s$assays) < 1L) {
    lines <- c(lines, "  (none)")
  } else {
    for (a in s$assays) {
      lines <- c(lines, sprintf("  %-22s %s (%s H5)", a$assay, a$modality, a$n_h5))
    }
  }
  lines <- c(
    lines,
    paste(rep("-", 48L), collapse = ""),
    "  Reductions",
    paste(rep("-", 48L), collapse = "")
  )
  if (length(s$reductions) < 1L) {
    lines <- c(lines, "  (none)")
  } else {
    for (r in s$reductions) lines <- c(lines, paste0("  ", r))
  }
  lines <- c(lines, paste(rep("=", 48L), collapse = ""))
  paste(lines, collapse = "\n")
}

#' Render plain text docs (architecture) for in-page HTML with UTF-8 trees
#' @keywords internal
.ic2_soda_plain_to_html <- function(path) {
  if (!file.exists(path)) return("<p class=\"muted\">File not found.</p>")
  raw <- tryCatch(
    readBin(path, what = "raw", n = as.integer(file.info(path)$size %||% 0)),
    error = function(e) raw(0)
  )
  txt <- if (length(raw)) {
    tryCatch(rawToChar(raw), error = function(e) {
      paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    })
  } else {
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }
  Encoding(txt) <- "UTF-8"
  body <- .ic2_soda_console_html(txt)
  extra <- c(
    "&lt;U+2502&gt;" = "&#9474;",
    "&lt;U+251C&gt;" = "&#9500;",
    "&lt;U+2514&gt;" = "&#9492;",
    "&lt;U+252C&gt;" = "&#9516;",
    "&lt;U+2534&gt;" = "&#9524;",
    "&lt;U+2524&gt;" = "&#9508;",
    "&lt;U+253C&gt;" = "&#9532;",
    "&lt;U+2192&gt;" = "&rarr;",
    "&lt;U+2026&gt;" = "&hellip;",
    "&lt;U+2013&gt;" = "&ndash;",
    "<U+2502>" = "&#9474;",
    "<U+251C>" = "&#9500;",
    "<U+2514>" = "&#9492;",
    "<U+2192>" = "&rarr;",
    "<U+2026>" = "&hellip;",
    "<U+2013>" = "&ndash;",
    "<U+2014>" = "&mdash;",
    "<U+00D7>" = "&times;",
    "<U+0394>" = "&Delta;",
    "&lt;U+2014&gt;" = "&mdash;",
    "&lt;U+00D7&gt;" = "&times;",
    "&lt;U+0394&gt;" = "&Delta;"
  )
  for (nm in names(extra)) body <- gsub(nm, extra[[nm]], body, fixed = TRUE)
  body <- tryCatch({
    body <- gsub(intToUtf8(0x2502L), "&#9474;", body, fixed = TRUE, useBytes = TRUE)
    body <- gsub(intToUtf8(0x251CL), "&#9500;", body, fixed = TRUE, useBytes = TRUE)
    body <- gsub(intToUtf8(0x2514L), "&#9492;", body, fixed = TRUE, useBytes = TRUE)
    body <- gsub(intToUtf8(0x2192L), "&rarr;", body, fixed = TRUE, useBytes = TRUE)
    body <- gsub(intToUtf8(0x2026L), "&hellip;", body, fixed = TRUE, useBytes = TRUE)
    body <- gsub(intToUtf8(0x2013L), "&ndash;", body, fixed = TRUE, useBytes = TRUE)
    body <- gsub(intToUtf8(0x2014L), "&mdash;", body, fixed = TRUE, useBytes = TRUE)
    body <- gsub(intToUtf8(0x00D7L), "&times;", body, fixed = TRUE, useBytes = TRUE)
    body <- gsub(intToUtf8(0x0394L), "&Delta;", body, fixed = TRUE, useBytes = TRUE)
    body
  }, error = function(e) body)
  paste0("<pre class=\"doc-pre\">", body, "</pre>")
}

#' Map a SODA root document to its UTF-8 HTML viewer path (relative)
#' @keywords internal
.ic2_soda_doc_view_href <- function(file) {
  paste0("#", .ic2_soda_doc_view_id(file))
}

#' Hash id for a SODA root document inside index.html
#' @keywords internal
.ic2_soda_doc_view_id <- function(file) {
  file <- as.character(file)[1L]
  map <- c(
    "soda.info" = "soda-info",
    "README.md" = "readme",
    "SODA_SPEC.md" = "spec",
    "architecture.txt" = "architecture",
    "iCellR2.log" = "log",
    "benchmark.txt" = "benchmark",
    "MANIFEST.txt" = "manifest"
  )
  if (file %in% names(map)) map[[file]] else gsub("[^A-Za-z0-9]+", "-", file)
}

#' Shared CSS for SODA index + document viewer pages
#' @keywords internal
.ic2_soda_html_css <- function(viewer = FALSE) {
  pre_max <- if (isTRUE(viewer)) "none" else "28rem"
  paste0(
    ":root{--nyu:#57068c;--nyu-deep:#3b0764;--nyu-ink:#2a0448;--ink:#1a1a2e;--muted:#64748b;--line:#e2e8f0;--blue:#2563eb;--teal:#0d9488;--amber:#d97706;}\n",
    "*{box-sizing:border-box}\n",
    "html{scroll-behavior:smooth;}\n",
    "body{margin:0;font-family:'DM Sans',system-ui,sans-serif;background:linear-gradient(180deg,#f8fbff 0%,#ffffff 140px);color:var(--ink);line-height:1.55;}\n",
    ".hero{background:linear-gradient(105deg,var(--nyu-ink) 0%,var(--nyu-deep) 48%,var(--nyu) 100%);color:#fff;padding:2.4rem 1.25rem 2rem;border-bottom:1px solid var(--nyu-ink);}\n",
    ".hero-inner{max-width:980px;margin:0 auto;}\n",
    ".hero.compact{padding:1.35rem 1.25rem 1.2rem;}\n",
    ".eyebrow{text-transform:uppercase;letter-spacing:.1em;font-size:.68rem;color:#e9d5ff;font-weight:600;margin:0 0 .55rem;opacity:.95;}\n",
    "h1{font-family:'DM Sans',system-ui,-apple-system,sans-serif;font-weight:600;font-size:clamp(1.5rem,2.6vw,1.95rem);margin:0 0 .5rem;letter-spacing:-.015em;line-height:1.25;}\n",
    ".hero.compact h1{font-size:clamp(1.2rem,2vw,1.45rem);margin:0 0 .25rem;}\n",
    ".lede{color:#f3e8ff;max-width:40rem;margin:0;font-size:.9rem;line-height:1.55;font-weight:400;opacity:.95;}\n",
    ".hero code{font-size:.86em;background:rgba(255,255,255,.12);color:#fff;padding:.1rem .35rem;border-radius:4px;}\n",
    ".stats{display:flex;flex-wrap:wrap;gap:.55rem;margin:1.15rem 0 0;}\n",
    ".stat{background:rgba(255,255,255,.1);border:1px solid rgba(255,255,255,.18);border-radius:999px;padding:.32rem .7rem;font-size:.78rem;color:#f3e8ff;}\n",
    ".stat b{color:#fff;font-weight:700;}\n",
    ".jump{position:sticky;top:0;z-index:20;background:rgba(248,251,255,.92);backdrop-filter:blur(8px);border-bottom:1px solid var(--line);}\n",
    ".jump-inner{max-width:980px;margin:0 auto;padding:.55rem 1.25rem;display:flex;flex-wrap:wrap;gap:.35rem .9rem;font-size:.8rem;}\n",
    ".jump a{color:#475569;text-decoration:none;font-weight:500;}\n",
    ".jump a:hover{color:var(--nyu);}\n",
    "main{max-width:980px;margin:0 auto;padding:1.35rem 1.25rem 3.5rem;}\n",
    "main.wide{max-width:1100px;}\n",
    "h2{font-family:'DM Sans',system-ui,sans-serif;font-size:1.05rem;margin:0 0 .55rem;font-weight:600;}\n",
    ".grid{display:grid;gap:1rem;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));margin:1rem 0 1.25rem;}\n",
    ".card{background:#fff;border:1px solid var(--line);border-radius:12px;padding:1rem 1.1rem;margin-top:1rem;}\n",
    ".card.accent{border-left:4px solid var(--blue);background:linear-gradient(90deg,#eff6ff 0%,#ffffff 32%);margin-top:0;}\n",
    ".card.teal{border-left:4px solid var(--teal);background:linear-gradient(90deg,#f0fdfa 0%,#ffffff 32%);}\n",
    ".card.amber{border-left:4px solid var(--amber);background:linear-gradient(90deg,#fffbeb 0%,#ffffff 32%);}\n",
    "ul{margin:.2rem 0 0;padding-left:1.1rem;}\n",
    "li{margin:.28rem 0;}\n",
    "a{color:var(--blue);text-decoration:none;}\n",
    "a:hover{text-decoration:underline;}\n",
    ".muted{color:var(--muted);}\n",
    ".tiny{font-size:.78rem;line-height:1.4;margin-top:.15rem;}\n",
    ".file{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.78rem;color:#94a3b8;margin-left:.35rem;}\n",
    "pre{background:#0f172a;color:#e2e8f0;border-radius:10px;padding:1rem 1.1rem;overflow:auto;font-size:.84rem;}\n",
    "code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}\n",
    "table{width:100%;border-collapse:collapse;font-size:.86rem;}\n",
    "th,td{text-align:left;padding:.55rem .45rem;border-bottom:1px solid var(--line);vertical-align:top;}\n",
    "th{font-size:.72rem;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);font-weight:600;}\n",
    ".about p{margin:0 0 .75rem;color:#334155;line-height:1.55;}\n",
    ".about .why-title,.about .get-title{font-size:.68rem;font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:var(--nyu);margin:1rem 0 .45rem;}\n",
    ".about ul{margin:.2rem 0 0;padding-left:1.15rem;color:#334155;}\n",
    ".about li{margin:.28rem 0;line-height:1.45;}\n",
    ".md-body{font-size:.92rem;}\n",
    ".md-body h2.md-title{font-family:'DM Sans',system-ui,sans-serif;font-size:1.25rem;font-weight:600;color:var(--nyu);margin:0 0 .75rem;}\n",
    ".md-body h3{font-family:'DM Sans',system-ui,sans-serif;font-size:1.02rem;color:var(--nyu-deep);margin:1.15rem 0 .45rem;padding-bottom:.25rem;border-bottom:1px solid #ede9fe;}\n",
    ".md-body h4{font-size:.92rem;color:#0f172a;margin:.9rem 0 .35rem;}\n",
    ".md-body p{margin:.4rem 0;color:#334155;}\n",
    ".md-body .md-list{margin:.35rem 0 .6rem;padding-left:1.15rem;}\n",
    ".md-body .md-list li{margin:.22rem 0;color:#334155;}\n",
    ".md-body .md-code{background:#0f172a;color:#e2e8f0;border-radius:10px;padding:1rem 1.1rem;overflow:auto;font-size:.84rem;margin:.65rem 0;}\n",
    ".md-body code{background:#f1f5f9;color:#0b3d91;padding:.08rem .35rem;border-radius:4px;font-size:.86em;}\n",
    ".md-body .md-code code{background:transparent;color:inherit;padding:0;}\n",
    ".md-body .md-meta{font-size:.8rem;color:var(--muted);font-style:italic;}\n",
    ".md-body strong{color:var(--nyu-ink);}\n",
    ".kv-grid{display:grid;gap:.45rem;margin-top:.35rem;}\n",
    ".kv{display:grid;grid-template-columns:minmax(120px,0.35fr) 1fr;gap:.75rem;align-items:baseline;padding:.45rem .65rem;border-radius:8px;background:#f8fafc;border:1px solid var(--line);}\n",
    ".kv .k{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.78rem;color:var(--nyu);font-weight:600;}\n",
    ".kv .v{font-size:.9rem;color:#0f172a;word-break:break-word;}\n",
    ".panel-head{display:flex;align-items:baseline;justify-content:space-between;gap:1rem;flex-wrap:wrap;margin-bottom:.35rem;}\n",
    ".panel-head .file{margin-left:0;}\n",
    ".doc-pre{background:#0f172a;color:#e2e8f0;border-radius:12px;padding:1rem 1.15rem;overflow:auto;font-size:.78rem;line-height:1.45;margin:.35rem 0 0;white-space:pre;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;border:1px solid #1e293b;max-height:", pre_max, ";}\n",
    "details.fold{border:1px solid var(--line);border-radius:12px;padding:.75rem 1rem;background:#fff;}\n",
    "details.fold > summary{cursor:pointer;font-weight:600;color:#0f172a;list-style:none;}\n",
    "details.fold > summary::-webkit-details-marker{display:none;}\n",
    "details.fold > summary::before{content:\"\\25B8 \";color:var(--nyu);}\n",
    "details.fold[open] > summary::before{content:\"\\25BE \";}\n",
    "details.fold .md-body,details.fold .doc-pre,details.fold table{margin-top:.75rem;}\n",
    ".code-wrap{position:relative;}\n",
    ".copy-btn{position:absolute;top:.55rem;right:.55rem;border:1px solid #334155;background:#1e293b;color:#e2e8f0;border-radius:6px;padding:.25rem .55rem;font-size:.72rem;cursor:pointer;font-family:'DM Sans',sans-serif;}\n",
    ".copy-btn:hover{background:#334155;}\n",
    ".back{display:inline-block;margin:0 0 .85rem;font-size:.84rem;font-weight:500;}\n",
    ".doc-view{display:none;}\n",
    "body.mode-doc .home-view{display:none;}\n",
    "body.mode-doc .doc-view{display:block;}\n",
    "body.mode-doc .jump{display:none;}\n",
    "body.mode-doc .hero.home-hero{display:none;}\n",
    "body.mode-doc .doc-pre{max-height:none;}\n",
    ".doc-nav{display:flex;flex-wrap:wrap;gap:.35rem .85rem;margin:0 0 1rem;font-size:.82rem;}\n",
    ".doc-nav a{color:#475569;font-weight:500;text-decoration:none;}\n",
    ".doc-nav a:hover,.doc-nav a.active{color:var(--nyu);}\n",
    ".doc-page{display:none;}\n",
    ".doc-page.active{display:block;}\n",
    "footer{margin-top:1.5rem;color:var(--muted);font-size:.8rem;}\n"
  )
}

#' Write index.html into a SODA root (static archive home page)
#' @keywords internal
.ic2_soda_write_html_index <- function(root, man = NULL, x = NULL, inventory = NULL) {
  root <- as.character(root)[1L]
  if (!nzchar(root) || !dir.exists(root)) return(invisible(FALSE))
  s <- .ic2_soda_summary(root, man = man, x = x)
  esc <- .ic2_html_escape

  if (is.null(inventory) && file.exists(file.path(root, "inventory.json"))) {
    inventory <- tryCatch(
      jsonlite::fromJSON(file.path(root, "inventory.json"), simplifyVector = FALSE),
      error = function(e) NULL
    )
  }

  ## Docs with dedicated viewer pages / home chrome (omit from "other files" table)
  panel_docs <- c(
    "soda.json", "soda.info", "README.md", "SODA_SPEC.md", "architecture.txt",
    "iCellR2.log", "benchmark.txt",
    "index.html", "docs.html", "inventory.json", "MANIFEST.txt"
  )

  link_items <- vapply(s$links, function(z) {
    f <- z$file %||% z$href
    if (identical(f, "index.html")) return(NA_character_)
    href <- .ic2_soda_doc_view_href(f)
    sprintf(
      "<li><a href=\"%s\"><strong>%s</strong></a> <span class=\"file\">%s</span><div class=\"muted tiny\">%s</div></li>",
      esc(href), esc(z$label), esc(f), esc(z$blurb %||% "")
    )
  }, character(1))
  link_li <- paste(link_items[!is.na(link_items)], collapse = "\n")

  has_images <- length(s$images) >= 1L
  img_section <- ""
  if (has_images) {
    show <- utils::head(s$images, 24L)
    img_li <- paste(vapply(show, function(f) {
      sprintf("<li><a href=\"images/%s\">%s</a></li>", esc(f), esc(f))
    }, character(1)), collapse = "\n")
    img_section <- paste0(
      "<section id=\"panel-images\" class=\"card amber\"><h2>Images / plots</h2><ul>\n",
      img_li, "\n</ul></section>\n"
    )
  }

  inv_files <- list()
  if (is.list(inventory) && length(inventory$files)) inv_files <- inventory$files

  h5_rows <- character(0)
  other_rows <- character(0)
  cache_rows <- character(0)
  for (e in inv_files) {
    path <- as.character(e$path %||% "")[1L]
    if (path %in% panel_docs || startsWith(path, "view/") || (grepl("\\.html$", path) && !identical(path, "index.html"))) next
    dims <- ""
    if (!is.null(e$h5) && length(e$h5$dims)) {
      dims <- paste(unlist(e$h5$dims), collapse = "  x  ")
    } else if (!is.null(e$table)) {
      rr <- e$table$rows; cc <- e$table$cols
      rr_ok <- !is.null(rr) && length(rr) && !is.na(rr[1L]) && !identical(as.character(rr[1L]), "NA")
      cc_ok <- !is.null(cc) && length(cc) && !is.na(cc[1L]) && !identical(as.character(cc[1L]), "NA")
      if (rr_ok) {
        dims <- sprintf("%s rows  x  %s cols", rr[1L], if (cc_ok) cc[1L] else "?")
      } else if (cc_ok) {
        dims <- sprintf("? rows  x  %s cols", cc[1L])
      } else if (length(e$table$keys)) {
        dims <- paste("keys:", paste(utils::head(unlist(e$table$keys), 8L), collapse = ", "))
      }
    }
    extra <- ""
    if (!is.null(e$h5)) {
      extra <- sprintf(
        "%s &middot; %s%s",
        esc(e$h5$format %||% e$h5$kind %||% "hdf5"),
        if (nzchar(dims)) paste0(esc(dims), " &middot; ") else "",
        if (!is.null(e$h5$nnz) && !is.na(e$h5$nnz)) paste0("nnz ", esc(e$h5$nnz)) else ""
      )
      extra <- sub(" &middot; $", "", extra)
    } else if (nzchar(dims)) {
      extra <- esc(dims)
    }
    size <- if (!is.null(e$size_bytes) && !is.na(e$size_bytes)) {
      sb <- as.numeric(e$size_bytes)
      if (sb >= 1e9) sprintf("%.1f GB", sb / 1e9)
      else if (sb >= 1e6) sprintf("%.1f MB", sb / 1e6)
      else if (sb >= 1e3) sprintf("%.1f KB", sb / 1e3)
      else paste0(sb, " B")
    } else "-"
    row <- sprintf(
      "<tr><td><strong>%s</strong><div class=\"muted tiny\">%s</div></td><td><code>%s</code></td><td>%s</td><td>%s</td></tr>",
      esc(e$label %||% path), esc(e$blurb %||% ""), esc(path),
      if (nzchar(extra)) extra else "&mdash;", esc(size)
    )
    if (!is.null(e$h5)) {
      ## Keep assay/reduction matrices in main H5 table; stash cache H5 under cache
      if (grepl("^cache/", path)) cache_rows <- c(cache_rows, row) else h5_rows <- c(h5_rows, row)
    } else if (grepl("^cache/", path)) {
      cache_rows <- c(cache_rows, row)
    } else {
      other_rows <- c(other_rows, row)
    }
  }
  if (!length(h5_rows)) {
    h5_rows <- "<tr><td colspan=\"4\" class=\"muted\">No HDF5 matrices inventoried yet. Run <code>refresh.soda()</code>.</td></tr>"
  }
  if (!length(other_rows)) {
    other_rows <- "<tr><td colspan=\"4\" class=\"muted\">No additional files.</td></tr>"
  }

  ## Active modalities = those with discovered assays
  active_mods <- unique(vapply(s$assays, function(a) a$modality, character(1)))
  if (!length(active_mods)) {
    active_mods <- setdiff(s$modalities, c("misc_modalities", ""))
    ## Prefer non-empty assay dirs only
    active_mods <- active_mods[vapply(active_mods, function(m) {
      ad <- file.path(root, "modalities", m, "assays")
      dir.exists(ad) && length(list.dirs(ad, recursive = FALSE)) > 0L
    }, logical(1))]
  }
  mod_txt <- if (length(active_mods)) paste(active_mods, collapse = ", ") else "&mdash;"

  open_snip <- sprintf(
    "library(iCellR2)\nobj <- open.soda(\"%s\")\nvalidate.soda(\"%s\")\nrefresh.soda(\"%s\")",
    s$basename, s$basename, s$basename
  )

  ## UTF-8 document views live inside index.html (#log, #readme, ...)
  doc_specs <- list(
    list(file = "soda.info", title = "Project metadata (text)",
         body = function() .ic2_soda_info_to_html(file.path(root, "soda.info"))),
    list(file = "README.md", title = "README",
         body = function() paste0("<div class=\"md-body\">\n", .ic2_md_to_html(readLines(file.path(root, "README.md"), warn = FALSE, encoding = "UTF-8")), "\n</div>")),
    list(file = "SODA_SPEC.md", title = "Format specification",
         body = function() paste0("<div class=\"md-body\">\n", .ic2_md_to_html(readLines(file.path(root, "SODA_SPEC.md"), warn = FALSE, encoding = "UTF-8")), "\n</div>")),
    list(file = "architecture.txt", title = "Folder layout",
         body = function() .ic2_soda_plain_to_html(file.path(root, "architecture.txt"))),
    list(file = "iCellR2.log", title = "Run log",
         body = function() .ic2_soda_plain_to_html(file.path(root, "iCellR2.log"))),
    list(file = "benchmark.txt", title = "Timing / benchmark",
         body = function() .ic2_soda_plain_to_html(file.path(root, "benchmark.txt"))),
    list(file = "MANIFEST.txt", title = "File inventory (light)",
         body = function() .ic2_soda_plain_to_html(file.path(root, "MANIFEST.txt")))
  )
  docs <- list()
  for (spec in doc_specs) {
    if (!file.exists(file.path(root, spec$file))) next
    docs[[length(docs) + 1L]] <- list(
      id = .ic2_soda_doc_view_id(spec$file),
      title = spec$title,
      source = spec$file,
      body = spec$body()
    )
  }
  doc_ids_js <- paste(vapply(docs, function(d) sprintf("\"%s\":true", d$id), character(1)), collapse = ",")
  doc_nav_html <- paste(vapply(docs, function(d) {
    sprintf("<a href=\"#%s\">%s</a>", esc(d$id), esc(d$title))
  }, character(1)), collapse = "\n")
  doc_pages_html <- paste(vapply(docs, function(d) {
    paste0(
      "<article class=\"doc-page\" id=\"", esc(d$id), "\">\n",
      "<div class=\"panel-head\"><h2>", esc(d$title), "</h2>",
      "<span class=\"file\">", esc(d$source), "</span></div>\n",
      d$body, "\n</article>\n"
    )
  }, character(1)), collapse = "\n")

  cache_block <- ""
  if (length(cache_rows)) {
    cache_block <- paste0(
      "<details class=\"fold\" style=\"margin-top:1rem\">\n",
      "<summary>Cache files <span class=\"muted tiny\">(", length(cache_rows), ")</span></summary>\n",
      "<table><thead><tr><th>What</th><th>Path</th><th>Summary</th><th>Size</th></tr></thead><tbody>\n",
      paste(cache_rows, collapse = "\n"),
      "\n</tbody></table>\n</details>\n"
    )
  }

  nav_imgs <- if (has_images) "<a href=\"#panel-images\">Images</a>" else ""
  docs_grid <- if (has_images) {
    paste0(
      "<div class=\"grid\">\n",
      "<section id=\"panel-docs\" class=\"card teal\"><h2>Archive documents</h2><ul>\n", link_li, "\n</ul></section>\n",
      img_section,
      "</div>\n"
    )
  } else {
    paste0(
      "<section id=\"panel-docs\" class=\"card teal\">\n",
      "<h2>Archive documents</h2><ul>\n", link_li, "\n</ul>\n</section>\n"
    )
  }

  html <- paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n",
    "<meta charset=\"utf-8\"/>\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>\n",
    "<title>", esc(s$project), " &mdash; SODA</title>\n",
    "<link rel=\"preconnect\" href=\"https://fonts.googleapis.com\"/>\n",
    "<link rel=\"preconnect\" href=\"https://fonts.gstatic.com\" crossorigin/>\n",
    "<link href=\"https://fonts.googleapis.com/css2?family=DM+Sans:ital,wght@0,400;0,500;0,600;0,700;1,400&display=swap\" rel=\"stylesheet\"/>\n",
    "<style>\n", .ic2_soda_html_css(viewer = FALSE), "</style>\n</head>\n<body class=\"mode-home\">\n",
    "<header class=\"hero home-hero\"><div class=\"hero-inner\">\n",
    "<p class=\"eyebrow\">SODA &middot; Single Omics Data Archive</p>\n",
    "<h1>", esc(s$project), "</h1>\n",
    "<p class=\"lede\">Portable on-disk single-cell project. Matrices stay in HDF5; open with <code>open.soda()</code>.</p>\n",
    "<div class=\"stats\">\n",
    "<div class=\"stat\">Cells <b>", esc(s$n_cells), "</b></div>\n",
    "<div class=\"stat\">Schema <b>", esc(s$schema_version), "</b></div>\n",
    "<div class=\"stat\">SODA <b>", esc(s$soda_version), "</b></div>\n",
    "<div class=\"stat\">iCellR2 <b>", esc(s$package_version), "</b></div>\n",
    "<div class=\"stat\">Modalities <b>", esc(mod_txt), "</b></div>\n",
    "</div></div></header>\n",
    "<nav class=\"jump\" aria-label=\"Section navigation\"><div class=\"jump-inner\">\n",
    "<a href=\"#panel-docs\">Documents</a>\n",
    if (has_images) "<a href=\"#panel-images\">Images</a>\n" else "",
    "<a href=\"#panel-h5\">Matrices</a>\n",
    "<a href=\"#panel-files\">Files</a>\n",
    "<a href=\"#panel-open\">Open in R</a>\n",
    "</div></nav>\n",
    "<div class=\"home-view\"><main>\n",
    "<section class=\"card accent about\">\n",
    "<p><strong>SODA</strong> (single omics data archive) is a folder of portable files (JSON &middot; HDF5 &middot; TSV &middot; Parquet): assays, reductions, metadata, logs, and benchmarks. Like BAM for alignments, it is a tool-neutral layout, not a language-private object.</p>\n",
    "<p class=\"why-title\">Why it matters</p>\n",
    "<ul>\n",
    "<li><code>.rds</code> / AnnData <code>h5ad</code> lock you into one ecosystem and often into RAM.</li>\n",
    "<li>Mid-pipeline handoff should be a path, not an object conversion.</li>\n",
    "<li>Progress, command history, and timings belong with the data so runs are auditable.</li>\n",
    "</ul>\n",
    "<p class=\"get-title\">What you get</p>\n",
    "<ul>\n",
    "<li>Start in iCellR2, continue in iCellPy (or reverse) on the same <code>*.soda/</code>.</li>\n",
    "<li>Browse UMAP / clusters in the Java SODA viewer without re-exporting.</li>\n",
    "<li><code>iCellR2.log</code>, <code>benchmark.txt</code>, and session summaries travel with the archive.</li>\n",
    "</ul>\n",
    "</section>\n",
    docs_grid,
    "<section id=\"panel-h5\" class=\"card\">\n",
    "<h2>HDF5 matrices &amp; embeddings</h2>\n",
    "<p class=\"muted tiny\">Structure peeked on refresh (shape / nnz / datasets).</p>\n",
    "<table><thead><tr><th>What</th><th>Path</th><th>Structure</th><th>Size</th></tr></thead><tbody>\n",
    paste(h5_rows, collapse = "\n"),
    "\n</tbody></table>\n</section>\n",
    "<section id=\"panel-files\" class=\"card\">\n",
    "<h2>Other archive files</h2>\n",
    "<table><thead><tr><th>What</th><th>Path</th><th>Summary</th><th>Size</th></tr></thead><tbody>\n",
    paste(other_rows, collapse = "\n"),
    "\n</tbody></table>\n",
    cache_block,
    "</section>\n",
    "<section id=\"panel-open\" class=\"card\">\n",
    "<h2>Open in R</h2>\n",
    "<div class=\"code-wrap\"><button type=\"button\" class=\"copy-btn\" data-copy>Copy</button>\n",
    "<pre id=\"open-snippet\"><code>", esc(open_snip), "</code></pre></div>\n</section>\n",
    "<footer>Generated ", esc(s$generated),
    if (nzchar(s$created) && !is.na(s$created)) paste0(" &middot; archive created ", esc(s$created)) else "",
    " &middot; working form is this folder (zip only for transfer)</footer>\n",
    "</main></div>\n",
    "<div class=\"doc-view\">\n",
    "<header class=\"hero compact\"><div class=\"hero-inner\">\n",
    "<p class=\"eyebrow\">SODA document</p>\n",
    "<h1 id=\"doc-hero-title\">Archive document</h1>\n",
    "<p class=\"lede\" id=\"doc-hero-lede\">UTF-8 view</p>\n",
    "</div></header>\n",
    "<main class=\"wide\">\n",
    "<a class=\"back\" href=\"#\">&larr; Archive home</a>\n",
    "<nav class=\"doc-nav\" aria-label=\"Documents\">\n", doc_nav_html, "\n</nav>\n",
    "<section class=\"card\">\n", doc_pages_html, "</section>\n",
    "</main></div>\n",
    "<script>\n",
    "(function(){\n",
    "  var docs={", doc_ids_js, "};\n",
    "  var pages=[].slice.call(document.querySelectorAll('.doc-page'));\n",
    "  var links=[].slice.call(document.querySelectorAll('.doc-nav a'));\n",
    "  var homeTitle=\"", esc(s$project), " &mdash; SODA\";\n",
    "  function showDoc(id){\n",
    "    pages.forEach(function(p){ p.classList.toggle('active', p.id===id); });\n",
    "    links.forEach(function(a){ a.classList.toggle('active', a.getAttribute('href')==='#'+id); });\n",
    "    var active=document.getElementById(id);\n",
    "    if(!active) return;\n",
    "    var h=active.querySelector('h2');\n",
    "    var f=active.querySelector('.file');\n",
    "    var ht=document.getElementById('doc-hero-title');\n",
    "    var hl=document.getElementById('doc-hero-lede');\n",
    "    if(ht) ht.textContent=h?h.textContent:'Document';\n",
    "    if(hl) hl.innerHTML=f?'<code>'+f.textContent+'</code>':'';\n",
    "    document.title=(h?h.textContent:'Document')+' - SODA';\n",
    "    window.scrollTo(0,0);\n",
    "  }\n",
    "  function route(){\n",
    "    var id=(location.hash||'').replace(/^#/,'');\n",
    "    if(id && docs[id]){\n",
    "      document.body.className='mode-doc';\n",
    "      showDoc(id);\n",
    "    } else {\n",
    "      document.body.className='mode-home';\n",
    "      document.title=homeTitle;\n",
    "      pages.forEach(function(p){ p.classList.remove('active'); });\n",
    "      links.forEach(function(a){ a.classList.remove('active'); });\n",
    "    }\n",
    "  }\n",
    "  window.addEventListener('hashchange', route);\n",
    "  route();\n",
    "  var btn=document.querySelector('[data-copy]');\n",
    "  if(btn){\n",
    "    btn.addEventListener('click', function(){\n",
    "      var el=document.getElementById('open-snippet');\n",
    "      var text=el ? el.innerText : '';\n",
    "      function done(){ btn.textContent='Copied'; setTimeout(function(){ btn.textContent='Copy'; }, 1200); }\n",
    "      if(navigator.clipboard && navigator.clipboard.writeText){\n",
    "        navigator.clipboard.writeText(text).then(done).catch(function(){\n",
    "          var ta=document.createElement('textarea'); ta.value=text; document.body.appendChild(ta); ta.select(); document.execCommand('copy'); document.body.removeChild(ta); done();\n",
    "        });\n",
    "      } else {\n",
    "        var ta=document.createElement('textarea'); ta.value=text; document.body.appendChild(ta); ta.select(); document.execCommand('copy'); document.body.removeChild(ta); done();\n",
    "      }\n",
    "    });\n",
    "  }\n",
    "})();\n",
    "</script>\n",
    "</body>\n</html>\n"
  )

  ## Drop legacy multi-file / docs.html chrome
  for (f in c("docs.html", "soda.info.html", "README.html", "SODA_SPEC.html",
              "architecture.html", "iCellR2.log.html", "benchmark.html", "MANIFEST.html")) {
    fp <- file.path(root, f)
    if (file.exists(fp)) try(unlink(fp), silent = TRUE)
  }
  if (dir.exists(file.path(root, "view"))) try(unlink(file.path(root, "view"), recursive = TRUE), silent = TRUE)

  html <- tryCatch({
    html <- gsub(intToUtf8(0x2014L), "&mdash;", html, fixed = TRUE, useBytes = TRUE)
    html <- gsub(intToUtf8(0x2013L), "&ndash;", html, fixed = TRUE, useBytes = TRUE)
    html <- gsub(intToUtf8(0x00B7L), "&middot;", html, fixed = TRUE, useBytes = TRUE)
    html <- gsub(intToUtf8(0x00D7L), "&times;", html, fixed = TRUE, useBytes = TRUE)
    html
  }, error = function(e) html)
  html <- .ic2_html_fix_utf8(html)
  writeLines(enc2utf8(html), file.path(root, "index.html"), useBytes = TRUE)
  invisible(TRUE)
}


#' Write README.md into a SODA root
#' @keywords internal
.ic2_soda_write_readme <- function(root, man = NULL, x = NULL) {
  root <- as.character(root)[1L]
  if (!nzchar(root) || !dir.exists(root)) return(invisible(FALSE))
  if (is.null(man) && file.exists(file.path(root, "soda.json"))) {
    man <- tryCatch(
      jsonlite::fromJSON(file.path(root, "soda.json"), simplifyVector = TRUE, simplifyDataFrame = FALSE),
      error = function(e) NULL
    )
  }
  proj <- man$project_name %||% basename(root)
  ver <- man$soda_version %||% "1.0.0"
  schema <- man$schema_version %||% "1.0.0"
  pkg <- man$package_version %||% NA_character_
  n_cells <- NA_character_
  if (!is.null(x) && !is.null(x@cell_data)) {
    n_cells <- as.character(nrow(x@cell_data))
  } else {
    cd <- file.path(root, "metadata", "cell_data.tsv")
    if (file.exists(cd)) {
      n_cells <- tryCatch({
        as.character(max(0L, length(readLines(cd, warn = FALSE)) - 1L))
      }, error = function(e) NA_character_)
    }
  }
  mods <- names(man$modalities %||% list())
  if (length(mods) < 1L && dir.exists(file.path(root, "modalities"))) {
    mods <- list.dirs(file.path(root, "modalities"), full.names = FALSE, recursive = FALSE)
  }
  lines <- c(
    paste0("# ", proj),
    "",
    "This directory is a **SODA** (Single Omics Data Archive)  -  a portable,",
    "on-disk single-cell project layout.",
    "",
    "## Quick open (R)",
    "",
    "```r",
    "library(iCellR2)",
    paste0("obj <- open.soda(\"", basename(root), "\")"),
    "validate.soda(\"", basename(root), "\")  # optional",
    "```",
    "",
    "## Layout",
    "",
    "- `index.html`  -  browse this archive in a browser",
    "- `soda.json`  -  manifest (schema + relative paths)",
    "- `metadata/`  -  cell/feature tables (TSV/JSON)",
    "- `modalities/<mod>/assays/`  -  sparse/dense HDF5 assays",
    "- `reductions/` or `modalities/<mod>/reductions/`  -  PCA/UMAP/...",
    "- `logs/`  -  run logs",
    "- `MANIFEST.txt`  -  core file listing",
    "- `inventory.json`  -  per-file roles, dimensions, HDF5 structure",
    "- `architecture.txt`  -  design notes",
    "",
    "## Versions",
    "",
    paste0("- SODA format: ", ver),
    paste0("- Schema: ", schema),
    paste0("- Written by iCellR2: ", pkg),
    paste0("- Cells (approx): ", n_cells),
    paste0("- Modalities: ", if (length(mods)) paste(mods, collapse = ", ") else "(none yet)"),
    "",
    "## Spec",
    "",
    "See the package file `SODA_SPEC.md` (or online iCellR2 docs) for the",
    "stable layout contract. Working form is a **folder**; zip only for transfer.",
    "",
    paste0("_Generated ", as.character(Sys.time()), "_"),
    ""
  )
  writeLines(lines, file.path(root, "README.md"), useBytes = TRUE)
  ## Also drop a copy of the spec when available
  spec_src <- system.file("soda", "SODA_SPEC.md", package = "iCellR2")
  if (!nzchar(spec_src) || !file.exists(spec_src)) {
    cand <- file.path(root, "..", "..", "iCellR2", "inst", "soda", "SODA_SPEC.md")
    if (file.exists(cand)) spec_src <- cand
  }
  if (nzchar(spec_src) && file.exists(spec_src)) {
    file.copy(spec_src, file.path(root, "SODA_SPEC.md"), overwrite = TRUE)
  }
  invisible(TRUE)
}

#' Write a light MANIFEST.txt of core SODA files
#' @keywords internal
.ic2_soda_write_file_manifest <- function(root) {
  root <- as.character(root)[1L]
  if (!dir.exists(root)) return(invisible(FALSE))
  patterns <- c(
    "soda.json", "soda.info", "index.html", "inventory.json", "README.md", "SODA_SPEC.md", "architecture.txt",
    "metadata/*", "modalities/*/assays/*/assay.json",
    "modalities/*/assays/*/*.h5", "reductions/*/*"
  )
  files <- character(0)
  for (pat in patterns) {
    hits <- Sys.glob(file.path(root, pat))
    files <- c(files, hits)
  }
  ## Also walk metadata shallow
  meta <- file.path(root, "metadata")
  if (dir.exists(meta)) {
    files <- c(files, list.files(meta, recursive = TRUE, full.names = TRUE))
  }
  files <- unique(normalizePath(files, winslash = "/", mustWork = FALSE))
  root_n <- normalizePath(root, winslash = "/", mustWork = FALSE)
  rel <- vapply(files, function(f) {
    .ic2_soda_path_rel(root_n, f) %||% basename(f)
  }, character(1))
  rel <- sort(unique(rel[nzchar(rel) & rel != "MANIFEST.txt"]))
  info <- file.info(file.path(root_n, rel))
  lines <- c(
    "# SODA MANIFEST",
    paste0("# root: ", basename(root_n)),
    paste0("# generated: ", as.character(Sys.time())),
    "# columns: path\\tsize_bytes",
    sprintf("%s\t%s", rel, as.character(ifelse(is.na(info$size), "", info$size)))
  )
  writeLines(lines, file.path(root, "MANIFEST.txt"), useBytes = TRUE)
  invisible(TRUE)
}


#' Copy package architecture.txt into a SODA root
#' @keywords internal
.ic2_soda_write_architecture <- function(root) {
  if (is.null(root) || !nzchar(as.character(root)[1L])) return(invisible(FALSE))
  root <- as.character(root)[1L]
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(root, "architecture.txt")
  src <- system.file("soda", "architecture.txt", package = "iCellR2")
  if (!nzchar(src) || !file.exists(src)) {
    roots <- unique(c(
      getwd(),
      file.path(getwd(), "..", "iCellR2"),
      file.path(getwd(), "iCellR2")
    ))
    for (r in roots) {
      p <- normalizePath(
        file.path(r, "inst", "soda", "architecture.txt"),
        mustWork = FALSE
      )
      if (file.exists(p)) {
        src <- p
        break
      }
    }
  }
  if (!nzchar(src) || !file.exists(src)) return(invisible(FALSE))
  ok <- tryCatch({
    file.copy(src, dest, overwrite = TRUE)
  }, error = function(e) FALSE)
  invisible(isTRUE(ok))
}

#' Rebuild assay handle list from a SODA assay directory
#' @keywords internal
.ic2_soda_load_assay_ref <- function(assay_dir, assay.name) {
  if (!dir.exists(assay_dir)) return(NULL)

  sample_dirs <- list.dirs(assay_dir, recursive = FALSE, full.names = TRUE)
  top_h5 <- list.files(assay_dir, pattern = "\\.(h5|hdf5)$", full.names = TRUE, ignore.case = TRUE)
  has_sample_h5 <- any(vapply(sample_dirs, function(sd) {
    length(list.files(sd, pattern = "\\.(h5|hdf5)$", full.names = FALSE, ignore.case = TRUE)) > 0L
  }, logical(1)))
  has_matrix <- length(top_h5) > 0L || isTRUE(has_sample_h5)

  ## Lazy norm: size factors only (no matrix written under this assay dir)
  sf_path <- file.path(assay_dir, "size_factors.tsv")
  assay_json <- file.path(assay_dir, "assay.json")
  if (file.exists(sf_path) && !has_matrix) {
    meta <- list()
    if (file.exists(assay_json)) {
      meta <- tryCatch(
        jsonlite::fromJSON(assay_json, simplifyVector = TRUE, simplifyDataFrame = FALSE),
        error = function(e) list()
      )
    }
    fu <- meta$feature_universe %||% character(0)
    if (is.list(fu)) fu <- unlist(fu, use.names = FALSE)
    return(list(
      aggregation_type = "lazy_norm",
      source_assay = as.character(meta$source_assay %||% "counts")[1L],
      size_factors_path = normalizePath(sf_path, winslash = "/", mustWork = FALSE),
      feature_universe = as.character(fu),
      n_cells = as.integer(meta$n_cells %||% NA_integer_),
      n_features = as.integer(
        meta$n_features %||% length(as.character(fu))
      )
    ))
  }

  handles <- list()
  paths <- list()

  ## Single-file assay (non-virtual)
  if (length(sample_dirs) < 1L && length(top_h5) >= 1L) {
    h <- .ic2_assay_handle(path = top_h5[[1L]], backend = "hdf5_sparse", group = "counts")
    return(list(
      aggregation_type = "single",
      handle = h,
      path = top_h5[[1L]]
    ))
  }

  if (length(sample_dirs) < 1L) {
    ## Maybe files are directly named sample*.h5 under assay_dir
    sample_dirs <- assay_dir
  }

  for (sd in sample_dirs) {
    sid <- basename(sd)
    h5s <- list.files(sd, pattern = "\\.(h5|hdf5)$", full.names = TRUE, ignore.case = TRUE)
    json_ptr <- file.path(sd, "tenx_pointer.json")
    rds_ptr <- file.path(sd, "tenx_pointer.rds")
    ptr_file <- if (file.exists(json_ptr)) json_ptr else if (file.exists(rds_ptr)) rds_ptr else NA_character_
    if (!is.na(ptr_file)) {
      meta <- if (grepl("\\.json$", ptr_file, ignore.case = TRUE)) {
        jsonlite::fromJSON(ptr_file, simplifyVector = TRUE, simplifyDataFrame = FALSE)
      } else {
        readRDS(ptr_file)
      }
      sh <- list(
        type = "iCellR2_assay_handle",
        backend = "tenx_h5",
        path = meta$source %||% meta$path,
        group = meta$group %||% "matrix",
        n_features = meta$n_features,
        n_cells = meta$n_cells,
        sample_id = sid,
        meta_path = meta$meta_path
      )
      class(sh) <- c("iCellR2_assay_handle", "list")
      handles[[length(handles) + 1L]] <- sh
      paths[[length(paths) + 1L]] <- sh$path
      next
    }
    if (length(h5s) < 1L) next
    ## Prefer sample-named h5
    pref <- file.path(sd, paste0(sid, ".h5"))
    hp <- if (file.exists(pref)) pref else h5s[[1L]]
    ## Detect tenx vs icellr2 csc by structure
    backend <- "hdf5_sparse"
    group <- "counts"
    try({
      hh <- hdf5r::H5File$new(hp, mode = "r")
      on.exit(try(hh$close_all(), silent = TRUE), add = TRUE)
      if (hh$exists("matrix")) {
        backend <- "tenx_h5"
        group <- "matrix"
      }
      try(hh$close_all(), silent = TRUE)
    }, silent = TRUE)
    meta_path <- NULL
    if (identical(backend, "tenx_h5")) {
      meta_cands <- c(
        file.path(dirname(hp), "tenx_meta.json"),
        file.path(dirname(hp), "tenx_meta.rds"),
        file.path(sd, "tenx_meta.json"),
        file.path(sd, "tenx_meta.rds")
      )
      meta_hit <- meta_cands[file.exists(meta_cands)]
      if (length(meta_hit)) meta_path <- meta_hit[[1L]]
    }
    sh <- list(
      type = "iCellR2_assay_handle",
      backend = backend,
      path = normalizePath(hp, winslash = "/", mustWork = FALSE),
      group = group,
      sample_id = sid,
      meta_path = meta_path
    )
    class(sh) <- c("iCellR2_assay_handle", "list")
    handles[[length(handles) + 1L]] <- sh
    paths[[length(paths) + 1L]] <- sh$path
  }

  if (length(handles) < 1L) return(NULL)
  if (length(handles) == 1L) {
    return(list(
      aggregation_type = "single",
      handle = handles[[1L]],
      path = paths[[1L]]
    ))
  }
  list(
    aggregation_type = "virtual",
    sample_handles = handles,
    sample_assay_paths = paths
  )
}

#' Load reduction refs from a reductions directory
#' @keywords internal
.ic2_soda_load_reductions <- function(red_root) {
  out <- list()
  if (!dir.exists(red_root)) return(out)
  ## Support both modalities/<mod>/reductions and legacy <soda>/reductions
  dirs <- list.dirs(red_root, recursive = FALSE, full.names = TRUE)
  for (d in dirs) {
    name <- basename(d)
    emb <- c(
      file.path(d, "cell_embeddings.h5"),
      file.path(d, "cell_embeddings.rds")
    )
    emb <- emb[file.exists(emb)][1L]
    loadp <- c(
      file.path(d, "feature_loadings.h5"),
      file.path(d, "feature_loadings.rds")
    )
    loadp <- loadp[file.exists(loadp)][1L]
    stdev <- c(file.path(d, "stdev.tsv"), file.path(d, "stdev.rds"))
    stdev <- stdev[file.exists(stdev)][1L]
    meta <- c(file.path(d, "reduction.json"),
              file.path(d, paste0(name, "_meta.rds")), file.path(d, "pca_meta.rds"),
              file.path(d, "umap_meta.rds"))
    meta <- meta[file.exists(meta)][1L]
    if (is.na(emb) && is.na(meta)) next
    ref <- list(
      cell_embeddings = if (!is.na(emb)) emb else NULL,
      feature_loadings = if (!is.na(loadp)) loadp else NULL,
      meta_path = if (!is.na(meta)) meta else NULL
    )
    if (!is.na(stdev) && grepl("\\.rds$", stdev)) {
      ref$stdev <- tryCatch(readRDS(stdev), error = function(e) NULL)
    } else if (!is.na(stdev) && grepl("\\.tsv$", stdev)) {
      ref$stdev <- tryCatch({
        df <- utils::read.delim(stdev, sep = "\t", stringsAsFactors = FALSE)
        as.numeric(df$stdev %||% df[[ncol(df)]])
      }, error = function(e) NULL)
    }
    if (!is.na(meta)) {
      mm <- tryCatch(.ic2_soda_read_meta(meta), error = function(e) NULL)
      if (is.list(mm)) {
        if (!is.null(mm$opt.pcs)) ref$opt.pcs <- mm$opt.pcs
        if (!is.null(mm$npcs)) ref$npcs <- mm$npcs
        if (!is.null(mm$assay.name) || !is.null(mm$assay)) {
          ref$assay <- mm$assay.name %||% mm$assay
        }
        if (!is.null(mm$stdev)) ref$stdev <- mm$stdev
        if (!is.null(mm$paths$embeddings) && is.null(ref$cell_embeddings)) {
          ref$cell_embeddings <- mm$paths$embeddings
        }
      }
    }
    out[[name]] <- ref
  }
  out
}

#' Open a SODA archive as an iCellR2 object
#'
#' Rebuilds a handle-only \code{iCellR2} object from an on-disk SODA folder.
#' Prefers schema 1.0.0 (JSON + TSV + HDF5) and falls back to legacy
#' \code{.rds} metadata.
#'
#' @param path Path to a \code{*.soda} directory (or stem without \code{.soda}).
#' @param modality Preferred modality to activate. Default first available / \code{"rna"}.
#' @param verbose Progress style.
#'
#' @return An S4 \code{iCellR2} object (metadata + on-disk handles only).
#'
#' @examples
#' \dontrun{
#' obj <- open.soda("project.soda")
#' obj <- open.soda("/path/to/project")
#' }
#'
#' @export
open.soda <- function(path, modality = NULL, verbose = TRUE) {
  if (missing(path) || is.null(path) || !nzchar(as.character(path)[1L])) {
    stop("'path' must be a SODA directory.")
  }
  path <- as.character(path)[1L]
  if (!dir.exists(path) && dir.exists(paste0(path, ".soda"))) {
    path <- paste0(path, ".soda")
  }
  if (!dir.exists(path)) stop("SODA directory not found: ", path)
  root <- normalizePath(path, winslash = "/", mustWork = TRUE)

  rep <- .ic2_new_reporter("open.soda", verbose, soda_root = root, persist = FALSE)
  .ic2_step(rep, paste0("[1/3] Reading manifest from ", root))

  man <- NULL
  soda_json <- file.path(root, "soda.json")
  soda_info <- file.path(root, "soda.info")
  soda_tsv_legacy <- file.path(root, "soda.tsv")
  man_rds <- file.path(root, "manifest.rds")
  if (file.exists(soda_json)) {
    man <- jsonlite::fromJSON(soda_json, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  } else if (file.exists(soda_info)) {
    man <- .ic2_read_kv_tsv(soda_info)
  } else if (file.exists(soda_tsv_legacy)) {
    man <- .ic2_read_kv_tsv(soda_tsv_legacy)
  } else if (file.exists(man_rds)) {
    man <- readRDS(man_rds)
  } else {
    man <- list(
      format = "SODA",
      soda_version = "0.2.0",
      schema_version = "0.2.0",
      root_dir = root,
      project_name = basename(root)
    )
  }

  meta.dir <- file.path(root, "metadata")
  .ic2_step(rep, "[2/3] Loading metadata tables")

  cell_data <- .ic2_soda_try_table(file.path(meta.dir, "cell_data"))
  sample_data <- .ic2_soda_try_table(file.path(meta.dir, "sample_data"))
  project_info <- .ic2_soda_try_json_or_rds(file.path(meta.dir, "project_info"))
  version <- .ic2_soda_try_json_or_rds(file.path(meta.dir, "version"))
  params <- .ic2_soda_try_json_or_rds(file.path(meta.dir, "params"))
  commands <- .ic2_soda_try_json_or_rds(file.path(meta.dir, "commands"))

  if (is.null(project_info)) {
    project_info <- list(
      project_id = make.names(basename(root)),
      project_name = man$project_name %||% basename(root),
      date_created = man$created %||% NA_character_,
      package_version = man$package_version %||% NA_character_
    )
  }
  if (is.null(version)) {
    version <- list(
      object_version = "0.2.0",
      schema_version = man$schema_version %||% "1.0.0",
      package_version = .ic2_package_version(),
      soda_version = man$soda_version %||% "1.0.0"
    )
  }
  if (is.null(params)) params <- list()
  if (is.null(commands)) {
    commands <- list(command_log = list())
  } else if (is.null(commands$command_log)) {
    ## legacy may store log at top level
    if (!is.null(commands$make.obj) || length(commands) > 0L) {
      commands <- list(command_log = commands)
    } else {
      commands$command_log <- list()
    }
  }
  if (is.null(cell_data)) cell_data <- data.frame(row.names = character(0))
  if (is.null(sample_data)) {
    sample_data <- data.frame(sample_id = character(0), stringsAsFactors = FALSE)
  }

  .ic2_step(rep, "[3/3] Rebuilding modality handles")
  modalities <- list()

  mod_names <- names(man$modalities %||% list())
  mod_disk <- character(0)
  mods_root <- file.path(root, "modalities")
  if (dir.exists(mods_root)) {
    mod_disk <- list.dirs(mods_root, recursive = FALSE, full.names = FALSE)
  }
  mod_names <- unique(c(mod_names, mod_disk))
  if (length(mod_names) < 1L) mod_names <- "rna"

  for (mod in mod_names) {
    mdir <- file.path(root, "modalities", mod)
    assays <- list()
    active_assay <- NULL
    assay_root <- file.path(mdir, "assays")
    assay_names <- if (dir.exists(assay_root)) {
      list.dirs(assay_root, recursive = FALSE, full.names = FALSE)
    } else {
      character(0)
    }
    ## Also from manifest
    man_assays <- names(man$modalities[[mod]]$assays %||% list())
    assay_names <- unique(c(assay_names, man_assays))
    for (an in assay_names) {
      aref <- .ic2_soda_load_assay_ref(file.path(assay_root, an), an)
      if (!is.null(aref)) {
        assays[[an]] <- aref
        if (is.null(active_assay)) active_assay <- an
      }
    }
    if (!is.null(man$modalities[[mod]]$active_assay) &&
        man$modalities[[mod]]$active_assay %in% names(assays)) {
      active_assay <- man$modalities[[mod]]$active_assay
    }
    ## Prefer normalized / lazy data when present (manifest may still say counts)
    if ("data" %in% names(assays) &&
        (is.null(active_assay) || identical(active_assay, "counts"))) {
      active_assay <- "data"
    } else if (is.null(active_assay) && length(assays)) {
      active_assay <- names(assays)[[1L]]
    }

    ## feature_data
    fd <- .ic2_soda_try_table(file.path(meta.dir, paste0(mod, "_feature_data")))
    if (is.null(fd)) {
      fd <- .ic2_soda_try_table(file.path(mdir, "feature_data"))
    }

    ## reductions: modality-local first, then legacy soda/reductions
    reds <- .ic2_soda_load_reductions(file.path(mdir, "reductions"))
    if (length(reds) < 1L) {
      reds <- .ic2_soda_load_reductions(file.path(root, "reductions"))
    }

    modalities[[mod]] <- list(
      assays = assays,
      active_assay = active_assay,
      active_reduction = man$modalities[[mod]]$active_reduction %||%
        if ("umap" %in% names(reds)) "umap" else if ("pca" %in% names(reds)) "pca" else NULL,
      reductions = reds,
      feature_data = fd,
      sketches = list(),
      graphs = list()
    )
  }

  if (is.null(modality)) {
    modality <- if ("rna" %in% names(modalities)) "rna" else names(modalities)[[1L]]
  }
  modality <- tolower(modality)
  if (!modality %in% names(modalities)) {
    stop("Modality '", modality, "' not found in SODA. Available: ",
         paste(names(modalities), collapse = ", "))
  }

  ## metadata paths registry
  cd_path <- .ic2_soda_find_table(file.path(meta.dir, "cell_data"))
  sd_path <- .ic2_soda_find_table(file.path(meta.dir, "sample_data"))
  metadata_paths <- list(
    cell_data = if (!is.na(cd_path)) cd_path else file.path(meta.dir, "cell_data.tsv"),
    sample_data = if (!is.na(sd_path)) sd_path else file.path(meta.dir, "sample_data.tsv"),
    project_info = file.path(meta.dir, "project_info.json"),
    commands = file.path(meta.dir, "commands.json"),
    params = file.path(meta.dir, "params.json"),
    version = file.path(meta.dir, "version.json"),
    manifest = if (file.exists(soda_json)) soda_json else man_rds
  )
  if (!is.null(man$metadata_paths) && is.list(man$metadata_paths)) {
    for (nm in names(man$metadata_paths)) {
      if (!is.null(man$metadata_paths[[nm]])) {
        metadata_paths[[nm]] <- .ic2_soda_path_abs(root, man$metadata_paths[[nm]])
      }
    }
  }

  init_assay <- modalities[[modality]]$active_assay
  n_cells <- nrow(cell_data)
  n_features <- {
    fd <- modalities[[modality]]$feature_data
    if (!is.null(fd)) nrow(fd) else NA_integer_
  }

  obj <- methods::new(
    "iCellR2",
    project_info = project_info,
    samples = list(sample_data = sample_data),
    cell_data = cell_data,
    modalities = modalities,
    multimodal = .ic2_multimodal_from_disk(
      root, .ic2_empty_multimodal()
    ),
    reductions = list(
      integrated_pca = NULL, integrated_umap = NULL,
      integrated_tsne = NULL, integrated_knetl = NULL, misc = list()
    ),
    graphs = list(
      integrated_knn = NULL, integrated_snn = NULL, wnn = NULL, misc = list()
    ),
    storage = list(
      root_dir = root,
      backend = "soda",
      assay_backend = "hdf5_sparse",
      modality_paths = {
        mp <- list()
        for (m in names(modalities)) {
          mp[[m]] <- file.path(root, "modalities", m)
        }
        mp
      },
      metadata_paths = metadata_paths,
      reduction_paths = list(),
      graph_paths = list(),
      image_paths = list(),
      molecule_paths = list(),
      segment_paths = list(),
      cache = list(dir = file.path(root, "cache")),
      compression = NA_character_,
      assay_paths = list()
    ),
    commands = commands,
    params = params,
    version = version,
    misc = list(
      n_features = n_features,
      n_cells = n_cells,
      n_samples = nrow(sample_data),
      initialized_modality = modality,
      initialized_assay = init_assay,
      handle_only = TRUE,
      opened_from = root
    ),
    validity = list(
      valid = TRUE,
      checked = as.character(Sys.time()),
      messages = "Opened from SODA (handle-only).",
      schema_notes = paste0("open.soda schema ", version$schema_version %||% "?")
    )
  )

  ## Mirror opt.pcs into misc when present
  pca <- modalities[[modality]]$reductions$pca
  if (!is.null(pca$opt.pcs)) obj@misc$opt.pcs <- pca$opt.pcs

  .ic2_finish(rep, extra_lines = paste0(
    "Opened: ", basename(root),
    "   |   modality: ", modality,
    "   |   cells: ", .ic2_fmt_num(n_cells),
    "   |   assays: ", paste(names(modalities[[modality]]$assays), collapse = ",")
  ))
  obj
}


#' Validate a SODA archive layout
#'
#' Checks that a SODA directory has the required manifest, metadata tables,
#' and that referenced paths resolve. Returns a list; set
#' \code{error = TRUE} to stop on failure.
#'
#' @param path SODA directory (or project name with \code{.soda} suffix).
#' @param error If \code{TRUE}, stop when invalid.
#' @param verbose Print a short report.
#' @return A list with \code{ok}, \code{errors}, \code{warnings}, \code{checks}.
#' @export
validate.soda <- function(path, error = FALSE, verbose = TRUE) {
  errors <- character(0)
  warnings <- character(0)
  checks <- list()

  if (missing(path) || is.null(path) || !nzchar(as.character(path)[1L])) {
    stop("'path' must be a SODA directory.")
  }
  path <- as.character(path)[1L]
  if (!dir.exists(path) && dir.exists(paste0(path, ".soda"))) {
    path <- paste0(path, ".soda")
  }
  if (!dir.exists(path)) {
    errors <- c(errors, paste0("SODA directory not found: ", path))
    out <- list(ok = FALSE, errors = errors, warnings = warnings, checks = checks)
    if (isTRUE(error)) stop(paste(errors, collapse = "\n"))
    return(out)
  }
  root <- normalizePath(path, winslash = "/", mustWork = TRUE)
  checks$root <- root

  soda_json <- file.path(root, "soda.json")
  if (!file.exists(soda_json)) {
    errors <- c(errors, "Missing soda.json")
  } else {
    man <- tryCatch(
      jsonlite::fromJSON(soda_json, simplifyVector = TRUE, simplifyDataFrame = FALSE),
      error = function(e) {
        errors <<- c(errors, paste0("soda.json parse error: ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(man)) {
      checks$format <- man$format %||% NA_character_
      checks$schema_version <- man$schema_version %||% NA_character_
      checks$soda_version <- man$soda_version %||% NA_character_
      if (!identical(as.character(man$format %||% ""), "SODA")) {
        warnings <- c(warnings, "soda.json format field is not 'SODA'")
      }
      if (is.null(man$schema_version) || !nzchar(as.character(man$schema_version)[1L])) {
        warnings <- c(warnings, "soda.json missing schema_version")
      }
      ## Absolute root_dir is discouraged
      rd <- as.character(man$root_dir %||% "")[1L]
      if (nzchar(rd) && grepl("^(~/|/|[A-Za-z]:)", rd) && !identical(rd, ".")) {
        warnings <- c(
          warnings,
          "soda.json root_dir is absolute; prefer '.' (re-save with a recent iCellR2)"
        )
      }
      ## Resolve metadata paths
      mp <- man$metadata_paths %||% list()
      for (nm in names(mp)) {
        p <- .ic2_soda_path_abs(root, mp[[nm]])
        if (is.null(p) || !file.exists(p)) {
          warnings <- c(warnings, paste0("metadata_paths$", nm, " missing: ", mp[[nm]]))
        }
      }
      ## Assay paths from manifest
      for (mod in names(man$modalities %||% list())) {
        assays <- man$modalities[[mod]]$assays %||% list()
        for (an in names(assays)) {
          ap <- assays[[an]]$path
          if (!is.null(ap) && nzchar(as.character(ap)[1L])) {
            ap_abs <- .ic2_soda_path_abs(root, ap)
            if (is.null(ap_abs) || !file.exists(ap_abs)) {
              ## virtual / multi-file assays may point at a dir
              adir <- file.path(root, "modalities", mod, "assays", an)
              if (!dir.exists(adir)) {
                warnings <- c(
                  warnings,
                  paste0("Assay path missing for ", mod, "/", an, ": ", ap)
                )
              }
            }
          }
        }
      }
    }
  }

  req_meta <- c("cell_data.tsv", "project_info.json", "version.json")
  for (f in req_meta) {
    p <- file.path(root, "metadata", f)
    ## cell_data may be .tsv or stem
    if (identical(f, "cell_data.tsv")) {
      alt <- .ic2_soda_find_table(file.path(root, "metadata", "cell_data"))
      if (is.na(alt) || !file.exists(alt)) {
        errors <- c(errors, "Missing metadata/cell_data table")
      }
    } else if (!file.exists(p)) {
      ## project_info / version may be .rds legacy
      stem <- sub("\\.json$", "", f)
      if (!file.exists(file.path(root, "metadata", paste0(stem, ".rds")))) {
        warnings <- c(warnings, paste0("Missing metadata/", f))
      }
    }
  }

  if (!dir.exists(file.path(root, "modalities"))) {
    warnings <- c(warnings, "Missing modalities/ directory")
  }
  if (!file.exists(file.path(root, "README.md"))) {
    warnings <- c(warnings, "Missing README.md (will be written on next manifest refresh)")
  }

  ok <- length(errors) < 1L
  checks$n_errors <- length(errors)
  checks$n_warnings <- length(warnings)

  if (isTRUE(verbose)) {
    cat("SODA validate: ", root, "\n", sep = "")
    cat("  status: ", if (ok) "OK" else "FAIL", "\n", sep = "")
    if (length(checks$schema_version)) {
      cat("  schema: ", checks$schema_version, "  soda: ", checks$soda_version, "\n", sep = "")
    }
    if (length(errors)) {
      cat("  errors:\n")
      for (e in errors) cat("   - ", e, "\n", sep = "")
    }
    if (length(warnings)) {
      cat("  warnings:\n")
      for (w in warnings) cat("   - ", w, "\n", sep = "")
    }
  }

  out <- list(ok = ok, errors = errors, warnings = warnings, checks = checks)
  if (isTRUE(error) && !ok) stop(paste(errors, collapse = "\n"))
  invisible(out)
}

#' Refresh SODA portable sidecars (manifest, README, MANIFEST.txt)
#'
#' @param x An \code{iCellR2} object with \code{x@storage$root_dir} set,
#'   or a SODA directory path (character).
#' @return Invisible path to \code{soda.json} or the validation list for a path.
#' @export
refresh.soda <- function(x) {
  if (is.character(x)) {
    root <- x[[1L]]
    if (!dir.exists(root) && dir.exists(paste0(root, ".soda"))) {
      root <- paste0(root, ".soda")
    }
    if (!dir.exists(root)) stop("SODA directory not found: ", root)
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    man <- NULL
    sj <- file.path(root, "soda.json")
    if (file.exists(sj)) {
      man <- jsonlite::fromJSON(sj, simplifyVector = TRUE, simplifyDataFrame = FALSE)
      ## Rewrite paths as relative if absolute
      man$root_dir <- "."
      man$portable <- TRUE
      if (!is.null(man$metadata_paths)) {
        man$metadata_paths <- lapply(man$metadata_paths, function(p) .ic2_soda_path_rel(root, p))
      }
      man$modalities <- .ic2_soda_paths_make_rel(root, man$modalities)
      man$created <- as.character(Sys.time())
      jsonlite::write_json(man, sj, auto_unbox = TRUE, pretty = TRUE, null = "null")
    }
    inv <- tryCatch(.ic2_soda_write_inventory(root, man = man), error = function(e) NULL)
    .ic2_soda_write_readme(root, man = man)
    .ic2_soda_write_html_index(root, man = man, inventory = inv)
    .ic2_soda_write_file_manifest(root)
    .ic2_soda_write_architecture(root)
    return(invisible(sj))
  }
  .ic2_require_object(x)
  invisible(.ic2_soda_write_manifest(x))
}

