#' Create an iCellR2 S4 project object with a SODA backend
#'
#' Creates a **handle-only** S4 \code{iCellR2} object plus a directory-based
#' SODA (Single Omics Data Archive) on disk. Assay matrices are **not** kept
#' in the S4 object.
#'
#' Assay registration depends on the input backend:
#' \itemize{
#'   \item \strong{Native 10x HDF5} (\code{tenx_h5}): by default
#'   (\code{copy.raw = TRUE}) the \code{.h5} is \strong{copied into SODA}
#'   so the archive is self-contained. Set \code{copy.raw = FALSE} to store
#'   only a pointer to the original path (no multi-GB copy).
#'   \item \strong{Legacy iCellR2 CSC HDF5}: converted to \code{tenx_h5}
#'   when copied into SODA (raw assays are always Cell Ranger-style 10x H5).
#'   \item \strong{In-memory matrices}: written as \code{tenx_h5} under SODA
#'   (fine for modest size; not for atlas-scale counts).
#' }
#'
#' @param x Sparse matrix, assay handle, \code{iCellR2_sample}, or
#'   \code{iCellR2_aggregation}.
#' @param out.dir Output directory stem; archive is created at
#'   \code{paste0(out.dir, ".soda")}.
#' @param project_name Project name.
#' @param modality Modality to initialize.
#' @param assay.name Assay name inside the modality. Default \code{"counts"}.
#' @param organism Optional organism.
#' @param genome_build Optional genome build.
#' @param overwrite Replace an existing SODA archive.
#' @param copy.raw If \code{TRUE} (default), copy \code{tenx_h5} count
#'   matrices into SODA so the archive stands alone. If \code{FALSE}, keep
#'   an external pointer (original file must remain at a stable path).
#' @param storage.backend Preferred backend when \emph{materializing} matrices:
#'   \code{"hdf5_sparse"} writes modern \strong{10x HDF5} (\code{/matrix})
#'   for raw counts, or \code{"rds"}.
#' @param verbose Progress style: \code{TRUE}/\code{"full"}, \code{"compact"},
#'   or \code{FALSE}/\code{"silent"}.
#'
#' @return An S4 object of class \code{iCellR2} (metadata + paths only).
#'
#' @importFrom methods as new
#' @export
make.obj <- function(
    x = NULL,
    out.dir,
    project_name = "iCellR2_project",
    modality = "rna",
    assay.name = "counts",
    organism = NA_character_,
    genome_build = NA_character_,
    overwrite = FALSE,
    copy.raw = TRUE,
    storage.backend = c("hdf5_sparse", "rds"),
    verbose = TRUE
) {
  storage.backend <- match.arg(storage.backend)
  copy.raw <- isTRUE(copy.raw)
  rep <- .ic2_new_reporter("make.obj", verbose)

  if (is.null(x)) stop("'x' must be provided.")
  if (missing(out.dir) || is.null(out.dir) || !is.character(out.dir) ||
      length(out.dir) != 1L || !nzchar(out.dir)) {
    stop("'out.dir' must be a non-empty character string.")
  }
  if (!is.character(project_name) || length(project_name) != 1L || !nzchar(project_name)) {
    stop("'project_name' must be a non-empty character string.")
  }
  if (!is.character(modality) || length(modality) != 1L || !nzchar(modality)) {
    stop("'modality' must be a non-empty character string.")
  }
  modality <- tolower(modality)
  allowed_modalities <- c(
    "rna", "atac", "adt", "vdj", "spatial", "visium", "xenium", "proteomics"
  )
  if (!modality %in% allowed_modalities) {
    stop("'modality' must be one of: ", paste(allowed_modalities, collapse = ", "), ".")
  }
  if (!is.character(assay.name) || length(assay.name) != 1L || !nzchar(assay.name)) {
    stop("'assay.name' must be a non-empty character string.")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.")
  }

  ## Guard: huge in-memory matrices are the wrong input for atlas-scale projects
  if (!inherits(x, "iCellR2_aggregation") &&
      !inherits(x, "iCellR2_sample") &&
      !.ic2_is_assay_handle(x)) {
    nc <- tryCatch(ncol(x), error = function(e) NA_integer_)
    if (is.finite(nc) && nc > 500000L) {
      stop(
        "In-memory input has ", .ic2_fmt_num(nc), " cells. ",
        "For atlas-scale data use load10x(..., output = \"virtual\") ",
        "(native tenx_h5 handle) and optionally data.aggregation(..., output = \"virtual\")."
      )
    }
    if (is.finite(nc) && nc > 100000L) {
      warning(
        "Large in-memory matrix (", .ic2_fmt_num(nc), " cells). ",
        "Prefer load10x(..., output = \"virtual\") to avoid densifying / copying counts.",
        call. = FALSE
      )
    }
  }

  ## Efficient cell_data from aggregation maps (avoids repeated rbind copies)
  .ic2_cell_data_from_maps <- function(cell_maps, modality) {
    n <- sum(vapply(cell_maps, function(cm) nrow(cm), integer(1)))
    cell_id <- character(n)
    barcode <- character(n)
    sample_id <- character(n)
    condition_id <- character(n)
    pos <- 1L
    for (cm in cell_maps) {
      k <- nrow(cm)
      if (k < 1L) next
      idx <- seq.int(pos, length.out = k)
      cell_id[idx] <- as.character(cm$global_cell)
      barcode[idx] <- as.character(cm$local_cell)
      sample_id[idx] <- as.character(cm$sample_id)
      condition_id[idx] <- as.character(cm$condition_id)
      pos <- pos + k
    }
    entity <- if (modality %in% c("spatial", "visium")) "spot" else "cell"
    out <- data.frame(
      cell_id = cell_id,
      barcode = barcode,
      sample_id = sample_id,
      donor_id = NA_character_,
      batch_id = NA_character_,
      condition_id = condition_id,
      entity_type = entity,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rownames(out) <- cell_id
    out
  }

  ## Split "sample_barcode" on the first underscore only (faster than strsplit)
  .ic2_split_sample_barcode <- function(cell_ids) {
    cell_ids <- as.character(cell_ids)
    sp <- regexpr("_", cell_ids, fixed = TRUE)
    if (length(sp) && all(sp > 0L)) {
      list(
        sample_id = substr(cell_ids, 1L, sp - 1L),
        barcode = substring(cell_ids, sp + 1L)
      )
    } else {
      list(
        sample_id = rep("sample1", length(cell_ids)),
        barcode = cell_ids
      )
    }
  }

  soda_root <- normalizePath(paste0(out.dir, ".soda"), winslash = "/", mustWork = FALSE)
  if (dir.exists(soda_root)) {
    if (!overwrite) {
      stop("SODA archive already exists at '", soda_root, "'. Use overwrite = TRUE.")
    }
    unlink(soda_root, recursive = TRUE, force = TRUE)
  }
  dir.create(soda_root, recursive = TRUE, showWarnings = FALSE)
  rep <- .ic2_attach_soda_log(rep, soda_root)
  ## Input: source assay/handle path when available; output: SODA archive
  src_path <- NULL
  if (inherits(x, "iCellR2_sample") && !is.null(x$handle$path)) {
    src_path <- x$handle$path
  } else if (inherits(x, "iCellR2_aggregation") && !is.null(x$sample_handles[[1L]]$path)) {
    src_path <- x$sample_handles[[1L]]$path
  } else if (.ic2_is_assay_handle(x) && !is.null(x$path)) {
    src_path <- x$path
  }
  rep <- .ic2_set_io_dirs(
    rep,
    input_dir = if (!is.null(src_path)) src_path else getwd(),
    output_dir = soda_root
  )
  if (!is.null(src_path)) {
    rep <- .ic2_add_input(rep, "Source", normalizePath(src_path, winslash = "/", mustWork = FALSE))
  }

  input_kind <- if (inherits(x, "iCellR2_aggregation")) {
    "virtual aggregation"
  } else if (inherits(x, "iCellR2_sample")) {
    "sample"
  } else if (.ic2_is_assay_handle(x)) {
    "assay handle"
  } else {
    "matrix"
  }
  rep <- .ic2_banner_start(rep, c(
    paste0("Project: ", project_name, "   |   modality: ", modality, "   |   assay: ", assay.name),
    paste0("Input: ", input_kind, "   |   SODA: ", soda_root)
  ))
  .ic2_step(rep, "[1/3] Creating SODA layout")
  rep <- .ic2_progress_bar(rep, 1L, 3L, prefix = "  make.obj")

  modality.dir <- file.path(soda_root, "modalities", modality)
  assay.dir <- file.path(modality.dir, "assays", assay.name)
  meta.dir <- file.path(soda_root, "metadata")
  reductions.dir <- file.path(soda_root, "reductions")
  graphs.dir <- file.path(soda_root, "graphs")
  images.dir <- file.path(soda_root, "images")
  molecules.dir <- file.path(soda_root, "molecules")
  segments.dir <- file.path(soda_root, "segments")
  multimodal.dir <- file.path(soda_root, "multimodal")
  sketch.dir <- file.path(modality.dir, "sketches")
  indexes.dir <- file.path(soda_root, "indexes")
  cache.dir <- file.path(soda_root, "cache")

  invisible(lapply(
    c(modality.dir, assay.dir, meta.dir, reductions.dir, graphs.dir,
      images.dir, molecules.dir, segments.dir, multimodal.dir,
      sketch.dir, indexes.dir, cache.dir),
    dir.create, recursive = TRUE, showWarnings = FALSE
  ))

  empty_modality <- function() {
    list(
      assays = list(),
      active_assay = NULL,
      feature_data = NULL,
      reductions = list(),
      active_reduction = NULL,
      graphs = list(),
      sketches = list(
        sketch_ids = NULL,
        sketch_models = list(),
        projections = list()
      ),
      markers = list(),
      models = list(),
      params = list(),
      misc = list()
    )
  }

  modalities <- list(
    rna = empty_modality(),
    atac = empty_modality(),
    adt = empty_modality(),
    vdj = list(
      contigs = NULL, clonotypes = NULL, chains = NULL, repertoire = NULL,
      feature_data = NULL, params = list(), misc = list()
    ),
    spatial = empty_modality(),
    visium = empty_modality(),
    xenium = empty_modality(),
    proteomics = empty_modality(),
    misc_modalities = list()
  )

  copy_handle_into_soda <- function(handle, dest_path) {
    handle <- .ic2_enrich_handle(handle, load_dimnames = FALSE)
    dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
    if (is.null(handle$path) || !nzchar(as.character(handle$path)[1L])) {
      stop("Assay handle has no path.")
    }
    if (!file.exists(handle$path)) {
      stop("Assay file not found: ", handle$path)
    }
    ## Native 10x H5
    if (identical(handle$backend, "tenx_h5")) {
      meta_dest <- file.path(dirname(dest_path), "tenx_meta.json")
      if (is.null(handle$meta_path) || !file.exists(handle$meta_path)) {
        stop("tenx_h5 handle missing meta_path.")
      }
      h5_dest <- if (grepl("\\.h5$", dest_path, ignore.case = TRUE)) {
        dest_path
      } else {
        paste0(tools::file_path_sans_ext(dest_path), ".h5")
      }
      meta0 <- .ic2_soda_read_tenx_meta(handle$meta_path)
      n_src <- as.integer(meta0$n_features_src %||% 0L)
      n_keep <- length(meta0$row_keep %||% integer(0))
      mixed <- !is.null(meta0$feature.type) && n_src > 0L && n_keep > 0L && n_keep < n_src
      if (isTRUE(mixed)) {
        .ic2_step(rep, paste0(
          "  Writing ", paste(meta0$feature.type, collapse = ", "),
          " only  ->  ", basename(h5_dest)
        ))
        cells <- meta0$cell_ids
        if (!is.null(handle$cell_ids) && length(handle$cell_ids) == length(cells)) {
          cells <- as.character(handle$cell_ids)
        }
        split_h <- .ic2_split_10x_feature_type(
          src_path = handle$path,
          dest_path = h5_dest,
          feature.type = meta0$feature.type,
          cell.ids = cells,
          sample_id = handle$sample_id,
          modality = handle$modality %||% modality
        )
        return(split_h)
      }
      if (isTRUE(copy.raw)) {
        ## Own a full copy inside SODA (archive is self-contained)
        if (!identical(
          normalizePath(handle$path, mustWork = FALSE),
          normalizePath(h5_dest, mustWork = FALSE)
        )) {
          .ic2_step(rep, paste0("  Copying raw 10x H5  ->  ", basename(h5_dest)))
          if (file.exists(h5_dest)) unlink(h5_dest)
          if (!file.copy(handle$path, h5_dest, overwrite = TRUE)) {
            stop("Failed to copy '", handle$path, "'  ->  '", h5_dest, "'.")
          }
        }
        handle$path <- normalizePath(h5_dest, winslash = "/", mustWork = TRUE)
        meta <- .ic2_soda_read_tenx_meta(handle$meta_path)
        if (is.null(meta)) stop("Could not read tenx meta: ", handle$meta_path)
        ## Keep aggregation-prefixed barcodes if present on the handle
        if (!is.null(handle$cell_ids) && length(handle$cell_ids)) {
          meta$cell_ids <- as.character(handle$cell_ids)
        }
        if (!is.null(handle$feature_ids) && length(handle$feature_ids)) {
          meta$feature_ids <- as.character(handle$feature_ids)
        }
        meta$source <- handle$path
        meta_dest <- .ic2_soda_write_tenx_meta(meta, meta_dest)
        .ic2_soda_unlink_rds(
          file.path(dirname(dest_path), "tenx_meta.rds"),
          file.path(dirname(dest_path), "tenx_pointer.rds")
        )
        handle$meta_path <- normalizePath(meta_dest, winslash = "/", mustWork = TRUE)
        return(handle)
      }
      ## Pointer-only (external file must remain at a stable path)
      meta <- .ic2_soda_read_tenx_meta(handle$meta_path)
      if (is.null(meta)) stop("Could not read tenx meta: ", handle$meta_path)
      meta_dest <- .ic2_soda_write_tenx_meta(meta, meta_dest)
      handle$meta_path <- normalizePath(meta_dest, winslash = "/", mustWork = TRUE)
      handle$path <- normalizePath(handle$path, winslash = "/", mustWork = TRUE)
      ptr_path <- file.path(dirname(dest_path), "tenx_pointer.json")
      .ic2_soda_write_json(
        list(
          type = "tenx_h5_soda_pointer",
          source = handle$path,
          meta_path = handle$meta_path,
          group = handle$group,
          n_features = handle$n_features,
          n_cells = handle$n_cells
        ),
        ptr_path
      )
      .ic2_soda_unlink_rds(
        file.path(dirname(dest_path), "tenx_meta.rds"),
        file.path(dirname(dest_path), "tenx_pointer.rds")
      )
      return(handle)
    }
    ## Legacy CSC  ->  convert to 10x H5 inside SODA (raw assays are always tenx_h5)
    if (identical(handle$backend, "hdf5_sparse")) {
      h5_dest <- if (grepl("\\.h5$", dest_path, ignore.case = TRUE)) {
        dest_path
      } else {
        paste0(tools::file_path_sans_ext(dest_path), ".h5")
      }
      .ic2_step(rep, paste0("  Converting CSC  ->  10x H5  ->  ", basename(h5_dest)))
      return(.ic2_convert_csc_handle_to_tenx(handle, h5_dest, overwrite = TRUE))
    }
    if (!identical(normalizePath(handle$path, mustWork = FALSE),
                   normalizePath(dest_path, mustWork = FALSE))) {
      file.copy(handle$path, dest_path, overwrite = TRUE)
    }
    handle$path <- normalizePath(dest_path, winslash = "/", mustWork = TRUE)
    handle
  }

  build_object <- function(
      project_info, sample_data, cell_data, feature_data,
      assay_ref, n_features, n_cells, n_samples, sample_ids,
      paths, commands, params, version, schema_notes
  ) {
    modalities[[modality]]$assays[[assay.name]] <- assay_ref
    modalities[[modality]]$active_assay <- assay.name
    ## feature_data kept as compact data.frame (ids only); large extras go on disk
    modalities[[modality]]$feature_data <- feature_data

    methods::new(
      "iCellR2",
      project_info = project_info,
      samples = list(sample_data = sample_data),
      cell_data = cell_data,
      modalities = modalities,
      multimodal = .ic2_empty_multimodal(),
      reductions = list(
        integrated_pca = NULL, integrated_umap = NULL,
        integrated_tsne = NULL, integrated_knetl = NULL, misc = list()
      ),
      graphs = list(
        integrated_knn = NULL, integrated_snn = NULL, wnn = NULL, misc = list()
      ),
      storage = list(
        root_dir = soda_root,
        backend = "soda",
        assay_backend = storage.backend,
        modality_paths = list(
          rna = file.path(soda_root, "modalities", "rna"),
          atac = file.path(soda_root, "modalities", "atac"),
          adt = file.path(soda_root, "modalities", "adt"),
          vdj = file.path(soda_root, "modalities", "vdj"),
          spatial = file.path(soda_root, "modalities", "spatial"),
          visium = file.path(soda_root, "modalities", "visium"),
          xenium = file.path(soda_root, "modalities", "xenium"),
          proteomics = file.path(soda_root, "modalities", "proteomics")
        ),
        metadata_paths = paths$metadata_paths,
        reduction_paths = list(),
        graph_paths = list(),
        image_paths = list(),
        molecule_paths = list(),
        segment_paths = list(),
        chunk_index = file.path(indexes.dir, "chunk_index.rds"),
        cache = list(dir = cache.dir),
        compression = NA_character_,
        assay_paths = paths$assay_paths
      ),
      commands = commands,
      params = params,
      version = version,
      misc = list(
        n_features = n_features,
        n_cells = n_cells,
        n_samples = n_samples,
        initialized_modality = modality,
        initialized_assay = assay.name,
        sample_ids = sample_ids,
        handle_only = TRUE
      ),
      validity = list(
        valid = TRUE,
        checked = as.character(Sys.time()),
        messages = "Object initialized successfully (handle-only).",
        schema_notes = schema_notes
      )
    )
  }

  save_meta <- function(obj, path) .ic2_soda_write_json(obj, path)

  ## Portable metadata writers (TSV tables + JSON lists)
  save_table <- function(df, stem) {
    .ic2_soda_write_table(df, stem)
  }
  save_json <- function(obj, path) {
    .ic2_soda_write_json(obj, path)
  }

  ## ---- Virtual aggregation (handles only) ----
  if (inherits(x, "iCellR2_aggregation")) {
    .ic2_step(rep, "[2/3] Registering virtual sample handles")
    modality <- x$modality
    handles <- if (!is.null(x$handles)) x$handles else x$samples
    sample_data <- x$sample_data
    feature_universe <- x$feature_universe
    feature_maps <- x$feature_maps
    cell_maps <- x$cell_maps
    rep <- .ic2_progress_bar(rep, 2L, 3L, prefix = "  make.obj")

    feature_data <- data.frame(
      feature_id = feature_universe,
      feature_name = feature_universe,
      feature_type = modality,
      stringsAsFactors = FALSE,
      row.names = feature_universe
    )
    .ic2_step(rep, "  Building cell_data from maps")
    cell_data <- .ic2_cell_data_from_maps(cell_maps, modality)

    feature_data_path <- file.path(meta.dir, paste0(modality, "_feature_data"))
    cell_data_path <- file.path(meta.dir, "cell_data")
    sample_data_path <- file.path(meta.dir, "sample_data")
    project_info_path <- file.path(meta.dir, "project_info.json")
    manifest_path <- file.path(soda_root, "soda.json")
    command_log_path <- file.path(meta.dir, "commands.json")
    params_path <- file.path(meta.dir, "params.json")
    version_path <- file.path(meta.dir, "version.json")
    feature_maps_path <- file.path(meta.dir, "feature_maps")
    cell_maps_path <- file.path(meta.dir, "cell_maps")

    feature_data_path <- save_table(feature_data, feature_data_path)
    cell_data_path <- save_table(cell_data, cell_data_path)
    sample_data_path <- save_table(sample_data, sample_data_path)
    feature_maps_path <- .ic2_soda_write_map_list(feature_maps, feature_maps_path)
    cell_maps_path <- .ic2_soda_write_map_list(cell_maps, cell_maps_path)

    sample_handles <- vector("list", length(handles))
    sample_assay_paths <- character(length(handles))
    external_sources <- character(0)

    for (i in seq_along(handles)) {
      sid <- sample_data$sample_id[i]
      ## Assays in SODA are always HDF5 (never .rds)
      ext <- ".h5"
      ## File named by sample.id for a self-describing SODA tree
      dest <- file.path(assay.dir, sid, paste0(sid, ext))
      h <- handles[[i]]
      if (.ic2_is_assay_handle(h)) {
        h <- copy_handle_into_soda(h, dest)
      } else {
        ## legacy: matrix in aggregation$samples
        h <- .ic2_materialize_handle(
          h, filepath = dest, group = "counts",
          overwrite = TRUE, backend = storage.backend
        )
      }
      sh <- list(
        type = "iCellR2_assay_handle",
        backend = h$backend,
        path = h$path,
        group = if (!is.null(h$group)) h$group else "counts",
        n_features = h$n_features,
        n_cells = h$n_cells,
        sample_id = sid
      )
      if (identical(h$backend, "tenx_h5")) {
        sh$meta_path <- h$meta_path
        ## Only flag as external when we did not copy into SODA
        if (!isTRUE(copy.raw)) {
          external_sources <- c(external_sources, h$path)
        }
      }
      sample_handles[[i]] <- sh
      class(sample_handles[[i]]) <- c("iCellR2_assay_handle", "list")
      sample_assay_paths[i] <- h$path
    }

    assay_ref <- list(
      aggregation_type = "virtual",
      sample_handles = sample_handles,
      sample_assay_paths = as.list(sample_assay_paths),
      feature_universe = feature_universe
    )

    project_info <- list(
      project_id = make.names(project_name, unique = TRUE),
      project_name = project_name,
      description = NA_character_,
      organism = organism,
      genome_build = genome_build,
      technologies = modality,
      date_created = as.character(Sys.time()),
      package_version = .ic2_package_version()
    )
    save_json(project_info, project_info_path)

    commands <- list(
      command_log = list(
        make.obj = list(
          function_name = "make.obj",
          timestamp = as.character(Sys.time()),
          arguments = list(
            out.dir = out.dir, project_name = project_name, modality = modality,
            assay.name = assay.name, storage.backend = storage.backend,
            overwrite = overwrite
          ),
          outputs = list(
            modality = modality, assay = assay.name,
            n_features = length(feature_universe),
            n_cells = nrow(cell_data), n_samples = nrow(sample_data)
          )
        )
      ),
      timestamps = as.character(Sys.time()),
      arguments = list(),
      outputs = list()
    )
    save_json(commands, command_log_path)

    params <- list(
      import = list(
        modality = modality, assay.name = assay.name,
        n_samples = nrow(sample_data), aggregation = "virtual",
        storage.backend = storage.backend
      ),
      qc = list(), normalization = list(), feature_selection = list(),
      dimensionality_reduction = list(), graph_construction = list(),
      clustering = list(), integration = list(), visualization = list(),
      sketch = list(), misc = list()
    )
    save_json(params, params_path)

    version <- list(
      object_version = "0.2.0",
      schema_version = "1.0.0",
      package_version = .ic2_package_version(),
      soda_version = "1.0.0"
    )
    save_json(version, version_path)

    manifest <- list(
      format = "SODA",
      format_name = "Single Omics Data Archive",
      soda_version = "1.0.0",
      schema_version = "1.0.0",
      root_dir = soda_root,
      modality = modality,
      assay = assay.name,
      aggregation = "virtual",
      handle_only = TRUE,
      external_assay_sources = unique(external_sources),
      metadata_paths = list(
        project_info = project_info_path,
        feature_data = feature_data_path,
        cell_data = cell_data_path,
        sample_data = sample_data_path,
        feature_maps = feature_maps_path,
        cell_maps = cell_maps_path,
        commands = command_log_path,
        params = params_path,
        version = version_path
      )
    )
    ## soda.json written after object build via .ic2_soda_write_manifest()
    unlink(file.path(soda_root, "manifest.rds"), force = TRUE)

    .ic2_step(rep, "[3/3] Building handle-only S4 object")
    rep <- .ic2_progress_bar(rep, 3L, 3L, prefix = "  make.obj")
    rep <- .ic2_add_output(rep, "SODA", soda_root)
    rep <- .ic2_add_output(rep, "Log", file.path(soda_root, "iCellR2.log"))
    extra <- paste0(
      "Features: ", .ic2_fmt_num(length(feature_universe)),
      "   |   Cells: ", .ic2_fmt_num(nrow(cell_data)),
      "   |   Samples: ", nrow(sample_data)
    )
    if (length(external_sources)) {
      extra <- c(
        extra,
        paste0(
          "External tenx_h5 sources (do not move/delete): ",
          length(unique(external_sources)), " file(s)"
        )
      )
      .ic2_step(rep, paste0(
        "  Note: counts remain in original 10x H5  -  keep paths stable or use copy.raw = TRUE."
      ))
    } else if (isTRUE(copy.raw)) {
      extra <- c(extra, "Raw assays: copied into SODA (self-contained)")
    }
    .ic2_finish(rep, extra_lines = extra)

    obj <- build_object(
      project_info = project_info,
      sample_data = sample_data,
      cell_data = cell_data,
      feature_data = feature_data,
      assay_ref = assay_ref,
      n_features = length(feature_universe),
      n_cells = nrow(cell_data),
      n_samples = nrow(sample_data),
      sample_ids = sample_data$sample_id,
      paths = list(
        metadata_paths = list(
          project_info = project_info_path,
          feature_data = feature_data_path,
          cell_data = cell_data_path,
          sample_data = sample_data_path,
          feature_maps = feature_maps_path,
          cell_maps = cell_maps_path,
          manifest = manifest_path,
          commands = command_log_path,
          params = params_path,
          version = version_path
        ),
        assay_paths = list(
          modality = modality,
          assay = assay.name,
          aggregation = "virtual",
          sample_paths = as.list(sample_assay_paths)
        )
      ),
      commands = commands,
      params = params,
      version = version,
      schema_notes = "Handle-only multimodal iCellR2 + SODA with virtual aggregation."
    )
    .ic2_soda_write_manifest(obj)
    return(obj)
  }

  ## ---- Single sample / matrix / handle ----
  .ic2_step(rep, "[2/3] Registering assay in SODA")
  if (inherits(x, "iCellR2_sample")) {
    modality <- x$modality
    if (!is.null(x$handle)) {
      x <- x$handle
    } else if (!is.null(x$assay)) {
      x <- x$assay
    } else {
      stop("iCellR2_sample has neither handle nor assay.")
    }
  }

  ## Assays in SODA are always HDF5 (never .rds)
  ext <- ".h5"
  sid0 <- NULL
  if (.ic2_is_assay_handle(x) && !is.null(x$sample_id)) {
    sid0 <- as.character(x$sample_id)[1L]
  }
  if (is.null(sid0) || !nzchar(sid0)) sid0 <- "sample1"
  assay_path <- file.path(assay.dir, sid0, paste0(sid0, ext))

  if (.ic2_is_assay_handle(x)) {
    handle <- copy_handle_into_soda(x, assay_path)
  } else {
    handle <- .ic2_materialize_handle(
      x, filepath = assay_path, group = "counts",
      overwrite = TRUE, backend = storage.backend
    )
  }

  ## Load dimnames once (from tenx meta or CSC); drop from handle after building tables
  handle <- .ic2_enrich_handle(handle, load_dimnames = TRUE)
  feature_ids <- as.character(handle$feature_ids)
  cell_ids <- as.character(handle$cell_ids)
  handle$feature_ids <- NULL
  handle$cell_ids <- NULL

  .ic2_step(rep, "  Building cell_data")
  split_ids <- .ic2_split_sample_barcode(cell_ids)
  sample_id <- split_ids$sample_id
  barcode <- split_ids$barcode
  entity_type <- if (modality %in% c("spatial", "visium")) "spot" else "cell"

  feature_data <- data.frame(
    feature_id = feature_ids,
    feature_name = feature_ids,
    feature_type = modality,
    stringsAsFactors = FALSE,
    row.names = feature_ids,
    check.names = FALSE
  )
  cell_data <- data.frame(
    cell_id = cell_ids,
    barcode = barcode,
    sample_id = sample_id,
    donor_id = NA_character_,
    batch_id = NA_character_,
    condition_id = sample_id,
    entity_type = entity_type,
    stringsAsFactors = FALSE,
    row.names = cell_ids,
    check.names = FALSE
  )
  unique_samples <- unique(sample_id)
  sample_data <- unique(data.frame(
    sample_id = unique_samples,
    donor_id = NA_character_,
    batch_id = NA_character_,
    condition_id = unique_samples,
    tissue = NA_character_,
    platform = modality,
    modality = modality,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  rownames(sample_data) <- sample_data$sample_id

  feature_data_path <- file.path(meta.dir, paste0(modality, "_feature_data"))
  cell_data_path <- file.path(meta.dir, "cell_data")
  sample_data_path <- file.path(meta.dir, "sample_data")
  project_info_path <- file.path(meta.dir, "project_info.json")
  manifest_path <- file.path(soda_root, "soda.json")
  command_log_path <- file.path(meta.dir, "commands.json")
  params_path <- file.path(meta.dir, "params.json")
  version_path <- file.path(meta.dir, "version.json")

  feature_data_path <- save_table(feature_data, feature_data_path)
  cell_data_path <- save_table(cell_data, cell_data_path)
  sample_data_path <- save_table(sample_data, sample_data_path)

  assay_ref <- list(
    aggregation_type = "single",
    handle = handle
  )

  project_info <- list(
    project_id = make.names(project_name, unique = TRUE),
    project_name = project_name,
    description = NA_character_,
    organism = organism,
    genome_build = genome_build,
    technologies = modality,
    date_created = as.character(Sys.time()),
    package_version = .ic2_package_version()
  )
  save_json(project_info, project_info_path)

  commands <- list(
    command_log = list(
      make.obj = list(
        function_name = "make.obj",
        timestamp = as.character(Sys.time()),
        arguments = list(
          out.dir = out.dir, project_name = project_name, modality = modality,
          assay.name = assay.name, storage.backend = storage.backend,
          overwrite = overwrite
        ),
        outputs = list(
          modality = modality, assay = assay.name,
          n_features = length(feature_ids),
          n_cells = length(cell_ids),
          n_samples = length(unique_samples)
        )
      )
    ),
    timestamps = as.character(Sys.time()),
    arguments = list(),
    outputs = list()
  )
  save_json(commands, command_log_path)

  params <- list(
    import = list(
      modality = modality, assay.name = assay.name,
      n_samples = length(unique_samples),
      storage.backend = storage.backend
    ),
    qc = list(), normalization = list(), feature_selection = list(),
    dimensionality_reduction = list(), graph_construction = list(),
    clustering = list(), integration = list(), visualization = list(),
    sketch = list(), misc = list()
  )
  save_json(params, params_path)

  version <- list(
    object_version = "0.2.0",
    schema_version = "1.0.0",
    package_version = .ic2_package_version(),
    soda_version = "1.0.0"
  )
  save_json(version, version_path)

  manifest <- list(
    format = "SODA",
    format_name = "Single Omics Data Archive",
    soda_version = "1.0.0",
    schema_version = "1.0.0",
    root_dir = soda_root,
    modality = modality,
    assay = assay.name,
    assay_path = handle$path,
    handle_only = TRUE,
    external_assay_sources = if (identical(handle$backend, "tenx_h5") && !isTRUE(copy.raw)) {
      handle$path
    } else {
      character(0)
    },
    metadata_paths = list(
      project_info = project_info_path,
      feature_data = feature_data_path,
      cell_data = cell_data_path,
      sample_data = sample_data_path,
      commands = command_log_path,
      params = params_path,
      version = version_path
    )
  )
  ## soda.json written after object build via .ic2_soda_write_manifest()
  unlink(file.path(soda_root, "manifest.rds"), force = TRUE)

  rep <- .ic2_progress_bar(rep, 2L, 3L, prefix = "  make.obj")
  .ic2_step(rep, "[3/3] Building handle-only S4 object")
  rep <- .ic2_progress_bar(rep, 3L, 3L, prefix = "  make.obj")
  rep <- .ic2_add_output(rep, "SODA", soda_root)
  rep <- .ic2_add_output(rep, "Log", file.path(soda_root, "iCellR2.log"))
  extra <- paste0(
    "Features: ", .ic2_fmt_num(length(feature_ids)),
    "   |   Cells: ", .ic2_fmt_num(length(cell_ids)),
    "   |   Samples: ", length(unique_samples)
  )
  if (identical(handle$backend, "tenx_h5") && !isTRUE(copy.raw)) {
    extra <- c(
      extra,
      paste0("External tenx_h5 source (keep path stable): ", handle$path)
    )
  } else if (identical(handle$backend, "tenx_h5") && isTRUE(copy.raw)) {
    extra <- c(extra, "Raw assay: copied into SODA (self-contained)")
  }
  .ic2_finish(rep, extra_lines = extra)

  obj <- build_object(
    project_info = project_info,
    sample_data = sample_data,
    cell_data = cell_data,
    feature_data = feature_data,
    assay_ref = assay_ref,
    n_features = length(feature_ids),
    n_cells = length(cell_ids),
    n_samples = length(unique_samples),
    sample_ids = unique_samples,
    paths = list(
      metadata_paths = list(
        project_info = project_info_path,
        feature_data = feature_data_path,
        cell_data = cell_data_path,
        sample_data = sample_data_path,
        manifest = manifest_path,
        commands = command_log_path,
        params = params_path,
        version = version_path
      ),
      assay_paths = list(
        modality = modality,
        assay = assay.name,
        filepath = handle$path,
        backend = handle$backend,
        group = handle$group
      )
    ),
    commands = commands,
    params = params,
    version = version,
    schema_notes = "Handle-only multimodal iCellR2 project with SODA + HDF5 sparse backend."
  )
  .ic2_soda_write_manifest(obj)
  obj
}
