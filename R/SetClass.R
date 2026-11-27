## S4 class definition for iCellR2
setClass(
  "iCellR2",
  slots = c(
    project_info = "list",
    samples = "list",
    cell_data = "data.frame",
    modalities = "list",
    multimodal = "list",
    reductions = "list",
    graphs = "list",
    storage = "list",
    commands = "list",
    params = "list",
    version = "list",
    misc = "list",
    validity = "list",
    my.freq = "data.frame"
  )
)

## Internal helpers for show()
.ic2_show_utf8 <- function() {
  if (isTRUE(getOption("iCellR2.show.utf8"))) return(TRUE)
  isTRUE(l10n_info()[["UTF-8"]])
}

.ic2_show_hline <- function(style = c("heavy", "light"), width = 48L) {
  style <- match.arg(style)
  if (.ic2_show_utf8()) {
    ch <- if (identical(style, "heavy")) "=" else "-"
  } else {
    ch <- if (identical(style, "heavy")) "=" else "-"
  }
  paste(rep(ch, width), collapse = "")
}

.ic2_show_num <- function(x) {
  if (is.null(x) || length(x) < 1L) return("NA")
  x <- as.integer(x)
  out <- format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  out[is.na(x)] <- "NA"
  out
}

.ic2_show_mark <- function(done, partial = FALSE) {
  if (.ic2_show_utf8()) {
    if (isTRUE(done)) return("[x]")
    if (isTRUE(partial)) return("[~]")
    return("[ ]")
  }
  if (isTRUE(done)) return("[x]")
  if (isTRUE(partial)) return("[~]")
  "[ ]"
}

.ic2_show_assay_path <- function(assay_ref) {
  if (is.null(assay_ref)) return(NA_character_)
  if (is.list(assay_ref) && !is.null(assay_ref$handle$path)) {
    return(as.character(assay_ref$handle$path))
  }
  if (is.list(assay_ref) && identical(assay_ref$aggregation_type, "virtual")) {
    paths <- assay_ref$sample_assay_paths
    if (is.null(paths) && !is.null(assay_ref$sample_handles)) {
      paths <- lapply(assay_ref$sample_handles, function(h) h$path)
    }
    if (!is.null(paths) && length(paths) > 0L) {
      if (length(paths) == 1L) return(as.character(paths[[1L]]))
      return(paste0(as.character(paths[[1L]]), "  (+", length(paths) - 1L, " more)"))
    }
  }
  if (inherits(assay_ref, "iCellR2_assay_handle") ||
      (is.list(assay_ref) && identical(assay_ref$type, "iCellR2_assay_handle"))) {
    return(as.character(assay_ref$path))
  }
  NA_character_
}

.ic2_show_pipeline_row <- function(mark, step, detail, note = "") {
  ## step padded to 12, detail to 22
  line <- sprintf("  %s %-12s %-22s %s", mark, step, detail, note)
  message(trimws(line, which = "right"))
}

## Custom show() method for iCellR2
setMethod(
  "show",
  "iCellR2",
  function(object) {
    hl <- .ic2_show_hline("heavy")
    sl <- .ic2_show_hline("light")

    message(hl)
    message("   _  ____     _ _ ____  ____")
    message("  (_)/ ___|___| | |  _ \\|___ \\")
    message("  | | |   / _ \\ | | |_) | __) |")
    message("  | | |__|  __/ | |  _ < / __/")
    message("  |_|\\____\\___|_|_|_| \\_\\_____|")
    message("")
    message("  Massive-scale single-cell analysis")
    message("  with a familiar R interface")
    message(hl)

    pkg.ver <- object@version$package_version %||% NA_character_
    schema.ver <- object@version$schema_version %||% NA_character_
    soda.ver <- object@version$soda_version %||% NA_character_
    proj.name <- object@project_info$project_name %||% NA_character_
    root.dir <- object@storage$root_dir %||% NA_character_
    backend <- object@storage$backend %||% NA_character_
    assay.backend <- object@storage$assay_backend %||% "hdf5_sparse"
    handle.only <- isTRUE(object@misc$handle_only)

    n.features <- object@misc$n_features %||% {
      if (nrow(object@cell_data) >= 0 && !is.null(object@misc$initialized_modality)) {
        fd <- object@modalities[[object@misc$initialized_modality]]$feature_data
        if (!is.null(fd)) nrow(fd) else 0L
      } else 0L
    }
    n.cells <- object@misc$n_cells %||% nrow(object@cell_data)
    n.samples <- object@misc$n_samples %||% {
      sd <- object@samples$sample_data
      if (!is.null(sd)) nrow(sd) else 0L
    }

    init.mod <- object@misc$initialized_modality %||% {
      mods <- names(object@modalities)
      mods <- mods[mods != "misc_modalities"]
      hit <- mods[vapply(mods, function(m) {
        !is.null(object@modalities[[m]]$active_assay)
      }, logical(1))]
      if (length(hit)) hit[[1L]] else NA_character_
    }
    active.assay <- NA_character_
    if (!is.na(init.mod) && !is.null(object@modalities[[init.mod]])) {
      active.assay <- object@modalities[[init.mod]]$active_assay %||%
        object@misc$initialized_assay %||% NA_character_
    }

    message(sprintf("  %-11s %s", "Project", proj.name))
    message(sprintf(
      "  %-11s %s   |   schema %s   |   SODA %s",
      "Version", pkg.ver, schema.ver, soda.ver
    ))
    message(sprintf(
      "  %-11s %s / %s%s",
      "Backend", backend, assay.backend,
      if (handle.only) "   |   handle-only" else ""
    ))
    message(sprintf("  %-11s %s", "Archive", root.dir))
    message("")
    message(sprintf("  %-11s %s%s", "Modality", init.mod,
                    if (!is.na(init.mod)) " (active)" else ""))
    message(sprintf("  %-11s %s%s", "Assay", active.assay,
                    if (!is.na(active.assay)) " (active)" else ""))
    message(sprintf("  %-11s %s", "Cells", .ic2_show_num(n.cells)))
    message(sprintf("  %-11s %s", "Features", .ic2_show_num(n.features)))

    sample_line <- .ic2_show_num(n.samples)
    if (nrow(object@cell_data) > 0L && "condition_id" %in% colnames(object@cell_data)) {
      conds <- table(object@cell_data$condition_id)
      if (length(conds) > 0L) {
        sample_line <- paste0(
          sample_line, "   |   ",
          paste(sprintf("%s (%s)", names(conds), .ic2_show_num(as.integer(conds))),
                collapse = ", ")
        )
      }
    }
    message(sprintf("  %-11s %s", "Samples", sample_line))

    ## ---- Pipeline status ----
    message(sl)
    message("  Pipeline")
    message(sl)

    mod <- init.mod
    assays <- list()
    if (!is.na(mod) && !is.null(object@modalities[[mod]]$assays)) {
      assays <- object@modalities[[mod]]$assays
    }
    assay.names <- names(assays)

    import.done <- length(assay.names) > 0L || !is.null(object@commands$command_log$make.obj)
    import.detail <- if (length(assay.names) > 0L) {
      ## prefer raw-ish name
      cand <- intersect(c("counts", "raw", "assay"), assay.names)
      if (length(cand)) cand[[1L]] else assay.names[[1L]]
    } else " - "
    import.note <- if (import.done) {
      n0 <- object@commands$command_log$make.obj$outputs$n_cells %||% n.cells
      paste0(.ic2_show_num(n0), " cells")
    } else "not run"

    qc.cols <- c("nGenes", "UMIs", "mito.percent")
    qc.done <- nrow(object@cell_data) > 0L && all(qc.cols %in% colnames(object@cell_data))
    qc.note <- if (qc.done) "mito / nGenes / UMIs" else "not run"

    filter.assays <- assay.names[grepl("_filtered$", assay.names)]
    filter.done <- length(filter.assays) > 0L ||
      any(grepl("^qc\\..*_filter$", names(object@params$qc %||% list()))) ||
      any(grepl("^cell\\.filter", names(object@commands$command_log %||% list())))
    filter.detail <- if (length(filter.assays)) filter.assays[[length(filter.assays)]] else " - "
    filter.note <- if (filter.done) {
      fs <- object@misc$my.filters
      if (!is.null(fs$n_before) && !is.null(fs$n_after)) {
        paste0(.ic2_show_num(fs$n_before), "  ->  ", .ic2_show_num(fs$n_after))
      } else {
        "completed"
      }
    } else "not run"

    down.assays <- assay.names[grepl("_downsampled$", assay.names)]
    down.done <- length(down.assays) > 0L || !is.null(object@params$misc$down.sample)
    down.detail <- if (length(down.assays)) down.assays[[length(down.assays)]] else " - "
    down.note <- if (down.done) {
      ds <- object@params$misc$down.sample
      if (!is.null(ds$group.by)) {
        paste0("balanced by ", ds$group.by)
      } else "completed"
    } else "not run"

    norm.params <- object@params$normalization[[mod]]
    norm.done <- !is.null(norm.params) ||
      ("data" %in% assay.names) ||
      any(grepl("^norm\\.data", names(object@commands$command_log %||% list())))
    norm.detail <- if (!is.null(norm.params$new.assay.name)) {
      norm.params$new.assay.name
    } else if ("data" %in% assay.names) {
      "data"
    } else " - "
    norm.note <- if (norm.done) {
      method <- norm.params$norm.method %||% "normalized"
      if (!is.null(norm.params$top.rank) && !is.na(norm.params$top.rank)) {
        paste0(method, " (top ", norm.params$top.rank, ")")
      } else {
        as.character(method)
      }
    } else "not run"

    sketch.done <- FALSE
    sketch.detail <- " - "
    if (!is.na(mod) && !is.null(object@modalities[[mod]]$sketches)) {
      sk <- object@modalities[[mod]]$sketches
      if (!is.null(sk$sketch_ids) && length(sk$sketch_ids) > 0L) {
        sketch.done <- TRUE
        sketch.detail <- paste0(length(sk$sketch_ids), " cells")
      } else if (!is.null(sk$sketch_models) && length(sk$sketch_models) > 0L) {
        sketch.done <- TRUE
        sketch.detail <- paste(names(sk$sketch_models), collapse = ", ")
      }
    }
    sketch.note <- if (sketch.done) "registered" else "not run"

    pca.done <- !is.null(object@reductions$integrated_pca) ||
      (!is.na(mod) && !is.null(object@modalities[[mod]]$reductions$pca))
    umap.done <- !is.null(object@reductions$integrated_umap) ||
      (!is.na(mod) && !is.null(object@modalities[[mod]]$reductions$umap))
    tsne.done <- !is.null(object@reductions$integrated_tsne) ||
      (!is.na(mod) && !is.null(object@modalities[[mod]]$reductions$tsne))
    knetl.done <- !is.null(object@reductions$integrated_knetl) ||
      (!is.na(mod) && !is.null(object@modalities[[mod]]$reductions$knetl))
    reduce.done <- pca.done || umap.done || tsne.done || knetl.done
    reduce.bits <- c(
      if (pca.done) "PCA",
      if (umap.done) "UMAP",
      if (knetl.done) "KNetL"
    )
    reduce.detail <- if (reduce.done) paste(reduce.bits, collapse = " / ") else "PCA / UMAP"
    reduce.note <- if (reduce.done) {
      bits <- character(0)
      opt.pcs <- object@misc$opt.pcs %||%
        object@params$pca[[mod]]$opt.pcs %||%
        object@modalities[[mod]]$reductions$pca$opt.pcs
      if (!is.null(opt.pcs) && length(opt.pcs) == 1L && !is.na(opt.pcs)) {
        bits <- c(bits, paste0("opt.pcs=", opt.pcs))
      }
      if (umap.done) {
        um <- object@modalities[[mod]]$reductions$umap
        if (isTRUE(um$sketch_mode)) bits <- c(bits, "UMAP sketch -> transform")
        else if (!is.null(um)) bits <- c(bits, "UMAP")
      }
      if (length(bits)) paste(bits, collapse = "  |  ") else "completed"
    } else {
      "not run"
    }

    clustering.done <- "cluster_id" %in% colnames(object@cell_data)
    n.clusters <- if (clustering.done) {
      length(unique(object@cell_data$cluster_id))
    } else 0L
    cluster.note <- if (clustering.done) {
      paste0(.ic2_show_num(n.clusters), " clusters")
    } else "not run"

    .ic2_show_pipeline_row(.ic2_show_mark(import.done), "Import", import.detail, import.note)
    .ic2_show_pipeline_row(.ic2_show_mark(qc.done), "QC", if (qc.done) "qc.stats" else " - ", qc.note)
    .ic2_show_pipeline_row(.ic2_show_mark(filter.done), "Filter", filter.detail, filter.note)
    .ic2_show_pipeline_row(.ic2_show_mark(down.done), "Downsample", down.detail, down.note)
    .ic2_show_pipeline_row(.ic2_show_mark(norm.done), "Normalize", norm.detail, norm.note)

    ## Gene stats / model / scale
    gstats.done <- FALSE
    gmodel.done <- FALSE
    gstats.note <- "not run"
    gmodel.note <- "not run"
    scale.detail <- " - "
    scale.note <- "not run"
    scale.done <- FALSE
    gmodel <- NULL

    if (!is.na(mod)) {
      gstats.done <- !is.null(object@modalities[[mod]]$gene.data) ||
        (!is.null(object@storage$metadata_paths$gene_stats) &&
           file.exists(object@storage$metadata_paths$gene_stats)) ||
        any(grepl("^gene\\.stats", names(object@commands$command_log %||% list())))
      if (gstats.done) {
        gd <- object@modalities[[mod]]$gene.data
        gstats.note <- if (!is.null(gd)) {
          paste0(.ic2_show_num(length(unique(gd$genes[gd$condition == "all"]))), " genes")
        } else "completed"
      }

      gmodel <- object@modalities[[mod]]$gene.model %||% object@misc$gene.model
      gmodel.done <- !is.null(gmodel) && length(gmodel) > 0L
      gmodel.note <- if (gmodel.done) {
        paste0(.ic2_show_num(length(gmodel)), " genes")
      } else "not run"

      scale.assays <- assay.names[grepl("^scale", assay.names) | assay.names %in% c("scale.data", "scaled")]
      scale.done <- length(scale.assays) > 0L ||
        !is.null(object@params$normalization[[paste0(mod, "_scale")]]) ||
        any(grepl("^data\\.scale", names(object@commands$command_log %||% list())))
      scale.detail <- if (length(scale.assays)) {
        scale.assays[[length(scale.assays)]]
      } else {
        object@params$normalization[[paste0(mod, "_scale")]]$new.assay.name %||% " - "
      }
      scale.note <- if (scale.done) {
        sp <- object@params$normalization[[paste0(mod, "_scale")]]
        sp$method %||% "completed"
      } else "not run"
    }

    .ic2_show_pipeline_row(.ic2_show_mark(gstats.done), "Gene stats", if (gstats.done) "gene.stats" else " - ", gstats.note)
    .ic2_show_pipeline_row(
      .ic2_show_mark(gmodel.done), "Gene model",
      if (gmodel.done) "make.gene.model" else " - ", gmodel.note
    )
    .ic2_show_pipeline_row(.ic2_show_mark(scale.done), "Scale", scale.detail, scale.note)

    .ic2_show_pipeline_row(.ic2_show_mark(sketch.done), "Sketch", sketch.detail, sketch.note)
    .ic2_show_pipeline_row(.ic2_show_mark(reduce.done), "Reduce", reduce.detail, reduce.note)
    .ic2_show_pipeline_row(
      .ic2_show_mark(clustering.done), "Cluster",
      if (clustering.done) "cluster_id" else " - ", cluster.note
    )

    ## ---- Assays on disk ----
    message(sl)
    message("  Assays on disk")
    message(sl)
    if (length(assay.names) < 1L) {
      message("  (none)")
    } else {
      for (an in assay.names) {
        path <- .ic2_show_assay_path(assays[[an]])
        active.tag <- if (!is.na(active.assay) && identical(an, active.assay)) "  <- active" else ""
        if (is.na(path) || !nzchar(path)) {
          message(sprintf("  %-22s (in-object / path unavailable)%s", an, active.tag))
        } else {
          ## shorten long paths a bit for display
          disp <- path
          if (nchar(disp) > 64L) {
            disp <- paste0("...", substr(disp, nchar(disp) - 60L, nchar(disp)))
          }
          message(sprintf("  %-22s %s%s", an, disp, active.tag))
        }
      }
    }

    ## ---- Modalities ----
    message(sl)
    message("  Modalities")
    message(sl)
    available.mods <- names(object@modalities)
    available.mods <- available.mods[available.mods != "misc_modalities"]
    for (m in available.mods) {
      mod.obj <- object@modalities[[m]]
      if (!is.null(mod.obj$assays) && length(mod.obj$assays) > 0L) {
        message(sprintf("  %-10s assays: %s", m, paste(names(mod.obj$assays), collapse = ", ")))
      } else {
        message(sprintf("  %-10s assays: none", m))
      }
    }

    mm <- tryCatch(object@multimodal, error = function(e) NULL)
    if (is.list(mm)) {
      mm_names <- names(mm)
      mm_names <- mm_names[vapply(mm_names, function(nm) {
        z <- mm[[nm]]
        !is.null(z) && length(z) > 0L && !identical(nm, "misc")
      }, logical(1))]
      if (length(mm_names) > 0L) {
        message(sl)
        message("  Multimodal")
        message(sl)
        message("  ", paste(mm_names, collapse = ", "))
      }
    }

    if (!is.null(object@validity$messages)) {
      message(sl)
      message("  Status      ", object@validity$messages)
    }

    message(hl)
  }
)
