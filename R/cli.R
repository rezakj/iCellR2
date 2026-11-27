#' iCellR2 command-line interface (bwa / samtools style)
#'
#' Invoked by \code{exec/iCellR2}. Operates on SODA paths  -  each command is a
#' short R process that reads/writes the archive on disk.
#'
#' @param args Character vector of CLI args (default: \code{commandArgs(TRUE)}).
#' @return Invisibly \code{NULL}; exits with status via \code{quit()} when run
#'   as a script (not while interactive).
#' @keywords internal
.ic2_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  args <- as.character(args %||% character(0))
  if (length(args) && identical(args[[1L]], "--")) args <- args[-1L]

  if (!length(args) || args[[1L]] %in% c("-h", "--help", "help")) {
    .ic2_cli_top_help()
    return(invisible(.ic2_cli_exit(0L)))
  }
  if (args[[1L]] %in% c("-v", "--version", "version")) {
    cat(.ic2_cli_version_line(), "\n", sep = "")
    return(invisible(.ic2_cli_exit(0L)))
  }

  cmd <- args[[1L]]
  rest <- if (length(args) > 1L) args[-1L] else character(0)
  handlers <- .ic2_cli_handlers()
  if (!cmd %in% names(handlers)) {
    cat("iCellR2: unknown command '", cmd, "'\n\n", sep = "")
    .ic2_cli_top_help()
    return(invisible(.ic2_cli_exit(1L)))
  }
  tryCatch(
    {
      handlers[[cmd]](rest)
      invisible(.ic2_cli_exit(0L))
    },
    error = function(e) {
      cat("iCellR2 ", cmd, ": ", conditionMessage(e), "\n", sep = "")
      invisible(.ic2_cli_exit(1L))
    }
  )
}

.ic2_cli_exit <- function(status = 0L) {
  status <- as.integer(status)[1L]
  if (interactive() || isTRUE(getOption("iCellR2.cli.noquit"))) {
    return(invisible(status))
  }
  quit(save = "no", status = status, runLast = FALSE)
}

.ic2_cli_version_line <- function() {
  ver <- tryCatch(.ic2_package_version(), error = function(e) "?.?.?")
  paste0("iCellR2 ", ver)
}

.ic2_cli_top_help <- function() {
  ver <- .ic2_cli_version_line()
  cat(
    "\n",
    "Program: iCellR2 (large-scale single-cell analysis  |  SODA archive)\n",
    "Version: ", sub("^iCellR2 ", "", ver), "\n",
    "\n",
    "Usage:   iCellR2 <command> [options]\n",
    "\n",
    "Command:\n",
    "  Import / project\n",
    "         load10x           load 10x Genomics counts (dir or .h5)\n",
    "         data.aggregation  combine samples into a virtual aggregation\n",
    "         make.obj          create handle-only project + SODA folder\n",
    "\n",
    "  QC / normalize / features\n",
    "         qc.stats          QC metrics into SODA cell metadata\n",
    "         stats.plot        QC stats plots (png/html)\n",
    "         cell.filter       filter cells by mito / genes / UMIs\n",
    "         down.sample       downsample cells by group\n",
    "         hto.anno          demultiplex hashtag oligos (HTO)\n",
    "         norm.data         normalize assay\n",
    "         gene.stats        gene-level stats\n",
    "         make.gene.model   feature / HVG model (+ optional plot)\n",
    "\n",
    "  Reduce / cluster / trajectory\n",
    "         run.pca           PCA\n",
    "         run.umap          UMAP\n",
    "         run.knetl         KNetL map (kNN graph)\n",
    "         run.cgumap        cluster-guided UMAP\n",
    "         iclust            graph clustering\n",
    "         run.impute        kNN smooth impute (streamed)\n",
    "         iba               batch alignment CPCA / CCCA\n",
    "         run.pseudotime    graph-DPT + soft branches\n",
    "\n",
    "  Markers / plots\n",
    "         findMarkers       marker genes (writes TSV)\n",
    "         i.score           gene-set scores (streamed)\n",
    "         add.adt           CITE-seq ADT counts (own modality)\n",
    "         norm.adt          per-cell CLR of ADT counts\n",
    "         adt.rna.merge     joint RNA+ADT PCA (multimodal/rna_adt)\n",
    "         add.vdj           V(D)J contigs onto the RNA object\n",
    "         vdj.stats         clonotype expansion (TCR and BCR)\n",
    "         clono.plot        highlight a clonotype on the UMAP\n",
    "         cc                cell cycle S / G2M / Phase (streamed)\n",
    "         run.tsne          tSNE on existing PCA\n",
    "         add.visium        Visium / Visium HD positions onto the RNA object\n",
    "         add.xenium        Xenium centroids (microns); transcripts stay on disk\n",
    "         spatial.plot      points, or a raster for Visium HD\n",
    "         spatial.graph     hex, square, kNN, or micron radius\n",
    "         spatial.genes     Moran's I on the spatial graph\n",
    "         spatial.transcripts  one Xenium gene, streamed\n",
    "         merge.peaks       one peak set shared by all samples\n",
    "         add.atac          ATAC peaks, with or without RNA\n",
    "         norm.atac         TF-IDF of the peak assay\n",
    "         gene.activity     sum peaks into genes\n",
    "         add.multiome      RNA + ATAC (rna.input and atac.input)\n",
    "         atac.rna.merge    joint RNA+ATAC PCs (multimodal/rna_atac)\n",
    "         top.markers       top DE genes from a markers table\n",
    "         run.diff.exp      two-group differential expression\n",
    "         volcano.ma.plot   volcano / MA from a DE table\n",
    "         clust.avg.exp     cluster mean expression\n",
    "         clust.stats.plot  mito / gene / UMI box or violin by cluster\n",
    "         clust.cond.info   condition composition (pie / bar / table)\n",
    "         cluster.plot      PCA / UMAP / cgUMAP / KNetL / tSNE plots\n",
    "         pseudotime.arrows stream / arrow overlays\n",
    "         heatmap.gg.plot   gene heatmap\n",
    "         bubble.gg.plot    gene bubble plot\n",
    "         gene.plot         feature scatter / violin / box / bar\n",
    "\n",
    "  SODA utilities\n",
    "         open              open a SODA folder (summary)\n",
    "         validate          validate a SODA folder\n",
    "         refresh           rebuild inventory / index.html / manifest\n",
    "         install-cli       symlink launcher onto PATH (~/bin)\n",
    "         version           print version\n",
    "         help              this message\n",
    "\n",
    "Note: The project is a SODA folder on disk, not a language-private object.\n",
    "      Typical start:\n",
    "        iCellR2 make.obj -input raw/PBMC -sample.id pbmc -out.dir iCellR2\n",
    "        iCellR2 qc.stats -soda iCellR2.soda\n",
    "        iCellR2 stats.plot -soda iCellR2.soda -out.name plot1_QC -format png\n",
    "\n",
    "      Run `iCellR2 <command> -h` for command-specific options.\n",
    "\n",
    sep = ""
  )
}

.ic2_cli_handlers <- function() {
  list(
    "load10x" = .ic2_cli_cmd_load10x,
    "data.aggregation" = .ic2_cli_cmd_data_aggregation,
    "make.obj" = .ic2_cli_cmd_make_obj,
    "qc.stats" = .ic2_cli_cmd_qc_stats,
    "stats.plot" = .ic2_cli_cmd_stats_plot,
    "cell.filter" = .ic2_cli_cmd_cell_filter,
    "down.sample" = .ic2_cli_cmd_down_sample,
    "hto.anno" = .ic2_cli_cmd_hto_anno,
    "hto.demux" = .ic2_cli_cmd_hto_demux,
    "norm.data" = .ic2_cli_cmd_norm_data,
    "gene.stats" = .ic2_cli_cmd_gene_stats,
    "make.gene.model" = .ic2_cli_cmd_make_gene_model,
    "run.pca" = .ic2_cli_cmd_run_pca,
    "run.umap" = .ic2_cli_cmd_run_umap,
    "run.knetl" = .ic2_cli_cmd_run_knetl,
    "run.cgumap" = .ic2_cli_cmd_run_cgumap,
    "iclust" = .ic2_cli_cmd_iclust,
    "run.impute" = .ic2_cli_cmd_run_impute,
    "iba" = .ic2_cli_cmd_iba,
    "run.pseudotime" = .ic2_cli_cmd_run_pseudotime,
    "findMarkers" = .ic2_cli_cmd_find_markers,
    "i.score" = .ic2_cli_cmd_i_score,
    "add.adt" = .ic2_cli_cmd_add_adt,
    "norm.adt" = .ic2_cli_cmd_norm_adt,
    "adt.rna.merge" = .ic2_cli_cmd_adt_rna_merge,
    "add.vdj" = .ic2_cli_cmd_add_vdj,
    "vdj.stats" = .ic2_cli_cmd_vdj_stats,
    "clono.plot" = .ic2_cli_cmd_clono_plot,
    "cc" = .ic2_cli_cmd_cc,
    "cell.cycle" = .ic2_cli_cmd_cc,
    "run.tsne" = .ic2_cli_cmd_run_tsne,
    "add.visium" = .ic2_cli_cmd_add_image,
    "add.10x.image" = .ic2_cli_cmd_add_image,
    "add.xenium" = .ic2_cli_cmd_add_xenium,
    "spatial.plot" = .ic2_cli_cmd_spatial_plot,
    "spatial.graph" = .ic2_cli_cmd_spatial_graph,
    "spatial.genes" = .ic2_cli_cmd_spatial_genes,
    "spatial.transcripts" = .ic2_cli_cmd_spatial_tx,
    "merge.peaks" = .ic2_cli_cmd_merge_peaks,
    "add.atac" = .ic2_cli_cmd_add_atac,
    "norm.atac" = .ic2_cli_cmd_norm_atac,
    "gene.activity" = .ic2_cli_cmd_gene_activity,
    "add.multiome" = .ic2_cli_cmd_add_multiome,
    "atac.rna.merge" = .ic2_cli_cmd_atac_rna_merge,
    "top.markers" = .ic2_cli_cmd_top_markers,
    "run.diff.exp" = .ic2_cli_cmd_run_diff_exp,
    "volcano.ma.plot" = .ic2_cli_cmd_volcano,
    "clust.avg.exp" = .ic2_cli_cmd_clust_avg,
    "clust.stats.plot" = .ic2_cli_cmd_clust_stats_plot,
    "clust.cond.info" = .ic2_cli_cmd_clust_cond_info,
    "cluster.plot" = .ic2_cli_cmd_cluster_plot,
    "pseudotime.arrows" = .ic2_cli_cmd_pseudotime_arrows,
    "heatmap.gg.plot" = .ic2_cli_cmd_heatmap,
    "bubble.gg.plot" = .ic2_cli_cmd_bubble,
    "gene.plot" = .ic2_cli_cmd_gene_plot,
    "open" = .ic2_cli_cmd_open,
    "validate" = .ic2_cli_cmd_validate,
    "refresh" = .ic2_cli_cmd_refresh,
    "install-cli" = .ic2_cli_cmd_install,
    "help" = function(args) { .ic2_cli_top_help(); invisible(NULL) },
    "version" = function(args) {
      cat(.ic2_cli_version_line(), "\n", sep = "")
      invisible(NULL)
    }
  )
}

## ---- parsing helpers -------------------------------------------------------

.ic2_cli_parse <- function(args, defaults = list()) {
  out <- defaults
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      out[["__help__"]] <- TRUE
      i <- i + 1L
      next
    }
    if (startsWith(a, "--") && grepl("=", a, fixed = TRUE)) {
      kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1L]]
      key <- kv[[1L]]
      val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
      out <- .ic2_cli_store(out, key, val)
      i <- i + 1L
      next
    }
    if (startsWith(a, "-") && nchar(a) > 1L) {
      key <- sub("^-+", "", a)
      if (i == length(args) || startsWith(args[[i + 1L]], "-")) {
        out <- .ic2_cli_store(out, key, TRUE)
        i <- i + 1L
      } else {
        out <- .ic2_cli_store(out, key, args[[i + 1L]])
        i <- i + 2L
      }
      next
    }
    stop("Unexpected argument: ", a, call. = FALSE)
  }
  out
}

.ic2_cli_store <- function(out, key, val) {
  if (!is.null(out[[key]])) {
    prev <- out[[key]]
    out[[key]] <- c(as.character(prev), as.character(val))
  } else {
    out[[key]] <- val
  }
  out
}

.ic2_cli_flag <- function(opts, name, default = NULL) {
  if (!is.null(opts[[name]])) opts[[name]] else default
}

.ic2_cli_flag1 <- function(opts, name, default = NULL) {
  v <- .ic2_cli_flag(opts, name, default)
  if (is.null(v)) return(default)
  as.character(v)[1L]
}

.ic2_cli_as_logical <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(isTRUE(x[[1L]]))
  x <- tolower(as.character(x)[1L])
  if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

.ic2_cli_as_num <- function(x, default = NULL) {
  if (is.null(x)) return(default)
  as.numeric(as.character(x)[1L])
}

.ic2_cli_as_int <- function(x, default = NULL) {
  if (is.null(x)) return(default)
  as.integer(as.numeric(as.character(x)[1L]))
}

.ic2_cli_as_dims <- function(x, default = 1:10) {
  if (is.null(x)) return(default)
  eval(parse(text = as.character(x)[1L]))
}

.ic2_cli_as_format <- function(x, default = "png") {
  if (is.null(x)) return(default)
  v <- unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE)
  trimws(v)
}

.ic2_cli_genes <- function(opts) {
  g <- .ic2_cli_flag(opts, "gene", NULL)
  gf <- .ic2_cli_flag1(opts, "gene.file", NULL)
  genes <- character(0)
  if (!is.null(g)) {
    genes <- c(genes, unlist(strsplit(paste(as.character(g), collapse = ","), ",", fixed = TRUE)))
  }
  if (!is.null(gf) && nzchar(gf) && file.exists(gf)) {
    genes <- c(genes, scan(gf, what = character(), quiet = TRUE, sep = "\n"))
  }
  genes <- unique(trimws(genes))
  genes[nzchar(genes)]
}

.ic2_cli_stage_dir <- function() {
  d <- file.path(getwd(), ".icellr2_stage")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

.ic2_cli_stage_file <- function(name = "last_sample.rds") {
  file.path(.ic2_cli_stage_dir(), name)
}

.ic2_cli_resolve_soda <- function(opts) {
  soda <- .ic2_cli_flag1(opts, "soda", NULL)
  if (is.null(soda) || !nzchar(soda)) stop("missing required -soda PATH", call. = FALSE)
  soda
}

.ic2_cli_open <- function(opts, verbose = NULL) {
  path <- .ic2_cli_resolve_soda(opts)
  if (is.null(verbose)) verbose <- .ic2_cli_as_logical(opts$verbose, FALSE)
  open.soda(path, verbose = verbose)
}

.ic2_cli_mutate_soda <- function(opts, fun, label) {
  path <- .ic2_cli_resolve_soda(opts)
  obj <- open.soda(path, verbose = FALSE)
  obj <- fun(obj)
  cat("[iCellR2] ", label, " OK\n", sep = "")
  cat("  SODA: ", normalizePath(path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(obj)
}

.ic2_cli_save_gg <- function(plot_obj, out, width = 8, height = 8, dpi = 300) {
  if (is.null(out) || !nzchar(as.character(out)[1L])) {
    print(plot_obj)
    return(invisible(NULL))
  }
  out <- as.character(out)[1L]
  ext <- tolower(tools::file_ext(out))
  if (!nzchar(ext)) {
    out <- paste0(out, ".png")
    ext <- "png"
  }
  if (ext %in% c("png", "pdf", "jpg", "jpeg", "tiff")) {
    if (ext == "png") {
      grDevices::png(out, width = width, height = height, units = "in", res = as.integer(dpi))
    } else if (ext == "pdf") {
      grDevices::pdf(out, width = width, height = height)
    } else if (ext %in% c("jpg", "jpeg")) {
      grDevices::jpeg(out, width = width, height = height, units = "in", res = as.integer(dpi))
    } else {
      grDevices::tiff(out, width = width, height = height, units = "in", res = as.integer(dpi))
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot_obj)
    cat("  wrote: ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
  } else {
    stop("Unsupported plot extension: ", ext, call. = FALSE)
  }
  invisible(out)
}

.ic2_cli_help <- function(...) {
  cat("\n", paste0(..., collapse = ""), "\n", sep = "")
  invisible(NULL)
}

## ---- commands: import ------------------------------------------------------

.ic2_cli_cmd_load10x <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    input.format = "auto", output = "virtual", modality = "rna",
    feature.column = "2", feature.field = "name", overwrite = FALSE, verbose = TRUE,
    bin.size = "8"
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 load10x [options]\n\n",
      "Options: -input PATH          10x directory or .h5 (required)\n",
      "         -sample.id ID        sample id (required)\n",
      "         -input.format FMT    auto|10x_dir|10x_h5|icellr2_h5\n",
      "         -output MODE         virtual|10x_hdf5|sparse|matrix|data.frame\n",
      "         -modality NAME       [rna]\n",
      "         -feature.column N    [2]\n",
      "         -feature.field NAME  name|id|feature_type\n",
      "         -h5.file PATH        optional output h5\n",
      "         -feature.type TYPE   optional multimodal filter\n",
      "         -bin.size N          Visium HD bin um: 8 (default), 16, or 2\n",
      "         -out FILE            stage RDS [.icellr2_stage/last_sample.rds]\n",
      "         -overwrite -verbose BOOL\n"
    ))
  }
  input <- .ic2_cli_flag1(opts, "input")
  sample.id <- .ic2_cli_flag1(opts, "sample.id")
  if (is.null(input) || !nzchar(input)) stop("missing -input PATH", call. = FALSE)
  if (is.null(sample.id) || !nzchar(sample.id)) stop("missing -sample.id ID", call. = FALSE)
  out <- .ic2_cli_flag1(opts, "out", .ic2_cli_stage_file())
  obj <- load10x(
    input = input,
    input.format = .ic2_cli_flag1(opts, "input.format", "auto"),
    output = .ic2_cli_flag1(opts, "output", "virtual"),
    sample.id = sample.id,
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    feature.column = .ic2_cli_as_int(opts$feature.column, 2L),
    feature.field = .ic2_cli_flag1(opts, "feature.field", "name"),
    h5.file = .ic2_cli_flag1(opts, "h5.file", NULL),
    feature.type = .ic2_cli_flag1(opts, "feature.type", NULL),
    bin.size = .ic2_cli_as_int(opts$bin.size, 8L),
    overwrite = .ic2_cli_as_logical(opts$overwrite, FALSE),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  saveRDS(obj, out)
  cat("[iCellR2] load10x OK\n  staged: ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(NULL)
}

.ic2_cli_cmd_data_aggregation <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", output = "virtual", overwrite = FALSE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 data.aggregation [options]\n\n",
      "Options: -sample PATH_OR_NAME=PATH   repeatable staged RDS or NAME=10x_path\n",
      "         -modality NAME              [rna]\n",
      "         -output MODE                virtual|sparse|matrix|data.frame|hdf5\n",
      "         -condition.names A,B,...    optional\n",
      "         -out FILE                   stage RDS [.icellr2_stage/agg.rds]\n",
      "         -overwrite -verbose BOOL\n\n",
      "Example:\n",
      "  iCellR2 load10x -input raw/WT -sample.id WT -out .icellr2_stage/WT.rds\n",
      "  iCellR2 load10x -input raw/KO -sample.id KO -out .icellr2_stage/KO.rds\n",
      "  iCellR2 data.aggregation -sample .icellr2_stage/WT.rds -sample .icellr2_stage/KO.rds\n"
    ))
  }
  samp <- .ic2_cli_flag(opts, "sample", NULL)
  if (is.null(samp)) stop("provide one or more -sample PATH (or NAME=PATH)", call. = FALSE)
  samp <- as.character(samp)
  samples <- list()
  for (s in samp) {
    if (grepl("=", s, fixed = TRUE)) {
      kv <- strsplit(s, "=", fixed = TRUE)[[1L]]
      nm <- kv[[1L]]
      path <- paste(kv[-1L], collapse = "=")
      if (grepl("\\.rds$", path, ignore.case = TRUE) && file.exists(path)) {
        samples[[nm]] <- readRDS(path)
      } else {
        samples[[nm]] <- load10x(
          input = path, sample.id = nm,
          modality = .ic2_cli_flag1(opts, "modality", "rna"),
          verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
        )
      }
    } else {
      if (!file.exists(s)) stop("sample not found: ", s, call. = FALSE)
      obj <- readRDS(s)
      nm <- if (!is.null(obj$sample.id)) as.character(obj$sample.id)[1L] else tools::file_path_sans_ext(basename(s))
      samples[[nm]] <- obj
    }
  }
  cond <- .ic2_cli_flag1(opts, "condition.names", NULL)
  conds <- if (!is.null(cond) && nzchar(cond)) trimws(strsplit(cond, ",", fixed = TRUE)[[1L]]) else NULL
  agg <- data.aggregation(
    samples = samples,
    condition.names = conds,
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    output = .ic2_cli_flag1(opts, "output", "virtual"),
    overwrite = .ic2_cli_as_logical(opts$overwrite, FALSE),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  out <- .ic2_cli_flag1(opts, "out", .ic2_cli_stage_file("agg.rds"))
  saveRDS(agg, out)
  cat("[iCellR2] data.aggregation OK\n  staged: ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(NULL)
}

.ic2_cli_cmd_make_obj <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    project_name = "iCellR2_project", modality = "rna", assay.name = "counts",
    overwrite = FALSE, copy.raw = TRUE, storage.backend = "hdf5_sparse",
    input.format = "auto", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 make.obj [options]\n\n",
      "Options: -out.dir STEM          ->  STEM.soda (required)\n",
      "         -project_name NAME\n",
      "         -modality NAME        [rna]\n",
      "         -assay.name NAME      [counts]\n",
      "         -organism NAME -genome_build NAME\n",
      "         -input PATH           load10x then make.obj\n",
      "         -sample.id ID         required with -input\n",
      "         -input.format FMT     [auto]\n",
      "         -x FILE               staged sample/aggregation RDS\n",
      "         -overwrite -copy.raw BOOL -storage.backend NAME -verbose BOOL\n"
    ))
  }
  out.dir <- .ic2_cli_flag1(opts, "out.dir")
  if (is.null(out.dir) || !nzchar(out.dir)) stop("missing -out.dir STEM", call. = FALSE)

  x <- NULL
  input <- .ic2_cli_flag1(opts, "input", NULL)
  xfile <- .ic2_cli_flag1(opts, "x", .ic2_cli_stage_file())
  if (!is.null(input) && nzchar(input)) {
    sample.id <- .ic2_cli_flag1(opts, "sample.id")
    if (is.null(sample.id) || !nzchar(sample.id)) {
      stop("-sample.id required with -input", call. = FALSE)
    }
    x <- load10x(
      input = input,
      input.format = .ic2_cli_flag1(opts, "input.format", "auto"),
      sample.id = sample.id,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      overwrite = .ic2_cli_as_logical(opts$overwrite, FALSE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  } else if (file.exists(xfile)) {
    x <- readRDS(xfile)
  } else {
    agg <- .ic2_cli_stage_file("agg.rds")
    if (file.exists(agg)) {
      x <- readRDS(agg)
    } else {
      stop("provide -input/-sample.id, or -x FILE, or run load10x / data.aggregation first", call. = FALSE)
    }
  }

  obj <- make.obj(
    x = x,
    out.dir = out.dir,
    project_name = .ic2_cli_flag1(opts, "project_name", "iCellR2_project"),
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    assay.name = .ic2_cli_flag1(opts, "assay.name", "counts"),
    organism = {
      o <- .ic2_cli_flag1(opts, "organism", NA_character_)
      if (is.null(o) || !nzchar(o) || identical(o, "NA")) NA_character_ else o
    },
    genome_build = {
      o <- .ic2_cli_flag1(opts, "genome_build", NA_character_)
      if (is.null(o) || !nzchar(o) || identical(o, "NA")) NA_character_ else o
    },
    overwrite = .ic2_cli_as_logical(opts$overwrite, FALSE),
    copy.raw = .ic2_cli_as_logical(opts$copy.raw, TRUE),
    storage.backend = .ic2_cli_flag1(opts, "storage.backend", "hdf5_sparse"),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  root <- obj@storage$root_dir %||% paste0(out.dir, ".soda")
  cat("[iCellR2] make.obj OK\n  SODA: ", normalizePath(root, winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(NULL)
}

## ---- QC / norm -------------------------------------------------------------

.ic2_cli_cmd_qc_stats <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", chunk.size = "200000", use.rcpp = TRUE, update.object = TRUE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 qc.stats [options]\n\n",
      "Options: -soda PATH -modality NAME -assay.name NAME\n",
      "         -chunk.size N -use.rcpp BOOL -update.object BOOL -verbose BOOL\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    qc.stats(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", NULL),
      chunk.size = .ic2_cli_as_int(opts$chunk.size, 200000L),
      use.rcpp = .ic2_cli_as_logical(opts$use.rcpp, TRUE),
      update.object = .ic2_cli_as_logical(opts$update.object, TRUE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "qc.stats")
}

.ic2_cli_cmd_stats_plot <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    out.name = "stats_plot", format = "png", min.genes = "200", min.umis = "500",
    max.mito = "0.1", group.by = "condition_id", width = "12", height = "9",
    dpi = "300", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 stats.plot [options]\n\n",
      "Options: -soda PATH -out.name PREFIX -format png|html|png,html\n",
      "         -min.genes N -min.umis N -max.mito N -group.by COL\n",
      "         -width N -height N -dpi N -verbose BOOL\n"
    ))
  }
  obj <- .ic2_cli_open(opts)
  stats.plot(
    obj,
    out.name = .ic2_cli_flag1(opts, "out.name", "stats_plot"),
    format = .ic2_cli_as_format(opts$format, "png"),
    min.genes = .ic2_cli_as_num(opts$min.genes, 200),
    min.umis = .ic2_cli_as_num(opts$min.umis, 500),
    max.mito = .ic2_cli_as_num(opts$max.mito, 0.1),
    group.by = .ic2_cli_flag1(opts, "group.by", "condition_id"),
    width = .ic2_cli_as_num(opts$width, 12),
    height = .ic2_cli_as_num(opts$height, 9),
    dpi = .ic2_cli_as_int(opts$dpi, 300L),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  cat("[iCellR2] stats.plot OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_cell_filter <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", min.mito = "0", max.mito = "1", min.genes = "0",
    max.genes = "Inf", min.umis = "0", max.umis = "Inf",
    update.assay = FALSE, set.active = TRUE, overwrite = TRUE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 cell.filter [options]\n\n",
      "Options: -soda PATH -modality NAME -assay.name NAME\n",
      "         -min.mito N -max.mito N -min.genes N -max.genes N\n",
      "         -min.umis N -max.umis N -update.assay BOOL -verbose BOOL\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    cell.filter(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", NULL),
      min.mito = .ic2_cli_as_num(opts$min.mito, 0),
      max.mito = .ic2_cli_as_num(opts$max.mito, 1),
      min.genes = .ic2_cli_as_num(opts$min.genes, 0),
      max.genes = {
        mg <- .ic2_cli_flag1(opts, "max.genes", "Inf")
        if (identical(mg, "Inf")) Inf else as.numeric(mg)
      },
      min.umis = .ic2_cli_as_num(opts$min.umis, 0),
      max.umis = {
        mu <- .ic2_cli_flag1(opts, "max.umis", "Inf")
        if (identical(mu, "Inf")) Inf else as.numeric(mu)
      },
      update.assay = .ic2_cli_as_logical(opts$update.assay, FALSE),
      set.active = .ic2_cli_as_logical(opts$set.active, TRUE),
      overwrite = .ic2_cli_as_logical(opts$overwrite, TRUE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "cell.filter")
}

.ic2_cli_cmd_down_sample <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", group.by = "condition_id", update.assay = FALSE,
    set.active = TRUE, overwrite = TRUE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 down.sample [options]\n\n",
      "Options: -soda PATH -group.by COL -n.cells N -modality NAME -verbose BOOL\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    down.sample(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      group.by = .ic2_cli_flag1(opts, "group.by", "condition_id"),
      n.cells = .ic2_cli_as_int(opts$n.cells, NULL),
      update.assay = .ic2_cli_as_logical(opts$update.assay, FALSE),
      set.active = .ic2_cli_as_logical(opts$set.active, TRUE),
      overwrite = .ic2_cli_as_logical(opts$overwrite, TRUE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "down.sample")
}

.ic2_cli_cmd_hto_anno <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", assay.name = NULL,
    hto.features = NULL,
    hto.pattern = "^CMO|CMO|TotalSeq|Hashtag|hashtag|HTO",
    feature.type = "Antibody Capture,Multiplexing Capture",
    cov.thr = "10", assignment.thr = "80",
    write.tsv = TRUE, plot = FALSE, plot.format = "png", max.plot.points = "25000",
    out.name = "hto_anno", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 hto.anno [options]\n\n",
      "Options: -soda PATH -cov.thr N -assignment.thr N -hto.pattern REGEX\n",
      "         -hto.features F1,F2 -feature.type NAME -plot BOOL\n",
      "         -plot.format png|pdf -max.plot.points N -out.name STEM\n",
      "         -write.tsv BOOL -verbose BOOL\n",
      "Notes:   Annotate only. Filter + split with: iCellR2 hto.demux\n",
      "         CellPlex: -hto.pattern ^CMO -feature.type NULL\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    feats <- .ic2_cli_flag1(opts, "hto.features", NULL)
    if (!is.null(feats) && nzchar(feats)) {
      feats <- strsplit(feats, ",", fixed = TRUE)[[1L]]
      feats <- trimws(feats)
    } else {
      feats <- NULL
    }
    ftype <- .ic2_cli_flag1(opts, "feature.type", "Antibody Capture,Multiplexing Capture")
    if (identical(ftype, "NULL") || identical(ftype, "none") || !nzchar(ftype)) {
      ftype <- NULL
    } else {
      ftype <- trimws(strsplit(ftype, ",", fixed = TRUE)[[1L]])
    }
    hto.anno(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", NULL),
      hto.features = feats,
      hto.pattern = .ic2_cli_flag1(opts, "hto.pattern", "^CMO|CMO|TotalSeq|Hashtag|hashtag|HTO"),
      feature.type = ftype,
      cov.thr = .ic2_cli_as_num(opts$cov.thr, 10),
      assignment.thr = .ic2_cli_as_num(opts$assignment.thr, 80),
      write.tsv = .ic2_cli_as_logical(opts$write.tsv, TRUE),
      out.name = .ic2_cli_flag1(opts, "out.name", "hto_anno"),
      plot = .ic2_cli_as_logical(opts$plot, FALSE),
      plot.format = .ic2_cli_flag1(opts, "plot.format", "png"),
      max.plot.points = .ic2_cli_as_int(opts$max.plot.points, 25000L),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "hto.anno")
}

.ic2_cli_cmd_hto_demux <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", assay.name = NULL,
    min.percent = "70", min.coverage = NULL,
    keep.good.only = FALSE, drop.low.cov = FALSE,
    assignments = NULL, subset.cells = FALSE,
    write = TRUE, out.dir = NULL, out.format = "10x_dir",
    prefix = "RNAs",
    adt.pattern = "TotalFF|TotalSeq|adt_|ADT",
    hto.pattern = NULL,
    write.tsv = TRUE, out.name = "hto_demux", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 hto.demux [options]\n\n",
      "Options: -soda PATH -min.percent N -min.coverage N\n",
      "         -keep.good.only BOOL -drop.low.cov BOOL -subset.cells BOOL\n",
      "         -write BOOL -out.dir DIR -out.format 10x_dir|10x_hdf5|h5\n",
      "         -prefix RNAs -assignments TAG1,TAG2 -adt.pattern REGEX\n",
      "         -write.tsv BOOL -out.name STEM -verbose BOOL\n",
      "Notes:   Run hto.anno first. Classic filter: -min.percent 70\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    asg <- .ic2_cli_flag1(opts, "assignments", NULL)
    if (!is.null(asg) && nzchar(asg)) {
      asg <- trimws(strsplit(asg, ",", fixed = TRUE)[[1L]])
    } else {
      asg <- NULL
    }
    hto.demux(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", NULL),
      min.percent = .ic2_cli_as_num(opts$min.percent, 70),
      min.coverage = .ic2_cli_as_num(opts$min.coverage, NULL),
      keep.good.only = .ic2_cli_as_logical(opts$keep.good.only, FALSE),
      drop.low.cov = .ic2_cli_as_logical(opts$drop.low.cov, FALSE),
      assignments = asg,
      subset.cells = .ic2_cli_as_logical(opts$subset.cells, FALSE),
      write = .ic2_cli_as_logical(opts$write, TRUE),
      out.dir = .ic2_cli_flag1(opts, "out.dir", NULL),
      out.format = .ic2_cli_flag1(opts, "out.format", "10x_dir"),
      prefix = .ic2_cli_flag1(opts, "prefix", "RNAs"),
      adt.pattern = .ic2_cli_flag1(opts, "adt.pattern", "TotalFF|TotalSeq|adt_|ADT"),
      hto.pattern = .ic2_cli_flag1(opts, "hto.pattern", NULL),
      write.tsv = .ic2_cli_as_logical(opts$write.tsv, TRUE),
      out.name = .ic2_cli_flag1(opts, "out.name", "hto_demux"),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "hto.demux")
}

.ic2_cli_cmd_norm_data <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", new.assay.name = "data", norm.method = "global.glsf",
    top.rank = "500", materialize = TRUE, overwrite = TRUE, set.active = TRUE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 norm.data [options]\n\n",
      "Options: -soda PATH -norm.method NAME -top.rank N -assay.name NAME\n",
      "         -new.assay.name NAME -materialize BOOL -modality NAME -verbose BOOL\n",
      "Methods: global.glsf|rpm|LogNormalize|ranked.glsf|deseq|spike.in|no.norm\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    norm.data(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", NULL),
      new.assay.name = .ic2_cli_flag1(opts, "new.assay.name", "data"),
      norm.method = .ic2_cli_flag1(opts, "norm.method", "global.glsf"),
      top.rank = .ic2_cli_as_int(opts$top.rank, 500L),
      materialize = .ic2_cli_as_logical(opts$materialize, TRUE),
      overwrite = .ic2_cli_as_logical(opts$overwrite, TRUE),
      set.active = .ic2_cli_as_logical(opts$set.active, TRUE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "norm.data")
}

.ic2_cli_cmd_gene_stats <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", assay.name = "data", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 gene.stats [options]\n\n",
      "Options: -soda PATH -modality NAME -assay.name NAME -verbose BOOL\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    gene.stats(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", "data"),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "gene.stats")
}

.ic2_cli_cmd_make_gene_model <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", gate.mode = "classic", dispersion.limit = "1.5",
    base.mean.rank = "1000", gene.num.max = "2000", my.out.put = "both",
    out.name = "plot", format = "html", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 make.gene.model [options]\n\n",
      "Options: -soda PATH -dispersion.limit N -base.mean.rank N -gene.num.max N\n",
      "         -gate.mode classic|mean -out.name PREFIX -format html|png|NULL\n",
      "         -my.out.put data|plot|both -modality NAME -verbose BOOL\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    fmt <- .ic2_cli_flag1(opts, "format", "html")
    if (identical(fmt, "NULL") || identical(fmt, "none")) fmt <- NULL
    make.gene.model(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      gate.mode = .ic2_cli_flag1(opts, "gate.mode", "classic"),
      dispersion.limit = .ic2_cli_as_num(opts$dispersion.limit, 1.5),
      base.mean.rank = .ic2_cli_as_int(opts$base.mean.rank, 1000L),
      gene.num.max = .ic2_cli_as_int(opts$gene.num.max, 2000L),
      my.out.put = .ic2_cli_flag1(opts, "my.out.put", "both"),
      out.name = .ic2_cli_flag1(opts, "out.name", "plot"),
      format = fmt,
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "make.gene.model")
}

## ---- reduce / cluster ------------------------------------------------------

.ic2_cli_cmd_run_pca <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", method = "gene.model", npcs = "50", transform = "auto",
    assay.name = NULL, use.sketch = TRUE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 run.pca [options]\n\n",
      "Options: -soda PATH -method gene.model|base.mean.rank|features\n",
      "         -npcs N -transform auto|none|log2p1|log1p -top.rank N\n",
      "         -use.sketch BOOL -modality NAME -assay.name NAME -verbose BOOL\n",
      "Notes:   ATAC: -modality atac -assay.name tfidf. Stored under\n",
      "         modalities/atac/reductions/pca, not reductions/pca.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    run.pca(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", NULL),
      method = .ic2_cli_flag1(opts, "method", "gene.model"),
      npcs = .ic2_cli_as_int(opts$npcs, 50L),
      transform = .ic2_cli_flag1(opts, "transform", "auto"),
      top.rank = .ic2_cli_as_int(opts$top.rank, 500L),
      use.sketch = .ic2_cli_as_logical(opts$use.sketch, TRUE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "run.pca")
}

.ic2_cli_cmd_run_umap <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", dims = "1:10", use = "auto", use.sketch = FALSE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 run.umap [options]\n\n",
      "Options: -soda PATH -dims SPEC -use auto|pca|aligned_pca\n",
      "         -use.sketch BOOL -modality NAME -verbose BOOL\n",
      "Notes:   use=auto reads aligned_pca when it exists, otherwise pca.\n",
      "         Aligned UMAP is written to reductions/aligned_umap.\n"
    ))
  }
  dims <- .ic2_cli_as_dims(opts$dims, 1:10)
  .ic2_cli_mutate_soda(opts, function(obj) {
    run.umap(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      use = .ic2_cli_flag1(opts, "use", "auto"),
      dims = dims,
      use.sketch = .ic2_cli_as_logical(opts$use.sketch, FALSE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "run.umap")
}

.ic2_cli_cmd_run_knetl <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", zoom = "300", dims = "1:20", use = "auto",
    data.type = "pca", joint = FALSE, add.3d = FALSE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 run.knetl [options]\n\n",
      "Options: -soda PATH -zoom N -dims SPEC -use auto|pca|aligned_pca\n",
      "         -data.type pca|umap -joint BOOL -add.3d BOOL -modality NAME\n",
      "Notes:   kNN neighborhood (no n x n distance) embedded with UMAP.\n",
      "         use=auto prefers aligned_pca. Writes reductions/knetl\n",
      "         or reductions/aligned_knetl. Default zoom=300, dims=1:20.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    run.knetl(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      zoom = .ic2_cli_as_int(opts$zoom, 300L),
      dims = .ic2_cli_as_dims(opts$dims, 1:20),
      use = .ic2_cli_flag1(opts, "use", "auto"),
      data.type = .ic2_cli_flag1(opts, "data.type", "pca"),
      joint = .ic2_cli_as_logical(opts$joint, FALSE),
      add.3d = .ic2_cli_as_logical(opts$add.3d, FALSE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "run.knetl")
}

.ic2_cli_cmd_run_cgumap <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", dims = "1:3", scale = TRUE, cluster.weight = "1",
    cluster.col = "cluster_id", center = TRUE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 run.cgumap [options]\n\n",
      "Options: -soda PATH -dims SPEC -cluster.col COL -cluster.weight N\n",
      "         -scale BOOL -center BOOL -modality NAME -verbose BOOL\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    run.cgumap(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      dims = .ic2_cli_as_dims(opts$dims, 1:3),
      scale = .ic2_cli_as_logical(opts$scale, TRUE),
      cluster.weight = .ic2_cli_as_num(opts$cluster.weight, 1),
      cluster.col = .ic2_cli_flag1(opts, "cluster.col", "cluster_id"),
      center = .ic2_cli_as_logical(opts$center, TRUE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "run.cgumap")
}

.ic2_cli_cmd_iclust <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    data.type = "pca", use = "auto", dims = "1:10", sensitivity = "30", use.sketch = FALSE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 iclust [options]\n\n",
      "Options: -soda PATH -data.type pca|umap -use auto|pca|aligned_pca\n",
      "         -dims SPEC -sensitivity N -use.sketch BOOL -verbose BOOL\n",
      "Notes:   use=auto prefers aligned_pca / aligned_umap when present.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    iclust(
      obj,
      data.type = .ic2_cli_flag1(opts, "data.type", "pca"),
      use = .ic2_cli_flag1(opts, "use", "auto"),
      dims = .ic2_cli_as_dims(opts$dims, 1:10),
      sensitivity = .ic2_cli_as_num(opts$sensitivity, 30),
      use.sketch = .ic2_cli_as_logical(opts$use.sketch, FALSE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "iclust")
}

.ic2_cli_cmd_run_impute <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", reduction = "pca", use = "auto", nn = "10", dims = "1:10",
    chunk.cells = "30000", overwrite = TRUE, set.active = FALSE, verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 run.impute [options]\n\n",
      "Options: -soda PATH -nn N -dims SPEC -reduction pca|umap -use auto|pca|aligned_pca\n",
      "         -assay.name NAME -new.assay.name NAME -chunk.cells N\n",
      "         -modality NAME -overwrite BOOL -set.active BOOL -verbose BOOL\n",
      "Notes:   Streamed iCellR.imp kNN mean  ->  assay \"imputed\" (requires run.pca).\n",
      "         Default dims=1:10, nn=10, chunk.cells=30000.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    dims <- .ic2_cli_as_dims(opts$dims, 1:10)
    run.impute(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", NULL),
      new.assay.name = .ic2_cli_flag1(opts, "new.assay.name", NULL),
      reduction = .ic2_cli_flag1(opts, "reduction", "pca"),
      use = .ic2_cli_flag1(opts, "use", "auto"),
      nn = .ic2_cli_as_int(opts$nn, 10L),
      dims = dims,
      chunk.cells = .ic2_cli_as_int(opts$chunk.cells, 30000L),
      overwrite = .ic2_cli_as_logical(opts$overwrite, TRUE),
      set.active = .ic2_cli_as_logical(opts$set.active, FALSE),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "run.impute")
}

.ic2_cli_cmd_iba <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", ba.method = "CPCA", k = "10", dims = "1:30",
    refit = FALSE, condition.by = "condition_id",
    method = "base.mean.rank", top.rank = "500",
    chunk.cells = "30000", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 iba [options]\n\n",
      "Options: -soda PATH -k N -dims SPEC\n",
      "         -ba.method CPCA|CCCA  (default CPCA)\n",
      "         -condition.by condition_id -refit BOOL\n",
      "Notes:   CPCA averages neighbor PCs; CCCA streams neighbor-mean\n",
      "         expression then PCA. Writes reductions/aligned_pca (pca unchanged).\n",
      "         Default k=10, dims=1:30. Reuses PCA on disk unless -refit TRUE.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    iba(
      obj,
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      ba.method = .ic2_cli_flag1(opts, "ba.method", "CPCA"),
      k = .ic2_cli_as_int(opts$k, 10L),
      dims = .ic2_cli_as_dims(opts$dims, 1:30),
      refit = .ic2_cli_as_logical(opts$refit, FALSE),
      condition.by = .ic2_cli_flag1(opts, "condition.by", "condition_id"),
      method = .ic2_cli_flag1(opts, "method", "base.mean.rank"),
      top.rank = .ic2_cli_as_int(opts$top.rank, 500L),
      chunk.cells = .ic2_cli_as_int(opts$chunk.cells, 30000L),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "iba")
}

.ic2_cli_cmd_run_pseudotime <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    soft.branches = TRUE, n.tips = "2", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 run.pseudotime [options]\n\n",
      "Options: -soda PATH -root.cluster ID (default: largest cluster)\n",
      "         -soft.branches BOOL -n.tips N -verbose BOOL\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    root <- .ic2_cli_flag1(opts, "root.cluster", NULL)
    if (is.null(root) || !nzchar(root)) {
      cl_tab <- sort(table(obj@cell_data$cluster_id), decreasing = TRUE)
      root <- names(cl_tab)[1L]
    }
    run.pseudotime(
      obj,
      root.cluster = root,
      soft.branches = .ic2_cli_as_logical(opts$soft.branches, TRUE),
      n.tips = .ic2_cli_as_int(opts$n.tips, 2L),
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "run.pseudotime")
}

## ---- markers / plots -------------------------------------------------------

.ic2_cli_cmd_add_adt <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(input = NULL, assay.name = "counts"))
  if (isTRUE(opts[["__help__"]]) || is.null(opts$input) || !nzchar(as.character(opts$input)[1L])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 add.adt -soda PATH -input adt.tsv [options]\n\n",
      "Options: -assay.name counts\n",
      "Notes:   TSV is antibodies x cells (rownames = antibodies).\n",
      "         Writes modalities/adt/assays/counts. RNA is unchanged.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    add.adt(
      obj,
      adt.data = as.character(opts$input)[1L],
      assay.name = .ic2_cli_flag1(opts, "assay.name", "counts")
    )
  }, "add.adt")
}

.ic2_cli_cmd_norm_adt <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(assay.name = "counts", new.assay.name = "data"))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 norm.adt -soda PATH [options]\n\n",
      "Options: -assay.name counts -new.assay.name data\n",
      "Notes:   Per-cell CLR. Writes modalities/adt/assays/data.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    norm.adt(
      obj,
      assay.name = .ic2_cli_flag1(opts, "assay.name", "counts"),
      new.assay.name = .ic2_cli_flag1(opts, "new.assay.name", "data")
    )
  }, "norm.adt")
}

.ic2_cli_cmd_adt_rna_merge <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    rna.dims = "30", adt.dims = "10", adt.assay = "data"
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 adt.rna.merge -soda PATH [options]\n\n",
      "Options: -rna.dims 30 -adt.dims 10 -adt.assay data\n",
      "Notes:   Needs run.pca and norm.adt. Writes multimodal/rna_adt.\n",
      "         Does not replace reductions/pca or the RNA assay.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    adt.rna.merge(
      obj,
      rna.dims = .ic2_cli_as_int(opts$rna.dims, 30L),
      adt.dims = .ic2_cli_as_int(opts$adt.dims, 10L),
      adt.assay = .ic2_cli_flag1(opts, "adt.assay", "data")
    )
  }, "adt.rna.merge")
}

.ic2_cli_cmd_add_vdj <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(input = NULL, cond.name = NULL, productive.only = TRUE))
  if (isTRUE(opts[["__help__"]]) || is.null(opts$input) || !nzchar(as.character(opts$input)[1L])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 add.vdj -soda PATH -input filtered_contig_annotations.csv\n\n",
      "Options: -cond.name SAMPLE -productive.only TRUE|FALSE\n",
      "Notes:   10x CSV/TSV or AIRR. Writes modalities/vdj and cell_data clonotype columns.\n",
      "         Barcodes are matched, not rewritten. Frequency is cells, not rows/2.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    add.vdj(
      obj,
      input = as.character(opts$input)[1L],
      cond.name = opts$cond.name,
      productive.only = .ic2_cli_as_logical(opts$productive.only, TRUE)
    )
  }, "add.vdj")
}

.ic2_cli_cmd_vdj_stats <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(plot = FALSE, out.name = NULL))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 vdj.stats -soda PATH [options]\n\n",
      "Options: -plot TRUE|FALSE -out.name stats.png\n",
      "Notes:   Writes modalities/vdj/stats.tsv. TCR and BCR bins.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    vdj.stats(
      obj,
      plot = .ic2_cli_as_logical(opts$plot, FALSE),
      out.name = opts$out.name
    )
  }, "vdj.stats")
}

.ic2_cli_cmd_clono_plot <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    clono = NULL, plot.type = "umap", max.cells = "200000", out.name = NULL, interactive = FALSE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 clono.plot -soda PATH [-clono ID] [options]\n\n",
      "Options: -plot.type umap|pca|knetl -max.cells 200000 -out.name plot.png\n",
      "         -interactive TRUE  writes HTML of the sampled plot\n",
      "Notes:   Omit -clono to color by expansion bin. Highlighted cells are all kept.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    clono.plot(
      obj,
      clono = if (is.null(opts$clono) || !nzchar(as.character(opts$clono)[1L])) NULL else as.character(opts$clono),
      plot.type = .ic2_cli_flag1(opts, "plot.type", "umap"),
      max.cells = .ic2_cli_as_int(opts$max.cells, 200000L),
      interactive = .ic2_cli_as_logical(opts$interactive, FALSE),
      out.name = opts$out.name
    )
    obj
  }, "clono.plot")
}

.ic2_cli_cmd_cc <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(assay.name = "data"))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 cc -soda PATH [options]\n\n",
      "Options: -assay.name data\n",
      "Notes:   Streams S and G2M genes. Writes cell_data$Phase.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    cc(obj, assay.name = .ic2_cli_flag1(opts, "assay.name", "data"))
  }, "cc")
}

.ic2_cli_cmd_run_tsne <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    use = "auto", dims = "30", perplexity = "30", theta = "0.5", max_iter = "750"
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 run.tsne -soda PATH [options]\n\n",
      "Options: -use auto|pca|aligned_pca -dims 30 -perplexity 30 -theta 0.5 -max_iter 750\n",
      "Notes:   tSNE of the stored PCA. Does not read the gene matrix.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    nd <- .ic2_cli_as_int(opts$dims, 30L)
    run.tsne(
      obj,
      use = .ic2_cli_flag1(opts, "use", "auto"),
      dims = seq_len(nd),
      perplexity = .ic2_cli_as_num(opts$perplexity, 30),
      theta = .ic2_cli_as_num(opts$theta, 0.5),
      max_iter = .ic2_cli_as_int(opts$max_iter, 750L)
    )
  }, "run.tsne")
}

.ic2_cli_cmd_add_image <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(input = NULL, cond.name = NULL, bin.size = "8"))
  if (isTRUE(opts[["__help__"]]) || is.null(opts$input) || !nzchar(as.character(opts$input)[1L])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 add.visium -soda PATH -input spatial_dir [options]\n\n",
      "Options: -cond.name SAMPLE -bin.size 8\n",
      "Notes:   CSV or parquet positions. HD outs uses binned_outputs/square_*um.\n",
      "         Coordinates are also stored in microns (spatial_x_um).\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    add.visium(
      obj, input = as.character(opts$input)[1L],
      condition.names = opts$cond.name,
      bin.size = .ic2_cli_as_int(opts$bin.size, 8L)
    )
  }, "add.visium")
}

.ic2_cli_cmd_add_xenium <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(input = NULL, cond.name = NULL))
  if (isTRUE(opts[["__help__"]]) || is.null(opts$input) || !nzchar(as.character(opts$input)[1L])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 add.xenium -soda PATH -input xenium_outs [options]\n\n",
      "Options: -cond.name SAMPLE\n",
      "Notes:   cells.parquet centroids in microns. transcripts.parquet is not loaded.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    add.xenium(obj, input = as.character(opts$input)[1L], condition.names = opts$cond.name)
  }, "add.xenium")
}

.ic2_cli_cmd_spatial_plot <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    col.by = "clusters", gene = NULL, max.spots = "200000", out.name = NULL,
    mark = "auto", xlim = NULL, ylim = NULL
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 spatial.plot -soda PATH [options]\n\n",
      "Options: -col.by clusters|gene|cc|nGenes|UMIs|mt -gene GENE -out.name plot.png\n",
      "         -mark auto|point|raster -max.spots 200000\n",
      "         -xlim MIN,MAX -ylim MIN,MAX   crop; uses the hires image\n"
    ))
  }
  pair <- function(z) {
    if (is.null(z) || !nzchar(as.character(z)[1L])) return(NULL)
    as.numeric(strsplit(as.character(z)[1L], ",", fixed = TRUE)[[1L]])
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    spatial.plot(
      obj,
      col.by = .ic2_cli_flag1(opts, "col.by", "clusters"),
      gene = opts$gene,
      max.spots = .ic2_cli_as_int(opts$max.spots, 200000L),
      mark = .ic2_cli_flag1(opts, "mark", "auto"),
      xlim = pair(opts$xlim),
      ylim = pair(opts$ylim),
      out.name = opts$out.name
    )
    obj
  }, "spatial.plot")
}

.ic2_cli_cmd_spatial_graph <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(k = "6", method = "auto", radius.um = "30"))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 spatial.graph -soda PATH [options]\n\n",
      "Options: -method auto|hex|square|knn|radius -k 6 -radius.um 30\n",
      "Notes:   auto: hex Visium, square Visium HD, radius Xenium.\n",
      "         Skips in_tissue == 0. Writes graphs/spatial_knn/edges.tsv.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    spatial.graph(
      obj,
      method = .ic2_cli_flag1(opts, "method", "auto"),
      k = .ic2_cli_as_int(opts$k, 6L),
      radius.um = .ic2_cli_as_num(opts$radius.um, 30)
    )
  }, "spatial.graph")
}

.ic2_cli_cmd_spatial_tx <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    gene = NULL, max.points = "200000", out.name = NULL, xlim = NULL, ylim = NULL
  ))
  if (isTRUE(opts[["__help__"]]) || is.null(opts$gene) || !nzchar(as.character(opts$gene)[1L])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 spatial.transcripts -soda PATH -gene CD3E [options]\n\n",
      "Options: -max.points 200000 -out.name tx.png -xlim MIN,MAX -ylim MIN,MAX\n",
      "Notes:   Streams one gene from transcripts.parquet. Needs the arrow package.\n"
    ))
  }
  pair <- function(z) {
    if (is.null(z) || !nzchar(as.character(z)[1L])) return(NULL)
    as.numeric(strsplit(as.character(z)[1L], ",", fixed = TRUE)[[1L]])
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    spatial.transcripts(
      obj,
      gene = as.character(opts$gene)[1L],
      xlim = pair(opts$xlim),
      ylim = pair(opts$ylim),
      max.points = .ic2_cli_as_int(opts$max.points, 200000L),
      out.name = opts$out.name
    )
    obj
  }, "spatial.transcripts")
}

.ic2_cli_cmd_merge_peaks <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    input = NULL, gap = "0", min.width = "1", out = NULL
  ))
  inputs <- .ic2_cli_flag(opts, "input", NULL)
  if (isTRUE(opts[["__help__"]]) || is.null(inputs)) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 merge.peaks -input peaks_a.bed -input peaks_b.bed [options]\n\n",
      "Options: -gap 0 -min.width 1 -out consensus_peaks.tsv\n",
      "Notes:   Overlapping or book-ended peaks on one chromosome become one id\n",
      "         (chr_start_end). -gap -1 merges overlaps only. A second file,\n",
      "         *_map.tsv, links each input peak to that id. BED is 0-based.\n"
    ))
  }
  if (length(inputs) == 1L && grepl(",", inputs, fixed = TRUE)) {
    inputs <- strsplit(inputs, ",", fixed = TRUE)[[1L]]
  }
  shared <- merge.peaks(
    as.list(inputs),
    gap = .ic2_cli_as_int(opts$gap, 0L),
    min.width = .ic2_cli_as_int(opts$min.width, 1L)
  )
  out <- .ic2_cli_flag1(opts, "out", "consensus_peaks.tsv")
  map_out <- sub("(\\.[^.]+)$", "_map\\1", out)
  if (identical(map_out, out)) map_out <- paste0(out, "_map.tsv")
  utils::write.table(shared, out, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(attr(shared, "map"), map_out, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", nrow(shared), " consensus peaks to ", out)
  message("Wrote peak map to ", map_out)
  invisible(shared)
}

.ic2_cli_cmd_add_atac <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    input = NULL, sample.ids = NULL, peaks = NULL, assay.name = "counts"
  ))
  inputs <- .ic2_cli_flag(opts, "input", NULL)
  if (isTRUE(opts[["__help__"]]) || is.null(inputs)) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 add.atac -input ATAC_OUTS -sample.id ID -out.dir NAME\n",
      "         iCellR2 add.atac -soda PATH -input ATAC_OUTS\n\n",
      "Options: -input DIR  (repeat for each sample) -sample.ids a,b\n",
      "         -peaks consensus_peaks.tsv -out.dir NAME -overwrite\n",
      "Notes:   No -soda: ATAC only, same on-disk load as load10x.\n",
      "         With -soda: attach peaks to an existing object. RNA is unchanged.\n",
      "         Differing peak sets are streamed onto one shared set.\n"
    ))
  }
  if (length(inputs) == 1L && grepl(",", inputs, fixed = TRUE)) {
    inputs <- strsplit(inputs, ",", fixed = TRUE)[[1L]]
  }
  sids <- .ic2_cli_flag(opts, "sample.ids", NULL)
  if (!is.null(sids) && length(sids) == 1L && grepl(",", sids, fixed = TRUE)) {
    sids <- strsplit(sids, ",", fixed = TRUE)[[1L]]
  }
  soda <- .ic2_cli_flag1(opts, "soda", NULL)
  if (!is.null(soda) && nzchar(soda)) {
    return(.ic2_cli_mutate_soda(opts, function(obj) {
      add.atac(
        obj,
        input = as.character(inputs),
        sample.ids = sids,
        peaks = .ic2_cli_flag1(opts, "peaks", NULL),
        assay.name = .ic2_cli_flag1(opts, "assay.name", "counts"),
        overwrite = .ic2_cli_as_logical(opts$overwrite, FALSE)
      )
    }, "add.atac"))
  }
  add.atac(
    x = NULL,
    input = as.character(inputs),
    sample.ids = sids,
    peaks = .ic2_cli_flag1(opts, "peaks", NULL),
    assay.name = .ic2_cli_flag1(opts, "assay.name", "counts"),
    out.dir = .ic2_cli_flag1(opts, "out.dir", "iCellR2"),
    overwrite = .ic2_cli_as_logical(opts$overwrite, FALSE),
    verbose = TRUE
  )
}

.ic2_cli_cmd_norm_atac <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(assay.name = "counts", new.assay.name = "tfidf"))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 norm.atac -soda PATH [options]\n\n",
      "Options: -assay.name counts -new.assay.name tfidf\n",
      "Notes:   Signac TF-IDF method 1. Then:\n",
      "         iCellR2 run.pca -soda PATH -modality atac -assay.name tfidf\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    norm.atac(
      obj,
      assay.name = .ic2_cli_flag1(opts, "assay.name", "counts"),
      new.assay.name = .ic2_cli_flag1(opts, "new.assay.name", "tfidf")
    )
  }, "norm.atac")
}

.ic2_cli_cmd_gene_activity <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    genes = NULL, upstream = "2000", downstream = "0", assay.name = "counts"
  ))
  if (isTRUE(opts[["__help__"]]) || is.null(opts$genes)) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 gene.activity -soda PATH -genes genes.gtf [options]\n\n",
      "Options: -upstream 2000 -downstream 0 -assay.name counts\n",
      "Notes:   Sums peaks over each gene body plus upstream. Writes\n",
      "         modalities/atac/assays/gene.activity.\n"
    ))
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    gene.activity(
      obj,
      genes = as.character(opts$genes)[1L],
      upstream = .ic2_cli_as_int(opts$upstream, 2000L),
      downstream = .ic2_cli_as_int(opts$downstream, 0L),
      assay.name = .ic2_cli_flag1(opts, "assay.name", "counts")
    )
  }, "gene.activity")
}

.ic2_cli_cmd_add_multiome <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    input = NULL, rna.input = NULL, atac.input = NULL,
    sample.ids = NULL, sample.id = NULL, peaks = NULL, out.dir = "iCellR2"
  ))
  split_in <- function(z) {
    if (is.null(z)) return(NULL)
    if (length(z) == 1L && grepl(",", z, fixed = TRUE)) {
      z <- strsplit(as.character(z), ",", fixed = TRUE)[[1L]]
    }
    as.character(z)
  }
  rna <- split_in(.ic2_cli_flag(opts, "rna.input", NULL))
  atac <- split_in(.ic2_cli_flag(opts, "atac.input", NULL))
  both <- split_in(.ic2_cli_flag(opts, "input", NULL))
  if (isTRUE(opts[["__help__"]]) || (is.null(rna) && is.null(atac) && is.null(both))) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 add.multiome -rna.input RNA -atac.input ATAC -sample.id ID -out.dir NAME\n",
      "         iCellR2 add.multiome -input ARC_OUTS -sample.id ID -out.dir NAME\n\n",
      "Options: -sample.ids a,b -peaks consensus_peaks.tsv -soda PATH -overwrite\n",
      "Notes:   RNA + ATAC only. Inputs are load10x formats and stay on disk.\n",
      "         -input is both sides when one folder contains both matrices.\n"
    ))
  }
  sids <- .ic2_cli_flag(opts, "sample.ids", NULL)
  if (is.null(sids)) sids <- .ic2_cli_flag(opts, "sample.id", NULL)
  sids <- split_in(sids)
  soda <- .ic2_cli_flag1(opts, "soda", NULL)
  if (!is.null(soda) && nzchar(soda)) {
    return(.ic2_cli_mutate_soda(opts, function(obj) {
      add.multiome(
        obj,
        rna.input = rna,
        atac.input = atac,
        input = both,
        sample.ids = sids,
        peaks = .ic2_cli_flag1(opts, "peaks", NULL)
      )
    }, "add.multiome"))
  }
  add.multiome(
    x = NULL,
    rna.input = rna,
    atac.input = atac,
    input = both,
    sample.ids = sids,
    peaks = .ic2_cli_flag1(opts, "peaks", NULL),
    out.dir = .ic2_cli_flag1(opts, "out.dir", "iCellR2"),
    overwrite = .ic2_cli_as_logical(opts$overwrite, FALSE)
  )
}

.ic2_cli_cmd_atac_rna_merge <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(rna.dims = "1:30", atac.dims = "2:30"))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 atac.rna.merge -soda PATH [options]\n\n",
      "Options: -rna.dims 1:30 -atac.dims 2:30\n",
      "Notes:   Scaled concat of the two PCAs. Writes multimodal/rna_atac.\n",
      "         Does not replace reductions/pca. ATAC PC1 is dropped by default.\n"
    ))
  }
  dims <- function(z, default) {
    z <- as.character(z %||% default)[1L]
    if (grepl(":", z, fixed = TRUE)) {
      ab <- as.integer(strsplit(z, ":", fixed = TRUE)[[1L]])
      return(seq.int(ab[[1L]], ab[[2L]]))
    }
    as.integer(strsplit(z, ",", fixed = TRUE)[[1L]])
  }
  .ic2_cli_mutate_soda(opts, function(obj) {
    atac.rna.merge(
      obj,
      rna.dims = dims(opts$rna.dims, "1:30"),
      atac.dims = dims(opts$atac.dims, "2:30")
    )
  }, "atac.rna.merge")
}

.ic2_cli_cmd_spatial_genes <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(genes = NULL, assay.name = "data"))
  if (isTRUE(opts[["__help__"]]) || is.null(opts$genes) || !nzchar(as.character(opts$genes)[1L])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 spatial.genes -soda PATH -genes CD3E,MS4A1 [options]\n\n",
      "Options: -assay.name data\n",
      "Notes:   -genes is a comma-separated list or a file with one gene per line.\n",
      "         Needs spatial.graph. Writes metadata/spatial/svg.tsv.\n"
    ))
  }
  raw <- as.character(opts$genes)[1L]
  genes <- if (file.exists(raw)) raw else strsplit(raw, ",", fixed = TRUE)[[1L]]
  .ic2_cli_mutate_soda(opts, function(obj) {
    spatial.genes(obj, genes = genes, assay.name = .ic2_cli_flag1(opts, "assay.name", "data"))
  }, "spatial.genes")
}

.ic2_cli_cmd_i_score <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", assay.name = "data", method = "tirosh",
    genes = NULL, assignment.thr = "50", cov.thr = "10", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]]) || is.null(opts$genes) || !nzchar(as.character(opts$genes)[1L])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 i.score -soda PATH -genes sets.tsv [options]\n\n",
      "Options: -assay.name data|counts|imputed -method tirosh|mean|sum|zscore\n",
      "         -assignment.thr N -cov.thr N -modality NAME\n",
      "Notes:   sets.tsv has columns set and gene (one gene per row).\n",
      "         Scores stream from the assay. Writes cell_data$score_assignment\n",
      "         and metadata/scoring/scores.tsv.\n"
    ))
  }
  tab <- utils::read.delim(as.character(opts$genes)[1L], sep = "\t", stringsAsFactors = FALSE)
  if (!all(c("set", "gene") %in% colnames(tab))) {
    stop("Gene file must have columns: set, gene")
  }
  sets <- split(as.character(tab$gene), as.character(tab$set))
  .ic2_cli_mutate_soda(opts, function(obj) {
    i.score(
      obj,
      scoring.List = sets,
      scoring.method = .ic2_cli_flag1(opts, "method", "tirosh"),
      assay.name = .ic2_cli_flag1(opts, "assay.name", "data"),
      modality = .ic2_cli_flag1(opts, "modality", "rna"),
      assignment.thr = .ic2_cli_as_num(opts$assignment.thr, 50),
      cov.thr = .ic2_cli_as_num(opts$cov.thr, 10),
      return.stats = FALSE,
      verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
    )
  }, "i.score")
}

.ic2_cli_cmd_find_markers <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", transform = "auto", fold.change = "2",
    pval.test = "t.test", positive = TRUE, out = "markers.tsv", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 findMarkers [options]\n\n",
      "Options: -soda PATH -out FILE [markers.tsv]\n",
      "         -modality NAME -assay.name NAME -transform MODE\n",
      "         -fold.change N -pval.test t.test|wilcox.test -positive BOOL -verbose BOOL\n"
    ))
  }
  obj <- .ic2_cli_open(opts)
  markers <- findMarkers(
    obj,
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    assay.name = .ic2_cli_flag1(opts, "assay.name", NULL),
    transform = .ic2_cli_flag1(opts, "transform", "auto"),
    fold.change = .ic2_cli_as_num(opts$fold.change, 2),
    pval.test = .ic2_cli_flag1(opts, "pval.test", "t.test"),
    positive = .ic2_cli_as_logical(opts$positive, TRUE),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  out <- .ic2_cli_flag1(opts, "out", "markers.tsv")
  utils::write.table(markers, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[iCellR2] findMarkers OK\n  wrote: ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(NULL)
}

.ic2_cli_cmd_top_markers <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    markers = "markers.tsv", topde = "10", min.base.mean = "0.2",
    filt.ambig = TRUE, out = "top_markers.txt"
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 top.markers [options]\n\n",
      "Options: -markers FILE -topde N -min.base.mean N -filt.ambig BOOL -out FILE\n"
    ))
  }
  mf <- .ic2_cli_flag1(opts, "markers", "markers.tsv")
  if (!file.exists(mf)) stop("markers file not found: ", mf, call. = FALSE)
  markers <- utils::read.delim(mf, stringsAsFactors = FALSE, check.names = FALSE)
  top <- top.markers(
    markers,
    topde = .ic2_cli_as_int(opts$topde, 10L),
    min.base.mean = .ic2_cli_as_num(opts$min.base.mean, 0.2),
    filt.ambig = .ic2_cli_as_logical(opts$filt.ambig, TRUE)
  )
  out <- .ic2_cli_flag1(opts, "out", "top_markers.txt")
  ## top.markers may return character genes or a table
  if (is.character(top)) {
    writeLines(top, out)
  } else if (is.data.frame(top)) {
    utils::write.table(top, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    writeLines(as.character(unlist(top)), out)
  }
  cat("[iCellR2] top.markers OK\n  wrote: ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(NULL)
}

.ic2_cli_cmd_run_diff_exp <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    compare = "clusters", modality = "rna", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 run.diff.exp [options]\n\n",
      "Options: -soda PATH -compare clusters|conditions\n",
      "         -group.1 ID -group.2 ID  (required)\n",
      "         -out FILE.tsv -modality NAME -verbose BOOL\n"
    ))
  }
  g1 <- .ic2_cli_flag1(opts, "group.1", .ic2_cli_flag1(opts, "group_1", NULL))
  g2 <- .ic2_cli_flag1(opts, "group.2", .ic2_cli_flag1(opts, "group_2", NULL))
  if (is.null(g1) || is.null(g2)) {
    stop("provide -group.1 and -group.2", call. = FALSE)
  }
  obj <- .ic2_cli_open(opts)
  de <- run.diff.exp(
    obj,
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    compare = .ic2_cli_flag1(opts, "compare", "clusters"),
    group.1 = g1,
    group.2 = g2,
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  out <- .ic2_cli_flag1(opts, "out", paste0("de_", g1, "_vs_", g2, ".tsv"))
  if (is.data.frame(de)) {
    utils::write.table(de, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
  } else if (is.list(de) && !is.null(de$results) && is.data.frame(de$results)) {
    utils::write.table(de$results, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    ## object may store DE in misc / SODA; still try to coerce
    tab <- tryCatch(as.data.frame(de), error = function(e) NULL)
    if (!is.null(tab)) {
      utils::write.table(tab, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  cat("[iCellR2] run.diff.exp OK\n  wrote: ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(NULL)
}

.ic2_cli_cmd_volcano <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    plot.type = "volcano", interactive = FALSE,
    width = "7", height = "6", dpi = "150", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 volcano.ma.plot [options]\n\n",
      "Options: -de FILE.tsv | -soda PATH (uses latest DE in SODA)\n",
      "         -plot.type volcano|ma -out FILE.png -interactive BOOL\n"
    ))
  }
  de_path <- .ic2_cli_flag1(opts, "de", .ic2_cli_flag1(opts, "input", NULL))
  if (!is.null(de_path) && file.exists(de_path)) {
    de <- utils::read.delim(de_path, check.names = FALSE)
  } else {
    obj <- .ic2_cli_open(opts)
    de <- obj
  }
  p <- volcano.ma.plot(
    de,
    plot.type = .ic2_cli_flag1(opts, "plot.type", "volcano"),
    interactive = .ic2_cli_as_logical(opts$interactive, FALSE),
    out.name = .ic2_cli_flag1(opts, "out.name", NULL)
  )
  if (!isTRUE(.ic2_cli_as_logical(opts$interactive, FALSE))) {
    .ic2_cli_save_gg(
      p, .ic2_cli_flag1(opts, "out", "volcano.png"),
      width = .ic2_cli_as_num(opts$width, 7),
      height = .ic2_cli_as_num(opts$height, 6),
      dpi = .ic2_cli_as_int(opts$dpi, 150L)
    )
  }
  cat("[iCellR2] volcano.ma.plot OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_clust_avg <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(modality = "rna", verbose = TRUE))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 clust.avg.exp [options]\n\n",
      "Options: -soda PATH -conds.to.avg ID|NULL -modality NAME -verbose BOOL\n"
    ))
  }
  obj <- .ic2_cli_open(opts)
  conds <- .ic2_cli_flag1(opts, "conds.to.avg", NULL)
  if (!is.null(conds) && identical(toupper(as.character(conds)), "NULL")) conds <- NULL
  clust.avg.exp(
    obj,
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    conds.to.avg = conds,
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  cat("[iCellR2] clust.avg.exp OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_clust_stats_plot <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    plot.type = "box.mito", interactive = FALSE,
    width = "8", height = "5", dpi = "300", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 clust.stats.plot [options]\n\n",
      "Options: -soda PATH -plot.type box.mito|box.gene|box.umi|...\n",
      "         -out FILE.png -interactive BOOL\n"
    ))
  }
  obj <- .ic2_cli_open(opts)
  p <- clust.stats.plot(
    obj,
    plot.type = .ic2_cli_flag1(opts, "plot.type", "box.mito"),
    interactive = .ic2_cli_as_logical(opts$interactive, FALSE),
    out.name = .ic2_cli_flag1(opts, "out.name", NULL),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  if (!isTRUE(.ic2_cli_as_logical(opts$interactive, FALSE))) {
    .ic2_cli_save_gg(
      p, .ic2_cli_flag1(opts, "out", NULL),
      width = .ic2_cli_as_num(opts$width, 8),
      height = .ic2_cli_as_num(opts$height, 5),
      dpi = .ic2_cli_as_int(opts$dpi, 300L)
    )
  }
  cat("[iCellR2] clust.stats.plot OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_clust_cond_info <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    plot.type = "pie", my.out.put = "plot",
    normalize.ncell = TRUE, normalize.by = "percentage",
    interactive = FALSE, width = "10", height = "8", dpi = "300", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 clust.cond.info [options]\n\n",
      "Options: -soda PATH -plot.type pie|bar|table\n",
      "         -my.out.put plot|data -normalize.ncell BOOL\n",
      "         -normalize.by percentage|counts -out FILE.png\n"
    ))
  }
  obj <- .ic2_cli_open(opts)
  res <- clust.cond.info(
    obj,
    plot.type = .ic2_cli_flag1(opts, "plot.type", "pie"),
    my.out.put = .ic2_cli_flag1(opts, "my.out.put", "plot"),
    normalize.ncell = .ic2_cli_as_logical(opts[["normalize.ncell"]], TRUE),
    normalize.by = .ic2_cli_flag1(opts, "normalize.by", "percentage")
  )
  out_mode <- .ic2_cli_flag1(opts, "my.out.put", "plot")
  if (identical(out_mode, "plot") && !is.null(res)) {
    .ic2_cli_save_gg(
      res, .ic2_cli_flag1(opts, "out", NULL),
      width = .ic2_cli_as_num(opts$width, 10),
      height = .ic2_cli_as_num(opts$height, 8),
      dpi = .ic2_cli_as_int(opts$dpi, 300L)
    )
  }
  cat("[iCellR2] clust.cond.info OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_cluster_plot <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", plot.type = "umap", col.by = "clusters",
    dims = "1:2", interactive = FALSE, width = "8", height = "8", dpi = "300",
    verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 cluster.plot [options]\n\n",
      "Options: -soda PATH\n",
      "         -plot.type umap|pca|tsne|knetl|cgumap|cgpca|ptumap|\n",
      "                   aligned_umap|aligned_tsne|aligned_knetl|...\n",
      "         -col.by clusters|Pseudotime|branch_id|...\n",
      "         -out FILE.png -width N -height N -dpi N\n",
      "         -interactive BOOL -out.name PREFIX (html) -verbose BOOL\n"
    ))
  }
  obj <- .ic2_cli_open(opts)
  p <- cluster.plot(
    obj,
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    plot.type = .ic2_cli_flag1(opts, "plot.type", "umap"),
    dims = .ic2_cli_as_dims(opts$dims, 1:2),
    col.by = .ic2_cli_flag1(opts, "col.by", "clusters"),
    interactive = .ic2_cli_as_logical(opts$interactive, FALSE),
    out.name = .ic2_cli_flag1(opts, "out.name", NULL),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  out <- .ic2_cli_flag1(opts, "out", NULL)
  if (!isTRUE(.ic2_cli_as_logical(opts$interactive, FALSE))) {
    .ic2_cli_save_gg(
      p, out,
      width = .ic2_cli_as_num(opts$width, 8),
      height = .ic2_cli_as_num(opts$height, 8),
      dpi = .ic2_cli_as_int(opts$dpi, 300L)
    )
  }
  cat("[iCellR2] cluster.plot OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_pseudotime_arrows <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", plot.type = "umap", col.by = "clusters",
    style = "stream", width = "8", height = "6", dpi = "300", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 pseudotime.arrows [options]\n\n",
      "Options: -soda PATH -plot.type umap|pca|cgumap -col.by COL\n",
      "         -style stream|arrows|both -out FILE.png -width N -height N -dpi N\n"
    ))
  }
  obj <- .ic2_cli_open(opts)
  p <- pseudotime.arrows(
    obj,
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    plot.type = .ic2_cli_flag1(opts, "plot.type", "umap"),
    col.by = .ic2_cli_flag1(opts, "col.by", "clusters"),
    style = .ic2_cli_flag1(opts, "style", "stream")
  )
  ## may return a list for style=both
  if (is.list(p) && !inherits(p, "ggplot") && !inherits(p, "gg")) {
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      p <- do.call(gridExtra::arrangeGrob, c(p, list(ncol = length(p))))
    } else {
      p <- p[[1L]]
    }
  }
  .ic2_cli_save_gg(
    p, .ic2_cli_flag1(opts, "out", NULL),
    width = .ic2_cli_as_num(opts$width, 8),
    height = .ic2_cli_as_num(opts$height, 6),
    dpi = .ic2_cli_as_int(opts$dpi, 300L)
  )
  cat("[iCellR2] pseudotime.arrows OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_heatmap <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", group.by = "cells", interactive = FALSE,
    width = "8", height = "10", dpi = "300", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 heatmap.gg.plot [options]\n\n",
      "Options: -soda PATH -gene A,B,C | -gene.file FILE\n",
      "         -group.by cells|clusters|conditions -out FILE.png\n",
      "         -width N -height N -dpi N -interactive BOOL -out.name PREFIX\n"
    ))
  }
  genes <- .ic2_cli_genes(opts)
  if (!length(genes)) stop("provide -gene A,B,C or -gene.file FILE", call. = FALSE)
  obj <- .ic2_cli_open(opts)
  p <- heatmap.gg.plot(
    obj,
    gene = genes,
    group.by = .ic2_cli_flag1(opts, "group.by", "cells"),
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    interactive = .ic2_cli_as_logical(opts$interactive, FALSE),
    out.name = .ic2_cli_flag1(opts, "out.name", NULL),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  if (!isTRUE(.ic2_cli_as_logical(opts$interactive, FALSE))) {
    .ic2_cli_save_gg(
      p, .ic2_cli_flag1(opts, "out", NULL),
      width = .ic2_cli_as_num(opts$width, 8),
      height = .ic2_cli_as_num(opts$height, 10),
      dpi = .ic2_cli_as_int(opts$dpi, 300L)
    )
  }
  cat("[iCellR2] heatmap.gg.plot OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_bubble <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", group.by = "clusters", interactive = FALSE,
    width = "8", height = "10", dpi = "300", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 bubble.gg.plot [options]\n\n",
      "Options: -soda PATH -gene A,B,C | -gene.file FILE\n",
      "         -group.by clusters|conditions -out FILE.png -interactive BOOL\n"
    ))
  }
  genes <- .ic2_cli_genes(opts)
  if (!length(genes)) stop("provide -gene A,B,C or -gene.file FILE", call. = FALSE)
  obj <- .ic2_cli_open(opts)
  p <- bubble.gg.plot(
    obj,
    gene = genes,
    group.by = .ic2_cli_flag1(opts, "group.by", "clusters"),
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    interactive = .ic2_cli_as_logical(opts$interactive, FALSE),
    out.name = .ic2_cli_flag1(opts, "out.name", NULL),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  if (!isTRUE(.ic2_cli_as_logical(opts$interactive, FALSE))) {
    .ic2_cli_save_gg(
      p, .ic2_cli_flag1(opts, "out", NULL),
      width = .ic2_cli_as_num(opts$width, 8),
      height = .ic2_cli_as_num(opts$height, 10),
      dpi = .ic2_cli_as_int(opts$dpi, 300L)
    )
  }
  cat("[iCellR2] bubble.gg.plot OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_gene_plot <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(
    modality = "rna", plot.type = "scatter", reduction = "umap",
    interactive = FALSE, width = "12", height = "10", dpi = "300", verbose = TRUE
  ))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 gene.plot [options]\n\n",
      "Options: -soda PATH -gene A,B,C | -gene.file FILE\n",
      "         -plot.type scatter|violin\n",
      "         -reduction umap|pca|knetl|aligned_umap|aligned_knetl|cgumap|...\n",
      "         -out FILE.png -interactive BOOL -out.name PREFIX\n"
    ))
  }
  genes <- .ic2_cli_genes(opts)
  if (!length(genes)) stop("provide -gene A,B,C or -gene.file FILE", call. = FALSE)
  obj <- .ic2_cli_open(opts)
  p <- gene.plot(
    obj,
    gene = genes,
    modality = .ic2_cli_flag1(opts, "modality", "rna"),
    plot.type = .ic2_cli_flag1(opts, "plot.type", "scatter"),
    reduction = .ic2_cli_flag1(opts, "reduction", "umap"),
    interactive = .ic2_cli_as_logical(opts$interactive, FALSE),
    out.name = .ic2_cli_flag1(opts, "out.name", NULL),
    verbose = .ic2_cli_as_logical(opts$verbose, TRUE)
  )
  if (!isTRUE(.ic2_cli_as_logical(opts$interactive, FALSE))) {
    .ic2_cli_save_gg(
      p, .ic2_cli_flag1(opts, "out", NULL),
      width = .ic2_cli_as_num(opts$width, 12),
      height = .ic2_cli_as_num(opts$height, 10),
      dpi = .ic2_cli_as_int(opts$dpi, 300L)
    )
  }
  cat("[iCellR2] gene.plot OK\n")
  invisible(NULL)
}

## ---- SODA utils ------------------------------------------------------------

.ic2_cli_cmd_open <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(verbose = TRUE))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help("Usage:   iCellR2 open -soda PATH [-verbose BOOL]\n"))
  }
  path <- .ic2_cli_resolve_soda(opts)
  open.soda(path, verbose = .ic2_cli_as_logical(opts$verbose, TRUE))
  cat("[iCellR2] open OK  |  ", basename(normalizePath(path)), "\n", sep = "")
  invisible(NULL)
}

.ic2_cli_cmd_validate <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(verbose = TRUE))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help("Usage:   iCellR2 validate -soda PATH [-verbose BOOL]\n"))
  }
  validate.soda(.ic2_cli_resolve_soda(opts), verbose = .ic2_cli_as_logical(opts$verbose, TRUE))
  cat("[iCellR2] validate OK\n")
  invisible(NULL)
}

.ic2_cli_cmd_refresh <- function(args) {
  opts <- .ic2_cli_parse(args)
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help("Usage:   iCellR2 refresh -soda PATH\n"))
  }
  path <- .ic2_cli_resolve_soda(opts)
  refresh.soda(path)
  cat("[iCellR2] refresh OK\n  index: ", file.path(normalizePath(path), "index.html"), "\n", sep = "")
  invisible(NULL)
}

.ic2_cli_cmd_install <- function(args) {
  opts <- .ic2_cli_parse(args, defaults = list(dir = path.expand("~/bin")))
  if (isTRUE(opts[["__help__"]])) {
    return(.ic2_cli_help(
      "Usage:   iCellR2 install-cli [-dir PATH]\n\n",
      "Symlinks exec/iCellR2 into ~/bin (add that dir to PATH).\n"
    ))
  }
  src <- system.file("exec", "iCellR2", package = "iCellR2")
  if (!nzchar(src) || !file.exists(src)) {
    cand <- Sys.getenv("ICELLR2_CLI_HOME", unset = "")
    if (nzchar(cand) && file.exists(file.path(cand, "iCellR2"))) {
      src <- file.path(cand, "iCellR2")
    } else {
      stop("Could not find exec/iCellR2. Reinstall iCellR2 with the exec/ folder.", call. = FALSE)
    }
  }
  dest_dir <- path.expand(.ic2_cli_flag1(opts, "dir", "~/bin"))
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(dest_dir, "iCellR2")
  if (file.exists(dest)) unlink(dest)
  file.symlink(normalizePath(src, winslash = "/", mustWork = TRUE), dest)
  Sys.chmod(src, mode = "0755")
  cat("[iCellR2] install-cli OK\n  linked: ", dest, " -> ", normalizePath(src), "\n", sep = "")
  cat("  ensure PATH contains: ", dest_dir, "\n", sep = "")
  invisible(NULL)
}

#' Install the iCellR2 bash launcher onto PATH
#'
#' Symlinks \code{exec/iCellR2} into \code{~/bin} (or \code{dir}).
#'
#' @param dir Install directory. Default \code{~/bin}.
#' @return Path to the symlink (invisible).
#' @export
install.cli <- function(dir = path.expand("~/bin")) {
  old <- getOption("iCellR2.cli.noquit")
  options(iCellR2.cli.noquit = TRUE)
  on.exit(options(iCellR2.cli.noquit = old), add = TRUE)
  .ic2_cli_cmd_install(c("-dir", dir))
  invisible(file.path(path.expand(dir), "iCellR2"))
}
