#' Write interactive gene-model HTML (mean-dispersion explorer)
#'
#' Self-contained dashboard matching
#' \code{demo_iCellR2/plot_design/gene_model_mean.html}.
#'
#' @keywords internal
.ic2_gene_model_write_html <- function(
    data,
    file,
    gate.mode = "classic",
    dispersion.limit = 1.5,
    mean.limit = 0.4,
    base.mean.rank = 1000,
    cell.count.limit = NA_real_,
    gene.num.max = 2000L,
    no.mito.model = TRUE,
    mark.mito = FALSE,
    selected = character(0),
    n_all = NULL
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("Package 'jsonlite' is required for gene-model HTML; skipping.")
    return(FALSE)
  }
  file <- path.expand(as.character(file)[1L])
  parent <- dirname(file)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }

  genes <- as.character(data$genes)
  me <- as.numeric(data$meanExp)
  nc <- as.numeric(data$numberOfCells)
  d <- as.numeric(data$SDs)
  me[!is.finite(me)] <- 0
  nc[!is.finite(nc)] <- 0
  d[!is.finite(d)] <- 0
  mito <- grepl("^mt[-.]", genes, ignore.case = TRUE)
  sel_set <- as.character(selected)

  ## Compact gene records for the browser
  gene_df <- data.frame(
    g = genes,
    nc = nc,
    lnc = log2(nc + 1),
    me = me,
    lm = log2(me + 1),
    d = d,
    sel = genes %in% sel_set,
    mito = mito,
    stringsAsFactors = FALSE
  )

  payload <- list(
    n_all = if (is.null(n_all)) length(genes) else as.integer(n_all),
    n_plot = length(genes),
    n_selected = length(sel_set),
    thresholds = list(
      gate_mode = gate.mode,
      dispersion_limit = as.numeric(dispersion.limit),
      mean_limit = as.numeric(mean.limit),
      base_mean_rank = as.numeric(base.mean.rank),
      cell_count_limit = as.numeric(cell.count.limit),
      top_rank_line = log2(pmax(as.numeric(cell.count.limit), 1)),
      gene_num_max = as.integer(gene.num.max),
      no_mito = isTRUE(no.mito.model),
      mark_mito = isTRUE(mark.mito)
    ),
    genes = gene_df
  )

  json <- jsonlite::toJSON(
    payload,
    dataframe = "rows",
    auto_unbox = TRUE,
    digits = 6,
    null = "null"
  )

  tmpl_path <- system.file("templates", "gene_model_dashboard.html", package = "iCellR2")
  if (!nzchar(tmpl_path) || !file.exists(tmpl_path)) {
    ## Dev fallbacks (devtools::load_all / running from package or demo dir)
    roots <- unique(c(
      getwd(),
      file.path(getwd(), "..", "iCellR2"),
      file.path(getwd(), "iCellR2")
    ))
    for (root in roots) {
      p <- normalizePath(
        file.path(root, "inst", "templates", "gene_model_dashboard.html"),
        mustWork = FALSE
      )
      if (file.exists(p)) {
        tmpl_path <- p
        break
      }
    }
  }
  if (!nzchar(tmpl_path) || !file.exists(tmpl_path)) {
    warning("gene_model_dashboard.html template not found; skipping HTML.")
    return(FALSE)
  }

  html <- paste(readLines(tmpl_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (!grepl("__IC2_GENE_MODEL_JSON__", html, fixed = TRUE)) {
    warning("gene_model_dashboard.html missing JSON placeholder; skipping HTML.")
    return(FALSE)
  }
  html <- sub("__IC2_GENE_MODEL_JSON__", json, html, fixed = TRUE)
  ## writeLines can mangle UTF-8 on some locales; writeBin is safer
  con <- file(file, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(paste0(html, "\n"))), con)
  TRUE
}
