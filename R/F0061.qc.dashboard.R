#' Interactive QC filter dashboard
#'
#' Writes an interactive HTML QC explorer (Plotly + embedded JS) so viewers
#' can try \code{min.genes} / \code{min.umis} / \code{max.mito} cuts before
#' calling \code{\link{cell.filter}}. Optionally also writes static
#' publication exports.
#'
#' @param x An \code{iCellR2} object (run \code{\link{qc.stats}} first).
#' @param out.name Output path prefix (no extension), or a directory + stem.
#'   Default \code{"qc_dashboard"}.
#' @param format Character vector of formats to write. Default \code{"html"}.
#'   Allowed: \code{"html"}, \code{"pdf"}, \code{"png"}, \code{"svg"}.
#'   HTML is always recommended; PDF/PNG/SVG are static multi-panel snapshots.
#' @param min.genes,min.umis,max.mito Default filter cutoffs shown in the
#'   dashboard (and drawn on static exports). Defaults match common
#'   \code{\link{cell.filter}} usage.
#' @param group.by Metadata column for condition panels / keep table.
#'   Default \code{"condition_id"}.
#' @param max.cells Max cells embedded in HTML / drawn on static plots.
#'   Default \code{500000}.
#' @param width,height Device size in inches for PDF/PNG/SVG.
#' @param dpi DPI for PNG (ignored for PDF/SVG vector). Default \code{300}.
#' @param seed RNG seed when downsampling. Default \code{1}.
#' @param verbose Progress reporter style.
#'
#' @return Invisibly, a named character vector of files written.
#'
#' @examples
#' \dontrun{
#' qc.dashboard(obj, out.name = "PBMC_qc")                 # HTML only
#' qc.dashboard(obj, out.name = "PBMC_qc", format = c("html", "pdf", "png"))
#' }
#'
#' @export
qc.dashboard <- function(
    x,
    out.name = "qc_dashboard",
    format = "html",
    min.genes = 200,
    min.umis = 500,
    max.mito = 0.1,
    group.by = "condition_id",
    max.cells = 500000L,
    width = 12,
    height = 9,
    dpi = 300,
    seed = 1L,
    verbose = TRUE
) {
  .ic2_require_object(x)
  rep <- .ic2_new_reporter("qc.dashboard", verbose, soda_root = x@storage$root_dir)

  format <- unique(tolower(as.character(format)))
  ok <- c("html", "pdf", "png", "svg")
  bad <- setdiff(format, ok)
  if (length(bad)) {
    stop("Unknown format(s): ", paste(bad, collapse = ", "),
         ". Use: ", paste(ok, collapse = ", "))
  }
  if (!length(format)) format <- "html"

  meta <- x@cell_data
  if (!is.data.frame(meta) || nrow(meta) < 1L) {
    stop("'x@cell_data' is empty; run qc.stats() first.")
  }
  need <- c("nGenes", "UMIs", "mito.percent")
  miss <- setdiff(need, colnames(meta))
  if (length(miss)) {
    stop("Missing QC columns in cell_data: ", paste(miss, collapse = ", "),
         ". Run qc.stats().")
  }

  if (!group.by %in% colnames(meta)) {
    meta[[group.by]] <- "all"
  }
  # One sample / missing labels: treat blanks and NAs as a single group.
  cond <- as.character(meta[[group.by]])
  cond[is.na(cond) | !nzchar(cond)] <- "all"

  df <- data.frame(
    cell_id = as.character(rownames(meta)),
    condition_id = cond,
    nGenes = as.numeric(meta$nGenes),
    UMIs = as.numeric(meta$UMIs),
    mito.percent = as.numeric(meta$mito.percent),
    stringsAsFactors = FALSE
  )
  if ("cluster_id" %in% colnames(meta)) {
    df$cluster_id <- meta$cluster_id
  }

  n_all <- nrow(df)
  max.cells <- as.integer(max.cells)
  if (is.finite(max.cells) && n_all > max.cells) {
    set.seed(as.integer(seed))
    df <- df[sample.int(n_all, max.cells), , drop = FALSE]
    .ic2_step(rep, paste0(
      "Downsampled ", .ic2_fmt_num(n_all), "  ->  ", .ic2_fmt_num(nrow(df)), " cells"
    ))
  }

  filters <- list(
    min_genes = as.numeric(min.genes),
    min_umis = as.numeric(min.umis),
    max_mito = as.numeric(max.mito)
  )

  out_stem <- .ic2_qc_dash_stem(out.name)
  written <- character()

  if ("html" %in% format) {
    .ic2_step(rep, "Writing interactive HTML")
    fhtml <- paste0(out_stem, ".html")
    .ic2_qc_dash_write_html(df, filters, fhtml)
    written <- c(written, html = fhtml)
  }

  static_formats <- intersect(format, c("pdf", "png", "svg"))
  if (length(static_formats)) {
    .ic2_step(rep, paste0("Writing static export(s): ", paste(static_formats, collapse = ", ")))
    for (fmt in static_formats) {
      fpath <- paste0(out_stem, ".", fmt)
      .ic2_qc_dash_write_static(
        df, filters,
        file = fpath,
        format = fmt,
        width = width,
        height = height,
        dpi = dpi
      )
      written <- c(written, stats::setNames(fpath, fmt))
    }
  }

  .ic2_finish(rep, extra_lines = paste0(
    "Wrote: ", paste(basename(written), collapse = ", ")
  ))
  invisible(written)
}

#' @keywords internal
.ic2_qc_dash_stem <- function(out.name) {
  out.name <- path.expand(as.character(out.name)[1L])
  if (dir.exists(out.name) || grepl("[/\\\\]$", out.name)) {
    dir.create(out.name, recursive = TRUE, showWarnings = FALSE)
    return(file.path(sub("[/\\\\]$", "", out.name), "qc_dashboard"))
  }
  parent <- dirname(out.name)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  # strip extension if user passed one
  sub("\\.(html|pdf|png|svg)$", "", out.name, ignore.case = TRUE)
}

#' @keywords internal
.ic2_qc_dash_pass <- function(df, filters) {
  df$nGenes >= filters$min_genes &
    df$UMIs >= filters$min_umis &
    df$mito.percent <= filters$max_mito
}

#' @keywords internal
.ic2_qc_dash_write_static <- function(df, filters, file, format, width, height, dpi) {
  # Publication QC figure (pre-filter): mito-coloured scatter + 3 sample boxplots.
  # Filter cutoffs are NOT drawn  -  those belong in the interactive HTML explorer.
  cond_levels <- .ic2_qc_dash_cond_levels(df$condition_id)
  df$condition_id <- factor(as.character(df$condition_id), levels = cond_levels)
  n_cond <- length(cond_levels)

  ink <- "#1F1F1F"
  mute <- "#6E6E6E"
  theme_qc <- .ic2_qc_dash_theme(ink, mute, n_cond = n_cond)

  lim_g <- range(df$nGenes, na.rm = TRUE)
  lim_u <- range(df$UMIs, na.rm = TRUE)
  pad <- function(r, p = 0.04) {
    d <- diff(r)
    if (!is.finite(d) || d <= 0) d <- 1
    c(r[1] - p * d, r[2] + p * d)
  }
  lim_g <- pad(lim_g)
  lim_u <- pad(lim_u)

  p_scatter <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$nGenes, .data$UMIs, colour = .data$mito.percent)
  ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.75, stroke = 0) +
    ggplot2::scale_colour_gradientn(
      colours = c("#E8E8E8", "#F4A582", "#B2182B"),
      name = "mito%"
    ) +
    ggplot2::coord_cartesian(xlim = lim_g, ylim = lim_u, expand = FALSE) +
    ggplot2::labs(
      title = "UMI  x  genes",
      subtitle = "Coloured by mitochondrial fraction",
      x = "nGenes", y = "UMIs"
    ) +
    theme_qc +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
    )

  box_n <- .ic2_qc_dash_box(
    df, yvar = "nGenes", ylab = "nGenes",
    title = "Genes / cell",
    cut = NULL, theme_qc = theme_qc, ink = ink, n_cond = n_cond
  )
  box_u <- .ic2_qc_dash_box(
    df, yvar = "UMIs", ylab = "UMIs",
    title = "UMIs / cell",
    cut = NULL, theme_qc = theme_qc, ink = ink, n_cond = n_cond
  )
  box_m <- .ic2_qc_dash_box(
    df, yvar = "mito.percent", ylab = "mito.percent",
    title = "Mito % / cell",
    cut = NULL, theme_qc = theme_qc, ink = ink, n_cond = n_cond
  )

  header <- sprintf(
    "QC metrics    |    %s cells    |    %s sample%s",
    format(nrow(df), big.mark = ","),
    n_cond,
    if (n_cond == 1L) "" else "s"
  )

  width <- max(as.numeric(width), 11 + 0.25 * max(0, n_cond - 4))
  height <- max(as.numeric(height), 8.5)

  if (requireNamespace("ggpubr", quietly = TRUE)) {
    # 2 x 2: scatter spans conceptually top-left; three boxes fill rest  - 
    # layout: scatter | genes ; UMIs | mito
    grob <- ggpubr::ggarrange(
      p_scatter, box_n,
      box_u, box_m,
      ncol = 2, nrow = 2,
      align = "hv"
    )
    grob <- ggpubr::annotate_figure(
      grob,
      top = ggpubr::text_grob(header, face = "bold", size = 11, color = ink, hjust = 0, x = 0.01)
    )
  } else {
    main <- grid::textGrob(
      header,
      gp = grid::gpar(fontsize = 11, col = ink, fontface = "bold"),
      x = 0.01, hjust = 0
    )
    grob <- gridExtra::arrangeGrob(
      p_scatter, box_n, box_u, box_m,
      ncol = 2, top = main
    )
  }

  format <- tolower(format)
  if (identical(format, "pdf")) {
    ggplot2::ggsave(file, grob, width = width, height = height, device = grDevices::pdf)
  } else if (identical(format, "svg")) {
    ggplot2::ggsave(file, grob, width = width, height = height, device = grDevices::svg)
  } else if (identical(format, "png")) {
    ggplot2::ggsave(file, grob, width = width, height = height, dpi = dpi, device = "png")
  } else {
    stop("Unsupported static format: ", format)
  }
  invisible(file)
}

#' Shared modern open-axis theme; tweaks x labels for many conditions.
#' @keywords internal
.ic2_qc_dash_theme <- function(ink = "#1F1F1F", mute = "#6E6E6E", n_cond = 4L) {
  x_angle <- if (n_cond > 8L) 60 else if (n_cond > 5L) 35 else 0
  x_hjust <- if (x_angle > 0) 1 else 0.5
  ggplot2::theme_classic(base_size = 11, base_family = "") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(
        face = "bold", colour = ink, size = 11.5, hjust = 0,
        margin = ggplot2::margin(b = 2)
      ),
      plot.subtitle = ggplot2::element_text(
        colour = mute, size = 9, hjust = 0,
        margin = ggplot2::margin(b = 6)
      ),
      axis.title = ggplot2::element_text(colour = ink, size = 10),
      axis.text = ggplot2::element_text(colour = mute, size = 8.5),
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = x_hjust, vjust = 1),
      axis.line = ggplot2::element_line(colour = ink, linewidth = 0.45, lineend = "square"),
      axis.ticks = ggplot2::element_line(colour = ink, linewidth = 0.35),
      axis.ticks.length = grid::unit(3, "pt"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(colour = mute, size = 8.5),
      legend.title = ggplot2::element_text(colour = ink, size = 9),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(6, 10, 6, 10)
    )
}

#' Order conditions: alphabetical, but keep numeric-looking ids numeric-sorted.
#' @keywords internal
.ic2_qc_dash_cond_levels <- function(x) {
  u <- unique(as.character(x))
  u <- u[!is.na(u) & nzchar(u)]
  if (!length(u)) return("all")
  if (all(grepl("^[-]?[0-9]+$", u))) {
    return(as.character(sort(as.integer(u))))
  }
  sort(u)
}

#' Soft sample colours (cycles for any n).
#' @keywords internal
.ic2_qc_dash_sample_cols <- function(n) {
  base <- c(
    "#4C78A8", "#F58518", "#54A24B", "#E45756",
    "#72B7B2", "#B279A2", "#FF9DA6", "#9D755D",
    "#BAB0AC", "#1F7A6C", "#6B5B95", "#D4A373"
  )
  if (n <= length(base)) base[seq_len(n)] else rep(base, length.out = n)
}

#' Violin + box + jitter by condition; widths adapt to n samples.
#' @param cut Optional horizontal reference line (HTML filter preview).
#' @keywords internal
.ic2_qc_dash_box <- function(df, yvar, ylab, title, cut = NULL, theme_qc, ink, n_cond = NULL) {
  if (is.null(n_cond)) n_cond <- nlevels(df$condition_id)
  n_cond <- max(1L, as.integer(n_cond))
  box_w <- max(0.18, min(0.42, 1.15 / n_cond))
  jit_w <- max(0.03, min(0.14, 0.45 / n_cond))
  pt_size <- if (n_cond > 8L) 0.35 else if (n_cond > 4L) 0.45 else 0.55
  cols <- .ic2_qc_dash_sample_cols(n_cond)
  names(cols) <- levels(df$condition_id)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$condition_id, y = .data[[yvar]], fill = .data$condition_id)
  ) +
    ggplot2::geom_violin(
      trim = TRUE,
      colour = ink,
      linewidth = 0.3,
      alpha = 0.35,
      scale = "width",
      width = 0.9
    ) +
    ggplot2::geom_boxplot(
      width = box_w,
      fill = "white",
      colour = ink,
      linewidth = 0.35,
      outlier.shape = NA,
      alpha = 0.95
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(colour = .data$condition_id),
      height = 0,
      width = jit_w,
      size = pt_size,
      alpha = 0.45,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      fun = mean, geom = "point", size = 1.8, colour = ink, shape = 18
    ) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = if (n_cond == 1L) "Single sample" else "By sample",
      x = NULL,
      y = ylab
    ) +
    theme_qc +
    ggplot2::theme(legend.position = "none")

  if (!is.null(cut) && is.finite(cut)) {
    p <- p + ggplot2::geom_hline(
      yintercept = cut, colour = "#C43C3C", linetype = 2, linewidth = 0.55
    )
  }
  p
}

#' @keywords internal
.ic2_qc_dash_write_html <- function(df, filters, file) {
  # ranges for sliders
  rng <- list(
    min_genes = list(
      min = 0,
      max = max(500, as.numeric(stats::quantile(df$nGenes, 0.995, names = FALSE))),
      step = 10
    ),
    min_umis = list(
      min = 0,
      max = max(1000, as.numeric(stats::quantile(df$UMIs, 0.995, names = FALSE))),
      step = 50
    ),
    max_mito = list(min = 0, max = 0.25, step = 0.005)
  )

  payload <- list(
    qc = df[, c("cell_id", "condition_id", "nGenes", "UMIs", "mito.percent")],
    condition_order = as.list(.ic2_qc_dash_cond_levels(df$condition_id)),
    qc_filters = list(
      defaults = list(
        min_genes = filters$min_genes,
        min_umis = filters$min_umis,
        max_mito = filters$max_mito
      ),
      ranges = rng,
      note = "Interactive QC filter explorer"
    )
  )

  json <- jsonlite::toJSON(payload, dataframe = "rows", auto_unbox = TRUE, digits = 6, null = "null")
  html <- .ic2_qc_dash_html_template(json)
  # self-contained single file
  writeLines(html, file, useBytes = TRUE)
  invisible(file)
}

#' @keywords internal
.ic2_qc_dash_html_template <- function(json_payload) {
  paste0(
'<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>QC filter explorer</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root, [data-theme="light"]{
  --ink:#14242a; --mute:#5a6b72; --line:#d7e0e3; --paper:#eef3f4; --card:#ffffff;
  --accent:#0f766e; --pass:#0f766e; --fail:#c5cdd1; --cut:#be123c;
  --chip:#ccfbf1; --shadow:0 10px 30px rgba(15,55,60,.08);
}
[data-theme="dark"]{
  --ink:#e8f1f2; --mute:#9bb0b6; --line:#2a3f44; --paper:#0b1619; --card:#122226;
  --accent:#2dd4bf; --pass:#2dd4bf; --fail:#3a4a4e; --cut:#fb7185;
  --chip:#134e4a; --shadow:0 10px 30px rgba(0,0,0,.35);
}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);
  font-family:"Avenir Next","Segoe UI","Helvetica Neue",Arial,sans-serif}
.wrap{max-width:1280px;margin:0 auto;padding:22px 18px 56px}
.top{display:flex;align-items:flex-start;justify-content:space-between;gap:12px;flex-wrap:wrap;margin-bottom:10px}
h1{margin:0;font-size:26px;letter-spacing:-.03em}
.lede{color:var(--mute);font-size:14px;margin:6px 0 0;max-width:72ch;line-height:1.5}
.btn{appearance:none;border:1px solid var(--line);background:var(--card);color:var(--ink);
  border-radius:10px;padding:8px 12px;font:inherit;font-size:12.5px;cursor:pointer;font-weight:600}
.btn:hover{border-color:var(--accent);color:var(--accent)}
.controls{position:sticky;top:10px;z-index:5;background:color-mix(in srgb, var(--card) 92%, transparent);
  backdrop-filter:blur(10px);border:1px solid var(--line);border-radius:16px;padding:14px 16px;
  margin:14px 0 16px;display:grid;gap:12px;box-shadow:var(--shadow)}
.keep{font-weight:600;color:var(--pass);background:var(--chip);border-radius:999px;padding:7px 12px;display:inline-block}
.sliders{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}
label{display:flex;justify-content:space-between;font-size:12px;color:var(--mute)}
label b{color:var(--ink)}
input[type=range]{width:100%;accent-color:var(--accent)}
.cmd{font-family:ui-monospace,Menlo,monospace;font-size:12px;background:color-mix(in srgb,var(--accent) 8%, var(--card));
  border:1px solid var(--line);padding:9px 11px;border-radius:10px;overflow:auto;color:var(--ink)}
.grid2{display:grid;grid-template-columns:1.4fr .9fr;gap:14px;margin-bottom:14px}
.grid3{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:16px;padding:12px;box-shadow:var(--shadow)}
.card-head{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:6px}
.card-head .title{font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:var(--mute)}
.dl-group{display:flex;gap:6px}
.dl{appearance:none;border:1px solid var(--line);background:var(--paper);color:var(--ink);
  border-radius:8px;padding:5px 9px;font:inherit;font-size:11.5px;cursor:pointer;font-weight:600}
.dl:hover{border-color:var(--accent);color:var(--accent)}
.fig-legend{margin:8px 2px 0;padding:0;color:var(--mute);font-size:12px;line-height:1.45;font-style:italic}
.plot{width:100%;height:380px}
.plot.box{height:340px}
h2{margin:20px 0 10px;font-size:15px;letter-spacing:.08em;text-transform:uppercase;color:var(--mute)}
table{width:100%;border-collapse:collapse;font-size:13px}
th,td{padding:8px 8px;border-bottom:1px solid var(--line);text-align:left}
th{color:var(--mute);font-size:11px;text-transform:uppercase;letter-spacing:.04em}
td.num{text-align:right;font-variant-numeric:tabular-nums}
@media(max-width:960px){.sliders,.grid2,.grid3{grid-template-columns:1fr}}
</style>
</head>
<body>
<div class="wrap">
  <div class="top">
    <div>
      <h1>QC filter explorer</h1>
      <p class="lede">
        Inspect QC before filtering. Sliders preview proposed cell.filter() cuts
        (pass/fail + keep table). Download PNG/SVG for any plot at the current settings.
      </p>
    </div>
    <button class="btn" id="theme-btn" type="button">Theme</button>
  </div>

  <div class="controls">
    <div><span class="keep" id="keep-pill"> - </span>
      <button class="dl" id="reset" type="button" style="margin-left:8px">Reset defaults</button></div>
    <div class="sliders">
      <div><label>min.genes <b id="v-g">200</b></label><input id="min-genes" type="range"/></div>
      <div><label>min.umis <b id="v-u">500</b></label><input id="min-umis" type="range"/></div>
      <div><label>max.mito <b id="v-m">0.10</b></label><input id="max-mito" type="range"/></div>
    </div>
    <div class="cmd" id="cmd"></div>
  </div>

  <h2>Core QC</h2>
  <div class="grid2">
    <div class="card">
      <div class="card-head"><span class="title">UMI  x  genes</span>
        <div class="dl-group">
          <button class="dl" type="button" data-dl="sc-mito" data-fmt="png">PNG</button>
          <button class="dl" type="button" data-dl="sc-mito" data-fmt="svg">SVG</button>
        </div>
      </div>
      <div class="plot" id="sc-mito"></div>
      <p class="fig-legend">UMI  x  genes coloured by mitochondrial fraction  -  joint view of complexity, depth, and mito load.</p>
    </div>
    <div class="card" id="sample-card">
      <div class="card-head"><span class="title" id="sample-title">Pass by sample</span></div>
      <div style="max-height:380px;overflow:auto">
        <table id="tbl"><thead><tr><th id="th-sample">Sample</th><th class="num">Cells</th><th class="num">Pass preview</th><th class="num">% pass</th></tr></thead><tbody></tbody></table>
      </div>
      <p class="fig-legend" id="sample-why">Cells and preview-pass counts by sample under the current slider cuts (data still unfiltered).</p>
    </div>
  </div>
  <div class="grid3">
    <div class="card">
      <div class="card-head"><span class="title">Genes</span>
        <div class="dl-group">
          <button class="dl" type="button" data-dl="box-g" data-fmt="png">PNG</button>
          <button class="dl" type="button" data-dl="box-g" data-fmt="svg">SVG</button>
        </div>
      </div>
      <div class="plot box" id="box-g"></div>
    </div>
    <div class="card">
      <div class="card-head"><span class="title">UMIs</span>
        <div class="dl-group">
          <button class="dl" type="button" data-dl="box-u" data-fmt="png">PNG</button>
          <button class="dl" type="button" data-dl="box-u" data-fmt="svg">SVG</button>
        </div>
      </div>
      <div class="plot box" id="box-u"></div>
    </div>
    <div class="card">
      <div class="card-head"><span class="title">Mito %</span>
        <div class="dl-group">
          <button class="dl" type="button" data-dl="box-m" data-fmt="png">PNG</button>
          <button class="dl" type="button" data-dl="box-m" data-fmt="svg">SVG</button>
        </div>
      </div>
      <div class="plot box" id="box-m"></div>
    </div>
  </div>

  <h2>Filter preview</h2>
  <div class="card">
    <div class="card-head"><span class="title">Pass / fail preview</span>
      <div class="dl-group">
        <button class="dl" type="button" data-dl="sc-pass" data-fmt="png">PNG</button>
        <button class="dl" type="button" data-dl="sc-pass" data-fmt="svg">SVG</button>
      </div>
    </div>
    <div class="plot" id="sc-pass"></div>
    <p class="fig-legend">UMI  x  genes coloured by pass/fail under the current slider cuts. Dotted lines: min.genes and min.umis.</p>
  </div>
</div>
<script>
window.SODA_DASH = ',
    json_payload,
';
</script>
<script>
(function(){
  const D=window.SODA_DASH, F=D.qc_filters, qc=D.qc;
  const root=document.documentElement;
  const saved=localStorage.getItem("ic2-qc-theme");
  if(saved) root.setAttribute("data-theme", saved);
  document.getElementById("theme-btn").onclick=()=>{
    const next=root.getAttribute("data-theme")==="dark"?"light":"dark";
    root.setAttribute("data-theme", next);
    localStorage.setItem("ic2-qc-theme", next);
    update();
  };

  function theme(){
    const dark=root.getAttribute("data-theme")==="dark";
    return {
      ink: dark?"#e8f1f2":"#14242a",
      mute: dark?"#9bb0b6":"#5a6b72",
      card: dark?"#122226":"#ffffff",
      pass: dark?"#2dd4bf":"#0f766e",
      fail: dark?"#3a4a4e":"#c5cdd1",
      cut: dark?"#fb7185":"#be123c",
      grid: dark?"rgba(255,255,255,.06)":"rgba(20,36,42,.05)"
    };
  }
  const cols=["#4C78A8","#F58518","#54A24B","#E45756","#72B7B2","#B279A2","#FF9DA6","#9D755D","#BAB0AC","#1F7A6C","#6B5B95","#D4A373"];
  // Always at least one group (single-sample / missing labels  ->  "all")
  let order = (D.condition_order||[]).map(String).filter(s=>s && s!=="null" && s!=="undefined");
  if(!order.length){
    const seen={};
    qc.forEach(r=>{ const c=String(r.condition_id||"all"); if(c) seen[c]=1; });
    order=Object.keys(seen);
  }
  if(!order.length) order=["all"];
  const nCond=order.length;
  const single=nCond===1;

  if(single){
    document.getElementById("sample-title").textContent="Pass summary";
    document.getElementById("th-sample").textContent="Group";
    document.getElementById("sample-why").textContent="Cells and preview-pass counts under the current slider cuts (data still unfiltered).";
  }

  function baseLayout(){
    const t=theme();
    return {paper_bgcolor:t.card, plot_bgcolor:t.card,
      font:{family:"Avenir Next, Segoe UI, Helvetica Neue, Arial, sans-serif", color:t.ink, size:12},
      margin:{l:52,r:14,t:36,b:44}, hovermode:"closest"};
  }
  const cfg={responsive:true, displayModeBar:true, displaylogo:false,
    modeBarButtonsToRemove:["lasso2d","select2d"],
    toImageButtonOptions:{format:"png", height:800, width:1200, scale:2}};

  function setup(id,key,def){const el=document.getElementById(id),r=F.ranges[key];el.min=r.min;el.max=r.max;el.step=r.step;el.value=def;return el;}
  const sG=setup("min-genes","min_genes",F.defaults.min_genes);
  const sU=setup("min-umis","min_umis",F.defaults.min_umis);
  const sM=setup("max-mito","max_mito",F.defaults.max_mito);
  function filt(){return{min_genes:+sG.value,min_umis:+sU.value,max_mito:+sM.value};}
  function ok(r,f){return r.nGenes>=f.min_genes&&r.UMIs>=f.min_umis&&r["mito.percent"]<=f.max_mito;}
  function ax(){const t=theme(); return {tickfont:{color:t.mute,size:10},showgrid:true,gridcolor:t.grid,zeroline:false,linecolor:t.ink,ticks:"outside"};}
  function axOpen(){const t=theme(); return {tickfont:{color:t.mute,size:10},showgrid:false,zeroline:false,linecolor:t.ink,ticks:"outside", mirror:false};}
  function hexAlpha(hex,a){
    const h=String(hex).replace("#","");
    const r=parseInt(h.slice(0,2),16), g=parseInt(h.slice(2,4),16), b=parseInt(h.slice(4,6),16);
    return "rgba("+r+","+g+","+b+","+a+")";
  }
  function mean(vals){return vals.length?vals.reduce((s,v)=>s+v,0)/vals.length:null;}

  // Match PDF: violin + white box + jitter + mean diamond (ggplot .ic2_qc_dash_box)
  function boxTraces(metric){
    const t=theme();
    const jit=Math.max(0.03, Math.min(0.14, 0.45/nCond));
    const boxW=Math.max(0.18, Math.min(0.42, 1.15/nCond));
    const pt=nCond>8?2.8:nCond>4?3.2:3.8;
    const boxFill=root.getAttribute("data-theme")==="dark"?"#1a2e32":"#ffffff";
    const violins=order.map((cond,i)=>{
      const vals=qc.filter(r=>String(r.condition_id)===String(cond)).map(r=>r[metric]);
      const col=cols[i%cols.length];
      return {
        type:"violin", name:cond,
        x:vals.map(()=>cond), y:vals,
        points:"all", jitter:jit, pointpos:0,
        scalemode:"width", spanmode:"hard", width:0.9,
        fillcolor:hexAlpha(col,0.35),
        line:{color:t.ink,width:0.9},
        marker:{size:pt, opacity:0.45, color:col, line:{width:0}},
        box:{visible:true, width:boxW, fillcolor:boxFill, line:{color:t.ink,width:1}},
        meanline:{visible:false},
        hoverinfo:"y+name",
        side:"both"
      };
    });
    const means={
      type:"scatter", mode:"markers", name:"mean",
      x:order,
      y:order.map(cond=>{
        const vals=qc.filter(r=>String(r.condition_id)===String(cond)).map(r=>r[metric]);
        return mean(vals);
      }),
      marker:{symbol:"diamond", size:9, color:t.ink, line:{width:0}},
      hovertemplate:"mean %{y:.3g}<extra></extra>",
      showlegend:false
    };
    return violins.concat([means]);
  }
  function drawBoxes(f){
    const t=theme();
    const xAngle=nCond>8?60:nCond>5?35:0;
    [["box-g","nGenes","Genes / cell",f.min_genes],
     ["box-u","UMIs","UMIs / cell",f.min_umis],
     ["box-m","mito.percent","Mito % / cell",f.max_mito]].forEach(([el,metric,title,cut])=>{
      Plotly.react(el, boxTraces(metric), Object.assign({}, baseLayout(), {
        title:{text:title+(single?"  |  single sample":"  |  by sample"), font:{size:14}, x:0, xanchor:"left"},
        showlegend:false,
        violingap:0.15,
        violinmode:"group",
        yaxis:Object.assign({title:metric}, axOpen()),
        xaxis:Object.assign({
          type:"category", categoryorder:"array", categoryarray:order,
          tickangle:xAngle
        }, axOpen()),
        shapes:[{type:"line", y0:cut, y1:cut, x0:0, x1:1, xref:"paper", line:{color:t.cut, width:1.3, dash:"dot"}}],
        margin:{l:52,r:12,t:40,b:xAngle?70:48}
      }), Object.assign({}, cfg, {toImageButtonOptions:{format:"png",filename:"QC_"+el,height:800,width:1000,scale:2}}));
    });
  }

  function update(){
    const f=filt(), t=theme();
    document.getElementById("v-g").textContent=Math.round(f.min_genes);
    document.getElementById("v-u").textContent=Math.round(f.min_umis);
    document.getElementById("v-m").textContent=f.max_mito.toFixed(3);
    document.getElementById("cmd").textContent="cell.filter(obj, min.genes="+Math.round(f.min_genes)+", min.umis="+Math.round(f.min_umis)+", max.mito="+(+f.max_mito.toFixed(3))+")";
    const keep=[], fail=[]; qc.forEach(r=>(ok(r,f)?keep:fail).push(r));
    document.getElementById("keep-pill").innerHTML="Preview pass <b>"+keep.length.toLocaleString()+"</b> / "+qc.length.toLocaleString()+" ("+((100*keep.length)/qc.length).toFixed(1)+"%)";
    const tb=document.querySelector("#tbl tbody"); tb.innerHTML="";
    order.forEach(cond=>{
      const all=qc.filter(r=>String(r.condition_id)===String(cond)), k=all.filter(r=>ok(r,f));
      const tr=document.createElement("tr");
      const label=single?(cond==="all"?"all cells":cond):cond;
      tr.innerHTML="<td>"+label+"</td><td class=\'num\'>"+all.length+"</td><td class=\'num\'>"+k.length+"</td><td class=\'num\'>"+(all.length?((100*k.length)/all.length).toFixed(1):" - ")+"%</td>";
      tb.appendChild(tr);
    });

    Plotly.react("sc-mito",[{
      type:"scattergl", mode:"markers",
      x:qc.map(r=>r.nGenes), y:qc.map(r=>r.UMIs),
      customdata:qc.map(r=>r["mito.percent"]),
      text:qc.map(r=>r.condition_id),
      marker:{size:4, opacity:.8, color:qc.map(r=>r["mito.percent"]),
        colorscale: (root.getAttribute("data-theme")==="dark") ?
          [[0,"#334155"],[0.5,"#fb923c"],[1,"#f43f5e"]] :
          [[0,"#e2e8f0"],[0.5,"#fb923c"],[1,"#be123c"]],
        cmin:0, cmax:Math.max(0.05, ...qc.map(r=>r["mito.percent"])),
        colorbar:{title:"mito%", thickness:12, len:.55}},
      hovertemplate:"%{text}<br>genes=%{x}<br>UMIs=%{y}<br>mito=%{customdata:.3f}<extra></extra>"
    }], Object.assign({}, baseLayout(), {
      title:{text:"UMI  x  genes  |  mito colour", font:{size:14}, x:0, xanchor:"left"},
      xaxis:Object.assign({title:"nGenes"}, ax()),
      yaxis:Object.assign({title:"UMIs"}, ax()),
      shapes:[
        {type:"line", x0:f.min_genes, x1:f.min_genes, y0:0, y1:1, yref:"paper", line:{color:t.cut, width:1.2, dash:"dot"}},
        {type:"line", y0:f.min_umis, y1:f.min_umis, x0:0, x1:1, xref:"paper", line:{color:t.cut, width:1.2, dash:"dot"}}
      ]
    }), Object.assign({}, cfg, {toImageButtonOptions:{format:"png",filename:"QC_sc-mito",height:800,width:1200,scale:2}}));

    Plotly.react("sc-pass",[
      {type:"scattergl", mode:"markers", name:"fail preview", x:fail.map(r=>r.nGenes), y:fail.map(r=>r.UMIs), marker:{size:4, color:t.fail, opacity:.7}},
      {type:"scattergl", mode:"markers", name:"pass preview", x:keep.map(r=>r.nGenes), y:keep.map(r=>r.UMIs), marker:{size:4, color:t.pass, opacity:.75}}
    ], Object.assign({}, baseLayout(), {
      title:{text:"UMI  x  genes  |  filter preview", font:{size:14}, x:0, xanchor:"left"},
      xaxis:Object.assign({title:"nGenes"}, ax()),
      yaxis:Object.assign({title:"UMIs"}, ax()),
      legend:{orientation:"h", y:1.12},
      shapes:[
        {type:"line", x0:f.min_genes, x1:f.min_genes, y0:0, y1:1, yref:"paper", line:{color:t.cut, width:1.3, dash:"dot"}},
        {type:"line", y0:f.min_umis, y1:f.min_umis, x0:0, x1:1, xref:"paper", line:{color:t.cut, width:1.3, dash:"dot"}}
      ]
    }), Object.assign({}, cfg, {toImageButtonOptions:{format:"png",filename:"QC_sc-pass",height:800,width:1200,scale:2}}));

    drawBoxes(f);
  }

  document.querySelectorAll("[data-dl]").forEach(btn=>{
    btn.addEventListener("click", ()=>{
      const id=btn.getAttribute("data-dl");
      const fmt=btn.getAttribute("data-fmt")||"png";
      Plotly.downloadImage(id, {
        format:fmt, width:1200, height:800, scale:2,
        filename:"QC_"+id
      });
    });
  });

  ["min-genes","min-umis","max-mito"].forEach(id=>document.getElementById(id).addEventListener("input", update));
  document.getElementById("reset").onclick=()=>{sG.value=F.defaults.min_genes; sU.value=F.defaults.min_umis; sM.value=F.defaults.max_mito; update();};
  update();
})();
</script>
</body>
</html>'
  )
}
