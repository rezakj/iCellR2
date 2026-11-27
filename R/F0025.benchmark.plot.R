#' Plot SODA session benchmarks
#'
#' Reads timing from \code{metadata/session_stats.json} (preferred) or
#' \code{benchmark.txt} and builds a stacked report:
#' KPI strip, Gantt timeline, wall-time bars, time-share, RAM delta, throughput.
#'
#' @param x An \code{iCellR2} object, a SODA directory path, a
#'   \code{benchmark.txt} path, or a \code{session_stats.json} path.
#' @param out.file Output path. Default \code{"benchmark.pdf"} under
#'   \code{<soda>/plots/} when a SODA root is known; otherwise the working
#'   directory. Use \code{NULL} to skip writing a file.
#' @param latest.only If \code{TRUE}, keep only the last run of each
#'   function name. Default \code{FALSE} (full pipeline order).
#' @param order \code{"pipeline"} (default, run order) or \code{"time"}
#'   (slowest first) for the wall / RAM bar panels.
#' @param width,height PDF size in inches. Defaults \code{10}  x  \code{14}.
#' @param verbose Progress style.
#'
#' @return Invisibly, a list with \code{data} (tidy table), \code{plots}
#'   (named ggplots), \code{out.file}, and \code{gtable} (stacked layout).
#'
#' @examples
#' \dontrun{
#' benchmark.plot(obj, out.file = "bench.pdf")
#' benchmark.plot("/path/to/proj.soda")
#' }
#'
#' @export
benchmark.plot <- function(
    x,
    out.file = "benchmark.pdf",
    latest.only = FALSE,
    order = c("pipeline", "time"),
    width = 10,
    height = 14,
    verbose = TRUE
) {
  order <- match.arg(order)
  rep <- .ic2_new_reporter("benchmark.plot", verbose, soda_root = NULL, persist = FALSE)

  resolved <- .ic2_resolve_benchmark_source(x)
  soda_root <- resolved$soda_root
  df <- .ic2_benchmark_to_df(resolved$stats)
  if (nrow(df) < 1L) {
    stop("No finished tools found in session stats / benchmark.txt.")
  }
  if (isTRUE(latest.only)) {
    df <- df[!duplicated(df$function_name, fromLast = TRUE), , drop = FALSE]
    df$step <- seq_len(nrow(df))
  }

  .ic2_step(rep, paste0(
    "Tools: ", nrow(df),
    "   |   total wall: ", .ic2_fmt_time(sum(df$wall, na.rm = TRUE))
  ))

  ## Resolve output path before building grobs  -  arrangeGrob/ggplotGrob open
  ## the default device; on headless / broken X11 that fails with X11cairo.
  write_file <- !is.null(out.file) && nzchar(as.character(out.file)[1L])
  if (write_file) {
    out.file <- as.character(out.file)[1L]
    if (!grepl("[/\\\\]", out.file) && !is.null(soda_root)) {
      plot_dir <- file.path(soda_root, "plots")
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
      out.file <- file.path(plot_dir, out.file)
    }
    dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)
  } else {
    out.file <- NULL
  }

  opened <- FALSE
  if (write_file) {
    ## Prefer cairo_pdf when available; fall back to pdf() (no X11 needed)
    if (isTRUE(capabilities("cairo"))) {
      tryCatch({
        grDevices::cairo_pdf(out.file, width = width, height = height)
        opened <- TRUE
      }, error = function(e) NULL)
    }
    if (!opened) {
      grDevices::pdf(out.file, width = width, height = height)
      opened <- TRUE
    }
  } else {
    ## Null PDF device so grob layout works without a display
    grDevices::pdf(NULL, width = width, height = height)
    opened <- TRUE
  }
  if (opened) on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  plots <- .ic2_benchmark_build_plots(df, order = order, stats = resolved$stats)
  g <- .ic2_benchmark_arrange(plots)

  if (write_file) {
    grid::grid.draw(g)
    .ic2_step(rep, paste0("Wrote ", out.file))
  }

  .ic2_finish(rep, extra_lines = if (write_file) {
    paste0("PDF: ", out.file)
  } else {
    "Returned plot objects (no file written)"
  })

  invisible(list(
    data = df,
    plots = plots,
    out.file = out.file,
    gtable = g,
    soda_root = soda_root
  ))
}

#' Resolve benchmark source from object / path
#' @keywords internal
.ic2_resolve_benchmark_source <- function(x) {
  soda_root <- NULL
  stats <- NULL

  if (inherits(x, "iCellR2")) {
    soda_root <- x@storage$root_dir
    if (is.null(soda_root) || !nzchar(as.character(soda_root)[1L])) {
      stop("iCellR2 object has no SODA root_dir.")
    }
    soda_root <- normalizePath(soda_root, winslash = "/", mustWork = FALSE)
    stats <- .ic2_load_session_stats(soda_root)
    return(list(soda_root = soda_root, stats = stats))
  }

  if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
    stop("'x' must be an iCellR2 object or a path to a SODA / benchmark file.")
  }
  path <- normalizePath(x, winslash = "/", mustWork = FALSE)

  if (dir.exists(path)) {
    ## SODA root or parent
    cand_json <- file.path(path, "metadata", "session_stats.json")
    cand_rds <- file.path(path, "metadata", "session_stats.rds")
    cand_bm <- file.path(path, "benchmark.txt")
    if (file.exists(cand_json) || file.exists(cand_rds) || file.exists(cand_bm) ||
        grepl("\\.soda$", path, ignore.case = TRUE)) {
      soda_root <- path
      stats <- .ic2_load_session_stats(soda_root)
      if (length(stats$tools %||% list()) < 1L && file.exists(cand_bm)) {
        stats <- .ic2_parse_benchmark_txt(cand_bm, soda_root = soda_root)
      }
      return(list(soda_root = soda_root, stats = stats))
    }
  }

  if (file.exists(path) && grepl("session_stats\\.(json|rds)$", path, ignore.case = TRUE)) {
    soda_root <- dirname(dirname(path))
    stats <- .ic2_load_session_stats(soda_root)
    if (!is.list(stats) || length(stats$tools %||% list()) < 1L) {
      if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        stats <- tryCatch(readRDS(path), error = function(e) NULL)
      } else {
        stats <- tryCatch(
          jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE),
          error = function(e) NULL
        )
      }
    }
    if (!is.list(stats)) stop("Could not read session_stats: ", path)
    soda_root <- stats$soda_root %||% soda_root
    return(list(soda_root = soda_root, stats = stats))
  }

  if (file.exists(path) && grepl("benchmark\\.txt$", path, ignore.case = TRUE)) {
    soda_root <- dirname(path)
    stats <- .ic2_parse_benchmark_txt(path, soda_root = soda_root)
    return(list(soda_root = soda_root, stats = stats))
  }

  stop("Could not resolve benchmark source from: ", path)
}

#' Session stats tools  ->  tidy data.frame
#' @keywords internal
.ic2_benchmark_to_df <- function(stats) {
  tools <- stats$tools %||% list()
  if (length(tools) < 1L) {
    return(data.frame(
      step = integer(0),
      function_name = character(0),
      wall = numeric(0),
      cpu = numeric(0),
      peak_delta = numeric(0),
      cells_per_sec = numeric(0),
      timestamp = character(0),
      stringsAsFactors = FALSE
    ))
  }
  df <- data.frame(
    step = seq_along(tools),
    function_name = vapply(tools, function(t) as.character(t$function_name %||% "?"), character(1)),
    wall = vapply(tools, function(t) {
      w <- as.numeric(t$wall %||% NA_real_)
      if (!is.finite(w)) 0 else max(0, w)
    }, numeric(1)),
    cpu = vapply(tools, function(t) {
      w <- as.numeric(t$cpu %||% NA_real_)
      if (!is.finite(w)) 0 else max(0, w)
    }, numeric(1)),
    peak_delta = vapply(tools, function(t) {
      w <- as.numeric(t$peak_delta %||% NA_real_)
      if (!is.finite(w)) 0 else max(0, w)
    }, numeric(1)),
    cells_per_sec = vapply(tools, function(t) {
      w <- as.numeric(t$cells_per_sec %||% NA_real_)
      if (!is.finite(w)) NA_real_ else w
    }, numeric(1)),
    timestamp = vapply(tools, function(t) as.character(t$timestamp %||% NA_character_), character(1)),
    stringsAsFactors = FALSE
  )
  df$label <- sprintf("%d. %s", df$step, df$function_name)
  df$wall_share <- if (sum(df$wall) > 0) df$wall / sum(df$wall) else rep(0, nrow(df))
  ## Gantt cumulative starts
  df$gantt_start <- c(0, cumsum(df$wall)[-nrow(df)])
  df$gantt_end <- cumsum(df$wall)
  df
}

#' Parse benchmark.txt into a session-stats-like list
#' @keywords internal
.ic2_parse_benchmark_txt <- function(path, soda_root = NULL) {
  lines <- readLines(path, warn = FALSE)
  ## Match table rows: "   1  | run.pca                | 1.2m       |     12.3 |     +40 |       12,345 | 2026-..."
  re <- paste0(
    "^\\s*(\\d+)\\s*\\|\\s*([^|]+?)\\s*\\|\\s*([^|]+?)\\s*\\|\\s*",
    "([0-9.]+)\\s*\\|\\s*\\+?([-0-9.]+)\\s*\\|\\s*([^|]+?)\\s*\\|\\s*(.*)$"
  )
  tools <- list()
  for (ln in lines) {
    m <- regexec(re, ln, perl = TRUE)
    g <- regmatches(ln, m)[[1]]
    if (length(g) < 8L) next
    wall_txt <- trimws(g[4])
    wall <- .ic2_parse_fmt_time(wall_txt)
    rate_txt <- trimws(g[7])
    rate <- if (
      identical(rate_txt, "\u2014") ||
        identical(rate_txt, " - ") ||
        identical(rate_txt, "-") ||
        !nzchar(rate_txt)
    ) {
      NA_real_
    } else {
      suppressWarnings(as.numeric(gsub(",", "", rate_txt)))
    }
    tools[[length(tools) + 1L]] <- list(
      function_name = trimws(g[3]),
      wall = wall,
      cpu = suppressWarnings(as.numeric(g[5])),
      peak_delta = suppressWarnings(as.numeric(g[6])),
      cells_per_sec = rate,
      timestamp = trimws(g[8])
    )
  }
  list(
    soda_root = soda_root %||% dirname(path),
    n_tools = length(tools),
    total_wall = sum(vapply(tools, function(t) t$wall, numeric(1)), na.rm = TRUE),
    total_cpu = sum(vapply(tools, function(t) t$cpu, numeric(1)), na.rm = TRUE),
    peak_ram_delta = max(0, vapply(tools, function(t) t$peak_delta, numeric(1)), na.rm = TRUE),
    last_step = if (length(tools)) tools[[length(tools)]]$function_name else NA_character_,
    last_finished = if (length(tools)) tools[[length(tools)]]$timestamp else NA_character_,
    tools = tools
  )
}

#' Parse human time strings from benchmark.txt (e.g. 1.2m, 3.4s, 1h 2m)
#' @keywords internal
.ic2_parse_fmt_time <- function(txt) {
  txt <- trimws(as.character(txt)[1L])
  if (!nzchar(txt) || identical(txt, " - ") || identical(txt, "NA")) return(0)
  ## Already numeric seconds?
  n <- suppressWarnings(as.numeric(txt))
  if (is.finite(n)) return(n)
  secs <- 0
  mh <- regmatches(txt, regexec("([0-9.]+)\\s*h", txt, perl = TRUE))[[1]]
  if (length(mh) >= 2L) secs <- secs + as.numeric(mh[2]) * 3600
  mm <- regmatches(txt, regexec("([0-9.]+)\\s*m", txt, perl = TRUE))[[1]]
  if (length(mm) >= 2L) secs <- secs + as.numeric(mm[2]) * 60
  ms <- regmatches(txt, regexec("([0-9.]+)\\s*s", txt, perl = TRUE))[[1]]
  if (length(ms) >= 2L) secs <- secs + as.numeric(ms[2])
  if (!is.finite(secs)) 0 else secs
}

#' Build the five ggplot panels
#' @keywords internal
.ic2_benchmark_build_plots <- function(df, order = "pipeline", stats = NULL) {
  df$label_f <- factor(df$label, levels = rev(df$label)) ## top = first step for Gantt
  df_bar <- df
  if (identical(order, "time")) {
    ord <- order(df_bar$wall, decreasing = TRUE)
    df_bar <- df_bar[ord, , drop = FALSE]
    df_bar$label_f <- factor(df_bar$label, levels = rev(df_bar$label))
  } else {
    df_bar$label_f <- factor(df_bar$label, levels = rev(df_bar$label))
  }

  pal <- grDevices::hcl.colors(nrow(df), palette = "Dark 3")
  names(pal) <- df$label

  theme_bm <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      axis.title = ggplot2::element_text(size = 10),
      legend.position = "none"
    )

  ## KPI text panel
  tot_wall <- sum(df$wall, na.rm = TRUE)
  tot_cpu <- sum(df$cpu, na.rm = TRUE)
  peak_ram <- if (!is.null(stats$peak_ram_delta)) {
    max(0, as.numeric(stats$peak_ram_delta))
  } else {
    max(df$peak_delta, na.rm = TRUE)
  }
  kpi <- data.frame(
    x = 0.5, y = 0.5,
    lab = paste0(
      "iCellR2 benchmark\n",
      "Tools: ", nrow(df),
      "    |    Total wall: ", .ic2_fmt_time(tot_wall),
      "    |    Total CPU: ", sprintf("%.1fs", tot_cpu),
      "    |    Peak RAM \u0394: +", sprintf("%.0f", peak_ram), " MB",
      if (!is.null(stats$last_step) && !is.na(stats$last_step)) {
        paste0("\nLast step: ", stats$last_step,
               if (!is.null(stats$last_finished)) paste0("  (", stats$last_finished, ")") else "")
      } else {
        ""
      }
    ),
    stringsAsFactors = FALSE
  )
  p_kpi <- ggplot2::ggplot(kpi, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab)) +
    ggplot2::geom_text(size = 3.6, lineheight = 1.15, hjust = 0.5, vjust = 0.5) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#F7F7F7", colour = "#DDDDDD")
    )

  ## Gantt
  p_gantt <- ggplot2::ggplot(df, ggplot2::aes(
    y = .data$label_f, xmin = .data$gantt_start, xmax = .data$gantt_end,
    colour = .data$label, fill = .data$label
  )) +
    ggplot2::geom_linerange(linewidth = 5, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(
      title = "Pipeline timeline (Gantt)",
      x = "Cumulative wall time (s)", y = NULL
    ) +
    theme_bm

  ## Wall bars
  df_bar$wall_lab <- vapply(df_bar$wall, .ic2_fmt_time, character(1))
  p_wall <- ggplot2::ggplot(df_bar, ggplot2::aes(
    x = .data$wall, y = .data$label_f, fill = .data$label
  )) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$wall_lab),
      hjust = -0.1, size = 3
    ) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title = paste0("Wall time by step (", order, ")"),
      x = "Wall time (s)", y = NULL
    ) +
    theme_bm

  ## Time share (donut-ish via stacked bar for clarity)
  df_share <- df
  df_share$label_f2 <- factor(df_share$label, levels = df_share$label)
  p_share <- ggplot2::ggplot(df_share, ggplot2::aes(
    x = "share", y = .data$wall_share, fill = .data$label_f2
  )) +
    ggplot2::geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::scale_fill_manual(values = pal, guide = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::labs(
      title = "Where time went (% total wall)",
      x = NULL, y = NULL, fill = NULL
    ) +
    theme_bm +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 7)
    )

  ## RAM
  p_ram <- ggplot2::ggplot(df_bar, ggplot2::aes(
    x = .data$peak_delta, y = .data$label_f, fill = .data$label
  )) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(
      title = "Peak RAM growth (\u0394 MB)",
      x = "RAM \u0394 (MB)", y = NULL
    ) +
    theme_bm

  ## Throughput
  df_rate <- df[!is.na(df$cells_per_sec) & df$cells_per_sec > 0, , drop = FALSE]
  if (nrow(df_rate) > 0L) {
    df_rate$label_f <- factor(df_rate$label, levels = rev(df_rate$label))
    p_rate <- ggplot2::ggplot(df_rate, ggplot2::aes(
      x = .data$cells_per_sec, y = .data$label_f, fill = .data$label
    )) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::scale_x_continuous(labels = function(x) {
        vapply(x, function(z) {
          if (!is.finite(z)) return("")
          .ic2_fmt_num(round(z))
        }, character(1))
      }) +
      ggplot2::labs(
        title = "Throughput (cells / s)",
        x = "Cells per second", y = NULL
      ) +
      theme_bm
  } else {
    p_rate <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0.5, y = 0.5,
        label = "No cells/s recorded for these steps",
        size = 3.5
      ) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::labs(title = "Throughput (cells / s)") +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12))
  }

  list(
    kpi = p_kpi,
    gantt = p_gantt,
    wall = p_wall,
    share = p_share,
    ram = p_ram,
    throughput = p_rate
  )
}

#' Stack benchmark panels with gridExtra
#' @keywords internal
.ic2_benchmark_arrange <- function(plots) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required for benchmark.plot().")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required for benchmark.plot().")
  }
  gridExtra::arrangeGrob(
    plots$kpi,
    plots$gantt,
    plots$wall,
    plots$share,
    plots$ram,
    plots$throughput,
    ncol = 1,
    heights = c(0.7, 1.6, 1.6, 1.1, 1.4, 1.4)
  )
}
