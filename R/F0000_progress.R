#' Progress, timing, and output messaging for iCellR2
#'
#' Console messaging is unchanged. When a SODA root is attached, the same
#' banners / steps / finish summaries (not progress-bar redraws) are also
#' appended to \code{<soda>/iCellR2.log}. A compact timing table is kept in
#' \code{<soda>/benchmark.txt} and rewritten after every finished
#' \strong{persistent} function. Tools whose names end in \code{.plot}
#' (e.g. \code{cluster.plot}) are console-only and do not touch the log or
#' benchmark files.
#'
#' @keywords internal
#' @name iCellR2-progress
NULL

.ic2_verbose_match <- function(verbose = TRUE) {
  if (isTRUE(verbose)) return("full")
  if (isFALSE(verbose)) return("silent")
  if (is.character(verbose) && length(verbose) == 1L && !is.na(verbose)) {
    v <- tolower(verbose)
    if (v %in% c("full", "compact", "silent")) return(v)
  }
  stop("'verbose' must be TRUE/FALSE or one of 'full', 'compact', 'silent'.")
}

.ic2_utf8 <- function() {
  isTRUE(l10n_info()[["UTF-8"]])
}

.ic2_chars <- function() {
  if (.ic2_utf8()) {
    list(
      border = "=",
      fill = "#",
      empty = ".",
      pipe = "|",
      check = "x"
    )
  } else {
    list(
      border = "=",
      fill = "#",
      empty = "-",
      pipe = "|",
      check = "OK"
    )
  }
}

.ic2_mem_mb <- function() {
  g <- tryCatch(gc(reset = FALSE), error = function(e) NULL)
  if (is.null(g) || !is.matrix(g)) return(NA_real_)
  if (ncol(g) >= 2L) {
    return(as.numeric(sum(g[, 2L], na.rm = TRUE)))
  }
  NA_real_
}

.ic2_fmt_time <- function(secs) {
  secs <- as.numeric(secs)
  if (!is.finite(secs) || secs < 0) return("NA")
  if (secs < 60) return(sprintf("%0.1fs", secs))
  m <- floor(secs / 60)
  s <- round(secs - 60 * m)
  if (m < 60) return(sprintf("%d:%02d", as.integer(m), as.integer(s)))
  h <- floor(m / 60)
  m <- m - 60 * h
  sprintf("%d:%02d:%02d", as.integer(h), as.integer(m), as.integer(s))
}

.ic2_fmt_num <- function(x) {
  x <- as.numeric(x)[1L]
  if (!is.finite(x)) return("NA")
  if (abs(x) <= .Machine$integer.max) {
    return(format(as.integer(round(x)), big.mark = ",", scientific = FALSE, trim = TRUE))
  }
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

.ic2_hline <- function(ch = NULL, width = 48L) {
  if (is.null(ch)) ch <- .ic2_chars()$border
  paste(rep(ch, width), collapse = "")
}

.ic2_soda_log_path <- function(soda_root) {
  if (is.null(soda_root) || !is.character(soda_root) || length(soda_root) != 1L) {
    return(NULL)
  }
  if (!nzchar(soda_root) || is.na(soda_root)) return(NULL)
  file.path(soda_root, "iCellR2.log")
}

.ic2_benchmark_path <- function(soda_root) {
  if (is.null(soda_root) || !is.character(soda_root) || length(soda_root) != 1L) {
    return(NULL)
  }
  if (!nzchar(soda_root) || is.na(soda_root)) return(NULL)
  file.path(soda_root, "benchmark.txt")
}

#' Human-readable per-function timing table for a SODA archive
#' @keywords internal
.ic2_write_benchmark_txt <- function(soda_root, st = NULL) {
  path <- .ic2_benchmark_path(soda_root)
  if (is.null(path)) return(invisible(NULL))
  if (is.null(st)) st <- .ic2_load_session_stats(soda_root)
  soda_norm <- normalizePath(soda_root, winslash = "/", mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  tools <- st$tools %||% list()
  rows <- character(0)
  path_rows <- character(0)
  if (length(tools) > 0L) {
    rows <- vapply(seq_along(tools), function(i) {
      t <- tools[[i]]
      rate <- if (!is.null(t$cells_per_sec) && is.finite(t$cells_per_sec)) {
        .ic2_fmt_num(round(as.numeric(t$cells_per_sec)))
      } else {
        "\u2014"
      }
      sprintf(
        "%4d  | %-22s | %-10s | %8.1f | %+8.0f | %12s | %s",
        i,
        substr(as.character(t$function_name %||% "?"), 1L, 22L),
        .ic2_fmt_time(t$wall),
        as.numeric(t$cpu %||% 0),
        max(0, as.numeric(t$peak_delta %||% 0)),
        rate,
        as.character(t$timestamp %||% "\u2014")
      )
    }, character(1))
    path_rows <- unlist(lapply(seq_along(tools), function(i) {
      t <- tools[[i]]
      lines <- c(
        sprintf(
          "# %4d  %s",
          i,
          as.character(t$function_name %||% "?")
        ),
        sprintf("#        input:   %s", t$input_dir %||% "\u2014"),
        sprintf("#        output:  %s", t$output_dir %||% "\u2014")
      )
      inns <- t$inputs %||% list()
      outs <- t$outputs %||% list()
      if (length(inns) > 0L) {
        lines <- c(lines, "#        inputs:")
        nms <- names(inns) %||% paste0("path", seq_along(inns))
        for (j in seq_along(inns)) {
          lines <- c(lines, sprintf("#          - %s: %s", nms[[j]], inns[[j]]))
        }
      }
      if (length(outs) > 0L) {
        lines <- c(lines, "#        outputs:")
        nms <- names(outs) %||% paste0("path", seq_along(outs))
        for (j in seq_along(outs)) {
          lines <- c(lines, sprintf("#          - %s: %s", nms[[j]], outs[[j]]))
        }
      }
      lines
    }), use.names = FALSE)
  } else {
    rows <- "   (no tools finished yet)"
    path_rows <- "#   (none)"
  }

  sep <- paste0(
    "------+------------------------+------------+----------+----------+--------------+---------------------------"
  )
  hdr <- paste0(
    "   #  | function               | wall       |  CPU (s) | RAM \u0394 MB |     cells/s | finished"
  )

  out <- c(
    "# iCellR2 benchmark.txt \u2014 timing for every function run in this SODA archive",
    paste0("# SODA:    ", soda_norm),
    paste0("# Updated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "#",
    "# Wall = elapsed clock time. CPU = user+system process time.",
    "# RAM \u0394 = peak R memory growth during the call (MB). cells/s when available.",
    "# This file is rewritten after each finished iCellR2 function.",
    "#",
    sep,
    hdr,
    sep,
    rows,
    sep,
    "#",
    "# INPUT / OUTPUT PATHS",
    "# --------------------",
    path_rows,
    "#",
    "# TOTALS",
    sprintf("#   Tools finished: %s", .ic2_fmt_num(st$n_tools %||% length(tools))),
    sprintf("#   Total wall:     %s", .ic2_fmt_time(st$total_wall %||% 0)),
    sprintf("#   Total CPU:      %.1fs", as.numeric(st$total_cpu %||% 0)),
    sprintf("#   Peak RAM \u0394:     +%.0f MB", max(0, as.numeric(st$peak_ram_delta %||% 0))),
    sprintf("#   Last step:      %s", st$last_step %||% "\u2014"),
    sprintf("#   Last finished:  %s", st$last_finished %||% "\u2014"),
    "#"
  )

  tryCatch(
    writeLines(out, con = path, useBytes = FALSE),
    error = function(e) invisible(NULL)
  )
  invisible(path)
}

.ic2_fmt_log_args <- function(arguments) {
  if (is.null(arguments) || length(arguments) == 0L) return("  (none)")
  nms <- names(arguments)
  if (is.null(nms)) nms <- paste0("arg", seq_along(arguments))
  lines <- vapply(seq_along(arguments), function(i) {
    nm <- nms[[i]]
    if (!nzchar(nm)) nm <- paste0("arg", i)
    v <- arguments[[i]]
    if (is.null(v)) {
      val <- "NULL"
    } else if (is.atomic(v)) {
      if (length(v) == 0L) {
        val <- paste0(typeof(v), "(0)")
      } else if (length(v) > 8L) {
        val <- sprintf("%s[%s]", typeof(v), .ic2_fmt_num(length(v)))
      } else if (is.character(v)) {
        val <- paste(encodeString(v, quote = "\""), collapse = ", ")
      } else {
        val <- paste(as.character(v), collapse = ", ")
      }
    } else {
      val <- paste0("<", paste(class(v), collapse = "/"), ">")
    }
    paste0("  ", nm, " = ", val)
  }, character(1))
  paste(lines, collapse = "\n")
}

.ic2_summary_marker <- "======== SESSION SUMMARY ========"

.ic2_session_stats_path <- function(soda_root) {
  file.path(soda_root, "metadata", "session_stats.json")
}

.ic2_load_session_stats <- function(soda_root) {
  path <- .ic2_session_stats_path(soda_root)
  legacy <- file.path(soda_root, "metadata", "session_stats.rds")
  if (file.exists(path)) {
    st <- tryCatch(
      jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE),
      error = function(e) NULL
    )
    if (is.list(st)) return(st)
  }
  if (file.exists(legacy)) {
    st <- tryCatch(readRDS(legacy), error = function(e) NULL)
    if (is.list(st)) {
      ## Migrate off RDS
      tryCatch(.ic2_save_session_stats(soda_root, st), error = function(e) NULL)
      return(st)
    }
  }
  list(
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    soda_root = normalizePath(soda_root, winslash = "/", mustWork = FALSE),
    n_tools = 0L,
    total_wall = 0,
    total_cpu = 0,
    peak_ram_delta = 0,
    last_step = NA_character_,
    last_finished = NA_character_,
    tools = list()
  )
}

.ic2_save_session_stats <- function(soda_root, stats) {
  path <- .ic2_session_stats_path(soda_root)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    jsonlite::write_json(stats, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    .ic2_soda_unlink_rds(file.path(soda_root, "metadata", "session_stats.rds"))
  }, error = function(e) invisible(NULL))
  invisible(path)
}

.ic2_record_session_tool <- function(
    soda_root,
    fn_name,
    wall = NA_real_,
    cpu = NA_real_,
    peak_delta = NA_real_,
    cells_per_sec = NULL,
    input_dir = NULL,
    output_dir = NULL,
    inputs = NULL,
    outputs = NULL
) {
  if (is.null(soda_root) || !nzchar(as.character(soda_root)[1L])) {
    return(invisible(NULL))
  }
  st <- .ic2_load_session_stats(soda_root)
  wall <- as.numeric(wall)
  cpu <- as.numeric(cpu)
  peak_delta <- as.numeric(peak_delta)
  if (!is.finite(wall)) wall <- 0
  if (!is.finite(cpu)) cpu <- 0
  if (!is.finite(peak_delta)) peak_delta <- 0

  st$n_tools <- as.integer(st$n_tools) + 1L
  st$total_wall <- as.numeric(st$total_wall) + wall
  st$total_cpu <- as.numeric(st$total_cpu) + cpu
  st$peak_ram_delta <- max(as.numeric(st$peak_ram_delta), peak_delta, na.rm = TRUE)
  st$last_step <- as.character(fn_name)
  st$last_finished <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  entry <- list(
    function_name = fn_name,
    wall = wall,
    cpu = cpu,
    peak_delta = peak_delta,
    cells_per_sec = cells_per_sec,
    input_dir = if (!is.null(input_dir) && nzchar(as.character(input_dir)[1L])) {
      as.character(input_dir)[1L]
    } else {
      NA_character_
    },
    output_dir = if (!is.null(output_dir) && nzchar(as.character(output_dir)[1L])) {
      as.character(output_dir)[1L]
    } else {
      NA_character_
    },
    inputs = .ic2_io_paths_as_list(inputs),
    outputs = .ic2_io_paths_as_list(outputs),
    timestamp = st$last_finished
  )
  st$tools[[length(st$tools) + 1L]] <- entry
  .ic2_save_session_stats(soda_root, st)
  .ic2_write_benchmark_txt(soda_root, st)
  invisible(st)
}

#' Normalize named path vectors for session stats
#' @keywords internal
.ic2_io_paths_as_list <- function(x) {
  if (is.null(x) || length(x) < 1L) return(list())
  if (is.list(x)) {
    nms <- names(x) %||% paste0("path", seq_along(x))
    out <- lapply(x, function(v) as.character(v)[1L])
    names(out) <- nms
    return(out)
  }
  x <- as.character(x)
  nms <- names(x)
  if (is.null(nms) || !any(nzchar(nms))) {
    nms <- paste0("path", seq_along(x))
  }
  out <- as.list(x)
  names(out) <- nms
  out
}

#' Labels that are archive bookkeeping, not primary I/O
#' @keywords internal
.ic2_io_meta_label <- function(label) {
  grepl(
    "^(Archive|Log|Benchmark|Active assay|QC cols)$",
    as.character(label)[1L],
    ignore.case = TRUE
  )
}

#' Pick a primary path from a named path vector
#' @keywords internal
.ic2_pick_primary_path <- function(paths, soda_root = NULL) {
  if (is.null(paths) || length(paths) < 1L) return(NULL)
  paths <- unlist(paths, use.names = TRUE)
  paths <- paths[nzchar(as.character(paths))]
  if (length(paths) < 1L) return(NULL)
  soda_norm <- if (!is.null(soda_root) && nzchar(as.character(soda_root)[1L])) {
    normalizePath(as.character(soda_root)[1L], winslash = "/", mustWork = FALSE)
  } else {
    NULL
  }
  nms <- names(paths) %||% rep("", length(paths))
  keep <- !vapply(nms, .ic2_io_meta_label, logical(1))
  if (!is.null(soda_norm)) {
    keep <- keep & !(normalizePath(paths, winslash = "/", mustWork = FALSE) == soda_norm)
  }
  cands <- if (any(keep)) paths[keep] else paths
  as.character(cands[[1L]])
}

#' If I/O dirs are still the SODA root, replace with concrete assay/result paths
#' @keywords internal
.ic2_finalize_io_dirs <- function(reporter) {
  soda <- reporter$log_env$soda_root %||% NULL
  soda_norm <- if (!is.null(soda) && nzchar(as.character(soda)[1L])) {
    normalizePath(as.character(soda)[1L], winslash = "/", mustWork = FALSE)
  } else {
    NULL
  }
  is_generic <- function(p) {
    if (is.null(p) || !nzchar(as.character(p)[1L])) return(TRUE)
    if (is.null(soda_norm)) return(FALSE)
    identical(
      normalizePath(as.character(p)[1L], winslash = "/", mustWork = FALSE),
      soda_norm
    )
  }
  prim_in <- .ic2_pick_primary_path(reporter$input_paths, soda_norm)
  prim_out <- .ic2_pick_primary_path(reporter$output_paths, soda_norm)
  if (is_generic(reporter$input_dir) && !is.null(prim_in)) {
    reporter$input_dir <- prim_in
  }
  if (is_generic(reporter$output_dir) && !is.null(prim_out)) {
    reporter$output_dir <- prim_out
  }
  reporter
}

#' Register assay handle path(s) as reporter inputs
#' @keywords internal
.ic2_register_assay_inputs <- function(
    reporter,
    x,
    modality,
    assay.name,
    label = "Assay"
) {
  handles <- tryCatch(
    .ic2_get_assay_handles(x, modality, assay.name),
    error = function(e) NULL
  )
  if (is.null(handles) || length(handles) < 1L) return(reporter)
  for (i in seq_along(handles)) {
    h <- tryCatch(.ic2_enrich_handle(handles[[i]], FALSE), error = function(e) handles[[i]])
    p <- h$path %||% NA_character_
    if (is.null(p) || !nzchar(as.character(p)[1L]) || is.na(p)) next
    lab <- if (length(handles) == 1L) {
      as.character(label)[1L]
    } else {
      paste0(as.character(label)[1L], " ", i)
    }
    reporter <- .ic2_add_input(reporter, lab, p)
  }
  reporter
}

.ic2_session_summary_lines <- function(soda_root) {
  st <- .ic2_load_session_stats(soda_root)
  tool_lines <- character(0)
  if (length(st$tools) > 0L) {
    ## Keep last few tools for a compact processing trail
    keep <- st$tools
    if (length(keep) > 8L) keep <- keep[(length(keep) - 7L):length(keep)]
    tool_lines <- unlist(lapply(keep, function(t) {
      lines <- c(
        sprintf(
          "  - %-18s  wall %s   |   CPU %.1fs   |   RAM %+0.0f MB",
          t$function_name,
          .ic2_fmt_time(t$wall),
          t$cpu,
          max(0, t$peak_delta)
        ),
        sprintf("      in:  %s", t$input_dir %||% " - "),
        sprintf("      out: %s", t$output_dir %||% " - ")
      )
      outs <- t$outputs %||% list()
      if (length(outs) > 0L) {
        nms <- names(outs) %||% paste0("path", seq_along(outs))
        ## Keep summary compact: first 3 concrete outputs
        show <- seq_len(min(3L, length(outs)))
        for (j in show) {
          lines <- c(lines, sprintf("         %s: %s", nms[[j]], outs[[j]]))
        }
        if (length(outs) > 3L) {
          lines <- c(lines, sprintf("         ... +%d more", length(outs) - 3L))
        }
      }
      lines
    }), use.names = FALSE)
    if (length(st$tools) > 8L) {
      tool_lines <- c(
        sprintf("  ... %d earlier tool(s)", length(st$tools) - 8L),
        tool_lines
      )
    }
  } else {
    tool_lines <- "  (none yet)"
  }

  c(
    .ic2_summary_marker,
    sprintf("Updated:       %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("SODA:          %s", st$soda_root %||% soda_root),
    sprintf("Benchmark:     %s", file.path(st$soda_root %||% soda_root, "benchmark.txt")),
    sprintf("Tools finished:%s", paste0(" ", .ic2_fmt_num(st$n_tools))),
    sprintf("Total wall:    %s", .ic2_fmt_time(st$total_wall)),
    sprintf("Total CPU:     %.1fs", as.numeric(st$total_cpu)),
    sprintf("Peak RAM \u0394:    +%.0f MB", max(0, as.numeric(st$peak_ram_delta))),
    sprintf("Last step:     %s", st$last_step %||% "\u2014"),
    sprintf("Last finished: %s", st$last_finished %||% "\u2014"),
    "Processing:",
    tool_lines,
    "=================================",
    ""
  )
}

.ic2_log_body_lines <- function(path) {
  if (!file.exists(path)) return(character(0))
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (length(lines) == 0L) return(character(0))
  idx <- which(grepl("^======== SESSION SUMMARY", lines))
  if (length(idx) > 0L) {
    cut <- idx[[1L]] - 1L
    if (cut < 1L) return(character(0))
    lines <- lines[seq_len(cut)]
  }
  ## Drop trailing blank lines so summary stays tight
  while (length(lines) > 0L && !nzchar(lines[[length(lines)]])) {
    lines <- lines[-length(lines)]
  }
  lines
}

.ic2_write_soda_log <- function(soda_root, body_lines) {
  path <- .ic2_soda_log_path(soda_root)
  if (is.null(path)) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  soda_norm <- normalizePath(soda_root, winslash = "/", mustWork = FALSE)

  if (length(body_lines) == 0L || !any(grepl("^# iCellR2\\.log", body_lines))) {
    header <- c(
      "# iCellR2.log  -  human-readable run log for this SODA archive",
      paste0("# SODA: ", soda_norm),
      paste0("# Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      "#"
    )
    body_lines <- c(header, body_lines)
  }

  summary <- .ic2_session_summary_lines(soda_root)
  ## Avoid duplicate blank lines before summary
  while (length(body_lines) > 0L && !nzchar(body_lines[[length(body_lines)]])) {
    body_lines <- body_lines[-length(body_lines)]
  }
  out <- c(body_lines, "", summary)
  tryCatch(
    writeLines(out, con = path, useBytes = FALSE),
    error = function(e) invisible(NULL)
  )
  invisible(path)
}

.ic2_append_soda_log <- function(soda_root, lines) {
  path <- .ic2_soda_log_path(soda_root)
  if (is.null(path)) return(invisible(NULL))
  lines <- as.character(lines)
  lines <- lines[!is.na(lines)]
  body <- .ic2_log_body_lines(path)
  if (length(lines) > 0L) {
    body <- c(body, lines)
  }
  .ic2_write_soda_log(soda_root, body)
}

.ic2_buf_log <- function(reporter, ...) {
  env <- reporter$log_env
  if (is.null(env) || !isTRUE(reporter$persist %||% TRUE)) {
    return(invisible(NULL))
  }
  if (is.null(env$path) || !nzchar(env$path)) {
    return(invisible(NULL))
  }
  line <- paste0(...)
  env$lines <- c(env$lines, line)
  invisible(NULL)
}

.ic2_flush_log <- function(reporter, record_tool = TRUE) {
  env <- reporter$log_env
  if (is.null(env) || is.null(env$path) || !nzchar(env$path)) {
    return(invisible(NULL))
  }
  if (!isTRUE(reporter$persist %||% TRUE)) {
    env$lines <- character(0)
    return(invisible(NULL))
  }
  soda_root <- env$soda_root
  if (isTRUE(record_tool)) {
    .ic2_record_session_tool(
      soda_root = soda_root,
      fn_name = reporter$fn,
      wall = reporter$stats$wall %||% NA_real_,
      cpu = reporter$stats$cpu %||% NA_real_,
      peak_delta = reporter$peak_delta,
      cells_per_sec = reporter$stats$cells_per_sec,
      input_dir = reporter$input_dir,
      output_dir = reporter$output_dir,
      inputs = reporter$input_paths,
      outputs = reporter$output_paths
    )
  }
  if (length(env$lines) > 0L) {
    .ic2_append_soda_log(soda_root, env$lines)
    env$lines <- character(0)
  } else {
    ## Still refresh trailing summary if stats changed
    .ic2_append_soda_log(soda_root, character(0))
  }
  invisible(env$path)
}

.ic2_attach_soda_log <- function(reporter, soda_root) {
  path <- .ic2_soda_log_path(soda_root)
  if (is.null(path)) return(invisible(reporter))
  soda_root <- normalizePath(soda_root, winslash = "/", mustWork = FALSE)
  path <- file.path(soda_root, "iCellR2.log")
  reporter$log_env$path <- path
  reporter$log_env$soda_root <- soda_root
  if (is.null(reporter$output_dir) || !nzchar(as.character(reporter$output_dir %||% "")[1L])) {
    reporter$output_dir <- soda_root
  }
  if (is.null(reporter$input_dir) || !nzchar(as.character(reporter$input_dir %||% "")[1L])) {
    reporter$input_dir <- soda_root
  }
  ## Ensure metadata + empty log / benchmark scaffold exist
  dir.create(file.path(soda_root, "metadata"), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    .ic2_write_soda_log(soda_root, character(0))
  }
  bm <- .ic2_benchmark_path(soda_root)
  if (!is.null(bm) && !file.exists(bm)) {
    .ic2_write_benchmark_txt(soda_root)
  }
  invisible(reporter)
}

.ic2_msg <- function(reporter, ..., to_log = TRUE) {
  line <- paste0(...)
  if (isTRUE(to_log)) .ic2_buf_log(reporter, line)
  if (!identical(reporter$mode, "silent")) {
    message(line)
  }
  invisible(NULL)
}

.ic2_new_reporter <- function(fn_name, verbose = TRUE, soda_root = NULL, persist = NULL) {
  mode <- .ic2_verbose_match(verbose)
  mem0 <- .ic2_mem_mb()
  ## Plot helpers should not spam iCellR2.log / benchmark.txt
  if (is.null(persist)) {
    persist <- !grepl("\\.plot$", as.character(fn_name)[1L])
  }
  persist <- isTRUE(persist)
  log_env <- new.env(parent = emptyenv())
  log_env$path <- NULL
  log_env$soda_root <- NULL
  log_env$lines <- character(0)
  reporter <- list(
    fn = fn_name,
    mode = mode,
    persist = persist,
    t0 = proc.time(),
    wall0 = Sys.time(),
    mem0 = mem0,
    peak_delta = 0,
    meta = list(),
    inputs = character(0),
    outputs = character(0),
    input_paths = character(0),
    output_paths = character(0),
    input_dir = NULL,
    output_dir = NULL,
    scale_notes = character(0),
    stats = list(),
    chars = .ic2_chars(),
    log_env = log_env,
    bar_open = FALSE,
    bar_last_t = -Inf,
    bar_last_pct = -1L,
    bar_last_width = 0L
  )
  if (!is.null(soda_root) && persist) {
    reporter <- .ic2_attach_soda_log(reporter, soda_root)
    ## Default I/O dirs for SODA-backed tools (overridable per function)
    soda_norm <- normalizePath(soda_root, winslash = "/", mustWork = FALSE)
    reporter$input_dir <- soda_norm
    reporter$output_dir <- soda_norm
  }
  reporter
}

.ic2_touch_mem <- function(reporter) {
  now <- .ic2_mem_mb()
  if (is.finite(now) && is.finite(reporter$mem0)) {
    delta <- now - reporter$mem0
    if (is.finite(delta) && delta > reporter$peak_delta) {
      reporter$peak_delta <- delta
    }
  }
  reporter
}

.ic2_banner_start <- function(reporter, lines = character(0)) {
  hl <- .ic2_hline(reporter$chars$border)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ## Always buffer for iCellR2.log; console only in full mode
  .ic2_buf_log(reporter, paste0("[", stamp, "]"))
  .ic2_buf_log(reporter, hl)
  .ic2_buf_log(reporter, paste0("  iCellR2  ::  ", reporter$fn))
  .ic2_buf_log(reporter, hl)
  for (ln in lines) .ic2_buf_log(reporter, paste0("  ", ln))
  if (length(lines) > 0L) .ic2_buf_log(reporter, "")

  if (identical(reporter$mode, "full")) {
    message(hl)
    message(paste0("  iCellR2  ::  ", reporter$fn))
    message(hl)
    for (ln in lines) message(paste0("  ", ln))
    if (length(lines) > 0L) message("")
  }
  invisible(reporter)
}

.ic2_progress_bar <- function(reporter, current, total, prefix = "", width = 26L,
                              min_interval = 2.0) {
  ## Console-only (too noisy for the log file)
  if (!identical(reporter$mode, "full")) return(invisible(reporter))
  total_n <- as.numeric(total)[1L]
  if (!is.finite(total_n) || total_n < 1) return(invisible(reporter))
  current_n <- as.numeric(current)[1L]
  if (!is.finite(current_n)) current_n <- 0
  current_n <- max(0, min(current_n, total_n))
  frac <- current_n / total_n
  done <- current_n >= total_n
  pct_i <- as.integer(floor(100 * frac + 1e-9))

  ## Throttle: on % advance, ~every 2s heartbeat, or finish (Rcpp fires often)
  now <- proc.time()[["elapsed"]]
  last_t <- reporter$bar_last_t
  if (is.null(last_t) || !is.finite(last_t)) last_t <- -Inf
  last_pct <- reporter$bar_last_pct
  if (is.null(last_pct)) last_pct <- -1L
  pct_advanced <- pct_i > last_pct
  if (done && last_pct >= 100L && !isTRUE(reporter$bar_open)) {
    return(invisible(reporter))
  }
  if (!done && !pct_advanced && (now - last_t) < min_interval) {
    return(invisible(reporter))
  }

  filled <- as.integer(round(frac * width))
  ch <- reporter$chars
  bar <- paste0(
    ch$pipe,
    paste(rep(ch$fill, filled), collapse = ""),
    paste(rep(ch$empty, max(0L, width - filled)), collapse = ""),
    ch$pipe
  )
  elapsed <- as.numeric(difftime(Sys.time(), reporter$wall0, units = "secs"))
  reporter <- .ic2_touch_mem(reporter)
  ram <- if (is.finite(reporter$peak_delta)) {
    sprintf("+%.0f MB", max(0, reporter$peak_delta))
  } else {
    "RAM n/a"
  }
  eta <- ""
  if (is.finite(elapsed) && elapsed > 1 && current_n > 0 && frac < 1) {
    remain <- elapsed * (1 - frac) / frac
    if (is.finite(remain) && remain >= 0) {
      eta <- paste0("  ETA ", .ic2_fmt_time(remain))
    }
  }
  txt <- sprintf(
    "  %s %s %3.0f%%  %s/%s   %s%s   RAM %s",
    prefix, bar, 100 * frac,
    .ic2_fmt_num(current_n), .ic2_fmt_num(total_n),
    .ic2_fmt_time(elapsed), eta, ram
  )

  ## Always one line per update. \\r is unreliable in RStudio / remote sessions
  ## and was gluing ticks into one unreadable line.
  cat(txt, "\n", sep = "")
  utils::flush.console()

  reporter$bar_last_t <- now
  reporter$bar_last_pct <- pct_i
  reporter$bar_last_width <- nchar(txt, type = "width")
  reporter$bar_open <- !done
  invisible(reporter)
}

.ic2_end_progress_bar <- function(reporter) {
  ## Newline mode already ends each tick with \\n; just clear open-bar state.
  reporter$bar_open <- FALSE
  reporter$bar_last_t <- -Inf
  reporter$bar_last_pct <- -1L
  reporter$bar_last_width <- 0L
  invisible(reporter)
}

.ic2_step <- function(reporter, text) {
  reporter <- .ic2_end_progress_bar(reporter)
  .ic2_buf_log(reporter, text)
  if (identical(reporter$mode, "full")) {
    message(text)
  }
  invisible(reporter)
}

.ic2_add_output <- function(reporter, label, path) {
  if (is.null(path) || !nzchar(as.character(path)[1L])) return(reporter)
  path <- as.character(path)[1L]
  if (file.exists(path) || dir.exists(path)) {
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  }
  lab <- as.character(label)[1L]
  reporter$outputs <- c(
    reporter$outputs,
    sprintf("  %-12s %s", paste0(lab, ":"), path)
  )
  if (is.null(reporter$output_paths)) reporter$output_paths <- character(0)
  reporter$output_paths <- c(reporter$output_paths, stats::setNames(path, lab))
  reporter
}

.ic2_add_input <- function(reporter, label, path) {
  if (is.null(path) || !nzchar(as.character(path)[1L])) return(reporter)
  path <- as.character(path)[1L]
  if (file.exists(path) || dir.exists(path)) {
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  }
  lab <- as.character(label)[1L]
  reporter$inputs <- c(
    reporter$inputs,
    sprintf("  %-12s %s", paste0(lab, ":"), path)
  )
  if (is.null(reporter$input_paths)) reporter$input_paths <- character(0)
  reporter$input_paths <- c(reporter$input_paths, stats::setNames(path, lab))
  reporter
}

#' Set primary input/output paths for logging / benchmark (files or dirs)
#' @keywords internal
.ic2_set_io_dirs <- function(reporter, input_dir = NULL, output_dir = NULL) {
  if (!is.null(input_dir) && length(input_dir) >= 1L && nzchar(as.character(input_dir)[1L])) {
    p <- as.character(input_dir)[1L]
    if (file.exists(p) || dir.exists(p)) {
      p <- normalizePath(p, winslash = "/", mustWork = FALSE)
    }
    reporter$input_dir <- p
  }
  if (!is.null(output_dir) && length(output_dir) >= 1L && nzchar(as.character(output_dir)[1L])) {
    p <- as.character(output_dir)[1L]
    if (file.exists(p) || dir.exists(p)) {
      p <- normalizePath(p, winslash = "/", mustWork = FALSE)
    }
    reporter$output_dir <- p
  }
  reporter
}

## Kept for call-site compatibility; notes are no longer printed.
.ic2_add_scale <- function(reporter, notes) {
  invisible(reporter)
}

.ic2_finish <- function(reporter, extra_lines = character(0)) {
  reporter <- .ic2_end_progress_bar(reporter)
  pt <- proc.time() - reporter$t0
  wall <- as.numeric(difftime(Sys.time(), reporter$wall0, units = "secs"))
  reporter <- .ic2_touch_mem(reporter)
  cpu <- pt[["user.self"]] + pt[["sys.self"]]
  ch <- reporter$chars
  ram <- if (is.finite(reporter$peak_delta)) {
    sprintf("+%.0f MB peak", max(0, reporter$peak_delta))
  } else {
    "RAM n/a"
  }

  done_line <- paste0(
    "  Done in ", .ic2_fmt_time(wall),
    "   |   CPU ", sprintf("%.1fs", cpu),
    "   |   ", ram
  )

  reporter$stats$wall <- wall
  reporter$stats$cpu <- cpu

  ## Fill missing I/O dirs from SODA root when attached
  soda_root <- reporter$log_env$soda_root
  if (!is.null(soda_root) && nzchar(as.character(soda_root)[1L])) {
    if (is.null(reporter$input_dir) || !nzchar(as.character(reporter$input_dir)[1L])) {
      reporter$input_dir <- soda_root
    }
    if (is.null(reporter$output_dir) || !nzchar(as.character(reporter$output_dir)[1L])) {
      reporter$output_dir <- soda_root
    }
    bm <- .ic2_benchmark_path(soda_root)
    if (!is.null(bm)) {
      reporter <- .ic2_add_output(reporter, "Benchmark", bm)
    }
  }

  ## Prefer concrete assay/result paths over SODA-root placeholders
  reporter <- .ic2_finalize_io_dirs(reporter)

  in_dir <- reporter$input_dir %||% " - "
  out_dir <- reporter$output_dir %||% " - "

  ## Record finish block in iCellR2.log / benchmark only for persistent tools
  if (isTRUE(reporter$persist %||% TRUE)) {
    .ic2_buf_log(reporter, .ic2_hline("-", 48L))
    .ic2_buf_log(reporter, sprintf("  Input dir:   %s", in_dir))
    .ic2_buf_log(reporter, sprintf("  Output dir:  %s", out_dir))
    if (length(reporter$inputs) > 0L) {
      .ic2_buf_log(reporter, "  Inputs")
      for (o in reporter$inputs) .ic2_buf_log(reporter, o)
    }
    if (length(reporter$outputs) > 0L) {
      .ic2_buf_log(reporter, "  Outputs")
      for (o in reporter$outputs) .ic2_buf_log(reporter, o)
    }
    .ic2_buf_log(reporter, "")
    for (ln in extra_lines) {
      if (nzchar(ln)) .ic2_buf_log(reporter, paste0("  ", ln))
    }
    .ic2_buf_log(reporter, done_line)
    .ic2_buf_log(reporter, .ic2_hline(ch$border, 48L))
    .ic2_flush_log(reporter, record_tool = TRUE)
  }

  if (identical(reporter$mode, "silent")) {
    return(invisible(reporter))
  }

  if (identical(reporter$mode, "compact")) {
    outs <- paste0(
      "\n  Input dir:  ", in_dir,
      "\n  Output dir: ", out_dir,
      if (length(reporter$outputs) > 0L) {
        paste0("\n", paste(reporter$outputs, collapse = "\n"))
      } else {
        ""
      }
    )
    message(sprintf(
      "%s %s  %s  |  CPU %.1fs  |  %s%s",
      reporter$fn, ch$check, .ic2_fmt_time(wall), cpu, ram, outs
    ))
    return(invisible(reporter))
  }

  message(.ic2_hline("-", 48L))
  message(sprintf("  Input dir:   %s", in_dir))
  message(sprintf("  Output dir:  %s", out_dir))
  if (length(reporter$inputs) > 0L) {
    message("  Inputs")
    for (o in reporter$inputs) message(o)
  }
  if (length(reporter$outputs) > 0L) {
    message("  Outputs")
    for (o in reporter$outputs) message(o)
  }
  message("")
  for (ln in extra_lines) {
    if (nzchar(ln)) message(paste0("  ", ln))
  }
  message(done_line)
  message(.ic2_hline(ch$border, 48L))
  invisible(reporter)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
