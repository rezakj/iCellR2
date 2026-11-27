#' Sum ATAC peaks into gene activity
#'
#' Each gene is the sum of peaks that overlap its body, extended upstream
#' (default 2000 bp) on the transcribed side. The sum is streamed in cell
#' chunks from the peak count assay and written to
#' \code{modalities/atac/assays/gene.activity}. Fragments are not read.
#' Pass a GTF (\code{gene} rows, \code{gene_name}) or a BED/table with
#' chromosome, start, end, and gene.
#'
#' @param x An iCellR2 object with an ATAC count assay.
#' @param genes GTF path, BED path, or a data.frame with columns
#'   \code{gene}, \code{chromosome}, \code{start}, \code{end}. Optional
#'   \code{strand} (\code{+} or \code{-}).
#' @param upstream Bases upstream of the TSS included in the gene. Default \code{2000}.
#' @param downstream Bases past the TES. Default \code{0}.
#' @param assay.name Peak assay to sum. Default \code{"counts"}.
#' @param new.assay.name Output assay. Default \code{"gene.activity"}.
#' @param chunk.size Cells per read.
#' @param overwrite Replace an existing gene-activity file.
#' @param verbose Progress style.
#' @return The iCellR2 object.
#'
#' @examples
#' \dontrun{
#' obj <- gene.activity(obj, genes = "genes.gtf", upstream = 2000)
#' }
#'
#' @export
gene.activity <- function(
    x = NULL,
    genes = NULL,
    upstream = 2000L,
    downstream = 0L,
    assay.name = "counts",
    new.assay.name = "gene.activity",
    chunk.size = 20000L,
    overwrite = TRUE,
    verbose = TRUE
) {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  if (is.null(genes)) stop("Provide genes (a GTF, BED, or data.frame).")
  if (is.null(x@modalities$atac$assays[[assay.name]])) {
    stop("No ATAC assay '", assay.name, "'. Run add.atac() first.")
  }
  upstream <- as.integer(upstream)[1L]
  downstream <- as.integer(downstream)[1L]
  chunk.size <- as.integer(chunk.size)
  rep <- .ic2_new_reporter("gene.activity", verbose)
  genes <- .ic2_atac_read_genes(genes)
  genes <- .ic2_atac_gene_windows(genes, upstream, downstream)

  handles <- .ic2_get_assay_handles(x, "atac", assay.name)
  handles <- lapply(handles, .ic2_enrich_handle, load_dimnames = TRUE)
  peak_ids <- as.character(handles[[1L]]$feature_ids)
  peaks <- .ic2_parse_peak_id(peak_ids)
  hit <- .ic2_atac_gene_peak_hits(genes, peaks)
  if (!length(hit$gene) || !length(hit$peak)) {
    stop("No peaks overlapped the gene windows.")
  }
  G <- Matrix::sparseMatrix(
    i = hit$gene, j = hit$peak, x = 1,
    dims = c(nrow(genes), length(peak_ids)),
    repr = "C"
  )
  .ic2_step(rep, paste0(
    nrow(genes), " genes   |   ", length(unique(hit$peak)), " peaks overlap a gene"
  ))

  cell_ids <- .ic2_adt_cell_ids(x)
  pieces <- list()
  for (h in handles) {
    starts <- seq.int(1L, h$n_cells, by = chunk.size)
    for (s in starts) {
      e <- min(s + chunk.size - 1L, h$n_cells)
      block <- .ic2_to_dgCMatrix(.ic2_read_cols(h, s:e))
      if (!identical(rownames(block), peak_ids)) {
        block <- block[peak_ids, , drop = FALSE]
      }
      scored <- G %*% block
      colnames(scored) <- h$cell_ids[s:e]
      pieces[[length(pieces) + 1L]] <- .ic2_to_dgCMatrix(scored)
    }
  }
  mat <- do.call(cbind, pieces)
  mat <- mat[, cell_ids, drop = FALSE]
  rownames(mat) <- genes$gene
  dest <- .ic2_soda_assay_path(x, "atac", new.assay.name, ext = ".h5")
  handle <- write_sparse_hdf5(mat, filepath = dest, group = "counts", overwrite = overwrite)
  handle$modality <- "atac"
  x <- .ic2_assign_assay(
    x, modality = "atac", assay.name = new.assay.name,
    handles = list(handle), virtual = FALSE, set.active = FALSE
  )
  root <- x@storage$root_dir
  gpath <- .ic2_soda_write_table(genes, file.path(root, "metadata", "atac_genes"))
  x@params$atac$gene.activity <- list(
    genes = gpath, upstream = upstream, downstream = downstream, source = assay.name
  )
  x <- .ic2_log_command(
    x, "gene.activity",
    list(upstream = upstream, downstream = downstream, n_genes = nrow(mat)),
    list(gene.activity = dest, genes = gpath)
  )
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  .ic2_finish(rep, extra_lines = paste0("Genes: ", nrow(mat)))
  x
}

#' @keywords internal
.ic2_atac_read_genes <- function(genes) {
  if (is.character(genes) && length(genes) == 1L && file.exists(genes)) {
    if (grepl("\\.(gtf|gff|gff3)(\\.gz)?$", genes, ignore.case = TRUE)) {
      return(.ic2_atac_read_gtf(genes))
    }
    tab <- utils::read.delim(genes, sep = "\t", comment.char = "#",
                             header = FALSE, stringsAsFactors = FALSE)
    if (ncol(tab) < 4L) stop("Gene BED needs chromosome, start, end, and gene.")
    ## BED is 0-based half-open. Store 1-based closed.
    return(data.frame(
      chromosome = tab[[1L]],
      start = as.integer(tab[[2L]]) + 1L,
      end = as.integer(tab[[3L]]),
      gene = as.character(tab[[4L]]),
      strand = if (ncol(tab) >= 6L) as.character(tab[[6L]]) else rep("+", nrow(tab)),
      stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(genes)) stop("'genes' must be a GTF, BED, or data.frame.")
  nm <- tolower(names(genes))
  names(genes) <- nm
  chrom <- intersect(c("chromosome", "chrom", "chr", "seqnames"), nm)[1L]
  start <- intersect(c("start", "chromstart"), nm)[1L]
  end <- intersect(c("end", "chromend", "stop"), nm)[1L]
  gene <- intersect(c("gene", "gene_name", "gene_id", "name"), nm)[1L]
  if (anyNA(c(chrom, start, end, gene))) {
    stop("Gene table needs gene, chromosome, start, and end.")
  }
  strand <- if ("strand" %in% nm) as.character(genes$strand) else rep("+", nrow(genes))
  data.frame(
    chromosome = as.character(genes[[chrom]]),
    start = as.integer(genes[[start]]),
    end = as.integer(genes[[end]]),
    gene = as.character(genes[[gene]]),
    strand = strand,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.ic2_atac_read_gtf <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else path
  lines <- readLines(con, warn = FALSE)
  if (inherits(con, "connection") && isOpen(con)) close(con)
  lines <- lines[!startsWith(lines, "#")]
  if (!length(lines)) stop("GTF is empty: ", path)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  keep <- lengths(parts) >= 9L & vapply(parts, function(z) z[[3L]] == "gene", logical(1))
  parts <- parts[keep]
  if (!length(parts)) stop("GTF has no gene rows: ", path)
  attr_of <- function(text, key) {
    m <- regexec(paste0(key, " \"([^\"]+)\""), text)
    r <- regmatches(text, m)[[1L]]
    if (length(r) == 2L) r[[2L]] else NA_character_
  }
  data.frame(
    chromosome = vapply(parts, `[[`, "", 1L),
    start = as.integer(vapply(parts, `[[`, "", 4L)),
    end = as.integer(vapply(parts, `[[`, "", 5L)),
    strand = vapply(parts, `[[`, "", 7L),
    gene = vapply(parts, function(z) {
      nm <- attr_of(z[[9L]], "gene_name")
      if (is.na(nm)) nm <- attr_of(z[[9L]], "gene_id")
      nm
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.ic2_atac_gene_windows <- function(genes, upstream, downstream) {
  ok <- !is.na(genes$gene) & nzchar(genes$gene) & !is.na(genes$start) & !is.na(genes$end)
  genes <- genes[ok, , drop = FALSE]
  if (!nrow(genes)) stop("No genes left after parsing.")
  genes <- genes[!duplicated(genes$gene), , drop = FALSE]
  minus <- genes$strand %in% c("-", "-1")
  tss_pad <- ifelse(minus, downstream, upstream)
  tes_pad <- ifelse(minus, upstream, downstream)
  genes$start <- pmax(1L, genes$start - tss_pad)
  genes$end <- genes$end + tes_pad
  genes
}

#' @keywords internal
.ic2_atac_gene_peak_hits <- function(genes, peaks) {
  gi <- integer(0)
  pi <- integer(0)
  bin <- 100000L
  for (chr in unique(genes$chromosome)) {
    g <- which(genes$chromosome == chr)
    p <- which(peaks$chromosome == chr & !is.na(peaks$start))
    if (!length(g) || !length(p)) next
    ps <- peaks$start[p]
    pe <- peaks$end[p]
    width <- max(1L, max(pe - ps))
    buckets <- split(seq_along(p), floor(ps / bin))
    for (ii in g) {
      gs <- genes$start[[ii]]
      ge <- genes$end[[ii]]
      b1 <- floor((gs - width) / bin)
      b2 <- floor(ge / bin)
      cand <- integer(0)
      for (bb in seq.int(b1, b2)) {
        hit <- buckets[[as.character(bb)]]
        if (!is.null(hit)) cand <- c(cand, hit)
      }
      if (!length(cand)) next
      ov <- pmin(ge, pe[cand]) - pmax(gs, ps[cand]) + 1L
      ok <- cand[ov > 0L]
      if (!length(ok)) next
      gi <- c(gi, rep.int(ii, length(ok)))
      pi <- c(pi, p[ok])
    }
  }
  list(gene = gi, peak = pi)
}
