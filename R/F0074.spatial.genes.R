#' Spatially variable genes from the neighbor graph
#'
#' Moran's I for each requested gene. The gene is read in cell chunks; the
#' edge list is applied once. Pass \code{genes}, or a TSV with one gene per line.
#'
#' @param x An iCellR2 object after \code{\link{spatial.graph}}.
#' @param genes Gene ids, or a file path.
#' @param assay.name Assay. Default \code{"data"}.
#' @return The object. The table is \code{metadata/spatial/svg.tsv}.
#' @examples
#' \dontrun{
#' obj <- spatial.genes(obj, genes = c("CD3E", "MS4A1"))
#' }
#' @export
spatial.genes <- function(x = NULL, genes = NULL, assay.name = "data") {
  if (!inherits(x, "iCellR2")) stop("x should be an object of class iCellR2")
  genes <- .ic2_spatial_gene_list(genes)
  root <- x@storage$root_dir
  epath <- x@graphs$spatial_knn %||% file.path(root, "graphs", "spatial_knn", "edges.tsv")
  if (!file.exists(epath)) stop("No spatial graph. Run spatial.graph() first.")
  edges <- utils::read.delim(epath, sep = "\t", stringsAsFactors = FALSE)
  rows <- lapply(genes, function(g) {
    vals <- .ic2_spatial_gene(x, g, assay.name)
    .ic2_spatial_moran(g, vals, edges)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$moran_i, decreasing = TRUE), , drop = FALSE]
  out <- file.path(root, "metadata", "spatial")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  path <- .ic2_vdj_write_tsv(tab, file.path(out, "svg.tsv"))
  x <- .ic2_log_command(x, "spatial.genes", list(n_genes = nrow(tab)), list(svg = path))
  try(.ic2_soda_write_manifest(x), silent = TRUE)
  x
}
.ic2_spatial_gene_list <- function(genes) {
  if (is.null(genes)) stop("Provide genes: a character vector or a file with one gene per line.")
  if (is.character(genes) && length(genes) == 1L && file.exists(genes)) {
    genes <- readLines(genes, warn = FALSE)
  }
  genes <- unique(as.character(genes))
  genes <- genes[nzchar(genes)]
  if (length(genes) < 1L) stop("No genes to score.")
  genes
}

#' @keywords internal
.ic2_spatial_moran <- function(gene, vals, edges) {
  ids <- names(vals)
  z <- as.numeric(vals)
  z[!is.finite(z)] <- 0
  mu <- mean(z)
  z <- z - mu
  den <- sum(z * z)
  i <- match(edges$from, ids)
  j <- match(edges$to, ids)
  ok <- !is.na(i) & !is.na(j)
  i <- i[ok]
  j <- j[ok]
  num <- if (length(i)) sum(z[i] * z[j]) else 0
  w <- length(i)
  n <- length(z)
  moran <- if (den > 0 && w > 0) (n / w) * num / den else NA_real_
  data.frame(gene = gene, moran_i = moran, n_edges = w, stringsAsFactors = FALSE)
}
