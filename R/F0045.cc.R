#' Cell-cycle scores (S and G2M) from the on-disk assay
#'
#' Streams the S and G2M gene sets through \code{\link{i.score}} (Tirosh).
#' Phase is G1 when both scores are at most 0, otherwise the higher positive
#' score. The gene matrix is not loaded. Results are \code{cell_data$score_S},
#' \code{score_G2M}, and \code{Phase}.
#'
#' @param object An iCellR2 object. Alias \code{x}.
#' @param x Alias for \code{object}.
#' @param s.genes,g2m.genes Gene ids. Default: the Tirosh human sets.
#' @param assay.name Assay to score. Default \code{"data"}.
#' @param modality Modality. Default \code{"rna"}.
#' @return The updated object.
#' @examples
#' \dontrun{
#' obj <- cc(obj)
#' table(obj@cell_data$Phase)
#' }
#' @export
cc <- function(
    object = NULL,
    x = NULL,
    s.genes = NULL,
    g2m.genes = NULL,
    assay.name = "data",
    modality = "rna"
) {
  if (is.null(object)) object <- x
  if (!inherits(object, "iCellR2")) stop("object should be an object of class iCellR2")
  if (is.null(s.genes)) s.genes <- .ic2_s_genes()
  if (is.null(g2m.genes)) g2m.genes <- .ic2_g2m_genes()
  object <- i.score(
    object,
    scoring.List = list(S = s.genes, G2M = g2m.genes),
    scoring.method = "tirosh",
    assay.name = assay.name,
    modality = modality
  )
  object@cell_data$Phase <- .ic2_cycle_phase(
    object@cell_data$score_S, object@cell_data$score_G2M
  )
  object <- .ic2_persist_cell_data(object)
  object <- .ic2_log_command(
    object, "cc",
    list(assay.name = assay.name),
    list(phase = "cell_data$Phase")
  )
  try(.ic2_soda_write_manifest(object), silent = TRUE)
  object
}

#' @keywords internal
.ic2_cycle_phase <- function(s_score, g2m_score) {
  s_score <- as.numeric(s_score)
  g2m_score <- as.numeric(g2m_score)
  phase <- rep("G1", length(s_score))
  phase[s_score > g2m_score & s_score > 0] <- "S"
  phase[g2m_score > s_score & g2m_score > 0] <- "G2M"
  phase[is.na(s_score) | is.na(g2m_score)] <- NA_character_
  phase
}

#' Tirosh human S-phase genes (same set used by Seurat cc.genes).
#' @keywords internal
.ic2_s_genes <- function() {
  c(
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2",
    "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"
  )
}

#' Tirosh human G2/M genes.
#' @keywords internal
.ic2_g2m_genes <- function() {
  c(
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
    "KIF20B", "HJURP", "CDCA3", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA",
    "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3",
    "CBX5", "CENPA"
  )
}
