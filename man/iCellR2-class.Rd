\name{iCellR2-class}
\docType{class}
\alias{iCellR2-class}
\alias{iCellR2}
\title{iCellR2 object}
\description{
S4 object returned by \code{\link{make.obj}} and reopened with
\code{\link{open.soda}}. Assays stay on disk under the SODA archive.
}
\section{Slots}{
\describe{
\item{\code{project_info}}{Project metadata.}
\item{\code{samples}}{Sample records.}
\item{\code{cell_data}}{Per-cell table.}
\item{\code{modalities}}{Assays and reductions for each modality.}
\item{\code{multimodal}}{Joint maps such as RNA plus ADT or RNA plus ATAC.}
\item{\code{reductions}}{RNA reductions, including \code{pca}.}
\item{\code{graphs}}{Neighbor graphs.}
\item{\code{storage}}{Archive paths.}
\item{\code{commands}}{Command log.}
\item{\code{params}}{Parameters from analysis steps.}
\item{\code{version}}{Object version.}
\item{\code{misc}}{Extra metadata.}
\item{\code{validity}}{Object checks.}
\item{\code{my.freq}}{Cluster-by-condition frequency table from
\code{\link{clust.cond.info}}.}
}
}
\seealso{\code{\link{make.obj}}, \code{\link{open.soda}}}
\keyword{classes}
