\name{iCellR2-internal}
\alias{cpp_mtx_info}
\alias{cpp_mtx_stream_sorted}
\alias{cpp_mtx_to_csc}
\alias{cpp_sparse_col_stats}
\alias{cpp_sparse_row_accum}
\alias{jaccard_coeff}
\alias{print.iCellR2_assay_handle}
\title{Internal iCellR2 helpers}
\description{
C++ helpers and the assay-handle printer. Not part of the public API.
}
\usage{
cpp_mtx_info(path)
cpp_mtx_to_csc(path, progress = NULL, progress_every = 2000000)
cpp_mtx_stream_sorted(
  path,
  flush_chunk,
  progress = NULL,
  chunk_nnz = 5000000,
  progress_every = 2000000
)
cpp_sparse_col_stats(mat)
cpp_sparse_row_accum(mat, sum, sumsq, nnz)
jaccard_coeff(idx)
\method{print}{iCellR2_assay_handle}(x, ...)
}
\arguments{
\item{path}{Path to a Matrix Market file.}
\item{progress}{Optional progress function.}
\item{progress_every}{How often to report progress.}
\item{flush_chunk}{Function called with each streamed chunk.}
\item{chunk_nnz}{Non-zeros per streamed chunk.}
\item{mat}{A sparse matrix.}
\item{sum, sumsq, nnz}{Numeric vectors updated in place.}
\item{idx}{Neighbor index matrix.}
\item{x}{An assay handle.}
\item{...}{Ignored.}
}
\keyword{internal}
