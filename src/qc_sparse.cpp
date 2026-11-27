#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_sparse_col_stats(const S4& mat) {
  IntegerVector p = mat.slot("p");
  NumericVector x = mat.slot("x");

  int ncol = p.size() - 1;
  NumericVector col_sums(ncol);
  IntegerVector nnz(ncol);

  for (int j = 0; j < ncol; ++j) {
    int start = p[j];
    int end = p[j + 1];

    double s = 0.0;
    int k = end - start;

    for (int idx = start; idx < end; ++idx) {
      s += x[idx];
    }

    col_sums[j] = s;
    nnz[j] = k;
  }

  return Rcpp::List::create(
    Rcpp::Named("col_sums") = col_sums,
    Rcpp::Named("nnz") = nnz
  );
}

// Accumulate row sum / sumsq / nnz from a dgCMatrix chunk (0-based i indices)
// [[Rcpp::export]]
void cpp_sparse_row_accum(
    const S4& mat,
    NumericVector sum,
    NumericVector sumsq,
    IntegerVector nnz
) {
  IntegerVector i = mat.slot("i");
  NumericVector x = mat.slot("x");
  int n = x.size();
  int nrow = sum.size();

  for (int k = 0; k < n; ++k) {
    int r = i[k];
    if (r < 0 || r >= nrow) continue;
    double v = x[k];
    sum[r] += v;
    sumsq[r] += v * v;
    nnz[r] += 1;
  }
}
