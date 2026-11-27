// MTX (Matrix Market coordinate) -> CSC for iCellR2
// Supports plain .mtx and gzipped .mtx.gz via zlib.

#include <Rcpp.h>
#include <zlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdio>

using namespace Rcpp;

namespace {

bool ends_with_ci(const std::string& s, const char* suffix) {
  const size_t n = std::strlen(suffix);
  if (s.size() < n) return false;
  for (size_t i = 0; i < n; ++i) {
    char a = s[s.size() - n + i];
    char b = suffix[i];
    if (a >= 'A' && a <= 'Z') a = static_cast<char>(a - 'A' + 'a');
    if (b >= 'A' && b <= 'Z') b = static_cast<char>(b - 'A' + 'a');
    if (a != b) return false;
  }
  return true;
}

class MtxReader {
public:
  explicit MtxReader(const std::string& path) : gz_(false), gzf_(Z_NULL), ok_(false) {
    if (ends_with_ci(path, ".gz")) {
      gz_ = true;
      gzf_ = gzopen(path.c_str(), "rb");
      ok_ = (gzf_ != Z_NULL);
    } else {
      ifs_.open(path.c_str(), std::ios::in);
      ok_ = ifs_.good();
    }
  }

  ~MtxReader() {
    if (gz_ && gzf_ != Z_NULL) gzclose(gzf_);
  }

  bool ok() const { return ok_; }

  bool getline(std::string& line) {
    line.clear();
    if (!ok_) return false;
    if (gz_) {
      char buf[1 << 16];
      if (gzgets(gzf_, buf, sizeof(buf)) == Z_NULL) return false;
      line.assign(buf);
      if (!line.empty() && line[line.size() - 1] == '\n') line.resize(line.size() - 1);
      if (!line.empty() && line[line.size() - 1] == '\r') line.resize(line.size() - 1);
      return true;
    }
    if (!std::getline(ifs_, line)) return false;
    if (!line.empty() && line[line.size() - 1] == '\r') line.resize(line.size() - 1);
    return true;
  }

private:
  bool gz_;
  gzFile gzf_;
  std::ifstream ifs_;
  bool ok_;
};

static void skip_header(MtxReader& r, long long& nr, long long& nc, long long& nnz) {
  std::string line;
  if (!r.getline(line) || line.find("MatrixMarket") == std::string::npos) {
    stop("Not a Matrix Market file (missing banner).");
  }
  if (line.find("coordinate") == std::string::npos) {
    stop("Only coordinate Matrix Market files are supported.");
  }
  while (r.getline(line)) {
    if (line.empty()) continue;
    if (!line.empty() && line[0] == '%') continue;
    long long a = 0, b = 0, c = 0;
    if (std::sscanf(line.c_str(), "%lld %lld %lld", &a, &b, &c) != 3) {
      stop("Invalid MTX size line: '%s'", line.c_str());
    }
    nr = a;
    nc = b;
    nnz = c;
    if (nr < 1 || nc < 1 || nnz < 0) {
      stop("Invalid MTX dimensions/nnz.");
    }
    return;
  }
  stop("Unexpected end of MTX header.");
}

static void call_progress(Nullable<Function>& progress, double cur, double tot) {
  if (progress.isNotNull()) {
    Function f(progress.get());
    f(cur, tot);
  }
}

static long long prog_step(double progress_every) {
  if (progress_every < 1.0) return 2000000LL;
  return static_cast<long long>(progress_every);
}

} // namespace

//' Read Matrix Market header only
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List cpp_mtx_info(std::string path) {
  MtxReader r(path);
  if (!r.ok()) stop("Cannot open MTX file: '%s'", path.c_str());
  long long nr = 0, nc = 0, nnz = 0;
  skip_header(r, nr, nc, nnz);
  return List::create(
    _["nrow"] = static_cast<double>(nr),
    _["ncol"] = static_cast<double>(nc),
    _["nnz"] = static_cast<double>(nnz)
  );
}

//' Parse MTX/MTX.GZ into full CSC vectors (high-RAM fast path)
//'
//' @param path Path to \code{.mtx} or \code{.mtx.gz}.
//' @param progress Optional R function \code{function(cur, tot)}.
//' @param progress_every Progress / interrupt cadence (entries).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List cpp_mtx_to_csc(
    std::string path,
    Nullable<Function> progress = R_NilValue,
    double progress_every = 2000000.0
) {
  MtxReader r1(path);
  if (!r1.ok()) stop("Cannot open MTX file: '%s'", path.c_str());

  long long nr = 0, nc = 0, nnz = 0;
  skip_header(r1, nr, nc, nnz);

  if (nnz > static_cast<long long>(R_XLEN_T_MAX)) {
    stop("MTX nnz exceeds R vector length limit.");
  }

  const R_xlen_t nnz_x = static_cast<R_xlen_t>(nnz);
  const R_xlen_t nc_x = static_cast<R_xlen_t>(nc);

  IntegerVector i_out(nnz_x);
  NumericVector x_out(nnz_x);
  NumericVector p_out(nc_x + 1);

  const long long every = prog_step(progress_every);
  call_progress(progress, 0.0, static_cast<double>(nnz));

  bool sorted = true;
  long long pos = 0;
  long long col_open = 0;
  p_out[0] = 0.0;

  std::string line;
  long long nread = 0;
  while (nread < nnz && r1.getline(line)) {
    if (line.empty() || line[0] == '%') continue;
    long long rr = 0, cc = 0;
    double xx = 0.0;
    if (std::sscanf(line.c_str(), "%lld %lld %lf", &rr, &cc, &xx) != 3) {
      stop("Invalid MTX triplet at entry %lld.", nread + 1);
    }
    if (rr < 1 || rr > nr || cc < 1 || cc > nc) {
      stop("MTX index out of range at entry %lld.", nread + 1);
    }
    const long long c0 = cc - 1;
    if (c0 < col_open) {
      sorted = false;
      break;
    }
    while (col_open < c0) {
      ++col_open;
      p_out[col_open] = static_cast<double>(pos);
    }
    i_out[pos] = static_cast<int>(rr - 1);
    x_out[pos] = xx;
    ++pos;
    ++nread;
    if (every > 0 && (nread % every) == 0) {
      checkUserInterrupt();
      call_progress(progress, static_cast<double>(nread), static_cast<double>(nnz));
    }
  }

  if (sorted) {
    if (nread != nnz) {
      stop("MTX ended early: read %lld of %lld entries.", nread, nnz);
    }
    while (col_open < nc) {
      ++col_open;
      p_out[col_open] = static_cast<double>(pos);
    }
    call_progress(progress, static_cast<double>(nnz), static_cast<double>(nnz));
    return List::create(
      _["i"] = i_out,
      _["p"] = p_out,
      _["x"] = x_out,
      _["nrow"] = static_cast<double>(nr),
      _["ncol"] = static_cast<double>(nc),
      _["nnz"] = static_cast<double>(nnz),
      _["sorted"] = true
    );
  }

  // Unsorted: two-pass count then fill (still full-RAM)
  MtxReader r2(path);
  if (!r2.ok()) stop("Cannot reopen MTX file for second pass: '%s'", path.c_str());
  long long nr2 = 0, nc2 = 0, nnz2 = 0;
  skip_header(r2, nr2, nc2, nnz2);

  std::vector<long long> counts(static_cast<size_t>(nc), 0);
  nread = 0;
  const double tot2 = static_cast<double>(nnz) * 2.0;
  while (nread < nnz && r2.getline(line)) {
    if (line.empty() || line[0] == '%') continue;
    long long rr = 0, cc = 0;
    double xx = 0.0;
    if (std::sscanf(line.c_str(), "%lld %lld %lf", &rr, &cc, &xx) != 3) {
      stop("Invalid MTX triplet on count pass at entry %lld.", nread + 1);
    }
    if (cc < 1 || cc > nc) stop("MTX column out of range.");
    counts[static_cast<size_t>(cc - 1)]++;
    ++nread;
    if (every > 0 && (nread % every) == 0) {
      checkUserInterrupt();
      call_progress(progress, static_cast<double>(nread), tot2);
    }
  }
  if (nread != nnz) stop("MTX ended early on count pass.");

  p_out[0] = 0.0;
  for (R_xlen_t j = 0; j < nc_x; ++j) {
    p_out[j + 1] = p_out[j] + static_cast<double>(counts[static_cast<size_t>(j)]);
  }

  std::vector<long long> next(static_cast<size_t>(nc));
  for (R_xlen_t j = 0; j < nc_x; ++j) {
    next[static_cast<size_t>(j)] = static_cast<long long>(p_out[j]);
  }

  MtxReader r3(path);
  if (!r3.ok()) stop("Cannot reopen MTX file for fill pass.");
  skip_header(r3, nr2, nc2, nnz2);
  nread = 0;
  while (nread < nnz && r3.getline(line)) {
    if (line.empty() || line[0] == '%') continue;
    long long rr = 0, cc = 0;
    double xx = 0.0;
    if (std::sscanf(line.c_str(), "%lld %lld %lf", &rr, &cc, &xx) != 3) {
      stop("Invalid MTX triplet on fill pass at entry %lld.", nread + 1);
    }
    const long long c0 = cc - 1;
    const long long dest = next[static_cast<size_t>(c0)]++;
    i_out[dest] = static_cast<int>(rr - 1);
    x_out[dest] = xx;
    ++nread;
    if (every > 0 && (nread % every) == 0) {
      checkUserInterrupt();
      call_progress(progress, static_cast<double>(nnz + nread), tot2);
    }
  }

  call_progress(progress, tot2, tot2);
  return List::create(
    _["i"] = i_out,
    _["p"] = p_out,
    _["x"] = x_out,
    _["nrow"] = static_cast<double>(nr),
    _["ncol"] = static_cast<double>(nc),
    _["nnz"] = static_cast<double>(nnz),
    _["sorted"] = false
  );
}

//' Stream column-sorted MTX into R via chunk flushes (low-RAM path)
//'
//' \code{flush_chunk(i, x)} is called with 0-based row indices. On success
//' returns CSC column pointer vector \code{p}. If the MTX is not column-sorted,
//' returns \code{sorted = FALSE} and does not call flush (caller should
//' fall back to the full-memory parser).
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List cpp_mtx_stream_sorted(
    std::string path,
    Function flush_chunk,
    Nullable<Function> progress = R_NilValue,
    double chunk_nnz = 5000000.0,
    double progress_every = 2000000.0
) {
  MtxReader r(path);
  if (!r.ok()) stop("Cannot open MTX file: '%s'", path.c_str());

  long long nr = 0, nc = 0, nnz = 0;
  skip_header(r, nr, nc, nnz);

  const long long chunk =
    (chunk_nnz < 1.0) ? 5000000LL : static_cast<long long>(chunk_nnz);
  const long long every = prog_step(progress_every);

  NumericVector p_out(static_cast<R_xlen_t>(nc) + 1);
  p_out[0] = 0.0;

  IntegerVector i_buf(static_cast<R_xlen_t>(chunk));
  NumericVector x_buf(static_cast<R_xlen_t>(chunk));
  long long buf_n = 0;
  long long out_nnz = 0;
  long long col_open = 0;
  long long nread = 0;
  bool sorted = true;

  auto flush = [&]() {
    if (buf_n < 1) return;
    IntegerVector i_view = i_buf[Range(0, buf_n - 1)];
    NumericVector x_view = x_buf[Range(0, buf_n - 1)];
    flush_chunk(i_view, x_view);
    out_nnz += buf_n;
    buf_n = 0;
  };

  call_progress(progress, 0.0, static_cast<double>(nnz));

  std::string line;
  while (nread < nnz && r.getline(line)) {
    if (line.empty() || line[0] == '%') continue;
    long long rr = 0, cc = 0;
    double xx = 0.0;
    if (std::sscanf(line.c_str(), "%lld %lld %lf", &rr, &cc, &xx) != 3) {
      stop("Invalid MTX triplet at entry %lld.", nread + 1);
    }
    if (rr < 1 || rr > nr || cc < 1 || cc > nc) {
      stop("MTX index out of range at entry %lld.", nread + 1);
    }
    const long long c0 = cc - 1;
    if (c0 < col_open) {
      sorted = false;
      break;
    }
    while (col_open < c0) {
      flush();
      ++col_open;
      p_out[col_open] = static_cast<double>(out_nnz);
    }
    i_buf[buf_n] = static_cast<int>(rr - 1);
    x_buf[buf_n] = xx;
    ++buf_n;
    ++nread;
    if (buf_n >= chunk) flush();
    if (every > 0 && (nread % every) == 0) {
      checkUserInterrupt();
      call_progress(progress, static_cast<double>(nread), static_cast<double>(nnz));
    }
  }

  if (!sorted) {
    return List::create(
      _["nrow"] = static_cast<double>(nr),
      _["ncol"] = static_cast<double>(nc),
      _["nnz"] = static_cast<double>(nnz),
      _["sorted"] = false,
      _["p"] = R_NilValue
    );
  }

  if (nread != nnz) {
    stop("MTX ended early: read %lld of %lld entries.", nread, nnz);
  }
  flush();
  while (col_open < nc) {
    ++col_open;
    p_out[col_open] = static_cast<double>(out_nnz);
  }
  call_progress(progress, static_cast<double>(nnz), static_cast<double>(nnz));

  return List::create(
    _["nrow"] = static_cast<double>(nr),
    _["ncol"] = static_cast<double>(nc),
    _["nnz"] = static_cast<double>(out_nnz),
    _["sorted"] = true,
    _["p"] = p_out
  );
}
