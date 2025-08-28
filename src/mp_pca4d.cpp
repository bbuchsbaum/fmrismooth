\
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
using namespace Rcpp;

inline size_t idx4(int x,int y,int z,int t,int nx,int ny,int nz,int nt){
  return (size_t)x + (size_t)nx * ((size_t)y + (size_t)ny * ((size_t)z + (size_t)nz * (size_t)t));
}
inline size_t idx3(int x,int y,int z,int nx,int ny,int nz){
  return (size_t)x + (size_t)nx * ((size_t)y + (size_t)ny * (size_t)z);
}

static double estimate_sigma_patch(const double* block, const int V, const int Tw) {
  std::vector<double> sdv(V, 0.0);
  for (int v=0; v<V; ++v) {
    const double* ts = block + (size_t)Tw * v;
    if (Tw < 3) { sdv[v] = 0.0; continue; }
    std::vector<double> diffs(Tw-1);
    for (int t=0; t<Tw-1; ++t) diffs[t] = ts[t+1] - ts[t];
    double med=0.0;
    { std::vector<double> tmp = diffs; size_t mid = tmp.size()/2;
      std::nth_element(tmp.begin(), tmp.begin()+mid, tmp.end()); med = tmp[mid]; }
    for (double &d: diffs) d = std::fabs(d - med);
    double mad=0.0;
    { std::vector<double> tmp = diffs; size_t mid = tmp.size()/2;
      std::nth_element(tmp.begin(), tmp.begin()+mid, tmp.end()); mad = tmp[mid]; }
    sdv[v] = 1.4826 * mad / std::sqrt(2.0);
  }
  std::sort(sdv.begin(), sdv.end());
  double med = sdv[sdv.size()/2];
  return med;
}

static bool svd_thin(std::vector<double>& A, int m, int n,
                     std::vector<double>& U, std::vector<double>& S, std::vector<double>& VT) {
  const int p = std::min(m,n);
  U.assign((size_t)m*p, 0.0); S.assign((size_t)p, 0.0); VT.assign((size_t)p*n, 0.0);
  char jobz = 'S'; int lda = m, ldu = m, ldvt = p; int lwork = -1, info = 0;
  std::vector<double> work(1); std::vector<int> iwork(8*p);
  int jobz_len = 1;  // Length of jobz character string
  F77_CALL(dgesdd)(&jobz, &m, &n, A.data(), &lda, S.data(),
                   U.data(), &ldu, VT.data(), &ldvt, work.data(), &lwork, iwork.data(), &info, jobz_len);
  if (info != 0) return false;
  lwork = (int)work[0]; work.assign((size_t)lwork, 0.0);
  F77_CALL(dgesdd)(&jobz, &m, &n, A.data(), &lda, S.data(),
                   U.data(), &ldu, VT.data(), &ldvt, work.data(), &lwork, iwork.data(), &info, jobz_len);
  return info == 0;
}

// [[Rcpp::export(name="mppca4d_core")]]
SEXP mppca4d_core(NumericVector vec4d, LogicalVector mask3,
                  IntegerVector patch, int tw, IntegerVector stride,
                  double sigma) {
  IntegerVector d = vec4d.attr("dim");
  if (d.size()!=4) stop("vec must be 4D array");
  const int nx=d[0], ny=d[1], nz=d[2], nt=d[3];
  IntegerVector dm = mask3.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");

  const int px = patch[0], py = patch[1], pz = patch[2];
  if (px<1 || py<1 || pz<1) stop("patch sizes must be >=1");
  tw = std::max(1, std::min(tw, nt));
  const int sx = stride[0], sy = stride[1], sz = stride[2], st = stride[3];

  // Accumulate reconstructed signals and weights; initialize accumulator to zero.
  NumericVector out(vec4d.size());
  std::vector<double> acc(out.size(), 0.0);
  std::vector<double> wts((size_t)nx*ny*nz*nt, 0.0);

  for (int z0=0; z0<nz; z0+=sz) {
    int z1 = std::min(z0 + pz, nz);
    for (int y0=0; y0<ny; y0+=sy) {
      int y1 = std::min(y0 + py, ny);
      for (int x0=0; x0<nx; x0+=sx) {
        int x1 = std::min(x0 + px, nx);
        bool anymask=false;
        for (int z=z0; z<z1 && !anymask; ++z)
          for (int y=y0; y<y1 && !anymask; ++y)
            for (int x=x0; x<x1; ++x)
              if (mask3[idx3(x,y,z,nx,ny,nz)]) anymask=true;
        if (!anymask) continue;

        for (int t0=0; t0<nt; t0+=st) {
          int t1 = std::min(t0 + tw, nt);
          const int Tw = t1 - t0;
          const int V  = (x1-x0)*(y1-y0)*(z1-z0);
          if (V==0 || Tw<2) continue;

          std::vector<double> Y((size_t)Tw*V, 0.0);
          std::vector<double> means(V, 0.0);
          int v = 0;
          for (int z=z0; z<z1; ++z) for (int y=y0; y<y1; ++y) for (int x=x0; x<x1; ++x, ++v) {
            if (!mask3[idx3(x,y,z,nx,ny,nz)]) continue;
            double mu=0.0;
            for (int t=t0; t<t1; ++t) mu += vec4d[idx4(x,y,z,t,nx,ny,nz,nt)];
            mu /= (double)Tw;
            means[v]=mu;
            for (int t=0; t<Tw; ++t) {
              Y[(size_t)t + (size_t)Tw*v] = vec4d[idx4(x,y,z,t0+t,nx,ny,nz,nt)] - mu;
            }
          }

          double sig = R_IsNA(sigma) ? estimate_sigma_patch(Y.data(), V, Tw) : sigma;

          std::vector<double> U, S, VT;
          if (Tw <= V) {
            if (!svd_thin(Y, Tw, V, U, S, VT)) continue;
            double gamma = (double)Tw / (double)V;
            double s_thresh = sig * (1.0 + std::sqrt(gamma)) * std::sqrt((double)V);
            int k=0; for (; k<(int)S.size(); ++k) if (S[k] <= s_thresh) break;
            if (k==0) continue;
            std::vector<double> Yhat((size_t)Tw*V, 0.0);
            for (int i=0;i<Tw;++i) for (int j=0;j<V;++j) {
              double sum=0.0; for (int r=0;r<k;++r)
                sum += U[(size_t)i + (size_t)Tw*r] * S[r] * VT[(size_t)r + (size_t)k*j];
              Yhat[(size_t)i + (size_t)Tw*j] = sum;
            }
            int v2 = 0;
            for (int z=z0; z<z1; ++z) for (int y=y0; y<y1; ++y) for (int x=x0; x<x1; ++x, ++v2) {
              if (!mask3[idx3(x,y,z,nx,ny,nz)]) continue;
              for (int t=0;t<Tw;++t) {
                const size_t idx = idx4(x,y,z,t0+t,nx,ny,nz,nt);
                acc[idx] += Yhat[(size_t)t + (size_t)Tw*v2] + means[v2]; wts[idx] += 1.0;
              }
            }
          } else {
            std::vector<double> At((size_t)V*Tw, 0.0);
            for (int vcol=0; vcol<V; ++vcol) for (int t=0;t<Tw;++t)
              At[(size_t)vcol + (size_t)V*t] = Y[(size_t)t + (size_t)Tw*vcol];
            if (!svd_thin(At, V, Tw, U, S, VT)) continue;
            double gamma = (double)V / (double)Tw;
            double s_thresh = sig * (1.0 + std::sqrt(gamma)) * std::sqrt((double)Tw);
            int k=0; for (; k<(int)S.size(); ++k) if (S[k] <= s_thresh) break;
            if (k==0) continue;
            std::vector<double> Yhat((size_t)Tw*V, 0.0);
            for (int vcol=0; vcol<V; ++vcol) for (int t=0;t<Tw;++t) {
              double sum=0.0; for (int r=0;r<k;++r)
                sum += U[(size_t)vcol + (size_t)V*r] * S[r] * VT[(size_t)r + (size_t)k*t];
              Yhat[(size_t)t + (size_t)Tw*vcol] = sum;
            }
            int v2 = 0;
            for (int z=z0; z<z1; ++z) for (int y=y0; y<y1; ++y) for (int x=x0; x<x1; ++x, ++v2) {
              if (!mask3[idx3(x,y,z,nx,ny,nz)]) continue;
              for (int t=0;t<Tw;++t) {
                const size_t idx = idx4(x,y,z,t0+t,nx,ny,nz,nt);
                acc[idx] += Yhat[(size_t)t + (size_t)Tw*v2] + means[v2]; wts[idx] += 1.0;
              }
            }
          }
        }
      }
    }
  }

  // Default to original signal where no patches contributed; otherwise average accumulated reconstructions.
  std::copy(vec4d.begin(), vec4d.end(), out.begin());
  for (size_t i=0;i<out.size();++i) if (wts[i] > 0.0) out[i] = acc[i] / wts[i];
  out.attr("dim") = IntegerVector::create(d[0],d[1],d[2],d[3]);
  return out;
}
