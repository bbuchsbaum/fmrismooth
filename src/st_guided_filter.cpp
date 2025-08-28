\
#include <Rcpp.h>
#include <algorithm>
#include <cstring>
using namespace Rcpp;

inline size_t idx3(const int x, const int y, const int z,
                   const int nx, const int ny, const int nz) {
  return (size_t)x + (size_t)nx * ((size_t)y + (size_t)ny * (size_t)z);
}

static inline void box1d_x_sum(const double* src, double* dst,
                               const int nx, const int radius) {
  double acc = 0.0;
  for (int i = 0; i <= std::min(radius, nx-1); ++i) acc += src[i];
  dst[0] = acc;
  for (int i = 1; i < nx; ++i) {
    int add = i + radius;
    int sub = i - radius - 1;
    if (add < nx) acc += src[add];
    if (sub >= 0) acc -= src[sub];
    dst[i] = acc;
  }
}

static inline void boxfilter3d_sum(const double* src, double* tmp, double* dst,
                                   const int nx, const int ny, const int nz,
                                   const int radius) {
  for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      const size_t base = idx3(0,y,z,nx,ny,nz);
      box1d_x_sum(src + base, tmp + base, nx, radius);
    }
  }
  for (int z = 0; z < nz; ++z) {
    for (int x = 0; x < nx; ++x) {
      double acc = 0.0;
      for (int y = 0; y <= std::min(radius, ny-1); ++y) acc += tmp[idx3(x,y,z,nx,ny,nz)];
      dst[idx3(x,0,z,nx,ny,nz)] = acc;
      for (int y = 1; y < ny; ++y) {
        int add = y + radius;
        int sub = y - radius - 1;
        if (add < ny) acc += tmp[idx3(x,add,z,nx,ny,nz)];
        if (sub >= 0) acc -= tmp[idx3(x,sub,z,nx,ny,nz)];
        dst[idx3(x,y,z,nx,ny,nz)] = acc;
      }
    }
  }
  for (int y = 0; y < ny; ++y) {
    for (int x = 0; x < nx; ++x) {
      double acc = 0.0;
      for (int z = 0; z <= std::min(radius, nz-1); ++z) acc += dst[idx3(x,y,z,nx,ny,nz)];
      tmp[idx3(x,y,0,nx,ny,nz)] = acc;
      for (int z = 1; z < nz; ++z) {
        int add = z + radius;
        int sub = z - radius - 1;
        if (add < nz) acc += dst[idx3(x,y,add,nx,ny,nz)];
        if (sub >= 0) acc -= dst[idx3(x,y,sub,nx,ny,nz)];
        tmp[idx3(x,y,z,nx,ny,nz)] = acc;
      }
    }
  }
  std::memcpy(dst, tmp, sizeof(double) * (size_t)nx*ny*nz);
}

static inline void wbox1d_sum(const double* x, const double* w,
                              double* sumx, double* sumw,
                              const int T, const int r) {
  double accx = 0.0, accw = 0.0;
  for (int t = 0; t <= std::min(r, T-1); ++t) { accx += w[t] * x[t]; accw += w[t]; }
  sumx[0] = accx; sumw[0] = accw;
  for (int t = 1; t < T; ++t) {
    int add = t + r;
    int sub = t - r - 1;
    if (add < T) { accx += w[add]*x[add]; accw += w[add]; }
    if (sub >= 0) { accx -= w[sub]*x[sub]; accw -= w[sub]; }
    sumx[t] = accx; sumw[t] = accw;
  }
}

// [[Rcpp::export(name = "guided3d_core")]]
SEXP guided3d_core(NumericVector vol, NumericVector guide, LogicalVector mask,
                   int radius, double eps) {
  IntegerVector dims = vol.attr("dim");
  if (dims.size() != 3) stop("vol must be 3D array");
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  IntegerVector d2 = guide.attr("dim");
  if (d2.size()!=3 || d2[0]!=nx || d2[1]!=ny || d2[2]!=nz) stop("guide dims mismatch");
  IntegerVector dm = mask.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");

  const size_t N = (size_t)nx*ny*nz;
  NumericVector out(N);
  out.attr("dim") = IntegerVector::create(nx,ny,nz);

  std::vector<double> I(N), p(N), M(N), II(N), Ip(N);
  std::vector<double> sumI(N), sumII(N), sumIp(N), sump(N), sumM(N);
  std::vector<double> tmp(N), dst(N);

  for (size_t i = 0; i < N; ++i) {
    const bool ok = mask[i];
    const double gi = guide[i];
    const double pi = vol[i];
    const double m  = ok ? 1.0 : 0.0;
    I[i]  = gi * m;
    p[i]  = pi * m;
    M[i]  = m;
    II[i] = gi*gi * m;
    Ip[i] = gi*pi * m;
  }

  boxfilter3d_sum(I.data(),  tmp.data(), sumI.data(),  nx,ny,nz, radius);
  boxfilter3d_sum(II.data(), tmp.data(), sumII.data(), nx,ny,nz, radius);
  boxfilter3d_sum(Ip.data(), tmp.data(), sumIp.data(), nx,ny,nz, radius);
  boxfilter3d_sum(p.data(),  tmp.data(), sump.data(),  nx,ny,nz, radius);
  boxfilter3d_sum(M.data(),  tmp.data(), sumM.data(),  nx,ny,nz, radius);

  std::vector<double> a(N), b(N), sumA(N), sumB(N);
  for (size_t i = 0; i < N; ++i) {
    const double cnt = sumM[i];
    if (cnt <= 0.0) { a[i] = 0.0; b[i] = vol[i]; continue; }
    const double meanI = sumI[i]  / cnt;
    const double meanII= sumII[i] / cnt;
    const double meanIp= sumIp[i] / cnt;
    const double meanp = sump[i]  / cnt;
    const double varI  = std::max(0.0, meanII - meanI*meanI);
    const double covIp = meanIp - meanI*meanp;
    const double ai = covIp / (varI + eps);
    const double bi = meanp - ai * meanI;
    a[i] = ai;
    b[i] = bi;
  }

  for (size_t i = 0; i < N; ++i) { a[i] *= M[i]; b[i] *= M[i]; }
  boxfilter3d_sum(a.data(), tmp.data(), sumA.data(), nx,ny,nz, radius);
  boxfilter3d_sum(b.data(), tmp.data(), sumB.data(), nx,ny,nz, radius);

  for (size_t i = 0; i < N; ++i) {
    const double cnt = sumM[i];
    if (cnt <= 0.0) { out[i] = vol[i]; continue; }
    const double meanA = sumA[i] / cnt;
    const double meanB = sumB[i] / cnt;
    out[i] = meanA * guide[i] + meanB;
  }
  return out;
}

// [[Rcpp::export(name = "temporal_guided1d_core")]]
SEXP temporal_guided1d_core(NumericVector vec4d,
                           SEXP gtemp,
                           NumericVector w,
                           LogicalVector mask3,
                           int radius_t,
                           double eps_t) {
  IntegerVector d = vec4d.attr("dim");
  if (d.size() != 4) stop("vec must be 4D array");
  const int nx=d[0], ny=d[1], nz=d[2], T=d[3];
  if (w.size() != T) stop("w must have length T");
  IntegerVector dm = mask3.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");

  std::vector<double> G(T), sumG(T), sumGG(T), sumW(T);
  bool use_global_guide = false;
  if (!Rf_isNull(gtemp)) {
    NumericVector gt(gtemp);
    if (gt.size() != T) stop("guide_temporal must be length T");
    for (int t=0;t<T;++t) G[t] = gt[t];
    use_global_guide = true;
  }

  NumericVector out(vec4d.size());
  out.attr("dim") = d;

  if (use_global_guide) {
    std::vector<double> GG(T);
    for (int t=0;t<T;++t) GG[t] = G[t]*G[t];
    wbox1d_sum(G.data(),  w.begin(), sumG.data(),  sumW.data(), T, radius_t);
    wbox1d_sum(GG.data(), w.begin(), sumGG.data(), sumW.data(), T, radius_t);
  }

  const size_t V = (size_t)nx*ny*nz;
  for (size_t v = 0; v < V; ++v) {
    if (!mask3[v]) {
      for (int t=0;t<T;++t) out[v + (size_t)V*t] = vec4d[v + (size_t)V*t];
      continue;
    }
    std::vector<double> P(T), sumP(T), sumWloc(T), GP(T), sumGP(T), A(T), B(T), sumA(T), sumB(T);
    for (int t=0;t<T;++t) P[t] = vec4d[v + (size_t)V*t];

    const double* I = nullptr;
    std::vector<double> Iown;
    if (!use_global_guide) {
      Iown = P; I = Iown.data();
      std::vector<double> II(T), sumI(T), sumII(T);
      for (int t=0;t<T;++t) II[t] = I[t]*I[t];
      wbox1d_sum(I,  w.begin(), sumI.data(),  sumWloc.data(), T, radius_t);
      wbox1d_sum(II.data(), w.begin(), sumII.data(), sumWloc.data(), T, radius_t);
      for (int t=0;t<T;++t) GP[t] = I[t]*P[t];
      wbox1d_sum(P.data(),  w.begin(), sumP.data(),  sumWloc.data(), T, radius_t);
      wbox1d_sum(GP.data(), w.begin(), sumGP.data(), sumWloc.data(), T, radius_t);

      for (int t=0;t<T;++t) {
        const double cnt = sumWloc[t];
        if (cnt <= 1e-12) { A[t]=0.0; B[t]=P[t]; continue; }
        const double meanI = sumI[t]/cnt;
        const double meanII= sumII[t]/cnt;
        const double meanP = sumP[t]/cnt;
        const double meanIP= sumGP[t]/cnt;
        const double varI  = std::max(0.0, meanII - meanI*meanI);
        const double covIP = meanIP - meanI*meanP;
        const double a = covIP / (varI + eps_t);
        const double b = meanP - a*meanI;
        A[t]=a; B[t]=b;
      }
      for (int t=0;t<T;++t) { A[t] *= w[t]; B[t] *= w[t]; }
      wbox1d_sum(A.data(), w.begin(), sumA.data(), sumWloc.data(), T, radius_t);
      wbox1d_sum(B.data(), w.begin(), sumB.data(), sumWloc.data(), T, radius_t);
      for (int t=0;t<T;++t) {
        const double cnt = sumWloc[t];
        const double meanA = (cnt>1e-12)? (sumA[t]/cnt) : 0.0;
        const double meanB = (cnt>1e-12)? (sumB[t]/cnt) : P[t];
        out[v + (size_t)V*t] = meanA*I[t] + meanB;
      }
    } else {
      I = G.data();
      for (int t=0;t<T;++t) GP[t] = I[t]*P[t];
      wbox1d_sum(P.data(),  w.begin(), sumP.data(),  sumWloc.data(), T, radius_t);
      wbox1d_sum(GP.data(), w.begin(), sumGP.data(), sumWloc.data(), T, radius_t);

      for (int t=0;t<T;++t) {
        const double cnt = sumW[t];
        if (cnt <= 1e-12) { A[t]=0.0; B[t]=P[t]; continue; }
        const double meanI = sumG[t]/cnt;
        const double meanII= sumGG[t]/cnt;
        const double meanP = sumP[t]/sumWloc[t];
        const double meanIP= sumGP[t]/sumWloc[t];
        const double varI  = std::max(0.0, meanII - meanI*meanI);
        const double covIP = meanIP - meanI*meanP;
        const double a = covIP / (varI + eps_t);
        const double b = meanP - a*meanI;
        A[t]=a; B[t]=b;
      }
      for (int t=0;t<T;++t) { A[t] *= w[t]; B[t] *= w[t]; }
      wbox1d_sum(A.data(), w.begin(), sumA.data(), sumWloc.data(), T, radius_t);
      wbox1d_sum(B.data(), w.begin(), sumB.data(), sumWloc.data(), T, radius_t);
      for (int t=0;t<T;++t) {
        const double cnt = sumWloc[t];
        const double meanA = (cnt>1e-12)? (sumA[t]/cnt) : 0.0;
        const double meanB = (cnt>1e-12)? (sumB[t]/cnt) : P[t];
        out[v + (size_t)V*t] = meanA*I[t] + meanB;
      }
    }
  }
  return out;
}
