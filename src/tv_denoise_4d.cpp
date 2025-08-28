\
#include <Rcpp.h>
using namespace Rcpp;

inline size_t idx4(int x,int y,int z,int t,int nx,int ny,int nz,int nt){
  return (size_t)x + (size_t)nx * ((size_t)y + (size_t)ny * ((size_t)z + (size_t)nz * (size_t)t));
}

// [[Rcpp::export(name = "tv_denoise_4d_cpp")]]
SEXP tv_denoise_4d_cpp(NumericVector vec4d, LogicalVector mask3,
                       double lambda_s, double lambda_t,
                       int iters, double tau, double sigma, double theta) {

  IntegerVector d = vec4d.attr("dim");
  if (d.size()!=4) stop("vec must be 4D array");
  const int nx=d[0], ny=d[1], nz=d[2], nt=d[3];
  IntegerVector dm = mask3.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");

  const size_t N = (size_t)nx*ny*nz*nt;
  NumericVector x(vec4d.size()), xbar(vec4d.size()), xprev(vec4d.size()), y(vec4d.size());
  std::copy(vec4d.begin(), vec4d.end(), y.begin());
  std::copy(y.begin(), y.end(), x.begin());
  std::copy(x.begin(), x.end(), xbar.begin());

  std::vector<double> px(N, 0.0), py(N, 0.0), pz(N, 0.0), pt(N, 0.0);

  auto m3 = [&](int x,int y,int z)->bool { return mask3[(size_t)x + (size_t)nx*((size_t)y + (size_t)ny*(size_t)z)]; };

  for (int k=0; k<iters; ++k) {
    for (int t=0; t<nt; ++t) {
      for (int z=0; z<nz; ++z) {
        for (int yv=0; yv<ny; ++yv) {
          for (int xw=0; xw<nx; ++xw) {
            const size_t i = idx4(xw,yv,z,t,nx,ny,nz,nt);
            if (!m3(xw,yv,z)) { px[i]=py[i]=pz[i]=pt[i]=0.0; continue; }
            double gx=0.0, gy=0.0, gz=0.0, gt=0.0;
            if (xw+1<nx && m3(xw+1,yv,z))  gx = xbar[idx4(xw+1,yv,z,t,nx,ny,nz,nt)] - xbar[i];
            if (yv+1<ny && m3(xw,yv+1,z))  gy = xbar[idx4(xw,yv+1,z,t,nx,ny,nz,nt)] - xbar[i];
            if (z+1<nz && m3(xw,yv,z+1))   gz = xbar[idx4(xw,yv,z+1,t,nx,ny,nz,nt)] - xbar[i];
            if (t+1<nt)                    gt = xbar[idx4(xw,yv,z,t+1,nx,ny,nz,nt)] - xbar[i];

            px[i] += sigma * gx;
            py[i] += sigma * gy;
            pz[i] += sigma * gz;
            pt[i] += sigma * gt;

            const double bx = lambda_s, by = lambda_s, bz = lambda_s, bt = lambda_t;
            if (std::fabs(px[i]) > bx) px[i] = (px[i]>0 ? bx : -bx);
            if (std::fabs(py[i]) > by) py[i] = (py[i]>0 ? by : -by);
            if (std::fabs(pz[i]) > bz) pz[i] = (pz[i]>0 ? bz : -bz);
            if (std::fabs(pt[i]) > bt) pt[i] = (pt[i]>0 ? bt : -bt);
          }
        }
      }
    }

    std::copy(x.begin(), x.end(), xprev.begin());
    for (int t=0; t<nt; ++t) {
      for (int z=0; z<nz; ++z) {
        for (int yv=0; yv<ny; ++yv) {
          for (int xw=0; xw<nx; ++xw) {
            const size_t i = idx4(xw,yv,z,t,nx,ny,nz,nt);
            if (!m3(xw,yv,z)) { x[i] = y[i]; continue; }
            double div = 0.0;
            if (m3(xw,yv,z))        div += px[i];
            if (xw-1>=0 && m3(xw-1,yv,z)) div -= px[idx4(xw-1,yv,z,t,nx,ny,nz,nt)];
            if (m3(xw,yv,z))        div += py[i];
            if (yv-1>=0 && m3(xw,yv-1,z)) div -= py[idx4(xw,yv-1,z,t,nx,ny,nz,nt)];
            if (m3(xw,yv,z))        div += pz[i];
            if (z-1>=0 && m3(xw,yv,z-1))  div -= pz[idx4(xw,yv,z-1,t,nx,ny,nz,nt)];
            if (m3(xw,yv,z))        div += pt[i];
            if (t-1>=0)             div -= pt[idx4(xw,yv,z,t-1,nx,ny,nz,nt)];

            const double xn = (x[i] + tau*(div) + tau*y[i]) / (1.0 + tau);
            x[i] = xn;
          }
        }
      }
    }
    for (size_t i=0;i<N;++i) xbar[i] = x[i] + theta*(x[i] - xprev[i]);
  }

  NumericVector out(N);
  std::copy(x.begin(), x.end(), out.begin());
  out.attr("dim") = d;
  return out;
}
