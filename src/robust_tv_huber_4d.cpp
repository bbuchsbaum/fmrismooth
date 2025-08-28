\
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

inline size_t idx4(int x,int y,int z,int t,int nx,int ny,int nz,int nt){
  return (size_t)x + (size_t)nx * ((size_t)y + (size_t)ny * ((size_t)z + (size_t)nz * (size_t)t));
}

// Compute divergence of dual fields p = (px,py,pz,pt) into div at xbar
static inline void divergence(const std::vector<double>& px,
                              const std::vector<double>& py,
                              const std::vector<double>& pz,
                              const std::vector<double>& pt,
                              const LogicalVector& mask3,
                              int nx,int ny,int nz,int nt,
                              std::vector<double>& div) {
  const size_t N3 = (size_t)nx*ny*nz;
  std::fill(div.begin(), div.end(), 0.0);
  auto inmask = [&](int x,int y,int z)->bool {
    return mask3[(size_t)x + (size_t)nx*((size_t)y + (size_t)ny*(size_t)z)];
  };
  for (int t=0; t<nt; ++t) {
    for (int z=0; z<nz; ++z) for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
      const size_t i = idx4(x,y,z,t,nx,ny,nz,nt);
      if (!inmask(x,y,z)) continue;
      double d = 0.0;
      // x
      d += px[i];
      if (x-1>=0 && inmask(x-1,y,z)) d -= px[idx4(x-1,y,z,t,nx,ny,nz,nt)];
      // y
      d += py[i];
      if (y-1>=0 && inmask(x,y-1,z)) d -= py[idx4(x,y-1,z,t,nx,ny,nz,nt)];
      // z
      d += pz[i];
      if (z-1>=0 && inmask(x,y,z-1)) d -= pz[idx4(x,y,z-1,t,nx,ny,nz,nt)];
      // t
      d += pt[i];
      if (t-1>=0 && inmask(x,y,z))   d -= pt[idx4(x,y,z,t-1,nx,ny,nz,nt)];
      div[i] = d;
    }
  }
}

// [[Rcpp::export(name="tv_robust_4d_cpp")]]
SEXP robust_tv_huber_4d_cpp(NumericVector vec4d, LogicalVector mask3,
                            double lambda_s, double lambda_t,
                            double delta, int outer_loops, int iters,
                            NumericVector wt_in) {
  IntegerVector d = vec4d.attr("dim");
  if (d.size()!=4) stop("vec must be 4D array");
  const int nx=d[0], ny=d[1], nz=d[2], nt=d[3];
  IntegerVector dm = mask3.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");
  if (wt_in.size()!=nt) stop("temporal weight vector must have length T");

  const size_t N = (size_t)nx*ny*nz*nt;
  NumericVector y(vec4d.size());
  std::copy(vec4d.begin(), vec4d.end(), y.begin());

  // Initialize x with y
  NumericVector x(vec4d.size());
  std::copy(y.begin(), y.end(), x.begin());
  NumericVector xbar(vec4d.size());
  std::copy(x.begin(), x.end(), xbar.begin());

  // Duals
  std::vector<double> px(N,0.0), py(N,0.0), pz(N,0.0), pt(N,0.0);
  std::vector<double> div(N,0.0);
  std::vector<double> w_data(N,1.0); // Huber weights for data term

  auto inmask = [&](int x,int y,int z)->bool {
    return mask3[(size_t)x + (size_t)nx*((size_t)y + (size_t)ny*(size_t)z)];
  };

  // Constant steps that satisfy tau*sigma*L^2 < 1 with L^2 <= 8
  const double tau   = 1.0 / std::sqrt(8.0);
  const double sigma = 1.0 / std::sqrt(8.0);
  const double theta = 1.0;

  // temporal weights vector
  std::vector<double> w_t(nt,1.0);
  for (int t=0;t<nt;++t) w_t[t] = wt_in[t];

  for (int outer=0; outer<outer_loops; ++outer) {
    // --- Update robust data weights w_data based on residual r = x - y ---
    for (int t=0; t<nt; ++t) {
      for (int z=0; z<nz; ++z) for (int yv=0; yv<ny; ++yv) for (int xv=0; xv<nx; ++xv) {
        const size_t i = idx4(xv,yv,z,t,nx,ny,nz,nt);
        if (!inmask(xv,yv,z)) { w_data[i] = 0.0; continue; }
        double r = x[i] - y[i];
        double a = std::fabs(r);
        if (a <= 1e-12) { w_data[i] = 1.0; }
        else {
          // Huber IRLS weight: psi(r)/r where psi = clamp(r, -delta, +delta)
          double psi = (a <= delta) ? r : (delta * (r > 0 ? 1.0 : -1.0));
          w_data[i] = std::fabs(psi / r); // in (0,1]
        }
      }
    }

    // --- Inner primal-dual iterations with fixed w_data ---
    for (int k=0; k<iters; ++k) {
      // Dual ascent: p += sigma * grad(xbar) ; project componentwise with bounds
      for (int t=0; t<nt; ++t) {
        for (int z=0; z<nz; ++z) for (int yv=0; yv<ny; ++yv) for (int xv=0; xv<nx; ++xv) {
          const size_t i = idx4(xv,yv,z,t,nx,ny,nz,nt);
          if (!inmask(xv,yv,z)) { px[i]=py[i]=pz[i]=pt[i]=0.0; continue; }
          // forward diffs (0 outside mask / boundary)
          double gx=0.0, gy=0.0, gz=0.0, gt=0.0;
          if (xv+1<nx && inmask(xv+1,yv,z))  gx = xbar[idx4(xv+1,yv,z,t,nx,ny,nz,nt)] - xbar[i];
          if (yv+1<ny && inmask(xv,yv+1,z))  gy = xbar[idx4(xv,yv+1,z,t,nx,ny,nz,nt)] - xbar[i];
          if (z+1<nz && inmask(xv,yv,z+1))   gz = xbar[idx4(xv,yv,z+1,t,nx,ny,nz,nt)] - xbar[i];
          if (t+1<nt && inmask(xv,yv,z))     gt = xbar[idx4(xv,yv,z,t+1,nx,ny,nz,nt)] - xbar[i];

          px[i] += sigma * gx;
          py[i] += sigma * gy;
          pz[i] += sigma * gz;
          pt[i] += sigma * gt;

          // Project with anisotropic bounds
          const double bx = lambda_s, by = lambda_s, bz = lambda_s;
          const double bt = lambda_t * w_t[t];
          // clamp
          if (std::fabs(px[i]) > bx) px[i] = (px[i]>0 ? bx : -bx);
          if (std::fabs(py[i]) > by) py[i] = (py[i]>0 ? by : -by);
          if (std::fabs(pz[i]) > bz) pz[i] = (pz[i]>0 ? bz : -bz);
          if (std::fabs(pt[i]) > bt) pt[i] = (pt[i]>0 ? bt : -bt);
        }
      }

      // Compute divergence of p
      divergence(px,py,pz,pt,mask3,nx,ny,nz,nt,div);

      // Primal descent + data proximal with weighted L2
      for (int t=0; t<nt; ++t) {
        for (int z=0; z<nz; ++z) for (int yv=0; yv<ny; ++yv) for (int xv=0; xv<nx; ++xv) {
          const size_t i = idx4(xv,yv,z,t,nx,ny,nz,nt);
          if (!inmask(xv,yv,z)) { x[i] = y[i]; continue; }
          const double wd = w_data[i];
          const double num = x[i] + tau * div[i] + tau * wd * y[i];
          const double den = 1.0 + tau * wd;
          const double xn = num / den;
          x[i] = xn;
        }
      }

      // Extrapolation
      for (size_t i=0;i<N;++i) {
        const double prev_xbar = xbar[i];
        xbar[i] = x[i] + theta * (x[i] - prev_xbar);
      }
    } // inner iters
  } // outer loops

  NumericVector out(N);
  std::copy(x.begin(), x.end(), out.begin());
  out.attr("dim") = d;
  return out;
}

// [[Rcpp::export(name="tv_tukey_4d_cpp")]]
SEXP tv_tukey_4d_cpp(NumericVector vec4d, LogicalVector mask3,
                     double lambda_s, double lambda_t,
                     double cthresh, double alpha,
                     NumericVector wt_edge, 
                     int iters, double tau, double sigma, double theta) {
  // For now, just use Huber loss as a placeholder
  // TODO: Implement proper Tukey bisquare loss
  // Convert cthresh to delta for Huber approximation
  double delta = cthresh / 3.5;  // Rough conversion
  
  // Note: This is calling Huber implementation - proper Tukey needs to be implemented
  return robust_tv_huber_4d_cpp(vec4d, mask3, lambda_s, lambda_t, 
                                 delta, 3, iters, 
                                 wt_edge.size() > 0 ? wt_edge : NumericVector(vec4d.attr("dim"))[3]);
}
