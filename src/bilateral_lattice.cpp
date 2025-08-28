\
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

struct Node { double v; double w; Node(): v(0.0), w(0.0) {} };

struct VecKey {
  std::vector<int> k; int D;
  bool operator==(VecKey const& other) const {
    if (D != other.D) return false;
    for (int i=0;i<D;++i) if (k[i] != other.k[i]) return false;
    return true;
  }
};
struct VecKeyHash {
  std::size_t operator()(VecKey const& a) const noexcept {
    std::size_t h = 1469598103934665603ULL;
    for (int i=0;i<a.D;++i) { std::size_t v = (std::size_t)(a.k[i] + 2147483647); h ^= v; h *= 1099511628211ULL; }
    return h;
  }
};

inline size_t idx3(int x,int y,int z,int nx,int ny,int nz){
  return (size_t)x + (size_t)nx * ((size_t)y + (size_t)ny * (size_t)z);
}
inline size_t idx4(int x,int y,int z,int t,int nx,int ny,int nz,int nt){
  return (size_t)x + (size_t)nx * ((size_t)y + (size_t)ny * ((size_t)z + (size_t)nz * (size_t)t));
}

static void splat_point(std::unordered_map<VecKey, Node, VecKeyHash>& grid,
                        const std::vector<int>& g0, const std::vector<double>& frac,
                        const int D, const double val) {
  const int corners = 1 << D;
  for (int mask=0; mask<corners; ++mask) {
    VecKey key; key.k.resize(D); key.D = D; double wt = 1.0;
    for (int j=0; j<D; ++j) { const int bit = (mask >> j) & 1; key.k[j] = g0[j] + bit; wt *= (bit ? frac[j] : (1.0 - frac[j])); }
    Node &n = grid[key]; n.v += wt * val; n.w += wt;
  }
}

static double slice_point(const std::unordered_map<VecKey, Node, VecKeyHash>& grid,
                          const std::vector<int>& g0, const std::vector<double>& frac,
                          const int D) {
  const int corners = 1 << D;
  double num = 0.0, den = 0.0;
  for (int mask=0; mask<corners; ++mask) {
    VecKey key; key.k.resize(D); key.D = D; double wt = 1.0;
    for (int j=0; j<D; ++j) { const int bit = (mask >> j) & 1; key.k[j] = g0[j] + bit; wt *= (bit ? frac[j] : (1.0 - frac[j])); }
    auto it = grid.find(key);
    if (it != grid.end()) { num += wt * it->second.v; den += wt * it->second.w; }
  }
  if (den <= 1e-12) return 0.0;
  return num / den;
}

static void blur_axis(std::unordered_map<VecKey, Node, VecKeyHash>& grid, int axis) {
  std::unordered_map<VecKey, Node, VecKeyHash> out;
  out.reserve(grid.size());
  for (auto const& kv : grid) {
    const VecKey& key = kv.first; const Node& n0 = kv.second;
    VecKey kminus = key; kminus.k[axis] -= 1; VecKey kplus  = key; kplus.k[axis]  += 1;
    auto itminus = grid.find(kminus); auto itplus  = grid.find(kplus);
    const Node& nm = (itminus!=grid.end() ? itminus->second : n0);
    const Node& np = (itplus !=grid.end() ? itplus->second  : n0);
    Node nnew; nnew.v = (nm.v + 2.0*n0.v + np.v) * 0.25; nnew.w = (nm.w + 2.0*n0.w + np.w) * 0.25;
    out[key] = nnew;
  }
  grid.swap(out);
}

// [[Rcpp::export(name="bilateral_lattice3d_cpp")]]
SEXP bilateral_lattice3d_cpp(NumericVector vol3d, LogicalVector mask3,
                             double sigma_sp, NumericVector sigma_r,
                             SEXP guide_, SEXP guides_, int blur_iters) {
  IntegerVector d = vol3d.attr("dim");
  if (d.size()!=3) stop("vol must be 3D array");
  const int nx=d[0], ny=d[1], nz=d[2];
  IntegerVector dm = mask3.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");

  NumericVector guide; bool has_guide = !Rf_isNull(guide_);
  if (has_guide) {
    guide = as<NumericVector>(guide_);
    IntegerVector dg = guide.attr("dim");
    if (dg.size()!=3 || dg[0]!=nx || dg[1]!=ny || dg[2]!=nz) stop("guide dims mismatch");
  }
  std::vector<NumericVector> glist;
  if (!Rf_isNull(guides_)) {
    List L(guides_); glist.reserve(L.size());
    for (int i=0;i<L.size();++i) {
      NumericVector g = as<NumericVector>(L[i]);
      IntegerVector dg = g.attr("dim");
      if (dg.size()!=3 || dg[0]!=nx || dg[1]!=ny || dg[2]!=nz) stop("guides dims mismatch");
      glist.push_back(g);
    }
  }
  const int extra = (has_guide?1:0) + (int)glist.size();
  const int D = 3 + extra;
  std::vector<double> sr = as<std::vector<double>>(sigma_r);
  if ((int)sr.size() < extra) sr.resize(extra, sr.empty()?1.0:sr.back());

  std::unordered_map<VecKey, Node, VecKeyHash> grid;
  grid.reserve((size_t)nx*ny*std::max(1, nz/4));

  for (int z=0; z<nz; ++z) for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
    if (!mask3[idx3(x,y,z,nx,ny,nz)]) continue;
    std::vector<double> p; p.reserve(D);
    p.push_back(x / sigma_sp); p.push_back(y / sigma_sp); p.push_back(z / sigma_sp);
    int gi = 0;
    if (has_guide) { p.push_back(guide[idx3(x,y,z,nx,ny,nz)] / sr[gi++]); }
    for (size_t k=0;k<glist.size();++k) p.push_back(glist[k][idx3(x,y,z,nx,ny,nz)] / sr[gi++]);
    std::vector<int> g0(D); std::vector<double> frac(D);
    for (int j=0;j<D;++j) { double u = p[j]; int b = (int)std::floor(u); g0[j]=b; frac[j]=u - (double)b; }
    const double val = vol3d[idx3(x,y,z,nx,ny,nz)];
    splat_point(grid, g0, frac, D, val);
  }
  for (int it=0; it<blur_iters; ++it) for (int ax=0; ax<D; ++ax) blur_axis(grid, ax);

  NumericVector out(no_init(vol3d.size()));
  for (int z=0; z<nz; ++z) for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
    if (!mask3[idx3(x,y,z,nx,ny,nz)]) { out[idx3(x,y,z,nx,ny,nz)] = vol3d[idx3(x,y,z,nx,ny,nz)]; continue; }
    std::vector<double> p; p.reserve(D);
    p.push_back(x / sigma_sp); p.push_back(y / sigma_sp); p.push_back(z / sigma_sp);
    int gi = 0;
    if (has_guide) { p.push_back(guide[idx3(x,y,z,nx,ny,nz)] / sr[gi++]); }
    for (size_t k=0;k<glist.size();++k) p.push_back(glist[k][idx3(x,y,z,nx,ny,nz)] / sr[gi++]);
    std::vector<int> g0(D); std::vector<double> frac(D);
    for (int j=0;j<D;++j) { double u = p[j]; int b = (int)std::floor(u); g0[j]=b; frac[j]=u - (double)b; }
    out[idx3(x,y,z,nx,ny,nz)] = slice_point(grid, g0, frac, D);
  }
  out.attr("dim") = d;
  return out;
}

// [[Rcpp::export(name="bilateral_lattice4d_cpp")]]
SEXP bilateral_lattice4d_cpp(NumericVector vec4d, LogicalVector mask3,
                             double sigma_sp, double sigma_t, NumericVector sigma_r,
                             SEXP guide_spatial_, SEXP guides_, SEXP design_, double sigma_d,
                             int blur_iters) {
  IntegerVector d = vec4d.attr("dim");
  if (d.size()!=4) stop("vec must be 4D array");
  const int nx=d[0], ny=d[1], nz=d[2], nt=d[3];
  IntegerVector dm = mask3.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");

  NumericVector guide; bool has_guide = !Rf_isNull(guide_spatial_);
  if (has_guide) {
    guide = as<NumericVector>(guide_spatial_);
    IntegerVector dg = guide.attr("dim");
    if (dg.size()!=3 || dg[0]!=nx || dg[1]!=ny || dg[2]!=nz) stop("guide dims mismatch");
  }
  std::vector<NumericVector> glist;
  if (!Rf_isNull(guides_)) {
    List L(guides_); glist.reserve(L.size());
    for (int i=0;i<L.size();++i) {
      NumericVector g = as<NumericVector>(L[i]);
      IntegerVector dg = g.attr("dim");
      if (dg.size()!=3 || dg[0]!=nx || dg[1]!=ny || dg[2]!=nz) stop("guides dims mismatch");
      glist.push_back(g);
    }
  }
  bool has_design = !Rf_isNull(design_);
  NumericVector design;
  if (has_design) { design = as<NumericVector>(design_); if ((int)design.size() != nt) stop("design length must equal T"); }

  const int extra = (has_guide?1:0) + (int)glist.size();
  const int D = 3 + (sigma_t>0.0 ? 1 : 0) + (has_design?1:0) + extra;
  std::vector<double> sr = as<std::vector<double>>(sigma_r);
  if ((int)sr.size() < extra) sr.resize(extra, sr.empty()?1.0:sr.back());

  std::unordered_map<VecKey, Node, VecKeyHash> grid;
  grid.reserve((size_t)nx*ny*std::max(1, nz/4));

  for (int t=0; t<nt; ++t) for (int z=0; z<nz; ++z) for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
    if (!mask3[idx3(x,y,z,nx,ny,nz)]) continue;
    std::vector<double> p; p.reserve(D);
    p.push_back(x / sigma_sp); p.push_back(y / sigma_sp); p.push_back(z / sigma_sp);
    if (sigma_t > 0.0) p.push_back(t / sigma_t);
    if (has_design) p.push_back(design[t] / sigma_d);
    int gi = 0;
    if (has_guide) { p.push_back(guide[idx3(x,y,z,nx,ny,nz)] / sr[gi++]); }
    for (size_t k=0;k<glist.size();++k) p.push_back(glist[k][idx3(x,y,z,nx,ny,nz)] / sr[gi++]);
    std::vector<int> g0(D); std::vector<double> frac(D);
    for (int j=0;j<D;++j) { double u = p[j]; int b = (int)std::floor(u); g0[j]=b; frac[j]=u - (double)b; }
    const double val = vec4d[idx4(x,y,z,t,nx,ny,nz,nt)];
    splat_point(grid, g0, frac, D, val);
  }
  for (int it=0; it<blur_iters; ++it) for (int ax=0; ax<D; ++ax) blur_axis(grid, ax);

  NumericVector out(no_init(vec4d.size()));
  for (int t=0; t<nt; ++t) for (int z=0; z<nz; ++z) for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
    if (!mask3[idx3(x,y,z,nx,ny,nz)]) { out[idx4(x,y,z,t,nx,ny,nz,nt)] = vec4d[idx4(x,y,z,t,nx,ny,nz,nt)]; continue; }
    std::vector<double> p; p.reserve(D);
    p.push_back(x / sigma_sp); p.push_back(y / sigma_sp); p.push_back(z / sigma_sp);
    if (sigma_t > 0.0) p.push_back(t / sigma_t);
    if (has_design) p.push_back(design[t] / sigma_d);
    int gi = 0;
    if (has_guide) { p.push_back(guide[idx3(x,y,z,nx,ny,nz)] / sr[gi++]); }
    for (size_t k=0;k<glist.size();++k) p.push_back(glist[k][idx3(x,y,z,nx,ny,nz)] / sr[gi++]);
    std::vector<int> g0(D); std::vector<double> frac(D);
    for (int j=0;j<D;++j) { double u = p[j]; int b = (int)std::floor(u); g0[j]=b; frac[j]=u - (double)b; }
    out[idx4(x,y,z,t,nx,ny,nz,nt)] = slice_point(grid, g0, frac, D);
  }
  out.attr("dim") = d;
  return out;
}
