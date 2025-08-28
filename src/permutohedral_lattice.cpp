
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstring>

using namespace Rcpp;

struct KeyHash {
  std::size_t operator()(const std::vector<int> &k) const noexcept {
    std::size_t h = 1469598103934665603ull;
    for (int v : k) {
      std::size_t x = static_cast<std::size_t>(v * 1315423911u);
      h ^= x; h *= 1099511628211ull;
    }
    return h;
  }
};
struct KeyEq {
  bool operator()(const std::vector<int>& a, const std::vector<int>& b) const noexcept {
    if (a.size() != b.size()) return false;
    for (size_t i=0;i<a.size();++i) if (a[i] != b[i]) return false;
    return true;
  }
};

struct Node { double val; double mass; Node():val(0.0),mass(0.0) {} };

static inline void barycentric_from_sorted(const std::vector<double> &r,
                                           const std::vector<int> &ord,
                                           std::vector<double> &w) {
  const int D = (int)r.size();
  w.assign(D+1, 0.0);
  w[0] = 1.0 - r[ord[0]];
  for (int k=1;k<D;++k) w[k] = r[ord[k-1]] - r[ord[k]];
  w[D] = r[ord[D-1]];
}

static inline void build_simplex_keys(const std::vector<int> &q,
                                      const std::vector<int> &ord,
                                      std::vector< std::vector<int> > &keys) {
  const int D = (int)q.size();
  keys.assign(D+1, std::vector<int>(D+1, 0));
  for (int i=0;i<D;++i) keys[0][i] = q[i];
  for (int k=1;k<=D;++k) {
    keys[k] = keys[k-1];
    keys[k][ ord[k-1] ] += 1;
  }
  for (int k=0;k<=D;++k) {
    int s = 0; for (int i=0;i<D;++i) s += keys[k][i];
    keys[k][D] = -s;
  }
}

static inline void neighbor_inc(int r, int Dp1, std::vector<int> &inc) {
  inc.assign(Dp1, 0);
  inc[r] = 1;
  inc[(r+1) % Dp1] = -1;
}

static void blur_axis(std::vector< std::vector<int> > &keys,
                      std::vector<Node> &nodes,
                      std::unordered_map< std::vector<int>, int, KeyHash, KeyEq > &lookup,
                      int r) {
  const int Dp1 = (int)keys[0].size();
  std::vector<int> inc; neighbor_inc(r, Dp1, inc);
  std::vector<int> dec(Dp1,0);
  for (int i=0;i<Dp1;++i) dec[i] = -inc[i];

  std::vector<Node> new_nodes(nodes.size());
  for (size_t idx=0; idx<nodes.size(); ++idx) {
    const std::vector<int> &k = keys[idx];
    std::vector<int> kp = k, km = k;
    for (int i=0;i<Dp1;++i) { kp[i] += inc[i]; km[i] += dec[i]; }
    double vp = 0.0, mp = 0.0, vm = 0.0, mm = 0.0;
    auto itp = lookup.find(kp);
    if (itp != lookup.end()) { const Node &np = nodes[itp->second]; vp = np.val; mp = np.mass; }
    auto itm = lookup.find(km);
    if (itm != lookup.end()) { const Node &nm = nodes[itm->second]; vm = nm.val; mm = nm.mass; }
    new_nodes[idx].val  = 0.5*nodes[idx].val + 0.25*(vp + vm);
    new_nodes[idx].mass = 0.5*nodes[idx].mass + 0.25*(mp + mm);
  }
  nodes.swap(new_nodes);
}

static inline void build_features_3d(int x,int y,int z,
                                     double sigma_sp,
                                     const double *guide,
                                     int nx,int ny,int nz,
                                     const std::vector<const double*> &guides,
                                     const std::vector<double> &sigma_r,
                                     std::vector<double> &f) {
  f.clear();
  f.push_back(x / sigma_sp);
  f.push_back(y / sigma_sp);
  f.push_back(z / sigma_sp);
  int gi = 0;
  if (guide) {
    double g = guide[x + nx*(y + (size_t)ny*z)];
    double sr = (sigma_r.size()>gi? sigma_r[gi]: 1.0); ++gi;
    f.push_back(g / sr);
  }
  for (size_t j=0;j<guides.size();++j) {
    double g = guides[j][ x + nx*(y + (size_t)ny*z) ];
    double sr = (sigma_r.size()>gi? sigma_r[gi]: 1.0); ++gi;
    f.push_back(g / sr);
  }
}

static inline void build_features_4d(int x,int y,int z,int t,
                                     double sigma_sp,double sigma_t,
                                     const double *guide,
                                     int nx,int ny,int nz,int nt,
                                     const std::vector<const double*> &guides,
                                     const std::vector<double> &sigma_r,
                                     const std::vector<const double*> &designs,
                                     const std::vector<double> &sigma_dv,
                                     std::vector<double> &f) {
  f.clear();
  f.push_back(x / sigma_sp);
  f.push_back(y / sigma_sp);
  f.push_back(z / sigma_sp);
  if (sigma_t > 0.0) f.push_back(t / sigma_t);
  if (!designs.empty()) {
    const int K = (int)designs.size();
    for (int k=0;k<K;++k) {
      double sd = (sigma_dv.size()>k? sigma_dv[k]: 1.0);
      f.push_back( designs[k][t] / sd );
    }
  }
  int gi = 0;
  if (guide) {
    double g = guide[x + nx*(y + (size_t)ny*z)];
    double sr = (sigma_r.size()>gi? sigma_r[gi]: 1.0); ++gi;
    f.push_back(g / sr);
  }
  for (size_t j=0;j<guides.size();++j) {
    double g = guides[j][ x + nx*(y + (size_t)ny*z) ];
    double sr = (sigma_r.size()>gi? sigma_r[gi]: 1.0); ++gi;
    f.push_back(g / sr);
  }
}

// [[Rcpp::export(name = "permutohedral_lattice3d_cpp")]]
SEXP permutohedral_lattice3d_cpp(SEXP vol3d_, SEXP mask3_,
                                 SEXP sigma_sp_, SEXP sigma_r_,
                                 SEXP guide_, SEXP guides_,
                                 SEXP passes_) {
  NumericVector vol3d(vol3d_);
  LogicalVector mask3(mask3_);
  double sigma_sp = as<double>(sigma_sp_);
  NumericVector sigma_r_nv(sigma_r_);
  int passes = as<int>(passes_);

  NumericVector guide_nv;
  bool has_guide = !Rf_isNull(guide_);
  if (has_guide) guide_nv = as<NumericVector>(guide_);

  std::vector<const double*> guides;
  if (!Rf_isNull(guides_)) {
    List gl(guides_);
    for (int i=0;i<gl.size();++i) {
      NumericVector g = as<NumericVector>(gl[i]);
      guides.push_back(g.begin());
    }
  }

  IntegerVector d = vol3d.attr("dim");
  if (d.size()!=3) stop("vol3d must be a 3D array");
  const int nx=d[0], ny=d[1], nz=d[2];
  IntegerVector dm = mask3.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");

  std::vector<double> sigma_r(sigma_r_nv.begin(), sigma_r_nv.end());
  std::vector<double> f; f.reserve(3 + (int)1 + (int)guides.size());

  const size_t N = (size_t)nx*ny*nz;
  std::vector< std::vector<int> > keys; keys.reserve(N*4);
  std::unordered_map< std::vector<int>, int, KeyHash, KeyEq > lookup; lookup.reserve(N*2);
  std::vector<Node> nodes; nodes.reserve(N*2);

  std::vector< std::array<int,8> > sample_idx; sample_idx.reserve(N);
  std::vector< std::array<double,8> > sample_w; sample_w.reserve(N);

  int Dmax = 0;
  size_t sample_counter = 0;
  for (int z=0; z<nz; ++z) for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
    if (!mask3[x + nx*(y + (size_t)ny*z)]) continue;
    build_features_3d(x,y,z, sigma_sp,
                      has_guide ? guide_nv.begin():nullptr,
                      nx,ny,nz,
                      guides, sigma_r, f);
    const int D = (int)f.size(); if (D>7) stop("Too many guide dims (supports up to 7 extras)");
    Dmax = std::max(Dmax, D);
    std::vector<int> q(D);
    std::vector<double> r(D);
    std::vector<int> ord(D);
    for (int i=0;i<D;++i) { double v=f[i]; int fl = (int)std::floor(v); q[i]=fl; r[i]=v-fl; ord[i]=i; }
    std::sort(ord.begin(), ord.end(), [&](int a,int b){return r[a]>r[b];});
    std::vector<double> w(D+1); barycentric_from_sorted(r, ord, w);
    std::vector< std::vector<int> > vkeys; build_simplex_keys(q, ord, vkeys);

    std::array<int,8> idxs; std::array<double,8> ws;
    for (int k=0;k<=D;++k) {
      auto it = lookup.find(vkeys[k]);
      int idx;
      if (it==lookup.end()) {
        idx = (int)nodes.size();
        lookup.emplace(vkeys[k], idx);
        keys.push_back(vkeys[k]);
        nodes.emplace_back();
      } else idx = it->second;
      nodes[idx].val  += w[k] * vol3d[x + (size_t)nx*(y + (size_t)ny*z)];
      nodes[idx].mass += w[k];
      idxs[k] = idx; ws[k] = w[k];
    }
    for (int k=D+1;k<8;++k) { idxs[k]=-1; ws[k]=0.0; }
    sample_idx.push_back(idxs);
    sample_w.push_back(ws);
    ++sample_counter;
  }

  for (int p=0;p<passes;++p) {
    for (int raxis=0; raxis<=Dmax; ++raxis) blur_axis(keys, nodes, lookup, raxis);
  }

  NumericVector out(N);
  std::fill(out.begin(), out.end(), 0.0);
  size_t n_ptr = 0;
  for (int z=0; z<nz; ++z) for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
    if (!mask3[x + nx*(y + (size_t)ny*z)]) { out[x + nx*(y + (size_t)ny*z)] = vol3d[x + nx*(y + (size_t)ny*z)]; continue; }
    const auto &idxs = sample_idx[n_ptr];
    const auto &ws   = sample_w[n_ptr];
    double num = 0.0, den = 0.0;
    for (int k=0;k<8;++k) {
      if (idxs[k] < 0) break;
      const Node &nd = nodes[idxs[k]];
      num += ws[k] * nd.val;
      den += ws[k] * nd.mass;
    }
    double v = (den > 1e-12) ? (num / den) : vol3d[x + nx*(y + (size_t)ny*z)];
    out[x + nx*(y + (size_t)ny*z)] = v;
    ++n_ptr;
  }

  out.attr("dim") = d;
  return out;
}

// [[Rcpp::export(name = "permutohedral_lattice4d_cpp")]]
SEXP permutohedral_lattice4d_cpp(SEXP vec4d_, SEXP mask3_,
                                 SEXP sigma_sp_, SEXP sigma_t_, SEXP sigma_r_,
                                 SEXP guide_spatial_, SEXP guides_,
                                 SEXP design_, SEXP sigma_d_,
                                 SEXP passes_) {
  NumericVector vec4d(vec4d_);
  LogicalVector mask3(mask3_);
  double sigma_sp = as<double>(sigma_sp_);
  double sigma_t  = as<double>(sigma_t_);
  NumericVector sigma_r_nv(sigma_r_);
  NumericVector guide_nv;
  bool has_guide = !Rf_isNull(guide_spatial_);
  if (has_guide) guide_nv = as<NumericVector>(guide_spatial_);
  NumericVector sigma_d_nv(sigma_d_);
  // design can be: NULL, numeric vector (length T), numeric matrix (nrow=T), or list of numeric vectors (each length T)
  std::vector<const double*> designs;
  std::vector<double> sigma_dv;
  int passes = as<int>(passes_);

  std::vector<const double*> guides;
  if (!Rf_isNull(guides_)) {
    List gl(guides_);
    for (int i=0;i<gl.size();++i) {
      NumericVector g = as<NumericVector>(gl[i]);
      guides.push_back(g.begin());
    }
  }

  IntegerVector d = vec4d.attr("dim");
  if (d.size()!=4) stop("vec4d must be a 4D array");
  const int nx=d[0], ny=d[1], nz=d[2], nt=d[3];
  IntegerVector dm = mask3.attr("dim");
  if (dm.size()!=3 || dm[0]!=nx || dm[1]!=ny || dm[2]!=nz) stop("mask dims mismatch");

  std::vector<double> sigma_r(sigma_r_nv.begin(), sigma_r_nv.end());
  const double *guide = has_guide ? guide_nv.begin() : nullptr;

  // Parse designs
  if (!Rf_isNull(design_)) {
    if (Rf_isMatrix(design_)) {
      NumericMatrix M(design_);
      if (M.nrow() != nt) stop("design matrix must have nrow = T (frames)");
      for (int j=0;j<M.ncol(); ++j) designs.push_back(&M(0,j));
    } else if (TYPEOF(design_) == VECSXP) {
      List dl(design_);
      for (int j=0;j<dl.size(); ++j) {
        NumericVector v = as<NumericVector>(dl[j]);
        if (v.size() != nt) stop("each design vector must have length T (frames)");
        designs.push_back(v.begin());
      }
    } else {
      NumericVector v(design_);
      if (v.size() != nt) stop("design must have length T (frames)");
      designs.push_back(v.begin());
    }
  }
  if (!designs.empty()) {
    if (sigma_d_nv.size() == 0) {
      sigma_dv.assign(designs.size(), 1.0);
    } else if (sigma_d_nv.size() == 1) {
      sigma_dv.assign(designs.size(), sigma_d_nv[0]);
    } else if ((size_t)sigma_d_nv.size() == designs.size()) {
      sigma_dv.assign(sigma_d_nv.begin(), sigma_d_nv.end());
    } else {
      stop("sigma_d must be length 1 or match number of design columns");
    }
  }

  std::vector<size_t> lin_index; lin_index.reserve((size_t)nx*ny*nz*nt);
  for (int t=0;t<nt;++t)
    for (int z=0; z<nz; ++z)
      for (int y=0; y<ny; ++y)
        for (int x=0; x<nx; ++x)
          if (mask3[x + (size_t)nx*(y + (size_t)ny*z)])
            lin_index.push_back(x + (size_t)nx*(y + (size_t)ny*z + (size_t)nz*t));

  const size_t N = lin_index.size();
  std::vector< std::vector<int> > keys; keys.reserve(N*6);
  std::unordered_map< std::vector<int>, int, KeyHash, KeyEq > lookup; lookup.reserve(N*2);
  std::vector<Node> nodes; nodes.reserve(N*2);
  std::vector< std::array<int,10> > sample_idx; sample_idx.reserve(N);
  std::vector< std::array<double,10> > sample_w; sample_w.reserve(N);
  int Dmax = 0;

  for (size_t idx=0; idx<N; ++idx) {
    size_t lin = lin_index[idx];
    int t = (int)(lin / ( (size_t)nx*ny*nz )); size_t rem = lin - (size_t)t*nx*ny*nz;
    int z = (int)(rem / ( (size_t)nx*ny )); rem -= (size_t)z*nx*ny;
    int y = (int)(rem / (size_t)nx); int x = (int)(rem - (size_t)y*nx);

    std::vector<double> f; f.reserve(6 + (int)guides.size());
    // (note: nt is unused here, kept for signature symmetry)
    std::vector<const double*> gtmp = guides;
    build_features_4d(x,y,z,t, sigma_sp, sigma_t,
                      guide, nx,ny,nz, /*nt*/0,
                      gtmp, sigma_r, designs, sigma_dv, f);
    const int D = (int)f.size(); if (D>9) stop("Too many features (supports up to 9)");
    Dmax = std::max(Dmax, D);

    std::vector<int> q(D);
    std::vector<double> r(D);
    std::vector<int> ord(D);
    for (int i=0;i<D;++i) { double v=f[i]; int fl = (int)std::floor(v); q[i]=fl; r[i]=v-fl; ord[i]=i; }
    std::sort(ord.begin(), ord.end(), [&](int a,int b){return r[a]>r[b];});
    std::vector<double> w(D+1); barycentric_from_sorted(r, ord, w);
    std::vector< std::vector<int> > vkeys; build_simplex_keys(q, ord, vkeys);

    std::array<int,10> idxs; std::array<double,10> ws;
    for (int k=0;k<=D;++k) {
      auto it = lookup.find(vkeys[k]);
      int li;
      if (it==lookup.end()) {
        li = (int)nodes.size();
        lookup.emplace(vkeys[k], li);
        keys.push_back(vkeys[k]);
        nodes.emplace_back();
      } else li = it->second;
      nodes[li].val  += w[k] * vec4d[lin];
      nodes[li].mass += w[k];
      idxs[k] = li; ws[k] = w[k];
    }
    for (int k=D+1;k<10;++k) { idxs[k] = -1; ws[k]=0.0; }
    sample_idx.push_back(idxs);
    sample_w.push_back(ws);
  }

  for (int p=0;p<passes;++p) {
    for (int raxis=0; raxis<=Dmax; ++raxis) blur_axis(keys, nodes, lookup, raxis);
  }

  NumericVector out(vec4d.size());
  std::copy(vec4d.begin(), vec4d.end(), out.begin());
  for (size_t si=0; si<N; ++si) {
    size_t lin = lin_index[si];
    const auto &idxs = sample_idx[si];
    const auto &ws   = sample_w[si];
    double num = 0.0, den = 0.0;
    for (int k=0;k<10;++k) {
      if (idxs[k] < 0) break;
      const Node &nd = nodes[idxs[k]];
      num += ws[k] * nd.val;
      den += ws[k] * nd.mass;
    }
    if (den > 1e-12) out[lin] = num / den;
  }

  // d already defined earlier, just reuse it
  out.attr("dim") = d;
  return out;
}
