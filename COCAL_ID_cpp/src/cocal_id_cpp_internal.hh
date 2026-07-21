#ifndef COCAL_ID_CPP_INTERNAL_HH
#define COCAL_ID_CPP_INTERNAL_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cocal_id_cpp {

using std::array;
using std::string;
using std::vector;

constexpr double pi = 3.141592653589793;
constexpr double g_cgs = 6.67428e-08;
constexpr double c_cgs = 2.99792458e+10;
constexpr double solmas = 1.98892e+33;
constexpr int nnpeos = 10;
constexpr int nnteos = 10001;
constexpr int nmpt = 3;
constexpr int lagrange_order = 4;
constexpr int lagrange_left_offset = 2;
constexpr int lagrange_upper_margin = lagrange_order - 1;
constexpr int angular_half_turn_divisor = 2;
constexpr int angular_quarter_turn_divisor = 4;
constexpr double cocal_zero_tolerance = 1.0e-14;
constexpr double grid_match_tolerance = 1.0e-10;
constexpr double bns_outer_patch_fraction = 0.7;
constexpr int bht_fill_points = 4;
constexpr int bht_nvars = 22;

enum BhtVar {
  BHT_ALP,
  BHT_BETAX,
  BHT_BETAY,
  BHT_BETAZ,
  BHT_GXX,
  BHT_GXY,
  BHT_GXZ,
  BHT_GYY,
  BHT_GYZ,
  BHT_GZZ,
  BHT_KXX,
  BHT_KXY,
  BHT_KXZ,
  BHT_KYY,
  BHT_KYZ,
  BHT_KZZ,
  BHT_RHO,
  BHT_PRESS,
  BHT_EPS,
  BHT_VX,
  BHT_VY,
  BHT_VZ
};

using BhtPoint = array<double, bht_nvars>;

[[noreturn]] void abort_cocal(const string &msg);
vector<string> data_tokens(string line);
vector<string> next_tokens(std::istream &is);
double parse_real(string tok, const string &context = {});
int parse_int(const string &tok, const string &context = {});

struct TextFile {
  explicit TextFile(const string &path) : path(path), in(path) {
    if (!in)
      abort_cocal("Could not open COCAL input file: " + path);
  }

  vector<string> line_tokens(size_t min_tokens = 1, const string &context = "record") {
    auto tokens = next_tokens(in);
    if (tokens.empty())
      abort_cocal("Unexpected EOF while reading " + path);
    if (tokens.size() < min_tokens) {
      abort_cocal("Malformed COCAL " + context + " in " + path + ": expected at least " +
                  std::to_string(min_tokens) + " tokens, got " + std::to_string(tokens.size()));
    }
    return tokens;
  }

  string path;
  std::ifstream in;
};

class TokenStream {
public:
  explicit TokenStream(const string &path) : path(path), in(path) {
    if (!in)
      abort_cocal("Could not open COCAL input file: " + path);
  }

  string next() {
    while (token_pos >= tokens.size()) {
      string line;
      if (!std::getline(in, line))
        abort_cocal("Unexpected EOF while reading " + path);
      tokens = data_tokens(line);
      token_pos = 0;
    }
    return tokens.at(token_pos++);
  }

  int next_int() { return parse_int(next(), path); }
  double next_real() { return parse_real(next(), path); }

private:
  string path;
  std::ifstream in;
  vector<string> tokens;
  size_t token_pos = 0;
};

template <typename T> class Array1D {
public:
  Array1D() = default;
  Array1D(int lo, int hi) { resize(lo, hi); }

  void resize(int lo, int hi) {
    lo_ = lo;
    hi_ = hi;
    data_.assign(static_cast<size_t>(hi - lo + 1), T{});
  }

  void clear() {
    lo_ = 0;
    hi_ = -1;
    vector<T>().swap(data_);
  }

  bool empty() const { return data_.empty(); }
  int lo() const { return lo_; }
  int hi() const { return hi_; }

  T &operator()(int i) { return element(static_cast<size_t>(i - lo_)); }
  const T &operator()(int i) const { return element(static_cast<size_t>(i - lo_)); }

  void fill(const T &value) { std::fill(data_.begin(), data_.end(), value); }

private:
  T &element(size_t index) {
#ifdef COCAL_ID_CPP_BOUNDS_CHECK
    return data_.at(index);
#else
    return data_[index];
#endif
  }

  const T &element(size_t index) const {
#ifdef COCAL_ID_CPP_BOUNDS_CHECK
    return data_.at(index);
#else
    return data_[index];
#endif
  }

  int lo_ = 0, hi_ = -1;
  vector<T> data_;
};

template <typename T> class Array2D {
public:
  Array2D() = default;
  Array2D(int lo1, int hi1, int lo2, int hi2) { resize(lo1, hi1, lo2, hi2); }

  void resize(int lo1, int hi1, int lo2, int hi2) {
    lo1_ = lo1;
    hi1_ = hi1;
    lo2_ = lo2;
    hi2_ = hi2;
    n1_ = hi1 - lo1 + 1;
    n2_ = hi2 - lo2 + 1;
    data_.assign(static_cast<size_t>(n1_ * n2_), T{});
  }

  void clear() {
    lo1_ = lo2_ = 0;
    hi1_ = hi2_ = -1;
    n1_ = n2_ = 0;
    vector<T>().swap(data_);
  }

  bool empty() const { return data_.empty(); }
  int lo1() const { return lo1_; }
  int hi1() const { return hi1_; }
  int lo2() const { return lo2_; }
  int hi2() const { return hi2_; }

  T &operator()(int i, int j) { return element(static_cast<size_t>((i - lo1_) + n1_ * (j - lo2_))); }
  const T &operator()(int i, int j) const {
    return element(static_cast<size_t>((i - lo1_) + n1_ * (j - lo2_)));
  }

  void fill(const T &value) { std::fill(data_.begin(), data_.end(), value); }

private:
  T &element(size_t index) {
#ifdef COCAL_ID_CPP_BOUNDS_CHECK
    return data_.at(index);
#else
    return data_[index];
#endif
  }

  const T &element(size_t index) const {
#ifdef COCAL_ID_CPP_BOUNDS_CHECK
    return data_.at(index);
#else
    return data_[index];
#endif
  }

  int lo1_ = 0, hi1_ = -1, lo2_ = 0, hi2_ = -1, n1_ = 0, n2_ = 0;
  vector<T> data_;
};

template <typename T> class Array3D {
public:
  Array3D() = default;
  Array3D(int lo1, int hi1, int lo2, int hi2, int lo3, int hi3) {
    resize(lo1, hi1, lo2, hi2, lo3, hi3);
  }

  void resize(int lo1, int hi1, int lo2, int hi2, int lo3, int hi3) {
    lo1_ = lo1;
    hi1_ = hi1;
    lo2_ = lo2;
    hi2_ = hi2;
    lo3_ = lo3;
    hi3_ = hi3;
    n1_ = hi1 - lo1 + 1;
    n2_ = hi2 - lo2 + 1;
    n3_ = hi3 - lo3 + 1;
    data_.assign(static_cast<size_t>(n1_ * n2_ * n3_), T{});
  }

  void clear() {
    lo1_ = lo2_ = lo3_ = 0;
    hi1_ = hi2_ = hi3_ = -1;
    n1_ = n2_ = n3_ = 0;
    vector<T>().swap(data_);
  }

  bool empty() const { return data_.empty(); }
  int lo1() const { return lo1_; }
  int hi1() const { return hi1_; }
  int lo2() const { return lo2_; }
  int hi2() const { return hi2_; }
  int lo3() const { return lo3_; }
  int hi3() const { return hi3_; }

  T &operator()(int i, int j, int k) {
    return element(static_cast<size_t>((i - lo1_) + n1_ * ((j - lo2_) + n2_ * (k - lo3_))));
  }
  const T &operator()(int i, int j, int k) const {
    return element(static_cast<size_t>((i - lo1_) + n1_ * ((j - lo2_) + n2_ * (k - lo3_))));
  }

  void fill(const T &value) { std::fill(data_.begin(), data_.end(), value); }

private:
  T &element(size_t index) {
#ifdef COCAL_ID_CPP_BOUNDS_CHECK
    return data_.at(index);
#else
    return data_[index];
#endif
  }

  const T &element(size_t index) const {
#ifdef COCAL_ID_CPP_BOUNDS_CHECK
    return data_.at(index);
#else
    return data_[index];
#endif
  }

  int lo1_ = 0, hi1_ = -1, lo2_ = 0, hi2_ = -1, lo3_ = 0, hi3_ = -1;
  int n1_ = 0, n2_ = 0, n3_ = 0;
  vector<T> data_;
};

template <typename T> void resize_like_grid(Array3D<T> &a, int nrg, int ntg, int npg) {
  a.resize(0, nrg, 0, ntg, 0, npg);
}

template <typename T> void resize_like_fluid(Array3D<T> &a, int nrf, int ntf, int npf) {
  a.resize(0, nrf, 0, ntf, 0, npf);
}

struct Eos {
  vector<double> abc;
  vector<double> abi;
  vector<double> rhoi_peos;
  vector<double> qi_peos;
  vector<double> hi_peos;
  vector<double> rhocgs_peos;
  vector<double> abccgs;

  vector<double> rhoi_teos;
  vector<double> qi_teos;
  vector<double> hi_teos;
  vector<double> prei;
  vector<double> enei;
  vector<double> rhocgs_teos;
  vector<double> precgs;
  vector<double> enecgs;

  double rhoini_cgs = 0.0;
  double rhoini_gcm1 = 0.0;
  double emdini_gcm1 = 0.0;
  double sgma = 0.0;
  double constqc = 0.0;
  double cbar = 0.0;
  double kappa_crust = 0.0;
  double gamma_crust = 0.0;
  int nphase = 0;
  int teos_cached_iphase = 1;
};

struct Grid {
  int nrg = 0, ntg = 0, npg = 0, nlg = 0;
  int nrf = 0, ntf = 0, npf = 0, nlf = 0;
  int nrf_deform = 0, nrgin = 0;
  int ntgpolp = 0, ntgpolm = 0, ntgeq = 0, ntgxy = 0;
  int npgxzp = 0, npgxzm = 0, npgyzp = 0, npgyzm = 0;
  int ntfpolp = 0, ntfpolm = 0, ntfeq = 0, ntfxy = 0;
  int npfxzp = 0, npfxzm = 0, npfyzp = 0, npfyzm = 0;
  int iter_max = 0, num_sol_seq = 0, deform_par = 0;
  string indata_type, outdata_type, NS_shape, EQ_point;
  string chrot, chgra, chope, sw_mass_iter, sw_art_deform;
  double rgin = 0.0, rgmid = 0.0, rgout = 0.0, ratio = 0.0;
  double conv_gra = 0.0, conv_den = 0.0, conv_vep = 0.0, conv_ini = 0.0;
  double conv0_gra = 0.0, conv0_den = 0.0, conv0_vep = 0.0;
  double eps_cocal = 0.0, mass_eps = 0.0;
  double pinx = 0.0, emdc = 0.0;
  double restmass_sph = 0.0, gravmass_sph = 0.0, MoverR_sph = 0.0, schwarz_radi_sph = 0.0;

  int nrg_1 = 0;
  string sw_eos;
  int sw_sepa = 0, sw_quant = 0, sw_spin = 0;
  double r_surf = 0.0, target_sepa = 0.0, target_qt = 0.0;
  double target_sx = 0.0, target_sy = 0.0, target_sz = 0.0;
  double omespx = 0.0, omespy = 0.0, omespz = 0.0;

  Array1D<double> rg, rginv, hrg, hrginv, drg, drginv;
  Array1D<double> thg, hthg, phig, hphig;
  double dthg = 0.0, dthginv = 0.0, dphig = 0.0, dphiginv = 0.0;

  Array1D<double> rgex, thgex, phigex;
  Array1D<int> irgex_r, itgex_th, ipgex_phi;
  Array2D<int> itgex_r, ipgex_r, ipgex_th;
};

struct BinaryExcision {
  int ex_nrg = 0, ex_ndis = 0;
  double ex_radius = 0.0, ex_rgin = 0.0, ex_rgmid = 0.0, ex_rgout = 0.0;
  double sepa = 0.0, dis = 0.0;
};

struct RnsData {
  double ome = 0.0, ber = 0.0, radi = 0.0;
  Grid grid;
  Eos eos;

  Array3D<double> utf, uxf, uyf, uzf;
  Array3D<double> emd, omef, psif, alphf, bvxuf, bvyuf, bvzuf;
  Array2D<double> rs;
  Array3D<double> psi, alph, bvxd, bvyd, bvzd, bvxu, bvyu, bvzu;
  Array3D<double> hxxd, hxyd, hxzd, hyyd, hyzd, hzzd;
  Array3D<double> hxxu, hxyu, hxzu, hyyu, hyzu, hzzu;
  Array3D<double> kxxa, kxya, kxza, kyya, kyza, kzza;
  Array3D<double> va, vaxd, vayd, vazd;
  Array3D<double> fxd, fyd, fzd, fxyd, fxzd, fyzd;
  bool have_read_data = false;
};

struct Patch {
  Grid grid;
  BinaryExcision ex;
  Eos eos;

  double ome = 0.0, ber = 0.0, radi = 0.0, r_surf = 0.0, confpow = 0.0;

  Array3D<double> emd, vep, vepxf, vepyf, vepzf;
  Array3D<double> wxspf, wyspf, wzspf;
  Array3D<double> psif, alphf, bvxdf, bvydf, bvzdf;
  Array2D<double> rs;
  Array3D<double> psi, alph, bvxd, bvyd, bvzd;
  Array3D<double> axx, axy, axz, ayy, ayz, azz;
};

struct BnsData {
  string path;
  string id_type;
  array<Patch, 3> patches;
  Eos active_eos;
  double rr3 = 0.0, dis_cm = 0.0;
  bool have_read_data = false;
};

struct BhtParameters {
  double r_surf = 0.0, mass_bh = 0.0;
  double ome_bh = 0.0, spin_bh = 0.0;
  double alph_bh = 0.0, psi_bh = 0.0;
  double th_spin_bh_deg = 0.0, phi_spin_bh_deg = 0.0;
  double mass_ratio = 0.0;
  string boundary_type, solution_type;
};

struct BhtData {
  Grid grid;
  Eos eos;
  BhtParameters parameters;
  double ome = 0.0, ber = 0.0, radi = 0.0;
  double rexc = 0.0, r_ah = 0.0, r_eh = 0.0, r_ring = 0.0;
  double center_x = 0.0, center_y = 0.0, center_z = 0.0;
  double extra_points_dr = 0.0, trace_width = 0.0, trace_power = 0.0;
  BhtPoint idpar_center{};
  BhtPoint smooth_center{};
  bool tabulated_eos = false;
  bool smooth_center_attempted = false;
  bool smooth_center_valid = false;
  bool have_read_data = false;

  Array3D<double> emd, omeg;
  Array3D<double> psi, alph, bvxd, bvyd, bvzd, bvxu, bvyu, bvzu;
  Array3D<double> hxxd, hxyd, hxzd, hyyd, hyzd, hzzd;
  Array3D<double> hxxu, hxyu, hxzu, hyyu, hyzu, hzzu;
  Array3D<double> kxxa, kxya, kxza, kyya, kyza, kzza;
};

struct MetricValues {
  double psica = 0.0, alphca = 0.0, bvxdca = 0.0, bvydca = 0.0, bvzdca = 0.0;
  double axx = 0.0, axy = 0.0, axz = 0.0, ayy = 0.0, ayz = 0.0, azz = 0.0;
};

struct FluidValues {
  double emdca = 0.0, psifca = 0.0, alphfca = 0.0, bvxdfca = 0.0, bvydfca = 0.0, bvzdfca = 0.0;
  double vepxfca = 0.0, vepyfca = 0.0, vepzfca = 0.0;
  double wxspfca = 0.0, wyspfca = 0.0, wzspfca = 0.0;
  bool inside = false;
};

struct BnsFluidPoint {
  double emdca = 0.0;
  double vxu = 0.0, vyu = 0.0, vzu = 0.0;
};

struct RnsPoint {
  double emdca = 0.0;
  double alphca = 0.0, bvxdca = 0.0, bvydca = 0.0, bvzdca = 0.0;
  double vxu = 0.0, vyu = 0.0, vzu = 0.0;
  double gxx1 = 0.0, gxy1 = 0.0, gxz1 = 0.0, gyy1 = 0.0, gyz1 = 0.0, gzz1 = 0.0;
  double kxx1 = 0.0, kxy1 = 0.0, kxz1 = 0.0, kyy1 = 0.0, kyz1 = 0.0, kzz1 = 0.0;
  double vaca = 0.0, vaxdca = 0.0, vaydca = 0.0, vazdca = 0.0;
  double fxydca = 0.0, fxzdca = 0.0, fyzdca = 0.0;
};

extern RnsData rns;
extern BnsData bns;
extern BhtData bht;
extern std::mutex rns_mutex;
extern std::mutex bns_mutex;
extern std::mutex bht_mutex;

Grid read_parameter_file(const string &path);
void read_surf_parameter_file(Grid &g, const string &path, int impt);
BinaryExcision read_binary_excision_file(const string &path, int impt);

void grid_r(Grid &gr);
void grid_r_bqs(Grid &gr);
void grid_r_bht(Grid &gr, double mass_bh, double spin_bh);
void grid_theta(Grid &gr);
void grid_phi(Grid &gr);
void grid_extended(Grid &gr);
void setup_coordinates(Grid &gr, bool bns_grid);
void calc_parameter_binary_excision(Patch &p);

void lagint_4th_weights(const array<double, lagrange_order> &x, double v, array<double, lagrange_order> &w);
double lagint_4th_apply(const array<double, lagrange_order> &w, const array<double, lagrange_order> &y);
double lagint_4th(const array<double, lagrange_order> &x, const array<double, lagrange_order> &y, double v);

struct Stencil4 {
  int ir[lagrange_order][lagrange_order][lagrange_order]{};
  int it[lagrange_order][lagrange_order][lagrange_order]{};
  int ip[lagrange_order][lagrange_order][lagrange_order]{};
  array<double, lagrange_order> wr{}, wth{}, wphi{};
};

void spherical_from_cart(double x, double y, double z, double &r, double &th, double &ph);
Stencil4 cgr_4th_setup_from_spherical(const Grid &gr, double rc, double thc, double phic, int nr);
Stencil4 gr2cgr_4th_setup(const Grid &gr, double x, double y, double z);
double cgr_4th_apply(const Array3D<double> &fnc, const Stencil4 &st);
double interpo_lag4th_2Dsurf(const Grid &gr, const Array2D<double> &fnc, double tv, double pv);
bool fl2cgr_4th_setup(const Grid &gr, const Array2D<double> &rs, double x, double y, double z, Stencil4 &st);

void read_3d_fields(const string &path, const vector<Array3D<double> *> &fields);
void read_2d_surface(const string &path, Array2D<double> &rs);
void read_rns_fluid(const string &path, Array3D<double> &emd, Array2D<double> &rs,
                    Array3D<double> &omef, double &ome, double &ber, double &radi);
void read_bht_fluid(const string &path, Array3D<double> &emd, Array3D<double> &omeg,
                    double &ome, double &ber, double &radi);
void read_bns_fluid_co(const string &path, Array3D<double> &emd, double &ome, double &ber, double &radi);
void read_bns_fluid_ir(const string &path, Array3D<double> &emd, Array3D<double> &vep,
                       double &ome, double &ber, double &radi);
void read_bns_fluid_sp(const string &path, Array3D<double> &emd, Array3D<double> &vep,
                       Array3D<double> &wxspf, Array3D<double> &wyspf, Array3D<double> &wzspf,
                       double &ome, double &ber, double &radi, double &confpow);

int peos_lookup(const Eos &eos, double qp, const vector<double> &qpar);
void peos_initialize(Eos &eos, const string &path);
void peos_q2hprho(Eos &eos, double &q, double &h, double &pre, double &rho, double &ened);
int teos_lookup(const Eos &eos, double qp, const vector<double> &qpar);
void teos_initialize(Eos &eos, const string &teos_path, const string &peos_path);
void teos_q2hprho(Eos &eos, double &q, double &h, double &pre, double &rho, double &ened);
void q2hprho(Eos &eos, bool tabulated, double &q, double &h, double &pre, double &rho, double &ened);

void interpo_gr2fl_export(const Grid &gr, const Array3D<double> &grv, Array3D<double> &flv,
                          const Array2D<double> &rs);
void interpo_gr2fl_metric_CF_export(const Grid &gr, const Array3D<double> &alph,
                                    const Array3D<double> &psi, const Array3D<double> &bvxd,
                                    const Array3D<double> &bvyd, const Array3D<double> &bvzd,
                                    Array3D<double> &alphf, Array3D<double> &psif,
                                    Array3D<double> &bvxdf, Array3D<double> &bvydf,
                                    Array3D<double> &bvzdf, const Array2D<double> &rs);
void invhij_WL_export(const Grid &gr, const Array3D<double> &hxxd, const Array3D<double> &hxyd,
                      const Array3D<double> &hxzd, const Array3D<double> &hyyd,
                      const Array3D<double> &hyzd, const Array3D<double> &hzzd,
                      Array3D<double> &hxxu, Array3D<double> &hxyu, Array3D<double> &hxzu,
                      Array3D<double> &hyyu, Array3D<double> &hyzu, Array3D<double> &hzzu);
void index_vec_down2up_export(const Grid &gr, const Array3D<double> &hxxu, const Array3D<double> &hxyu,
                              const Array3D<double> &hxzu, const Array3D<double> &hyyu,
                              const Array3D<double> &hyzu, const Array3D<double> &hzzu,
                              Array3D<double> &vecxu, Array3D<double> &vecyu, Array3D<double> &veczu,
                              const Array3D<double> &vecxd, const Array3D<double> &vecyd,
                              const Array3D<double> &veczd);

void allocate_rns_arrays(RnsData &d, bool mrns, bool read_potential);
void allocate_bht_arrays(BhtData &d);
void allocate_patch_arrays(Patch &p, const string &id_type, bool star_patch);
void sign_flip(Array3D<double> &a);
RnsPoint interpolate_rns_point(const RnsData &d, const string &rnstype, bool read_magnetic_potential,
                               double xcoc, double ycoc, double zcoc);
MetricValues interpolate_metric_patch(const Patch &p, double x, double y, double z);
FluidValues interpolate_fluid_patch(const Patch &p, const string &id_type, double x, double y, double z);
BnsFluidPoint interpolate_bns_fluid_patch(const Patch &p, const string &id_type, double x, double y, double z,
                                          double xcoc, double ycoc, double zcoc, double ecc_cor_velx);

} // namespace cocal_id_cpp

#endif // COCAL_ID_CPP_INTERNAL_HH
