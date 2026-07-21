#include "cocal_id_cpp_internal.hh"

namespace cocal_id_cpp {

namespace {
constexpr int bns_sp_trailing_metadata_values = 3;

string dims2(int n1, int n2) {
  return std::to_string(n1 + 1) + "x" + std::to_string(n2 + 1);
}

string dims3(int nr, int nt, int np) {
  return std::to_string(nr + 1) + "x" + std::to_string(nt + 1) + "x" + std::to_string(np + 1);
}

void require_dimensions(const string &path, const Array2D<double> &field, int n1, int n2) {
  if (field.empty() || field.lo1() != 0 || field.hi1() != n1 || field.lo2() != 0 || field.hi2() != n2) {
    abort_cocal("COCAL 2D field dimensions do not match " + path + ": file has " + dims2(n1, n2));
  }
}

void require_dimensions(const string &path, const Array3D<double> &field, int nr, int nt, int np) {
  if (field.empty() || field.lo1() != 0 || field.hi1() != nr || field.lo2() != 0 || field.hi2() != nt ||
      field.lo3() != 0 || field.hi3() != np) {
    abort_cocal("COCAL 3D field dimensions do not match " + path + ": file has " + dims3(nr, nt, np));
  }
}

void require_3d_fields(const string &path, const vector<Array3D<double> *> &fields, int nr, int nt, int np) {
  if (fields.empty())
    abort_cocal("Internal error: no destination fields for " + path);
  for (const Array3D<double> *field : fields) {
    if (field == nullptr)
      abort_cocal("Internal error: null destination field for " + path);
    require_dimensions(path, *field, nr, nt, np);
  }
}

double read_single_real(TextFile &f, const string &context) {
  return parse_real(f.line_tokens(1, context).at(0), f.path + " " + context);
}
} // namespace

Grid read_parameter_file(const string &path) {
  TextFile f(path);
  Grid g;
  auto t = f.line_tokens(4, "grid dimensions");
  g.nrg = parse_int(t.at(0));
  g.ntg = parse_int(t.at(1));
  g.npg = parse_int(t.at(2));
  g.nlg = parse_int(t.at(3));
  t = f.line_tokens(4, "fluid dimensions");
  g.nrf = parse_int(t.at(0));
  g.ntf = parse_int(t.at(1));
  g.npf = parse_int(t.at(2));
  g.nlf = parse_int(t.at(3));
  t = f.line_tokens(4, "deformation settings");
  g.nrf_deform = parse_int(t.at(0));
  g.nrgin = parse_int(t.at(1));
  g.NS_shape = t.at(2);
  g.EQ_point = t.at(3);
  t = f.line_tokens(3, "radial domain");
  g.rgin = parse_real(t.at(0));
  g.rgmid = parse_real(t.at(1));
  g.rgout = parse_real(t.at(2));
  t = f.line_tokens(3, "iteration switches");
  g.iter_max = parse_int(t.at(0));
  g.sw_mass_iter = t.at(1);
  g.sw_art_deform = t.at(2);
  t = f.line_tokens(2, "gravity convergence");
  g.conv0_gra = parse_real(t.at(0));
  g.conv_ini = parse_real(t.at(1));
  t = f.line_tokens(2, "fluid convergence");
  g.conv0_den = parse_real(t.at(0));
  g.conv0_vep = parse_real(t.at(1));
  t = f.line_tokens(3, "data type switches");
  g.indata_type = t.at(0);
  g.outdata_type = t.at(1);
  if (t.at(2).size() >= 3) {
    g.chrot = t.at(2).substr(0, 1);
    g.chgra = t.at(2).substr(1, 1);
    g.chope = t.at(2).substr(2, 1);
  } else if (t.size() >= 5) {
    g.chrot = t.at(2);
    g.chgra = t.at(3);
    g.chope = t.at(4);
  } else {
    abort_cocal("Malformed COCAL rotation/gravity/operator switches in " + path);
  }
  t = f.line_tokens(2, "mass convergence");
  g.eps_cocal = parse_real(t.at(0));
  g.mass_eps = parse_real(t.at(1));
  t = f.line_tokens(2, "sequence settings");
  g.num_sol_seq = parse_int(t.at(0));
  g.deform_par = parse_int(t.at(1));
  t = f.line_tokens(2, "central values");
  g.emdc = parse_real(t.at(0));
  g.pinx = parse_real(t.at(1));
  t = f.line_tokens(2, "spherical masses");
  g.restmass_sph = parse_real(t.at(0));
  g.gravmass_sph = parse_real(t.at(1));
  t = f.line_tokens(2, "spherical compactness");
  g.MoverR_sph = parse_real(t.at(0));
  g.schwarz_radi_sph = parse_real(t.at(1));

  g.ratio = static_cast<double>(g.nrf_deform) / static_cast<double>(g.nrf);
  g.ntgpolp = 0;
  g.ntgpolm = g.ntg;
  g.ntgeq = g.ntg / angular_half_turn_divisor;
  g.ntgxy = g.ntg / angular_half_turn_divisor;
  g.npgxzp = 0;
  g.npgxzm = g.npg / angular_half_turn_divisor;
  g.npgyzp = g.npg / angular_quarter_turn_divisor;
  g.npgyzm = 3 * (g.npg / angular_quarter_turn_divisor);
  g.ntfpolp = 0;
  g.ntfpolm = g.ntf;
  g.ntfeq = g.ntf / angular_half_turn_divisor;
  g.ntfxy = g.ntf / angular_half_turn_divisor;
  g.npfxzp = 0;
  g.npfxzm = g.npf / angular_half_turn_divisor;
  g.npfyzp = g.npf / angular_quarter_turn_divisor;
  g.npfyzm = 3 * (g.npf / angular_quarter_turn_divisor);
  return g;
}

void read_surf_parameter_file(Grid &g, const string &path, int impt) {
  if (impt > 2)
    return;
  TextFile f(path + "/rnspar_surf_mpt" + std::to_string(impt) + ".dat");
  auto t = f.line_tokens(5, "surface parameter header");
  g.nrg_1 = parse_int(t.at(0));
  g.sw_eos = t.at(1);
  g.sw_sepa = parse_int(t.at(2));
  g.sw_quant = parse_int(t.at(3));
  g.sw_spin = parse_int(t.at(4));
  g.r_surf = read_single_real(f, "surface radius");
  g.target_sepa = read_single_real(f, "target separation");
  g.target_qt = read_single_real(f, "target mass ratio");
  g.target_sx = read_single_real(f, "target spin x");
  g.target_sy = read_single_real(f, "target spin y");
  g.target_sz = read_single_real(f, "target spin z");
  g.omespx = read_single_real(f, "spin angular velocity x");
  g.omespy = read_single_real(f, "spin angular velocity y");
  g.omespz = read_single_real(f, "spin angular velocity z");

  if (g.indata_type == "3D") {
    TextFile rs(path + "/r_surf_mpt" + std::to_string(impt) + ".dat");
    g.r_surf = read_single_real(rs, "3D surface radius");
  }
}

BinaryExcision read_binary_excision_file(const string &path, int impt) {
  TextFile f(path + "/bin_ex_par_mpt" + std::to_string(impt) + ".dat");
  BinaryExcision ex;
  const auto t = f.line_tokens(2, "binary excision header");
  ex.ex_nrg = parse_int(t.at(0));
  ex.ex_ndis = parse_int(t.at(1));
  return ex;
}

void read_3d_fields(const string &path, const vector<Array3D<double> *> &fields) {
  TokenStream ts(path);
  const int nr = ts.next_int();
  const int nt = ts.next_int();
  const int np = ts.next_int();
  require_3d_fields(path, fields, nr, nt, np);
  for (int ip = 0; ip <= np; ++ip)
    for (int it = 0; it <= nt; ++it)
      for (int ir = 0; ir <= nr; ++ir)
        for (auto *field : fields)
          (*field)(ir, it, ip) = ts.next_real();
}

void read_2d_surface(const string &path, Array2D<double> &rs) {
  TokenStream ts(path);
  const int nt = ts.next_int();
  const int np = ts.next_int();
  require_dimensions(path, rs, nt, np);
  for (int ip = 0; ip <= np; ++ip)
    for (int it = 0; it <= nt; ++it)
      rs(it, ip) = ts.next_real();
}

void read_rns_fluid(const string &path, Array3D<double> &emd, Array2D<double> &rs,
                    Array3D<double> &omef, double &ome, double &ber, double &radi) {
  TokenStream ts(path);
  const int nr = ts.next_int();
  const int nt = ts.next_int();
  const int np = ts.next_int();
  require_dimensions(path, emd, nr, nt, np);
  require_dimensions(path, omef, nr, nt, np);
  require_dimensions(path, rs, nt, np);
  for (int ip = 0; ip <= np; ++ip)
    for (int it = 0; it <= nt; ++it)
      for (int ir = 0; ir <= nr; ++ir) {
        emd(ir, it, ip) = ts.next_real();
        rs(it, ip) = ts.next_real();
        omef(ir, it, ip) = ts.next_real();
      }
  ome = ts.next_real();
  ber = ts.next_real();
  radi = ts.next_real();
}

void read_bht_fluid(const string &path, Array3D<double> &emd, Array3D<double> &omeg,
                    double &ome, double &ber, double &radi) {
  TokenStream ts(path);
  const int nr = ts.next_int();
  const int nt = ts.next_int();
  const int np = ts.next_int();
  require_dimensions(path, emd, nr, nt, np);
  require_dimensions(path, omeg, nr, nt, np);
  for (int ip = 0; ip <= np; ++ip)
    for (int it = 0; it <= nt; ++it)
      for (int ir = 0; ir <= nr; ++ir) {
        emd(ir, it, ip) = ts.next_real();
        omeg(ir, it, ip) = ts.next_real();
      }
  ome = ts.next_real();
  ber = ts.next_real();
  radi = ts.next_real();
}

void read_bns_fluid_co(const string &path, Array3D<double> &emd, double &ome, double &ber, double &radi) {
  TextFile f(path);
  auto t = f.line_tokens(3, "BNS CO fluid dimensions");
  const int nr = parse_int(t.at(0));
  const int nt = parse_int(t.at(1));
  const int np = parse_int(t.at(2));
  require_dimensions(path, emd, nr, nt, np);
  for (int ip = 0; ip <= np; ++ip)
    for (int it = 0; it <= nt; ++it)
      for (int ir = 0; ir <= nr; ++ir) {
        t = f.line_tokens(1, "BNS CO fluid row");
        emd(ir, it, ip) = parse_real(t.at(0));
      }
  t = f.line_tokens(3, "BNS CO orbital parameters");
  ome = parse_real(t.at(0));
  ber = parse_real(t.at(1));
  radi = parse_real(t.at(2));
}

void read_bns_fluid_ir(const string &path, Array3D<double> &emd, Array3D<double> &vep,
                       double &ome, double &ber, double &radi) {
  TextFile f(path);
  auto t = f.line_tokens(3, "BNS IR fluid dimensions");
  const int nr = parse_int(t.at(0));
  const int nt = parse_int(t.at(1));
  const int np = parse_int(t.at(2));
  require_dimensions(path, emd, nr, nt, np);
  require_dimensions(path, vep, nr, nt, np);
  for (int ip = 0; ip <= np; ++ip)
    for (int it = 0; it <= nt; ++it)
      for (int ir = 0; ir <= nr; ++ir) {
        t = f.line_tokens(2, "BNS IR fluid row");
        emd(ir, it, ip) = parse_real(t.at(0));
        vep(ir, it, ip) = parse_real(t.at(1));
      }
  t = f.line_tokens(3, "BNS IR orbital parameters");
  ome = parse_real(t.at(0));
  ber = parse_real(t.at(1));
  radi = parse_real(t.at(2));
}

void read_bns_fluid_sp(const string &path, Array3D<double> &emd, Array3D<double> &vep,
                       Array3D<double> &wxspf, Array3D<double> &wyspf, Array3D<double> &wzspf,
                       double &ome, double &ber, double &radi, double &confpow) {
  TokenStream ts(path);
  const int nr = ts.next_int();
  const int nt = ts.next_int();
  const int np = ts.next_int();
  require_dimensions(path, emd, nr, nt, np);
  require_dimensions(path, vep, nr, nt, np);
  require_dimensions(path, wxspf, nr, nt, np);
  require_dimensions(path, wyspf, nr, nt, np);
  require_dimensions(path, wzspf, nr, nt, np);
  for (int ip = 0; ip <= np; ++ip)
    for (int it = 0; it <= nt; ++it)
      for (int ir = 0; ir <= nr; ++ir) {
        emd(ir, it, ip) = ts.next_real();
        vep(ir, it, ip) = ts.next_real();
        wxspf(ir, it, ip) = ts.next_real();
        wyspf(ir, it, ip) = ts.next_real();
        wzspf(ir, it, ip) = ts.next_real();
      }
  ome = ts.next_real();
  ber = ts.next_real();
  radi = ts.next_real();
  confpow = ts.next_real();
  for (int i = 0; i < bns_sp_trailing_metadata_values; ++i)
    (void)ts.next_real();
}

} // namespace cocal_id_cpp
