#include "cocal_id_cpp_internal.hh"

namespace cocal_id_cpp {

namespace {
constexpr double eos_q_floor = 1.0e-60;

void require_peos_phase_count(int nphase, const string &path) {
  if (nphase < 1 || nphase > nnpeos)
    abort_cocal("Unsupported PEOS phase count in " + path + ": " + std::to_string(nphase));
}

void require_teos_table_count(int nphase, const string &path) {
  if (nphase < lagrange_upper_margin || nphase > nnteos)
    abort_cocal("Unsupported TEOS table size in " + path + ": " + std::to_string(nphase + 1));
}

void resize_peos_tables(Eos &eos) {
  const size_t n = static_cast<size_t>(nnpeos + 1);
  eos.abc.assign(n, 0.0);
  eos.abi.assign(n, 0.0);
  eos.rhoi_peos.assign(n, 0.0);
  eos.qi_peos.assign(n, 0.0);
  eos.hi_peos.assign(n, 0.0);
  eos.rhocgs_peos.assign(n, 0.0);
  eos.abccgs.assign(n, 0.0);
}

void resize_teos_tables(Eos &eos) {
  const size_t n = static_cast<size_t>(nnteos + 1);
  eos.rhoi_teos.assign(n, 0.0);
  eos.qi_teos.assign(n, 0.0);
  eos.hi_teos.assign(n, 0.0);
  eos.prei.assign(n, 0.0);
  eos.enei.assign(n, 0.0);
  eos.rhocgs_teos.assign(n, 0.0);
  eos.precgs.assign(n, 0.0);
  eos.enecgs.assign(n, 0.0);
}
} // namespace

int peos_lookup(const Eos &eos, double qp, const vector<double> &qpar) {
  int iphase = 1;
  for (int ii = 1; ii <= eos.nphase; ++ii) {
    const double det = (qp - qpar.at(ii)) * (qp - qpar.at(ii - 1));
    if (det <= 0.0) {
      iphase = ii;
      break;
    }
  }
  return iphase;
}

void peos_initialize(Eos &eos, const string &path) {
  eos = Eos{};
  resize_peos_tables(eos);

  TextFile f(path);
  auto t = f.line_tokens(2, "PEOS header");
  eos.nphase = parse_int(t.at(0));
  require_peos_phase_count(eos.nphase, path);
  eos.rhoini_cgs = parse_real(t.at(1));
  t = f.line_tokens(2, "PEOS reference state");
  const double rho_0 = parse_real(t.at(0));
  const double pre_0 = parse_real(t.at(1));
  for (int ii = eos.nphase; ii >= 0; --ii) {
    t = f.line_tokens(2, "PEOS phase row");
    eos.rhocgs_peos.at(ii) = parse_real(t.at(0));
    eos.abi.at(ii) = parse_real(t.at(1));
  }

  const double facrho = std::pow(g_cgs / (c_cgs * c_cgs), 3) * solmas * solmas;
  const double facpre = std::pow(g_cgs, 3) * solmas * solmas / std::pow(c_cgs, 8);
  for (int ii = 0; ii <= eos.nphase; ++ii)
    eos.rhoi_peos.at(ii) = facrho * eos.rhocgs_peos.at(ii);

  const int iphase = peos_lookup(eos, rho_0, eos.rhocgs_peos);
  eos.abc.at(iphase) = facpre / std::pow(facrho, eos.abi.at(iphase)) *
                       (pre_0 / std::pow(rho_0, eos.abi.at(iphase)));
  eos.abccgs.at(iphase) = pre_0 / std::pow(rho_0, eos.abi.at(iphase));

  for (int ii = iphase - 1; ii >= 0; --ii) {
    eos.abc.at(ii) =
        std::pow(eos.rhoi_peos.at(ii), eos.abi.at(ii + 1) - eos.abi.at(ii)) * eos.abc.at(ii + 1);
    eos.abccgs.at(ii) =
        std::pow(eos.rhocgs_peos.at(ii), eos.abi.at(ii + 1) - eos.abi.at(ii)) * eos.abccgs.at(ii + 1);
  }
  for (int ii = iphase + 1; ii <= eos.nphase; ++ii) {
    eos.abc.at(ii) =
        std::pow(eos.rhoi_peos.at(ii - 1), eos.abi.at(ii - 1) - eos.abi.at(ii)) * eos.abc.at(ii - 1);
    eos.abccgs.at(ii) =
        std::pow(eos.rhocgs_peos.at(ii - 1), eos.abi.at(ii - 1) - eos.abi.at(ii)) *
        eos.abccgs.at(ii - 1);
  }
  for (int ii = 0; ii <= eos.nphase; ++ii)
    eos.qi_peos.at(ii) = eos.abc.at(ii) * std::pow(eos.rhoi_peos.at(ii), eos.abi.at(ii) - 1.0);

  eos.hi_peos.at(0) = 1.0;
  for (int ii = 1; ii <= eos.nphase; ++ii) {
    const double fac2 = eos.abi.at(ii) / (eos.abi.at(ii) - 1.0);
    eos.hi_peos.at(ii) = eos.hi_peos.at(ii - 1) + fac2 * (eos.qi_peos.at(ii) - eos.qi_peos.at(ii - 1));
  }

  eos.rhoini_gcm1 = facrho * eos.rhoini_cgs;
  const int iphase_ini = peos_lookup(eos, eos.rhoini_gcm1, eos.rhoi_peos);
  eos.emdini_gcm1 =
      eos.abc.at(iphase_ini) * std::pow(eos.rhoini_gcm1, eos.abi.at(iphase_ini) - 1.0);
}

void peos_q2hprho(Eos &eos, double &q, double &h, double &pre, double &rho, double &ened) {
  const int iphase = peos_lookup(eos, q, eos.qi_peos);
  const double hin = eos.hi_peos.at(iphase);
  const double qin = eos.qi_peos.at(iphase);
  const double abin = eos.abi.at(iphase);
  const double abct = eos.abc.at(iphase);
  const double fac1 = 1.0 / (abin - 1.0);
  const double fac2 = abin / (abin - 1.0);
  const double fack = std::pow(abct, -fac1);
  if (q <= eos_q_floor)
    q = eos_q_floor;
  h = hin + fac2 * (q - qin);
  if (h <= 1.0)
    h = 1.0;
  pre = fack * std::pow(q, fac2);
  rho = fack * std::pow(q, fac1);
  ened = rho * h - pre;
}

int teos_lookup(const Eos &eos, double qp, const vector<double> &qpar) {
  int iphase = 1;
  for (int ii = 1; ii <= eos.nphase; ++ii) {
    const double det = (qp - qpar.at(ii)) * (qp - qpar.at(ii - 1));
    if (det <= 0.0) {
      iphase = ii;
      break;
    }
  }
  return iphase;
}

void teos_initialize(Eos &eos, const string &teos_path, const string &peos_path) {
  eos = Eos{};
  resize_peos_tables(eos);
  resize_teos_tables(eos);

  {
    TextFile f(teos_path);
    eos.nphase = parse_int(f.line_tokens(1, "TEOS table length").at(0)) - 1;
    require_teos_table_count(eos.nphase, teos_path);
    for (int ii = 0; ii <= eos.nphase; ++ii) {
      const auto t = f.line_tokens(4, "TEOS table row");
      eos.enei.at(ii) = parse_real(t.at(0));
      eos.prei.at(ii) = parse_real(t.at(1));
      eos.hi_teos.at(ii) = parse_real(t.at(2));
      eos.rhoi_teos.at(ii) = parse_real(t.at(3));
      eos.qi_teos.at(ii) = eos.prei.at(ii) / eos.rhoi_teos.at(ii);
    }
  }
  {
    TextFile f(peos_path);
    auto t = f.line_tokens(2, "PEOS header");
    eos.rhoini_cgs = parse_real(t.at(1));
    t = f.line_tokens(2, "PEOS crust row");
    eos.kappa_crust = parse_real(t.at(0));
    eos.gamma_crust = parse_real(t.at(1));
  }
  const double facrho = std::pow(g_cgs / (c_cgs * c_cgs), 3) * solmas * solmas;
  eos.rhoini_gcm1 = facrho * eos.rhoini_cgs;
  const int iphase = teos_lookup(eos, eos.rhoini_gcm1, eos.rhoi_teos);
  const int i0 = std::min(std::max(iphase - lagrange_left_offset, 0), eos.nphase - lagrange_upper_margin);
  array<double, lagrange_order> x4{}, f4{};
  for (int i = 0; i < lagrange_order; ++i) {
    x4[i] = eos.rhoi_teos.at(i0 + i);
    f4[i] = eos.qi_teos.at(i0 + i);
  }
  eos.emdini_gcm1 = lagint_4th(x4, f4, eos.rhoini_gcm1);
  eos.teos_cached_iphase = 1;
}

void teos_q2hprho(Eos &eos, double &q, double &h, double &pre, double &rho, double &ened) {
  if (q <= eos.qi_teos.at(0)) {
    const double hin = eos.hi_teos.at(0);
    const double qin = eos.qi_teos.at(0);
    const double abin = eos.gamma_crust;
    const double abct = eos.kappa_crust;
    const double fac1 = 1.0 / (abin - 1.0);
    const double fac2 = abin / (abin - 1.0);
    const double fack = std::pow(abct, -fac1);
    if (q <= eos_q_floor)
      q = eos_q_floor;
    h = hin + fac2 * (q - qin);
    if (h <= 1.0)
      h = 1.0;
    pre = fack * std::pow(q, fac2);
    rho = fack * std::pow(q, fac1);
    ened = rho * h - pre;
    eos.teos_cached_iphase = 1;
    return;
  }

  int iphase = eos.teos_cached_iphase;
  if (iphase < 1 || iphase > eos.nphase)
    iphase = 1;
  const int ii = iphase;
  if (q < eos.qi_teos.at(iphase)) {
    for (int i = ii; i >= 1; --i) {
      if (q >= eos.qi_teos.at(i - 1) && q <= eos.qi_teos.at(i)) {
        iphase = i;
        break;
      }
    }
  } else {
    for (int i = ii; i <= eos.nphase; ++i) {
      if (q >= eos.qi_teos.at(i - 1) && q <= eos.qi_teos.at(i)) {
        iphase = i;
        break;
      }
    }
  }
  eos.teos_cached_iphase = iphase;

  const int i0 = std::min(std::max(iphase - lagrange_left_offset, 0), eos.nphase - lagrange_upper_margin);
  array<double, lagrange_order> x4{}, f4{};
  for (int i = 0; i < lagrange_order; ++i)
    x4[i] = eos.qi_teos.at(i0 + i);
  for (int i = 0; i < lagrange_order; ++i)
    f4[i] = eos.hi_teos.at(i0 + i);
  h = lagint_4th(x4, f4, q);
  if (h <= 1.0)
    h = 1.0;
  for (int i = 0; i < lagrange_order; ++i)
    f4[i] = eos.prei.at(i0 + i);
  pre = lagint_4th(x4, f4, q);
  for (int i = 0; i < lagrange_order; ++i)
    f4[i] = eos.rhoi_teos.at(i0 + i);
  rho = lagint_4th(x4, f4, q);
  for (int i = 0; i < lagrange_order; ++i)
    f4[i] = eos.enei.at(i0 + i);
  ened = lagint_4th(x4, f4, q);
}

void q2hprho(Eos &eos, bool tabulated, double &q, double &h, double &pre, double &rho, double &ened) {
  if (tabulated)
    teos_q2hprho(eos, q, h, pre, rho, ened);
  else
    peos_q2hprho(eos, q, h, pre, rho, ened);
}

} // namespace cocal_id_cpp
