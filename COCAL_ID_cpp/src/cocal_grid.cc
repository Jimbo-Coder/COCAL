#include "cocal_id_cpp_internal.hh"

namespace cocal_id_cpp {

void grid_r(Grid &gr) {
  gr.rg.resize(0, gr.nrg + 2);
  gr.rginv.resize(0, gr.nrg + 2);
  gr.hrg.resize(1, gr.nrg + 2);
  gr.hrginv.resize(1, gr.nrg + 2);
  gr.drg.resize(1, gr.nrg + 2);
  gr.drginv.resize(1, gr.nrg + 2);

  double drdr = 1.0 / static_cast<double>(gr.nrf);
  double drdrinv = 1.0 / drdr;
  const double rvdom = gr.rgout - gr.rgmid;
  const int nrout = gr.nrg - gr.nrgin;

  for (int ir = 1; ir <= gr.nrgin + 1; ++ir) {
    gr.drg(ir) = drdr;
    gr.drginv(ir) = drdrinv;
  }

  double ratrf = 0.0;
  if (gr.rgin >= 1.0e-13 && gr.rgin < 1.0) {
    const double rfdom = 1.0 - gr.rgin;
    ratrf = 0.1;
    for (;;) {
      const double ratrfb = ratrf;
      const double alge = -(std::pow(ratrf, gr.nrf) - 1.0 - rfdom * drdrinv * (ratrf - 1.0));
      const double dalge = static_cast<double>(gr.nrf) * std::pow(ratrf, gr.nrf - 1) - rfdom * drdrinv;
      ratrf += alge / dalge;
      const double error = 2.0 * (ratrf - ratrfb) / (ratrf + ratrfb);
      if (std::abs(error) <= cocal_zero_tolerance)
        break;
    }
    if (ratrf > 1.0)
      abort_cocal("ratrf not good");
    gr.drg(gr.nrf) = drdr;
    gr.drginv(gr.nrf) = 1.0 / gr.drg(gr.nrf);
    for (int ir = gr.nrf - 1; ir >= 1; --ir) {
      gr.drg(ir) = gr.drg(ir + 1) * ratrf;
      gr.drginv(ir) = 1.0 / gr.drg(ir);
    }
    for (int ir = gr.nrf + 1; ir <= gr.nrgin; ++ir) {
      gr.drg(ir) = gr.drg(gr.nrf);
      gr.drginv(ir) = 1.0 / gr.drg(ir);
    }
  }

  for (int itry = 0; itry <= 3; ++itry) {
    double ratrr = (itry == 0) ? 2.0e-01 : (itry == 1) ? 2.0 : (itry == 2) ? 10.0 : 50.0;
    for (;;) {
      const double ratrrb = ratrr;
      const double alge = -(std::pow(ratrr, nrout + 1) - ratrr - rvdom * drdrinv * (ratrr - 1.0));
      const double dalge = static_cast<double>(nrout + 1) * std::pow(ratrr, nrout) - 1.0 - rvdom * drdrinv;
      ratrr += alge / dalge;
      const double error = 2.0 * (ratrr - ratrrb) / (ratrr + ratrrb);
      if (std::abs(error) <= cocal_zero_tolerance)
        break;
    }
    for (int ir = gr.nrgin + 1; ir <= gr.nrg; ++ir) {
      gr.drg(ir) = ratrr * gr.drg(ir - 1);
      gr.drginv(ir) = 1.0 / gr.drg(ir);
    }
    gr.rg(0) = gr.rgin;
    gr.rginv(0) = (gr.rg(0) <= cocal_zero_tolerance) ? 0.0 : 1.0 / gr.rg(0);
    for (int ir = 1; ir <= gr.nrg; ++ir) {
      gr.rg(ir) = gr.rg(ir - 1) + gr.drg(ir);
      gr.hrg(ir) = 0.5 * (gr.rg(ir) + gr.rg(ir - 1));
      gr.rginv(ir) = 1.0 / gr.rg(ir);
      gr.hrginv(ir) = 1.0 / gr.hrg(ir);
    }
    const double rdet = (gr.rg(gr.nrg) - gr.rgout) / gr.rgout;
    if (std::abs(rdet) < grid_match_tolerance)
      return;
    if (itry == 3)
      abort_cocal("bad coordinate in grid_r");
  }
}

void grid_r_bht(Grid &gr, double mass_bh, double spin_bh) {
  gr.rg.resize(0, gr.nrg + 2);
  gr.rginv.resize(0, gr.nrg + 2);
  gr.hrg.resize(1, gr.nrg + 2);
  gr.hrginv.resize(1, gr.nrg + 2);
  gr.drg.resize(1, gr.nrg + 2);
  gr.drginv.resize(1, gr.nrg + 2);

  constexpr int horizon_intervals = 8;
  if (mass_bh <= 0.0)
    abort_cocal("Invalid BHT parameters: mass_bh must be positive.");
  if (std::abs(spin_bh) > mass_bh)
    abort_cocal("Invalid BHT parameters: abs(spin_bh) exceeds mass_bh.");
  if (gr.rgmid < gr.r_surf || gr.nrgin > gr.nrf || gr.nrgin <= horizon_intervals)
    abort_cocal("Invalid BHT radial-grid parameters.");

  const double horizon_radius = mass_bh + std::sqrt(mass_bh * mass_bh - spin_bh * spin_bh);
  const double uniform_dr = (gr.rgmid - gr.r_surf) / static_cast<double>(gr.nrf - gr.nrgin);
  for (int ir = 1; ir <= gr.nrf + 1; ++ir) {
    gr.drg(ir) = uniform_dr;
    gr.drginv(ir) = 1.0 / uniform_dr;
  }

  const double inner_dr = (horizon_radius - gr.rgin) / horizon_intervals;
  for (int ir = 1; ir <= horizon_intervals; ++ir) {
    gr.drg(ir) = inner_dr;
    gr.drginv(ir) = 1.0 / inner_dr;
  }
  gr.rg(0) = gr.rgin;
  for (int ir = 1; ir <= horizon_intervals; ++ir)
    gr.rg(ir) = gr.rg(ir - 1) + gr.drg(ir);

  const int surface_intervals = gr.nrgin - horizon_intervals;
  const double surface_span = gr.r_surf - horizon_radius;
  bool surface_grid_ok = false;
  for (double ratio : {1.01, 1.1, 2.0, 10.0}) {
    bool converged = false;
    for (int iteration = 0; iteration < 1000; ++iteration) {
      const double previous = ratio;
      const double residual =
          -(std::pow(ratio, surface_intervals) - 1.0 - surface_span / inner_dr * (ratio - 1.0));
      const double derivative = static_cast<double>(surface_intervals) *
                                    std::pow(ratio, surface_intervals - 1) -
                                surface_span / inner_dr;
      ratio += residual / derivative;
      const double error = 2.0 * (ratio - previous) / (ratio + previous);
      if (std::abs(error) <= 1.0e-15) {
        converged = true;
        break;
      }
    }
    if (!converged || !std::isfinite(ratio))
      continue;

    gr.drg(horizon_intervals + 1) = inner_dr;
    gr.drginv(horizon_intervals + 1) = 1.0 / inner_dr;
    for (int ir = horizon_intervals + 2; ir <= gr.nrgin; ++ir) {
      gr.drg(ir) = ratio * gr.drg(ir - 1);
      gr.drginv(ir) = 1.0 / gr.drg(ir);
    }
    for (int ir = horizon_intervals + 1; ir <= gr.nrgin; ++ir)
      gr.rg(ir) = gr.rg(ir - 1) + gr.drg(ir);
    for (int ir = gr.nrgin + 1; ir <= gr.nrf; ++ir)
      gr.rg(ir) = gr.rg(gr.nrgin) + static_cast<double>(ir - gr.nrgin) * uniform_dr;

    if (std::abs(gr.rg(gr.nrf) - gr.rgmid) < 2.0e-14) {
      surface_grid_ok = true;
      break;
    }
  }
  if (!surface_grid_ok)
    abort_cocal("Could not construct BHT radial grid through the torus surface.");

  const int outer_intervals = gr.nrg - gr.nrf;
  const double outer_span = gr.rgout - gr.rgmid;
  bool outer_grid_ok = false;
  for (double ratio : {1.05, 2.0, 10.0, 50.0}) {
    bool converged = false;
    const double first_outer_dr = gr.drg(gr.nrf);
    for (int iteration = 0; iteration < 1000; ++iteration) {
      const double previous = ratio;
      const double denominator = ratio + outer_span / first_outer_dr * (ratio - 1.0);
      const double residual =
          -(static_cast<double>(outer_intervals + 1) * std::log(ratio) - std::log(denominator));
      const double derivative = static_cast<double>(outer_intervals + 1) / ratio -
                                (1.0 + outer_span / first_outer_dr) / denominator;
      ratio += residual / derivative;
      const double error = 2.0 * (ratio - previous) / (ratio + previous);
      if (std::abs(error) <= 1.0e-16) {
        converged = true;
        break;
      }
    }
    if (!converged || !std::isfinite(ratio))
      continue;

    for (int ir = gr.nrf + 1; ir <= gr.nrg; ++ir) {
      gr.drg(ir) = ratio * gr.drg(ir - 1);
      gr.drginv(ir) = 1.0 / gr.drg(ir);
    }
    gr.rg(0) = gr.rgin;
    gr.rginv(0) = 1.0 / gr.rg(0);
    for (int ir = 1; ir <= gr.nrg; ++ir) {
      gr.rg(ir) = gr.rg(ir - 1) + gr.drg(ir);
      gr.hrg(ir) = 0.5 * (gr.rg(ir) + gr.rg(ir - 1));
      gr.rginv(ir) = 1.0 / gr.rg(ir);
      gr.hrginv(ir) = 1.0 / gr.hrg(ir);
    }
    if (std::abs((gr.rg(gr.nrg) - gr.rgout) / gr.rgout) < 5.0e-14) {
      outer_grid_ok = true;
      break;
    }
  }
  if (!outer_grid_ok)
    abort_cocal("Could not construct BHT outer radial grid.");
}

void grid_r_bqs(Grid &gr) {
  gr.rg.resize(0, gr.nrg + 2);
  gr.rginv.resize(0, gr.nrg + 2);
  gr.hrg.resize(1, gr.nrg + 2);
  gr.hrginv.resize(1, gr.nrg + 2);
  gr.drg.resize(1, gr.nrg + 2);
  gr.drginv.resize(1, gr.nrg + 2);

  if (gr.r_surf != 1.0)
    gr.r_surf = 1.0;

  if (gr.rgin > 0.0) {
    double drdr = 1.0 / static_cast<double>(gr.nrf);
    double drdrinv = 1.0 / drdr;
    for (int ir = 1; ir <= gr.nrgin + 1; ++ir) {
      gr.drg(ir) = drdr;
      gr.drginv(ir) = drdrinv;
    }
    const double rvdom = gr.rgout - gr.rgmid;
    const int nrout = gr.nrg - gr.nrgin;
    for (int itry = 0; itry <= 3; ++itry) {
      double ratrr = (itry == 0) ? 2.0e-01 : (itry == 1) ? 2.0 : (itry == 2) ? 10.0 : 50.0;
      for (;;) {
        const double ratrrb = ratrr;
        const double alge = -(std::pow(ratrr, nrout + 1) - ratrr - rvdom * drdrinv * (ratrr - 1.0));
        const double dalge = static_cast<double>(nrout + 1) * std::pow(ratrr, nrout) - 1.0 - rvdom * drdrinv;
        ratrr += alge / dalge;
        const double error = 2.0 * (ratrr - ratrrb) / (ratrr + ratrrb);
        if (std::abs(error) <= cocal_zero_tolerance)
          break;
      }
      for (int ir = gr.nrgin + 1; ir <= gr.nrg; ++ir) {
        gr.drg(ir) = ratrr * gr.drg(ir - 1);
        gr.drginv(ir) = 1.0 / gr.drg(ir);
      }
      gr.rg(0) = gr.rgin;
      gr.rginv(0) = (gr.rg(0) <= cocal_zero_tolerance) ? 0.0 : 1.0 / gr.rg(0);
      for (int ir = 1; ir <= gr.nrg; ++ir) {
        gr.rg(ir) = gr.rg(ir - 1) + gr.drg(ir);
        gr.hrg(ir) = 0.5 * (gr.rg(ir) + gr.rg(ir - 1));
        gr.rginv(ir) = 1.0 / gr.rg(ir);
        gr.hrginv(ir) = 1.0 / gr.hrg(ir);
      }
      const double rdet = std::abs(gr.rg(gr.nrg) - gr.rgout) / gr.rgout;
      if (rdet < grid_match_tolerance)
        return;
      if (itry == 3)
        abort_cocal("bad coordinate in grid_r_bqs");
    }
  } else {
    double drdr = 1.0 / static_cast<double>(gr.nrf);
    const int nr1 = gr.nrf + static_cast<int>((1.25 - 1.0) / drdr);
    for (int ir = 1; ir <= nr1; ++ir) {
      gr.drg(ir) = drdr;
      gr.drginv(ir) = 1.0 / drdr;
    }
    const double dr2 = 2.0 * drdr;
    int nr2 = nr1 + static_cast<int>((5.0 - 1.25) / dr2);
    for (int ir = nr1 + 1; ir <= nr2; ++ir) {
      gr.drg(ir) = dr2;
      gr.drginv(ir) = 1.0 / dr2;
    }
    gr.rg(0) = gr.rgin;
    for (int ir = 1; ir <= nr2; ++ir)
      gr.rg(ir) = gr.rg(ir - 1) + gr.drg(ir);

    for (int ir = 1; ir <= nr2; ++ir) {
      if ((gr.rg(ir) - grid_match_tolerance) <= gr.rgmid &&
          gr.rgmid < (gr.rg(ir + 1) - grid_match_tolerance)) {
        if (std::abs(gr.rg(ir) - gr.rgmid) <= cocal_zero_tolerance) {
          gr.nrgin = ir;
        } else {
          gr.drg(ir + 1) = gr.rgmid - gr.rg(ir);
          gr.drg(ir + 2) = gr.rg(ir + 1) - gr.rgmid;
          gr.nrgin = ir + 1;
          for (int jr = ir + 3; jr <= nr2; ++jr)
            gr.drg(jr) = dr2;
        }
        break;
      }
    }
    gr.rg(0) = gr.rgin;
    for (int ir = 1; ir <= nr2; ++ir)
      gr.rg(ir) = gr.rg(ir - 1) + gr.drg(ir);

    for (int itry = 0; itry <= 3; ++itry) {
      double ratrr = (itry == 0) ? 2.0e-01 : (itry == 1) ? 2.0 : (itry == 2) ? 10.0 : 50.0;
      drdr = gr.drg(nr2);
      const double drdrinv = 1.0 / dr2;
      const double rvdom = gr.rgout - gr.rg(nr2);
      const int nrout = gr.nrg - nr2;
      for (;;) {
        const double ratrrb = ratrr;
        const double alge = -(std::pow(ratrr, nrout + 1) - ratrr - rvdom * drdrinv * (ratrr - 1.0));
        const double dalge = static_cast<double>(nrout + 1) * std::pow(ratrr, nrout) - 1.0 - rvdom * drdrinv;
        ratrr += alge / dalge;
        const double error = 2.0 * (ratrr - ratrrb) / (ratrr + ratrrb);
        if (std::abs(error) <= cocal_zero_tolerance)
          break;
      }
      for (int ir = nr2 + 1; ir <= gr.nrg; ++ir) {
        gr.drg(ir) = ratrr * gr.drg(ir - 1);
        gr.drginv(ir) = 1.0 / gr.drg(ir);
      }
      gr.rg(0) = gr.rgin;
      gr.rginv(0) = 0.0;
      for (int ir = 1; ir <= gr.nrg; ++ir) {
        gr.rg(ir) = gr.rg(ir - 1) + gr.drg(ir);
        gr.hrg(ir) = 0.5 * (gr.rg(ir) + gr.rg(ir - 1));
        gr.rginv(ir) = 1.0 / gr.rg(ir);
        gr.hrginv(ir) = 1.0 / gr.hrg(ir);
      }
      const double rdet = std::abs(gr.rg(gr.nrg) - gr.rgout) / gr.rgout;
      if (rdet < grid_match_tolerance)
        return;
      if (itry == 3)
        abort_cocal("bad coordinate in grid_r_bqs");
    }
  }
}

void grid_theta(Grid &gr) {
  gr.thg.resize(0, gr.ntg + 2);
  gr.hthg.resize(1, gr.ntg + 2);
  gr.dthg = pi / static_cast<double>(gr.ntg);
  gr.dthginv = 1.0 / gr.dthg;
  gr.thg(0) = 0.0;
  for (int it = 1; it <= gr.ntg; ++it) {
    gr.thg(it) = static_cast<double>(it) * gr.dthg;
    gr.hthg(it) = 0.5 * (gr.thg(it) + gr.thg(it - 1));
  }
}

void grid_phi(Grid &gr) {
  gr.phig.resize(0, gr.npg + 2);
  gr.hphig.resize(1, gr.npg + 2);
  gr.dphig = 2.0 * pi / static_cast<double>(gr.npg);
  gr.dphiginv = 1.0 / gr.dphig;
  gr.phig(0) = 0.0;
  for (int ip = 1; ip <= gr.npg; ++ip) {
    gr.phig(ip) = static_cast<double>(ip) * gr.dphig;
    gr.hphig(ip) = 0.5 * (gr.phig(ip) + gr.phig(ip - 1));
  }
}

void grid_extended(Grid &gr) {
  gr.rgex.resize(-2, gr.nrg + 2);
  gr.thgex.resize(-2, gr.ntg + 2);
  gr.phigex.resize(-2, gr.npg + 2);
  gr.irgex_r.resize(-2, gr.nrg + 2);
  gr.itgex_th.resize(-2, gr.ntg + 2);
  gr.ipgex_phi.resize(-2, gr.npg + 2);
  gr.itgex_r.resize(0, gr.ntg, -2, gr.nrg + 2);
  gr.ipgex_r.resize(0, gr.npg, -2, gr.nrg + 2);
  gr.ipgex_th.resize(0, gr.npg, -2, gr.ntg + 2);

  for (int ir = 0; ir <= gr.nrg; ++ir)
    gr.rgex(ir) = gr.rg(ir);
  for (int it = 0; it <= gr.ntg; ++it)
    gr.thgex(it) = gr.thg(it);
  for (int ip = 0; ip <= gr.npg; ++ip)
    gr.phigex(ip) = gr.phig(ip);

  gr.rgex(-1) = gr.rg(0) - (gr.rg(1) - gr.rg(0));
  gr.rgex(-2) = gr.rg(0) - (gr.rg(2) - gr.rg(0));
  gr.rgex(gr.nrg + 1) = gr.rg(gr.nrg) + (gr.rg(gr.nrg) - gr.rg(gr.nrg - 1));
  gr.rgex(gr.nrg + 2) = gr.rg(gr.nrg) + 2.0 * (gr.rg(gr.nrg) - gr.rg(gr.nrg - 1));
  gr.thgex(-1) = -gr.thg(1);
  gr.thgex(-2) = -gr.thg(2);
  gr.thgex(gr.ntg + 1) = 2.0 * gr.thg(gr.ntg) - gr.thg(gr.ntg - 1);
  gr.thgex(gr.ntg + 2) = 2.0 * gr.thg(gr.ntg) - gr.thg(gr.ntg - 2);
  gr.phigex(-1) = gr.phig(gr.npg - 1) - gr.phig(gr.npg);
  gr.phigex(-2) = gr.phig(gr.npg - 2) - gr.phig(gr.npg);
  gr.phigex(gr.npg + 1) = gr.phig(gr.npg) + gr.phig(1);
  gr.phigex(gr.npg + 2) = gr.phig(gr.npg) + gr.phig(2);

  for (int ir = -2; ir <= gr.nrg; ++ir) {
    if (ir >= 0 && ir <= gr.nrg)
      gr.irgex_r(ir) = ir;
    if (ir <= -1)
      gr.irgex_r(ir) = std::abs(ir);
  }
  for (int it = -2; it <= gr.ntg + 2; ++it) {
    if (it >= 0 && it <= gr.ntg)
      gr.itgex_th(it) = it;
    if (it <= -1)
      gr.itgex_th(it) = std::abs(it);
    if (it >= gr.ntg + 1)
      gr.itgex_th(it) = 2 * gr.ntg - it;
  }
  for (int ip = -2; ip <= gr.npg + 2; ++ip) {
    if (ip >= 0 && ip <= gr.npg)
      gr.ipgex_phi(ip) = ip;
    if (ip <= -1)
      gr.ipgex_phi(ip) = gr.npg + ip;
    if (ip >= gr.npg + 1)
      gr.ipgex_phi(ip) = ip - gr.npg;
  }
  for (int ir = -2; ir <= gr.nrg; ++ir) {
    for (int it = 0; it <= gr.ntg; ++it)
      gr.itgex_r(it, ir) = (ir >= 0) ? it : gr.ntg - it;
    for (int ip = 0; ip <= gr.npg; ++ip)
      gr.ipgex_r(ip, ir) = (ir >= 0) ? ip : ((ip + gr.npg / 2) % gr.npg);
  }
  for (int it = -2; it <= gr.ntg + 2; ++it) {
    for (int ip = 0; ip <= gr.npg; ++ip) {
      if (it >= 0 && it <= gr.ntg)
        gr.ipgex_th(ip, it) = ip;
      else
        gr.ipgex_th(ip, it) = (ip + gr.npg / 2) % gr.npg;
    }
  }
}

void setup_coordinates(Grid &gr, bool bns_grid) {
  if (bns_grid)
    grid_r_bqs(gr);
  else
    grid_r(gr);
  grid_theta(gr);
  grid_phi(gr);
  grid_extended(gr);
}

void calc_parameter_binary_excision(Patch &p) {
  auto &ex = p.ex;
  auto &gr = p.grid;
  if (ex.ex_nrg == 0) {
    ex.sepa = 0.0;
    ex.dis = 0.0;
    ex.ex_radius = 0.0;
    ex.ex_rgin = 0.0;
    ex.ex_rgmid = 0.0;
    ex.ex_rgout = 0.0;
    return;
  }
  ex.sepa = 2.0 * gr.rg(ex.ex_nrg + ex.ex_ndis);
  ex.dis = gr.rg(ex.ex_nrg + ex.ex_ndis);
  ex.ex_radius = gr.rg(ex.ex_nrg);
  ex.ex_rgin = ex.dis + ex.dis - gr.rg(ex.ex_nrg);
  ex.ex_rgmid = ex.sepa;
  ex.ex_rgout = ex.sepa + gr.rg(ex.ex_nrg);
}

} // namespace cocal_id_cpp
