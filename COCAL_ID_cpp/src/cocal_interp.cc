#include "cocal_id_cpp_internal.hh"

namespace cocal_id_cpp {

void lagint_4th_weights(const array<double, lagrange_order> &x, double v, array<double, lagrange_order> &w) {
  const double dx12 = x[0] - x[1];
  const double dx13 = x[0] - x[2];
  const double dx14 = x[0] - x[3];
  const double dx23 = x[1] - x[2];
  const double dx24 = x[1] - x[3];
  const double dx34 = x[2] - x[3];
  const double dx21 = -dx12;
  const double dx31 = -dx13;
  const double dx32 = -dx23;
  const double dx41 = -dx14;
  const double dx42 = -dx24;
  const double dx43 = -dx34;
  const double xv1 = v - x[0];
  const double xv2 = v - x[1];
  const double xv3 = v - x[2];
  const double xv4 = v - x[3];
  w[0] = xv2 * xv3 * xv4 / (dx12 * dx13 * dx14);
  w[1] = xv1 * xv3 * xv4 / (dx21 * dx23 * dx24);
  w[2] = xv1 * xv2 * xv4 / (dx31 * dx32 * dx34);
  w[3] = xv1 * xv2 * xv3 / (dx41 * dx42 * dx43);
}

double lagint_4th_apply(const array<double, lagrange_order> &w, const array<double, lagrange_order> &y) {
  return w[0] * y[0] + w[1] * y[1] + w[2] * y[2] + w[3] * y[3];
}

double lagint_4th(const array<double, lagrange_order> &x, const array<double, lagrange_order> &y, double v) {
  array<double, lagrange_order> w{};
  lagint_4th_weights(x, v, w);
  return lagint_4th_apply(w, y);
}

void spherical_from_cart(double x, double y, double z, double &r, double &th, double &ph) {
  r = std::sqrt(std::abs(x * x + y * y + z * z));
  const double varpi = std::sqrt(std::abs(x * x + y * y));
  th = std::fmod(2.0 * pi + std::atan2(varpi, z), 2.0 * pi);
  ph = std::fmod(2.0 * pi + std::atan2(y, x), 2.0 * pi);
}

Stencil4 cgr_4th_setup_from_spherical(const Grid &gr, double rc, double thc, double phic, int nr) {
  Stencil4 st;
  int ir0 = 0, it0 = 0, ip0 = 0;
  for (int ir = 0; ir <= nr + 1; ++ir) {
    if (rc < gr.rgex(ir) && rc >= gr.rgex(ir - 1))
      ir0 = std::min(ir - lagrange_left_offset, nr - lagrange_upper_margin);
  }
  for (int it = 0; it <= gr.ntg + 1; ++it) {
    if (thc < gr.thgex(it) && thc >= gr.thgex(it - 1))
      it0 = it - lagrange_left_offset;
  }
  for (int ip = 0; ip <= gr.npg + 1; ++ip) {
    if (phic < gr.phigex(ip) && phic >= gr.phigex(ip - 1))
      ip0 = ip - lagrange_left_offset;
  }

  array<double, lagrange_order> r4{}, th4{}, phi4{};
  for (int ii = 0; ii < lagrange_order; ++ii) {
    r4[ii] = gr.rgex(ir0 + ii);
    th4[ii] = gr.thgex(it0 + ii);
    phi4[ii] = gr.phigex(ip0 + ii);
  }
  lagint_4th_weights(r4, rc, st.wr);
  lagint_4th_weights(th4, thc, st.wth);
  lagint_4th_weights(phi4, phic, st.wphi);

  for (int kk = 0; kk < lagrange_order; ++kk) {
    const int ipg0 = ip0 + kk;
    for (int jj = 0; jj < lagrange_order; ++jj) {
      const int itg0 = it0 + jj;
      for (int ii = 0; ii < lagrange_order; ++ii) {
        const int irg0 = ir0 + ii;
        st.ir[ii][jj][kk] = gr.irgex_r(irg0);
        st.it[ii][jj][kk] = gr.itgex_r(gr.itgex_th(itg0), irg0);
        st.ip[ii][jj][kk] = gr.ipgex_r(gr.ipgex_th(gr.ipgex_phi(ipg0), itg0), irg0);
      }
    }
  }
  return st;
}

Stencil4 gr2cgr_4th_setup(const Grid &gr, double x, double y, double z) {
  double r = 0.0, th = 0.0, ph = 0.0;
  spherical_from_cart(x, y, z, r, th, ph);
  return cgr_4th_setup_from_spherical(gr, r, th, ph, gr.nrg);
}

double cgr_4th_apply(const Array3D<double> &fnc, const Stencil4 &st) {
  array<double, lagrange_order> fr4{}, ft4{}, fp4{};
  for (int kk = 0; kk < lagrange_order; ++kk) {
    for (int jj = 0; jj < lagrange_order; ++jj) {
      for (int ii = 0; ii < lagrange_order; ++ii)
        fr4[ii] = fnc(st.ir[ii][jj][kk], st.it[ii][jj][kk], st.ip[ii][jj][kk]);
      ft4[jj] = lagint_4th_apply(st.wr, fr4);
    }
    fp4[kk] = lagint_4th_apply(st.wth, ft4);
  }
  return lagint_4th_apply(st.wphi, fp4);
}

double interpo_lag4th_2Dsurf(const Grid &gr, const Array2D<double> &fnc, double tv, double pv) {
  int it0 = 0, ip0 = 0;
  for (int it = 0; it <= gr.ntg + 1; ++it) {
    if (tv < gr.thgex(it) && tv >= gr.thgex(it - 1))
      it0 = it - lagrange_left_offset;
  }
  for (int ip = 0; ip <= gr.npg + 1; ++ip) {
    if (pv < gr.phigex(ip) && pv >= gr.phigex(ip - 1))
      ip0 = ip - lagrange_left_offset;
  }
  array<double, lagrange_order> th4{}, phi4{}, ft4{}, fp4{};
  for (int ii = 0; ii < lagrange_order; ++ii) {
    th4[ii] = gr.thgex(it0 + ii);
    phi4[ii] = gr.phigex(ip0 + ii);
  }
  for (int kk = 0; kk < lagrange_order; ++kk) {
    const int ipg0 = ip0 + kk;
    for (int jj = 0; jj < lagrange_order; ++jj) {
      const int itg0 = it0 + jj;
      const int itgex = gr.itgex_th(itg0);
      const int ipgex = gr.ipgex_th(gr.ipgex_phi(ipg0), itg0);
      ft4[jj] = fnc(itgex, ipgex);
    }
    fp4[kk] = lagint_4th(th4, ft4, tv);
  }
  return lagint_4th(phi4, fp4, pv);
}

bool fl2cgr_4th_setup(const Grid &gr, const Array2D<double> &rs, double x, double y, double z, Stencil4 &st) {
  double rc_gr = 0.0, th = 0.0, ph = 0.0;
  spherical_from_cart(x, y, z, rc_gr, th, ph);
  const double rsca = interpo_lag4th_2Dsurf(gr, rs, th, ph);
  const double rc = rc_gr / rsca;
  if (rc > gr.rg(gr.nrf))
    return false;
  st = cgr_4th_setup_from_spherical(gr, rc, th, ph, gr.nrf);
  return true;
}

void interpo_gr2fl_export(const Grid &gr, const Array3D<double> &grv, Array3D<double> &flv,
                          const Array2D<double> &rs) {
  flv.fill(0.0);
  for (int ipf = 0; ipf <= gr.npf; ++ipf) {
    for (int itf = 0; itf <= gr.ntf; ++itf) {
      for (int irf = 0; irf <= gr.nrf; ++irf) {
        const double rrff = rs(itf, ipf) * gr.rg(irf);
        int ir0 = 0;
        for (int irg = 0; irg <= gr.nrg - 1; ++irg) {
          if (rrff <= gr.rg(irg)) {
            ir0 = std::min(std::max(0, irg - lagrange_left_offset), gr.nrg - lagrange_upper_margin);
            break;
          }
        }
        array<double, lagrange_order> x{}, f{};
        for (int i = 0; i < lagrange_order; ++i) {
          x[i] = gr.rg(ir0 + i);
          f[i] = grv(ir0 + i, itf, ipf);
        }
        flv(irf, itf, ipf) = lagint_4th(x, f, rrff);
      }
    }
  }
}

void interpo_gr2fl_metric_CF_export(const Grid &gr, const Array3D<double> &alph,
                                    const Array3D<double> &psi, const Array3D<double> &bvxd,
                                    const Array3D<double> &bvyd, const Array3D<double> &bvzd,
                                    Array3D<double> &alphf, Array3D<double> &psif,
                                    Array3D<double> &bvxdf, Array3D<double> &bvydf,
                                    Array3D<double> &bvzdf, const Array2D<double> &rs) {
  interpo_gr2fl_export(gr, alph, alphf, rs);
  interpo_gr2fl_export(gr, psi, psif, rs);
  interpo_gr2fl_export(gr, bvxd, bvxdf, rs);
  interpo_gr2fl_export(gr, bvyd, bvydf, rs);
  interpo_gr2fl_export(gr, bvzd, bvzdf, rs);
}

void invhij_WL_export(const Grid &gr, const Array3D<double> &hxxd, const Array3D<double> &hxyd,
                      const Array3D<double> &hxzd, const Array3D<double> &hyyd,
                      const Array3D<double> &hyzd, const Array3D<double> &hzzd,
                      Array3D<double> &hxxu, Array3D<double> &hxyu, Array3D<double> &hxzu,
                      Array3D<double> &hyyu, Array3D<double> &hyzu, Array3D<double> &hzzu) {
  for (int ip = 0; ip <= gr.npg; ++ip) {
    for (int it = 0; it <= gr.ntg; ++it) {
      for (int ir = 0; ir <= gr.nrg; ++ir) {
        const double hxx = hxxd(ir, it, ip);
        const double hxy = hxyd(ir, it, ip);
        const double hxz = hxzd(ir, it, ip);
        const double hyy = hyyd(ir, it, ip);
        const double hyz = hyzd(ir, it, ip);
        const double hzz = hzzd(ir, it, ip);
        const double hod1 = hxx + hyy + hzz;
        const double hod2 = hxx * hyy + hxx * hzz + hyy * hzz - hxy * hxy - hxz * hxz - hyz * hyz;
        const double hod3 = hxx * hyy * hzz + hxy * hyz * hxz + hxz * hxy * hyz -
                            hxx * hyz * hyz - hxy * hxy * hzz - hxz * hyy * hxz;
        const double detgmi = 1.0 / (1.0 + hod1 + hod2 + hod3);

        const double gmxxu = (1.0 + hyy + hzz + hyy * hzz - hyz * hyz) * detgmi;
        const double gmxyu = (-hxy + hxz * hyz - hxy * hzz) * detgmi;
        const double gmxzu = (-hxz + hxy * hyz - hxz * hyy) * detgmi;
        const double gmyyu = (1.0 + hxx + hzz + hxx * hzz - hxz * hxz) * detgmi;
        const double gmyzu = (-hyz + hxz * hxy - hxx * hyz) * detgmi;
        const double gmzzu = (1.0 + hxx + hyy + hxx * hyy - hxy * hxy) * detgmi;

        hxxu(ir, it, ip) = gmxxu - 1.0;
        hxyu(ir, it, ip) = gmxyu;
        hxzu(ir, it, ip) = gmxzu;
        hyyu(ir, it, ip) = gmyyu - 1.0;
        hyzu(ir, it, ip) = gmyzu;
        hzzu(ir, it, ip) = gmzzu - 1.0;
      }
    }
  }
}

void index_vec_down2up_export(const Grid &gr, const Array3D<double> &hxxu, const Array3D<double> &hxyu,
                              const Array3D<double> &hxzu, const Array3D<double> &hyyu,
                              const Array3D<double> &hyzu, const Array3D<double> &hzzu,
                              Array3D<double> &vecxu, Array3D<double> &vecyu, Array3D<double> &veczu,
                              const Array3D<double> &vecxd, const Array3D<double> &vecyd,
                              const Array3D<double> &veczd) {
  for (int ip = 0; ip <= gr.npg; ++ip) {
    for (int it = 0; it <= gr.ntg; ++it) {
      for (int ir = 0; ir <= gr.nrg; ++ir) {
        const double gmxxu = 1.0 + hxxu(ir, it, ip);
        const double gmxyu = hxyu(ir, it, ip);
        const double gmxzu = hxzu(ir, it, ip);
        const double gmyyu = 1.0 + hyyu(ir, it, ip);
        const double gmyzu = hyzu(ir, it, ip);
        const double gmzzu = 1.0 + hzzu(ir, it, ip);
        vecxu(ir, it, ip) = gmxxu * vecxd(ir, it, ip) + gmxyu * vecyd(ir, it, ip) +
                            gmxzu * veczd(ir, it, ip);
        vecyu(ir, it, ip) = gmxyu * vecxd(ir, it, ip) + gmyyu * vecyd(ir, it, ip) +
                            gmyzu * veczd(ir, it, ip);
        veczu(ir, it, ip) = gmxzu * vecxd(ir, it, ip) + gmyzu * vecyd(ir, it, ip) +
                            gmzzu * veczd(ir, it, ip);
      }
    }
  }
}

} // namespace cocal_id_cpp
