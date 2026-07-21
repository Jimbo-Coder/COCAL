#include "cocal_id_cpp_internal.hh"

namespace cocal_id_cpp {

void allocate_bht_arrays(BhtData &d) {
  const auto &g = d.grid;
  for (Array3D<double> *a : {&d.emd,  &d.omeg, &d.psi,  &d.alph, &d.bvxd, &d.bvyd, &d.bvzd,
                             &d.bvxu, &d.bvyu, &d.bvzu, &d.hxxd, &d.hxyd, &d.hxzd, &d.hyyd,
                             &d.hyzd, &d.hzzd, &d.hxxu, &d.hxyu, &d.hxzu, &d.hyyu, &d.hyzu,
                             &d.hzzu, &d.kxxa, &d.kxya, &d.kxza, &d.kyya, &d.kyza, &d.kzza}) {
    resize_like_grid(*a, g.nrg, g.ntg, g.npg);
    a->fill(0.0);
  }
}

void allocate_rns_arrays(RnsData &d, bool mrns, bool read_potential) {
  const auto &g = d.grid;
  resize_like_fluid(d.emd, g.nrf, g.ntf, g.npf);
  resize_like_fluid(d.omef, g.nrf, g.ntf, g.npf);
  resize_like_fluid(d.psif, g.nrf, g.ntf, g.npf);
  resize_like_fluid(d.alphf, g.nrf, g.ntf, g.npf);
  resize_like_fluid(d.bvxuf, g.nrf, g.ntf, g.npf);
  resize_like_fluid(d.bvyuf, g.nrf, g.ntf, g.npf);
  resize_like_fluid(d.bvzuf, g.nrf, g.ntf, g.npf);
  d.rs.resize(0, g.ntf, 0, g.npf);

  for (Array3D<double> *a : {&d.psi, &d.alph, &d.bvxd, &d.bvyd, &d.bvzd, &d.bvxu, &d.bvyu, &d.bvzu,
                             &d.hxxd, &d.hxyd, &d.hxzd, &d.hyyd, &d.hyzd, &d.hzzd, &d.hxxu, &d.hxyu,
                             &d.hxzu, &d.hyyu, &d.hyzu, &d.hzzu, &d.kxxa, &d.kxya, &d.kxza, &d.kyya,
                             &d.kyza, &d.kzza}) {
    resize_like_grid(*a, g.nrg, g.ntg, g.npg);
    a->fill(0.0);
  }
  for (Array3D<double> *a : {&d.emd, &d.omef, &d.psif, &d.alphf, &d.bvxuf, &d.bvyuf, &d.bvzuf})
    a->fill(0.0);
  d.rs.fill(0.0);

  if (mrns) {
    for (Array3D<double> *a : {&d.utf, &d.uxf, &d.uyf, &d.uzf}) {
      resize_like_fluid(*a, g.nrf, g.ntf, g.npf);
      a->fill(0.0);
    }
    if (read_potential) {
      for (Array3D<double> *a : {&d.va, &d.vaxd, &d.vayd, &d.vazd}) {
        resize_like_grid(*a, g.nrg, g.ntg, g.npg);
        a->fill(0.0);
      }
    } else {
      for (Array3D<double> *a : {&d.fxd, &d.fyd, &d.fzd, &d.fxyd, &d.fxzd, &d.fyzd}) {
        resize_like_grid(*a, g.nrg, g.ntg, g.npg);
        a->fill(0.0);
      }
    }
  }
}

void allocate_patch_arrays(Patch &p, const string &id_type, bool star_patch) {
  const auto &g = p.grid;
  for (Array3D<double> *a : {&p.psi, &p.alph, &p.bvxd, &p.bvyd, &p.bvzd, &p.axx, &p.axy, &p.axz,
                             &p.ayy, &p.ayz, &p.azz}) {
    resize_like_grid(*a, g.nrg, g.ntg, g.npg);
    a->fill(0.0);
  }
  if (!star_patch)
    return;

  resize_like_fluid(p.emd, g.nrf, g.ntf, g.npf);
  p.emd.fill(0.0);
  p.rs.resize(0, g.ntf, 0, g.npf);
  p.rs.fill(0.0);
  for (Array3D<double> *a : {&p.psif, &p.alphf, &p.bvxdf, &p.bvydf, &p.bvzdf}) {
    resize_like_fluid(*a, g.nrf, g.ntf, g.npf);
    a->fill(0.0);
  }
  if (id_type == "IR" || id_type == "SP") {
    for (Array3D<double> *a : {&p.vep, &p.vepxf, &p.vepyf, &p.vepzf}) {
      resize_like_fluid(*a, g.nrf, g.ntf, g.npf);
      a->fill(0.0);
    }
    if (id_type == "SP") {
      for (Array3D<double> *a : {&p.wxspf, &p.wyspf, &p.wzspf}) {
        resize_like_fluid(*a, g.nrf, g.ntf, g.npf);
        a->fill(0.0);
      }
    }
  }
}

void sign_flip(Array3D<double> &a) {
  if (a.empty())
    return;
  for (int k = a.lo3(); k <= a.hi3(); ++k)
    for (int j = a.lo2(); j <= a.hi2(); ++j)
      for (int i = a.lo1(); i <= a.hi1(); ++i)
        a(i, j, k) = -a(i, j, k);
}

RnsPoint interpolate_rns_point(const RnsData &d, const string &rnstype, bool read_magnetic_potential,
                               double xcoc, double ycoc, double zcoc) {
  const bool mrns = rnstype == "MRNS_WL";
  const auto gr_st = gr2cgr_4th_setup(d.grid, xcoc, ycoc, zcoc);

  const double psica = cgr_4th_apply(d.psi, gr_st);
  const double psica4 = std::pow(psica, 4);
  const double hxxdca = cgr_4th_apply(d.hxxd, gr_st);
  const double hxydca = cgr_4th_apply(d.hxyd, gr_st);
  const double hxzdca = cgr_4th_apply(d.hxzd, gr_st);
  const double hyydca = cgr_4th_apply(d.hyyd, gr_st);
  const double hyzdca = cgr_4th_apply(d.hyzd, gr_st);
  const double hzzdca = cgr_4th_apply(d.hzzd, gr_st);

  RnsPoint point;
  point.alphca = cgr_4th_apply(d.alph, gr_st);
  point.bvxdca = cgr_4th_apply(d.bvxd, gr_st);
  point.bvydca = cgr_4th_apply(d.bvyd, gr_st);
  point.bvzdca = cgr_4th_apply(d.bvzd, gr_st);
  point.gxx1 = psica4 * (1.0 + hxxdca);
  point.gxy1 = psica4 * hxydca;
  point.gxz1 = psica4 * hxzdca;
  point.gyy1 = psica4 * (1.0 + hyydca);
  point.gyz1 = psica4 * hyzdca;
  point.gzz1 = psica4 * (1.0 + hzzdca);
  point.kxx1 = psica4 * cgr_4th_apply(d.kxxa, gr_st) / d.radi;
  point.kxy1 = psica4 * cgr_4th_apply(d.kxya, gr_st) / d.radi;
  point.kxz1 = psica4 * cgr_4th_apply(d.kxza, gr_st) / d.radi;
  point.kyy1 = psica4 * cgr_4th_apply(d.kyya, gr_st) / d.radi;
  point.kyz1 = psica4 * cgr_4th_apply(d.kyza, gr_st) / d.radi;
  point.kzz1 = psica4 * cgr_4th_apply(d.kzza, gr_st) / d.radi;

  Stencil4 fl_st;
  double omefca = 0.0, alphfca = 0.0;
  double bvxufca = 0.0, bvyufca = 0.0, bvzufca = 0.0;
  double utfca = 0.0, uxfca = 0.0, uyfca = 0.0, uzfca = 0.0;
  if (fl2cgr_4th_setup(d.grid, d.rs, xcoc, ycoc, zcoc, fl_st)) {
    point.emdca = cgr_4th_apply(d.emd, fl_st);
    omefca = cgr_4th_apply(d.omef, fl_st);
    bvxufca = cgr_4th_apply(d.bvxuf, fl_st);
    bvyufca = cgr_4th_apply(d.bvyuf, fl_st);
    bvzufca = cgr_4th_apply(d.bvzuf, fl_st);
    alphfca = cgr_4th_apply(d.alphf, fl_st);
    if (mrns) {
      utfca = cgr_4th_apply(d.utf, fl_st);
      uxfca = cgr_4th_apply(d.uxf, fl_st);
      uyfca = cgr_4th_apply(d.uyf, fl_st);
      uzfca = cgr_4th_apply(d.uzf, fl_st);
    }
  }

  if (mrns) {
    if (read_magnetic_potential) {
      point.vaca = cgr_4th_apply(d.va, gr_st);
      point.vaxdca = cgr_4th_apply(d.vaxd, gr_st);
      point.vaydca = cgr_4th_apply(d.vayd, gr_st);
      point.vazdca = cgr_4th_apply(d.vazd, gr_st);
    } else {
      point.fxydca = cgr_4th_apply(d.fxyd, gr_st);
      point.fxzdca = cgr_4th_apply(d.fxzd, gr_st);
      point.fyzdca = cgr_4th_apply(d.fyzd, gr_st);
    }
  }

  const double bxcor = bvxufca + omefca * (-ycoc);
  const double bycor = bvyufca + omefca * xcoc;
  const double bzcor = bvzufca;

  if (std::abs(point.emdca) > cocal_zero_tolerance) {
    if (rnstype == "RNS_WL" || rnstype == "RNS_CF") {
      point.vxu = bxcor / alphfca;
      point.vyu = bycor / alphfca;
      point.vzu = bzcor / alphfca;
    } else if (mrns) {
      point.vxu = (uxfca / utfca + bvxufca) / alphfca;
      point.vyu = (uyfca / utfca + bvyufca) / alphfca;
      point.vzu = (uzfca / utfca + bvzufca) / alphfca;
    }
  } else {
    point.emdca = 0.0;
  }

  if (point.emdca != point.emdca)
    point.emdca = 0.0;
  return point;
}

MetricValues interpolate_metric_patch(const Patch &p, double x, double y, double z) {
  const auto st = gr2cgr_4th_setup(p.grid, x, y, z);
  MetricValues v;
  v.psica = cgr_4th_apply(p.psi, st);
  v.alphca = cgr_4th_apply(p.alph, st);
  v.bvxdca = cgr_4th_apply(p.bvxd, st);
  v.bvydca = cgr_4th_apply(p.bvyd, st);
  v.bvzdca = cgr_4th_apply(p.bvzd, st);
  v.axx = cgr_4th_apply(p.axx, st);
  v.axy = cgr_4th_apply(p.axy, st);
  v.axz = cgr_4th_apply(p.axz, st);
  v.ayy = cgr_4th_apply(p.ayy, st);
  v.ayz = cgr_4th_apply(p.ayz, st);
  v.azz = cgr_4th_apply(p.azz, st);
  return v;
}

FluidValues interpolate_fluid_patch(const Patch &p, const string &id_type, double x, double y, double z) {
  FluidValues v;
  double rc = 0.0, th = 0.0, ph = 0.0;
  spherical_from_cart(x, y, z, rc, th, ph);
  const double rsca = interpo_lag4th_2Dsurf(p.grid, p.rs, th, ph);
  const double rcf = rc / rsca;
  if (rcf > p.grid.rg(p.grid.nrf))
    return v;

  v.inside = true;
  const auto st = cgr_4th_setup_from_spherical(p.grid, rcf, th, ph, p.grid.nrf);
  v.emdca = cgr_4th_apply(p.emd, st);
  v.psifca = cgr_4th_apply(p.psif, st);
  v.alphfca = cgr_4th_apply(p.alphf, st);
  v.bvxdfca = cgr_4th_apply(p.bvxdf, st);
  v.bvydfca = cgr_4th_apply(p.bvydf, st);
  v.bvzdfca = cgr_4th_apply(p.bvzdf, st);
  if (id_type == "IR" || id_type == "SP") {
    v.vepxfca = cgr_4th_apply(p.vepxf, st);
    v.vepyfca = cgr_4th_apply(p.vepyf, st);
    v.vepzfca = cgr_4th_apply(p.vepzf, st);
    if (id_type == "SP") {
      v.wxspfca = cgr_4th_apply(p.wxspf, st);
      v.wyspfca = cgr_4th_apply(p.wyspf, st);
      v.wzspfca = cgr_4th_apply(p.wzspf, st);
    }
  }
  return v;
}

BnsFluidPoint interpolate_bns_fluid_patch(const Patch &p, const string &id_type, double x, double y, double z,
                                          double xcoc, double ycoc, double zcoc, double ecc_cor_velx) {
  BnsFluidPoint point;
  const FluidValues fv = interpolate_fluid_patch(p, id_type, x, y, z);
  point.emdca = fv.emdca;
  if (!fv.inside)
    return point;

  const double psif4ca = std::pow(fv.psifca, 4);
  const double bxcor = fv.bvxdfca + p.ome * (-ycoc) + ecc_cor_velx * xcoc;
  const double bycor = fv.bvydfca + p.ome * xcoc + ecc_cor_velx * ycoc;
  const double bzcor = fv.bvzdfca + ecc_cor_velx * zcoc;

  if (id_type == "CO") {
    point.vxu = bxcor / fv.alphfca;
    point.vyu = bycor / fv.alphfca;
    point.vzu = bzcor / fv.alphfca;
  } else if (id_type == "IR") {
    const double lam = p.ber + bxcor * fv.vepxfca + bycor * fv.vepyfca + bzcor * fv.vepzfca;
    const double huta = lam / fv.alphfca;
    point.vxu = (fv.vepxfca / psif4ca) / huta;
    point.vyu = (fv.vepyfca / psif4ca) / huta;
    point.vzu = (fv.vepzfca / psif4ca) / huta;
  } else if (id_type == "SP") {
    const double psifcacp = std::pow(fv.psifca, p.confpow);
    const double lam = p.ber + bxcor * fv.vepxfca + bycor * fv.vepyfca + bzcor * fv.vepzfca;
    const double wdvep =
        (fv.wxspfca * fv.vepxfca + fv.wyspfca * fv.vepyfca + fv.wzspfca * fv.vepzfca) * psifcacp;
    const double w2 = psif4ca *
                      (fv.wxspfca * fv.wxspfca + fv.wyspfca * fv.wyspfca + fv.wzspfca * fv.wzspfca) *
                      psifcacp * psifcacp;
    const double huta = (lam + std::sqrt(lam * lam + 4.0 * fv.alphfca * fv.alphfca * (wdvep + w2))) /
                        (2.0 * fv.alphfca);
    point.vxu = (fv.vepxfca / psif4ca + psifcacp * fv.wxspfca) / huta;
    point.vyu = (fv.vepyfca / psif4ca + psifcacp * fv.wyspfca) / huta;
    point.vzu = (fv.vepzfca / psif4ca + psifcacp * fv.wzspfca) / huta;
  }
  return point;
}

} // namespace cocal_id_cpp
