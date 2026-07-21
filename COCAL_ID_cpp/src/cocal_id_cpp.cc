#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "cocal_id_cpp_internal.hh"

using namespace cocal_id_cpp;

namespace {
int gf_index(const cGH *cctkGH, int i, int j, int k) { return CCTK_GFINDEX3D(cctkGH, i, j, k); }
int vgf_index(const cGH *cctkGH, int i, int j, int k, int n) {
  return CCTK_VECTGFINDEX3D(cctkGH, i, j, k, n);
}

bool is_mrns_wl(const string &rnstype) { return rnstype == "MRNS_WL"; }

bool has_waveless_metric(const string &rnstype) {
  return rnstype == "RNS_WL" || rnstype == "MRNS_WL";
}

double interpolate_rns_field(const Array3D<double> &field, double x, double y, double z) {
  return cgr_4th_apply(field, gr2cgr_4th_setup(rns.grid, x, y, z));
}

void read_single_star_data(const string &dir, const string &rnstype, bool read_tabulated_eos,
                           bool read_magnetic_potential) {
  const bool mrns = is_mrns_wl(rnstype);

  rns = RnsData{};
  rns.grid = read_parameter_file(dir + "/rnspar.dat");
  if (read_tabulated_eos)
    teos_initialize(rns.eos, dir + "/teos_parameter.dat", dir + "/peos_parameter.dat");
  else
    peos_initialize(rns.eos, dir + "/peos_parameter.dat");
  setup_coordinates(rns.grid, false);
  allocate_rns_arrays(rns, mrns, read_magnetic_potential);

  read_3d_fields(dir + "/rnsgra_3D.las", {&rns.psi, &rns.alph, &rns.bvxd, &rns.bvyd, &rns.bvzd});
  read_rns_fluid(dir + "/rnsflu_3D.las", rns.emd, rns.rs, rns.omef, rns.ome, rns.ber, rns.radi);
  read_3d_fields(dir + "/rnsgra_Kij_3D.las", {&rns.kxxa, &rns.kxya, &rns.kxza, &rns.kyya, &rns.kyza, &rns.kzza});

  if (has_waveless_metric(rnstype)) {
    read_3d_fields(dir + "/rnsgra_hij_3D.las", {&rns.hxxd, &rns.hxyd, &rns.hxzd, &rns.hyyd, &rns.hyzd, &rns.hzzd});
  } else if (rnstype != "RNS_CF") {
    abort_cocal("Invalid COCAL_ID_rnstype");
  }

  if (mrns) {
    if (read_magnetic_potential)
      read_3d_fields(dir + "/rnsEMF_3D.las", {&rns.va, &rns.vaxd, &rns.vayd, &rns.vazd});
    else
      read_3d_fields(dir + "/rnsEMF_faraday_3D.las",
                     {&rns.fxd, &rns.fyd, &rns.fzd, &rns.fxyd, &rns.fxzd, &rns.fyzd});
    read_3d_fields(dir + "/rns4ve_3D.las", {&rns.utf, &rns.uxf, &rns.uyf, &rns.uzf});
  }

  invhij_WL_export(rns.grid, rns.hxxd, rns.hxyd, rns.hxzd, rns.hyyd, rns.hyzd, rns.hzzd, rns.hxxu,
                   rns.hxyu, rns.hxzu, rns.hyyu, rns.hyzu, rns.hzzu);
  index_vec_down2up_export(rns.grid, rns.hxxu, rns.hxyu, rns.hxzu, rns.hyyu, rns.hyzu, rns.hzzu,
                           rns.bvxu, rns.bvyu, rns.bvzu, rns.bvxd, rns.bvyd, rns.bvzd);
  interpo_gr2fl_metric_CF_export(rns.grid, rns.alph, rns.psi, rns.bvxu, rns.bvyu, rns.bvzu,
                                 rns.alphf, rns.psif, rns.bvxuf, rns.bvyuf, rns.bvzuf, rns.rs);
  rns.have_read_data = true;
}
} // namespace

extern "C" void COCAL_ID_cpp_read_rns_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  std::lock_guard<std::mutex> lock(rns_mutex);
  if (rns.have_read_data)
    return;

  read_single_star_data(COCAL_ID_PathToID, COCAL_ID_rnstype, COCAL_ID_read_tabulated_eos != 0,
                        COCAL_ID_read_magnetic_potential != 0);
}

extern "C" void COCAL_ID_cpp_deallocate_rns(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  std::lock_guard<std::mutex> lock(rns_mutex);
  rns = RnsData{};
}

extern "C" void COCAL_ID_cpp_rns(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!rns.have_read_data)
    abort_cocal("COCAL_ID_cpp_rns called before COCAL_ID_cpp_read_rns_data completed.");

  const string rnstype = COCAL_ID_rnstype;
  const bool mrns = rnstype == "MRNS_WL";
  const bool bool_metric_curv = CCTK_Equals(initial_data, "CocalRNS");
  const bool bool_lapse = CCTK_Equals(initial_lapse, "Cocal");
  const bool bool_shift = CCTK_Equals(initial_shift, "Cocal");
  const bool bool_hydro = CCTK_Equals(initial_hydro, "Cocal");
  const bool bool_Bvec = CCTK_Equals(initial_Bvec, "Cocal");
  const bool bool_Avec_Aphi = CCTK_Equals(initial_Avec, "Cocal") && CCTK_Equals(initial_Aphi, "Cocal");

  if (CCTK_Equals(initial_Avec, "Cocal") && !CCTK_Equals(initial_Aphi, "Cocal"))
    abort_cocal("Invalid initial data configuration: initial_Avec is set to Cocal but initial_Aphi is not.");
  if (!CCTK_Equals(initial_Avec, "Cocal") && CCTK_Equals(initial_Aphi, "Cocal"))
    abort_cocal("Invalid initial data configuration: initial_Aphi is set to Cocal but initial_Avec is not.");
  if (bool_Avec_Aphi && !COCAL_ID_read_magnetic_potential)
    abort_cocal("Invalid initial data configuration: Cocal initial_Avec/Aphi require COCAL_ID_read_magnetic_potential = yes.");
  if (bool_Bvec && COCAL_ID_read_magnetic_potential)
    abort_cocal("Invalid initial data configuration: Cocal initial_Bvec requires COCAL_ID_read_magnetic_potential = no.");

  const int nx = cctk_lsh[0];
  const int ny = cctk_lsh[1];
  const int nz = cctk_lsh[2];

  for (int kk = 0; kk < nz; ++kk) {
    for (int jj = 0; jj < ny; ++jj) {
      for (int ii = 0; ii < nx; ++ii) {
        const int idx = gf_index(cctkGH, ii, jj, kk);
        const double xcoc = x[idx] / rns.radi;
        const double ycoc = y[idx] / rns.radi;
        const double zcoc = z[idx] / rns.radi;

        RnsPoint point =
            interpolate_rns_point(rns, rnstype, COCAL_ID_read_magnetic_potential != 0, xcoc, ycoc, zcoc);
        double emdca = point.emdca;
        double hca = 0.0, preca = 0.0, rhoca = 0.0, eneca = 0.0;
        q2hprho(rns.eos, COCAL_ID_read_tabulated_eos != 0, emdca, hca, preca, rhoca, eneca);

        if (mrns) {
          if (COCAL_ID_read_magnetic_potential) {
            if (bool_Avec_Aphi) {
              double vaca = point.vaca;
              double vaxdca = point.vaxdca;
              double vaydca = point.vaydca;
              double vazdca = point.vazdca;
              if (COCAL_ID_offset_AvecAphi) {
                const double half_dx = 0.5 * CCTK_DELTA_SPACE(0) / rns.radi;
                const double half_dy = 0.5 * CCTK_DELTA_SPACE(1) / rns.radi;
                const double half_dz = 0.5 * CCTK_DELTA_SPACE(2) / rns.radi;
                vaxdca = interpolate_rns_field(rns.vaxd, xcoc, ycoc + half_dy,
                                               zcoc + half_dz);
                vaydca = interpolate_rns_field(rns.vayd, xcoc + half_dx, ycoc,
                                               zcoc + half_dz);
                vazdca = interpolate_rns_field(rns.vazd, xcoc + half_dx,
                                               ycoc + half_dy, zcoc);
                vaca = interpolate_rns_field(rns.va, xcoc + half_dx,
                                             ycoc + half_dy, zcoc + half_dz);
              }
              Aphi[idx] = vaca;
              Avec[vgf_index(cctkGH, ii, jj, kk, 0)] = vaxdca;
              Avec[vgf_index(cctkGH, ii, jj, kk, 1)] = vaydca;
              Avec[vgf_index(cctkGH, ii, jj, kk, 2)] = vazdca;
            }
          } else if (bool_Bvec) {
            Bvec[vgf_index(cctkGH, ii, jj, kk, 0)] = point.fyzdca / rns.radi;
            Bvec[vgf_index(cctkGH, ii, jj, kk, 1)] = -point.fxzdca / rns.radi;
            Bvec[vgf_index(cctkGH, ii, jj, kk, 2)] = point.fxydca / rns.radi;
          }
        }

        if (bool_metric_curv) {
          gxx[idx] = point.gxx1;
          gxy[idx] = point.gxy1;
          gxz[idx] = point.gxz1;
          gyy[idx] = point.gyy1;
          gyz[idx] = point.gyz1;
          gzz[idx] = point.gzz1;
          kxx[idx] = point.kxx1;
          kxy[idx] = point.kxy1;
          kxz[idx] = point.kxz1;
          kyy[idx] = point.kyy1;
          kyz[idx] = point.kyz1;
          kzz[idx] = point.kzz1;
          if (point.gxx1 != point.gxx1)
            abort_cocal("NaN in gxx");
        }
        if (bool_lapse)
          alp[idx] = point.alphca;
        if (bool_shift) {
          betax[idx] = point.bvxdca;
          betay[idx] = point.bvydca;
          betaz[idx] = point.bvzdca;
        }
        if (bool_hydro) {
          rho[idx] = rhoca;
          press[idx] = preca;
          eps[idx] = eneca / rhoca - 1.0;
          vel[vgf_index(cctkGH, ii, jj, kk, 0)] = point.vxu;
          vel[vgf_index(cctkGH, ii, jj, kk, 1)] = point.vyu;
          vel[vgf_index(cctkGH, ii, jj, kk, 2)] = point.vzu;
        }
      }
    }
  }
}

extern "C" void COCAL_ID_cpp_read_bns_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  std::lock_guard<std::mutex> lock(bns_mutex);
  if (bns.have_read_data)
    return;

  bns = BnsData{};
  bns.path = COCAL_ID_PathToID;

  {
    TextFile f(bns.path + "/rnspar_mpt1.dat");
    (void)f.line_tokens();
    (void)f.line_tokens();
    const auto t = f.line_tokens(3, "BNS ID type");
    bns.id_type = t.at(2);
  }

  for (int impt = 1; impt <= nmpt; ++impt) {
    Patch &p = bns.patches.at(static_cast<size_t>(impt - 1));
    p.grid = read_parameter_file(bns.path + "/rnspar_mpt" + std::to_string(impt) + ".dat");
    read_surf_parameter_file(p.grid, bns.path, impt);
    p.ex = read_binary_excision_file(bns.path, impt);
    if (COCAL_ID_read_tabulated_eos) {
      if (impt <= 2)
        teos_initialize(p.eos, bns.path + "/teos_parameter_mpt" + std::to_string(impt) + ".dat",
                        bns.path + "/peos_parameter_mpt" + std::to_string(impt) + ".dat");
      else
        p.eos = bns.patches[1].eos;
    } else {
      if (impt <= 2)
        peos_initialize(p.eos, bns.path + "/peos_parameter_mpt" + std::to_string(impt) + ".dat");
      else
        p.eos = bns.patches[1].eos;
    }
    setup_coordinates(p.grid, true);
    calc_parameter_binary_excision(p);
  }
  bns.active_eos = bns.patches[1].eos;

  for (int impt = 1; impt <= nmpt; ++impt) {
    Patch &p = bns.patches.at(static_cast<size_t>(impt - 1));
    const bool star_patch = impt <= 2;
    allocate_patch_arrays(p, bns.id_type, star_patch);

    read_3d_fields(bns.path + "/bnsgra_3D_mpt" + std::to_string(impt) + ".las",
                   {&p.psi, &p.alph, &p.bvxd, &p.bvyd, &p.bvzd});
    read_3d_fields(bns.path + "/bnsgra_Kij_3D_mpt" + std::to_string(impt) + ".las",
                   {&p.axx, &p.axy, &p.axz, &p.ayy, &p.ayz, &p.azz});
    if (!star_patch)
      continue;

    if (bns.id_type == "CO") {
      read_bns_fluid_co(bns.path + "/bnsflu_3D_mpt" + std::to_string(impt) + ".las", p.emd, p.ome,
                        p.ber, p.radi);
    } else if (bns.id_type == "IR") {
      read_bns_fluid_ir(bns.path + "/bnsflu_3D_mpt" + std::to_string(impt) + ".las", p.emd, p.vep,
                        p.ome, p.ber, p.radi);
    } else if (bns.id_type == "SP") {
      read_bns_fluid_sp(bns.path + "/bnsflu_3D_mpt" + std::to_string(impt) + ".las", p.emd, p.vep,
                        p.wxspf, p.wyspf, p.wzspf, p.ome, p.ber, p.radi, p.confpow);
    }
    p.r_surf = p.grid.r_surf;
    read_2d_surface(bns.path + "/bnssur_3D_mpt" + std::to_string(impt) + ".las", p.rs);
    interpo_gr2fl_metric_CF_export(p.grid, p.alph, p.psi, p.bvxd, p.bvyd, p.bvzd, p.alphf, p.psif,
                                   p.bvxdf, p.bvydf, p.bvzdf, p.rs);
    if (bns.id_type == "IR" || bns.id_type == "SP")
      read_3d_fields(bns.path + "/bnsdvep_3D_mpt" + std::to_string(impt) + ".las",
                     {&p.vepxf, &p.vepyf, &p.vepzf});

    if (impt == 2) {
      if (bns.id_type == "IR" || bns.id_type == "SP") {
        sign_flip(p.vepxf);
        sign_flip(p.vepyf);
        if (bns.id_type == "SP") {
          sign_flip(p.wxspf);
          sign_flip(p.wyspf);
        }
      }
      sign_flip(p.bvxdf);
      sign_flip(p.bvydf);
      sign_flip(p.bvxd);
      sign_flip(p.bvyd);
      sign_flip(p.axz);
      sign_flip(p.ayz);
    }
  }

  bns.rr3 = bns_outer_patch_fraction * (bns.patches[0].grid.rgout - bns.patches[0].grid.rgmid);
  bns.dis_cm = bns.patches[0].ex.dis;
  bns.have_read_data = true;
}

extern "C" void COCAL_ID_cpp_deallocate_bns(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  std::lock_guard<std::mutex> lock(bns_mutex);
  bns = BnsData{};
}

extern "C" void COCAL_ID_cpp_bns(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!bns.have_read_data)
    abort_cocal("COCAL_ID_cpp_bns called before COCAL_ID_cpp_read_bns_data completed.");

  const bool bool_metric_curv = CCTK_Equals(initial_data, "CocalBNS");
  const bool bool_lapse = CCTK_Equals(initial_lapse, "Cocal");
  const bool bool_shift = CCTK_Equals(initial_shift, "Cocal");
  const bool bool_hydro = CCTK_Equals(initial_hydro, "Cocal");
  const int nx = cctk_lsh[0];
  const int ny = cctk_lsh[1];
  const int nz = cctk_lsh[2];

  const Patch &p1 = bns.patches[0];
  const Patch &p2 = bns.patches[1];
  const Patch &p3 = bns.patches[2];

  for (int kk = 0; kk < nz; ++kk) {
    for (int jj = 0; jj < ny; ++jj) {
      for (int ii = 0; ii < nx; ++ii) {
        const int idx = gf_index(cctkGH, ii, jj, kk);
        const double xcoc = x[idx] / p1.radi;
        const double ycoc = y[idx] / p1.radi;
        const double zcoc = z[idx] / p1.radi;

        MetricValues mv;
        const double rrcm = std::sqrt(xcoc * xcoc + ycoc * ycoc + zcoc * zcoc);
        double vxu = 0.0, vyu = 0.0, vzu = 0.0;
        double emdca = 0.0;

        if (rrcm > bns.rr3) {
          mv = interpolate_metric_patch(p3, xcoc, ycoc, zcoc);
        } else if (xcoc <= 0.0) {
          const double xp = xcoc + bns.dis_cm;
          mv = interpolate_metric_patch(p1, xp, ycoc, zcoc);
          const BnsFluidPoint fluid =
              interpolate_bns_fluid_patch(p1, bns.id_type, xp, ycoc, zcoc, xcoc, ycoc, zcoc,
                                          COCAL_ID_ecc_cor_velx);
          emdca = fluid.emdca;
          vxu = fluid.vxu;
          vyu = fluid.vyu;
          vzu = fluid.vzu;
        } else {
          const double xp = -(xcoc - bns.dis_cm);
          const double yp = -ycoc;
          mv = interpolate_metric_patch(p2, xp, yp, zcoc);
          const BnsFluidPoint fluid =
              interpolate_bns_fluid_patch(p2, bns.id_type, xp, yp, zcoc, xcoc, ycoc, zcoc,
                                          COCAL_ID_ecc_cor_velx);
          emdca = fluid.emdca;
          vxu = fluid.vxu;
          vyu = fluid.vyu;
          vzu = fluid.vzu;
        }

        const double psi4ca = std::pow(mv.psica, 4);
        const double gxx1 = psi4ca;
        const double gyy1 = psi4ca;
        const double gzz1 = psi4ca;
        const double kxx1 = psi4ca * mv.axx / p1.radi;
        const double kxy1 = psi4ca * mv.axy / p1.radi;
        const double kxz1 = psi4ca * mv.axz / p1.radi;
        const double kyy1 = psi4ca * mv.ayy / p1.radi;
        const double kyz1 = psi4ca * mv.ayz / p1.radi;
        const double kzz1 = psi4ca * mv.azz / p1.radi;

        double hca = 0.0, preca = 0.0, rhoca = 0.0, eneca = 0.0;
        q2hprho(bns.active_eos, COCAL_ID_read_tabulated_eos != 0, emdca, hca, preca, rhoca, eneca);

        if (bool_metric_curv) {
          gxx[idx] = gxx1;
          gxy[idx] = 0.0;
          gxz[idx] = 0.0;
          gyy[idx] = gyy1;
          gyz[idx] = 0.0;
          gzz[idx] = gzz1;
          kxx[idx] = kxx1;
          kxy[idx] = kxy1;
          kxz[idx] = kxz1;
          kyy[idx] = kyy1;
          kyz[idx] = kyz1;
          kzz[idx] = kzz1;
          if (gxx1 != gxx1)
            abort_cocal("NaN in gxx");
        }
        if (bool_lapse)
          alp[idx] = mv.alphca;
        if (bool_shift) {
          betax[idx] = mv.bvxdca;
          betay[idx] = mv.bvydca;
          betaz[idx] = mv.bvzdca;
        }
        if (bool_hydro) {
          rho[idx] = rhoca;
          press[idx] = preca;
          eps[idx] = eneca / rhoca - 1.0;
          vel[vgf_index(cctkGH, ii, jj, kk, 0)] = vxu;
          vel[vgf_index(cctkGH, ii, jj, kk, 1)] = vyu;
          vel[vgf_index(cctkGH, ii, jj, kk, 2)] = vzu;
        }
      }
    }
  }
}
