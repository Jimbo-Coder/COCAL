#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "cocal_id_cpp_internal.hh"

// This is the C++ counterpart of COCAL_ID's BHT reader and fill path. The raw
// field conversion follows PCOCAL's coc2pri_bht_wl; the automatic center
// estimate, trace handling, and Hermite fill follow COCAL_ID's
// cocal_id_bht_fill.f90 and have no separate upstream PCOCAL routine.
namespace cocal_id_cpp {
namespace {

using Symmetric3 = array<double, 6>;

BhtParameters read_bht_parameters(const string &path) {
  TextFile file(path);
  BhtParameters p;
  auto row = file.line_tokens(2, "BHT radius and mass");
  p.r_surf = parse_real(row.at(0), path);
  p.mass_bh = parse_real(row.at(1), path);
  row = file.line_tokens(2, "BHT angular velocity and spin");
  p.ome_bh = parse_real(row.at(0), path);
  p.spin_bh = parse_real(row.at(1), path);
  row = file.line_tokens(2, "BHT lapse and conformal factor");
  p.alph_bh = parse_real(row.at(0), path);
  p.psi_bh = parse_real(row.at(1), path);
  row = file.line_tokens(2, "BHT spin direction");
  p.th_spin_bh_deg = parse_real(row.at(0), path);
  p.phi_spin_bh_deg = parse_real(row.at(1), path);
  row = file.line_tokens(3, "BHT boundary and solution types");
  p.boundary_type = row.at(0);
  p.solution_type = row.at(1);
  p.mass_ratio = parse_real(row.at(2), path);
  return p;
}

double read_idpar_real(TextFile &file, const string &context) {
  return parse_real(file.line_tokens(1, context).at(0), file.path);
}

void read_idpar_center(BhtData &d, const string &path) {
  TextFile file(path);
  d.center_x = read_idpar_real(file, "black-hole center x");
  d.center_y = read_idpar_real(file, "black-hole center y");
  d.center_z = read_idpar_real(file, "black-hole center z");
  d.r_ah = read_idpar_real(file, "fill radius");
  (void)file.line_tokens(1, "legacy extra-point count");
  d.extra_points_dr = read_idpar_real(file, "extra-point spacing");

  d.idpar_center.fill(0.0);
  d.idpar_center[BHT_ALP] = read_idpar_real(file, "central lapse");
  (void)read_idpar_real(file, "legacy central conformal factor");
  d.idpar_center[BHT_BETAX] = read_idpar_real(file, "central shift x");
  d.idpar_center[BHT_BETAY] = read_idpar_real(file, "central shift y");
  d.idpar_center[BHT_BETAZ] = read_idpar_real(file, "central shift z");
  d.idpar_center[BHT_GXX] = read_idpar_real(file, "central gxx");
  d.idpar_center[BHT_GXY] = read_idpar_real(file, "central gxy");
  d.idpar_center[BHT_GXZ] = read_idpar_real(file, "central gxz");
  d.idpar_center[BHT_GYY] = read_idpar_real(file, "central gyy");
  d.idpar_center[BHT_GYZ] = read_idpar_real(file, "central gyz");
  d.idpar_center[BHT_GZZ] = read_idpar_real(file, "central gzz");
  d.idpar_center[BHT_KXX] = read_idpar_real(file, "central Kxx");
  d.idpar_center[BHT_KXY] = read_idpar_real(file, "central Kxy");
  d.idpar_center[BHT_KXZ] = read_idpar_real(file, "central Kxz");
  d.idpar_center[BHT_KYZ] = read_idpar_real(file, "central Kyz");
  d.idpar_center[BHT_KZZ] = read_idpar_real(file, "central Kzz");
  d.idpar_center[BHT_RHO] = read_idpar_real(file, "central density");
  d.idpar_center[BHT_VX] = read_idpar_real(file, "central velocity x");
  d.idpar_center[BHT_VY] = read_idpar_real(file, "central velocity y");
  d.idpar_center[BHT_VZ] = read_idpar_real(file, "central velocity z");
  d.trace_width = read_idpar_real(file, "trace damping width");
  d.trace_power = read_idpar_real(file, "trace damping power");
}

void set_smooth_scales(BhtData &d) {
  double local_dr = 0.0;
  for (int ir = 0; ir < d.grid.nrg; ++ir) {
    if (d.grid.rg(ir) <= d.r_ah && d.grid.rg(ir + 1) > d.r_ah) {
      local_dr = d.grid.rg(ir + 1) - d.grid.rg(ir);
      break;
    }
  }
  const double safe_radius = std::max(d.r_ah, 1.0e-12);
  if (local_dr <= 0.0)
    local_dr = 0.05 * safe_radius;
  d.extra_points_dr = std::clamp(std::max(2.0 * local_dr, 0.20 * safe_radius),
                                 0.0, 0.25 * safe_radius);

  const double gap = std::max(0.0, d.r_ah - d.rexc);
  d.trace_width = gap > 0.0 ? std::clamp(1.25 * gap / safe_radius, 0.04, 0.16) : 0.08;
  d.trace_power = 2.0;
}

double metric_determinant(const BhtPoint &q) {
  return -q[BHT_GXZ] * q[BHT_GXZ] * q[BHT_GYY] +
         2.0 * q[BHT_GXY] * q[BHT_GXZ] * q[BHT_GYZ] -
         q[BHT_GXX] * q[BHT_GYZ] * q[BHT_GYZ] -
         q[BHT_GXY] * q[BHT_GXY] * q[BHT_GZZ] +
         q[BHT_GXX] * q[BHT_GYY] * q[BHT_GZZ];
}

bool metric_inverse(const BhtPoint &q, Symmetric3 &inverse) {
  const double det = metric_determinant(q);
  if (!(det > 0.0) || !std::isfinite(det))
    return false;
  inverse[0] = (-q[BHT_GYZ] * q[BHT_GYZ] + q[BHT_GYY] * q[BHT_GZZ]) / det;
  inverse[1] = (q[BHT_GXZ] * q[BHT_GYZ] - q[BHT_GXY] * q[BHT_GZZ]) / det;
  inverse[2] = (-q[BHT_GXZ] * q[BHT_GYY] + q[BHT_GXY] * q[BHT_GYZ]) / det;
  inverse[3] = (-q[BHT_GXZ] * q[BHT_GXZ] + q[BHT_GXX] * q[BHT_GZZ]) / det;
  inverse[4] = (q[BHT_GXY] * q[BHT_GXZ] - q[BHT_GXX] * q[BHT_GYZ]) / det;
  inverse[5] = (-q[BHT_GXY] * q[BHT_GXY] + q[BHT_GXX] * q[BHT_GYY]) / det;
  return std::all_of(inverse.begin(), inverse.end(), [](double value) { return std::isfinite(value); });
}

void zero_hydro(BhtPoint &q) {
  for (int variable = BHT_RHO; variable <= BHT_VZ; ++variable)
    q[static_cast<size_t>(variable)] = 0.0;
}

bool validate_state(const BhtPoint &q) {
  return std::all_of(q.begin(), q.end(), [](double value) { return std::isfinite(value); }) &&
         q[BHT_ALP] > 0.0 && metric_determinant(q) > 0.0;
}

bool split_curvature(const BhtPoint &q, Symmetric3 &tracefree, double &trace) {
  Symmetric3 inverse{};
  if (!metric_inverse(q, inverse))
    return false;
  trace = inverse[0] * q[BHT_KXX] + inverse[3] * q[BHT_KYY] + inverse[5] * q[BHT_KZZ] +
          2.0 * (inverse[1] * q[BHT_KXY] + inverse[2] * q[BHT_KXZ] +
                 inverse[4] * q[BHT_KYZ]);
  tracefree = {q[BHT_KXX] - q[BHT_GXX] * trace / 3.0,
               q[BHT_KXY] - q[BHT_GXY] * trace / 3.0,
               q[BHT_KXZ] - q[BHT_GXZ] * trace / 3.0,
               q[BHT_KYY] - q[BHT_GYY] * trace / 3.0,
               q[BHT_KYZ] - q[BHT_GYZ] * trace / 3.0,
               q[BHT_KZZ] - q[BHT_GZZ] * trace / 3.0};
  return std::isfinite(trace);
}

void store_tracefree(BhtPoint &q, const Symmetric3 &tracefree) {
  q[BHT_KXX] = tracefree[0];
  q[BHT_KXY] = tracefree[1];
  q[BHT_KXZ] = tracefree[2];
  q[BHT_KYY] = tracefree[3];
  q[BHT_KYZ] = tracefree[4];
  q[BHT_KZZ] = tracefree[5];
}

bool project_tracefree(BhtPoint &q) {
  Symmetric3 inverse{};
  if (!metric_inverse(q, inverse) || std::abs(inverse[3]) <= 1.0e-14)
    return false;
  q[BHT_KYY] = -(inverse[0] * q[BHT_KXX] + inverse[5] * q[BHT_KZZ] +
                   2.0 * (inverse[1] * q[BHT_KXY] + inverse[2] * q[BHT_KXZ] +
                          inverse[4] * q[BHT_KYZ])) /
                 inverse[3];
  return std::isfinite(q[BHT_KYY]);
}

void add_trace(BhtPoint &q, double trace) {
  q[BHT_KXX] += q[BHT_GXX] * trace / 3.0;
  q[BHT_KXY] += q[BHT_GXY] * trace / 3.0;
  q[BHT_KXZ] += q[BHT_GXZ] * trace / 3.0;
  q[BHT_KYY] += q[BHT_GYY] * trace / 3.0;
  q[BHT_KYZ] += q[BHT_GYZ] * trace / 3.0;
  q[BHT_KZZ] += q[BHT_GZZ] * trace / 3.0;
}

bool evaluate_raw(BhtData &d, double x, double y, double z, bool allow_inside, BhtPoint &q) {
  q.fill(0.0);
  const double xcoc = x / d.radi;
  const double ycoc = y / d.radi;
  const double zcoc = z / d.radi;
  const double radius = std::sqrt(xcoc * xcoc + ycoc * ycoc + zcoc * zcoc);
  if (!allow_inside && radius < d.rexc)
    return false;

  const Stencil4 stencil = gr2cgr_4th_setup(d.grid, xcoc, ycoc, zcoc);
  const double psi = cgr_4th_apply(d.psi, stencil);
  const double psi4 = std::pow(psi, 4);
  const double hxx = cgr_4th_apply(d.hxxd, stencil);
  const double hxy = cgr_4th_apply(d.hxyd, stencil);
  const double hxz = cgr_4th_apply(d.hxzd, stencil);
  const double hyy = cgr_4th_apply(d.hyyd, stencil);
  const double hyz = cgr_4th_apply(d.hyzd, stencil);
  const double hzz = cgr_4th_apply(d.hzzd, stencil);

  q[BHT_GXX] = psi4 * (1.0 + hxx);
  q[BHT_GXY] = psi4 * hxy;
  q[BHT_GXZ] = psi4 * hxz;
  q[BHT_GYY] = psi4 * (1.0 + hyy);
  q[BHT_GYZ] = psi4 * hyz;
  q[BHT_GZZ] = psi4 * (1.0 + hzz);
  q[BHT_KXX] = psi4 * cgr_4th_apply(d.kxxa, stencil) / d.radi;
  q[BHT_KXY] = psi4 * cgr_4th_apply(d.kxya, stencil) / d.radi;
  q[BHT_KXZ] = psi4 * cgr_4th_apply(d.kxza, stencil) / d.radi;
  q[BHT_KYY] = psi4 * cgr_4th_apply(d.kyya, stencil) / d.radi;
  q[BHT_KYZ] = psi4 * cgr_4th_apply(d.kyza, stencil) / d.radi;
  q[BHT_KZZ] = psi4 * cgr_4th_apply(d.kzza, stencil) / d.radi;

  q[BHT_ALP] = cgr_4th_apply(d.alph, stencil);
  q[BHT_BETAX] = psi4 * cgr_4th_apply(d.bvxd, stencil);
  q[BHT_BETAY] = psi4 * cgr_4th_apply(d.bvyd, stencil);
  q[BHT_BETAZ] = psi4 * cgr_4th_apply(d.bvzd, stencil);

  double enthalpy_variable = cgr_4th_apply(d.emd, stencil);
  if (std::abs(enthalpy_variable) > cocal_zero_tolerance) {
    const double alpha = q[BHT_ALP];
    if (std::abs(alpha) <= cocal_zero_tolerance)
      return false;
    const double rotation = cgr_4th_apply(d.omeg, stencil);
    q[BHT_VX] = (cgr_4th_apply(d.bvxu, stencil) - rotation * ycoc) / alpha;
    q[BHT_VY] = (cgr_4th_apply(d.bvyu, stencil) + rotation * xcoc) / alpha;
    q[BHT_VZ] = cgr_4th_apply(d.bvzu, stencil) / alpha;
    double h = 0.0, pressure = 0.0, density = 0.0, energy_density = 0.0;
    q2hprho(d.eos, d.tabulated_eos, enthalpy_variable, h, pressure, density, energy_density);
    q[BHT_RHO] = density;
    q[BHT_PRESS] = pressure;
    if (density > 0.0)
      q[BHT_EPS] = energy_density / density - 1.0;
  }
  return true;
}

bool compute_smooth_center(BhtData &d) {
  constexpr array<array<double, 3>, 14> directions{{
      {{1.0, 0.0, 0.0}},   {{-1.0, 0.0, 0.0}},  {{0.0, 1.0, 0.0}},
      {{0.0, -1.0, 0.0}},  {{0.0, 0.0, 1.0}},   {{0.0, 0.0, -1.0}},
      {{1.0, 1.0, 1.0}},   {{1.0, 1.0, -1.0}},  {{1.0, -1.0, 1.0}},
      {{1.0, -1.0, -1.0}}, {{-1.0, 1.0, 1.0}},  {{-1.0, 1.0, -1.0}},
      {{-1.0, -1.0, 1.0}}, {{-1.0, -1.0, -1.0}},
  }};

  d.smooth_center.fill(0.0);
  int valid_directions = 0;
  for (const auto &raw_direction : directions) {
    const double norm = std::sqrt(raw_direction[0] * raw_direction[0] +
                                  raw_direction[1] * raw_direction[1] +
                                  raw_direction[2] * raw_direction[2]);
    const array<double, 3> direction{raw_direction[0] / norm, raw_direction[1] / norm,
                                     raw_direction[2] / norm};
    array<double, bht_fill_points> radii_squared{};
    array<BhtPoint, bht_fill_points> samples{};
    bool direction_valid = true;
    for (int point = 0; point < bht_fill_points; ++point) {
      const double radius = d.r_ah + d.extra_points_dr * (static_cast<double>(point) + 1.0e-8);
      radii_squared[point] = radius * radius;
      if (!evaluate_raw(d, d.center_x + d.radi * radius * direction[0],
                        d.center_y + d.radi * radius * direction[1],
                        d.center_z + d.radi * radius * direction[2], true, samples[point])) {
        direction_valid = false;
        break;
      }
      Symmetric3 tracefree{};
      double trace = 0.0;
      if (!split_curvature(samples[point], tracefree, trace)) {
        direction_valid = false;
        break;
      }
      store_tracefree(samples[point], tracefree);
    }
    if (!direction_valid)
      continue;

    array<double, bht_fill_points> weights{};
    lagint_4th_weights(radii_squared, 0.0, weights);
    for (int variable = 0; variable < bht_nvars; ++variable) {
      array<double, bht_fill_points> values{};
      for (int point = 0; point < bht_fill_points; ++point)
        values[point] = samples[point][variable];
      d.smooth_center[variable] += lagint_4th_apply(weights, values);
    }
    ++valid_directions;
  }

  d.smooth_center_attempted = true;
  d.smooth_center_valid = valid_directions > 0;
  if (!d.smooth_center_valid)
    return false;
  for (double &value : d.smooth_center)
    value /= static_cast<double>(valid_directions);
  zero_hydro(d.smooth_center);
  return true;
}

bool fill_point(BhtData &d, double x, double y, double z, bool use_idpar, BhtPoint &q) {
  const double xrel = x - d.center_x;
  const double yrel = y - d.center_y;
  const double zrel = z - d.center_z;
  const double physical_radius = std::sqrt(xrel * xrel + yrel * yrel + zrel * zrel);
  const double radius = physical_radius / d.radi;
  if (radius >= d.r_ah)
    return evaluate_raw(d, x, y, z, false, q) && validate_state(q);

  array<double, 3> direction{1.0, 0.0, 0.0};
  if (physical_radius > 1.0e-14)
    direction = {xrel / physical_radius, yrel / physical_radius, zrel / physical_radius};

  double boundary_trace = 0.0;
  if (use_idpar) {
    array<double, bht_fill_points> radii{};
    array<BhtPoint, bht_fill_points> samples{};
    radii[0] = 0.0;
    samples[0] = d.idpar_center;
    for (int point = 1; point < bht_fill_points; ++point) {
      const double sample_radius =
          d.r_ah + d.extra_points_dr * (static_cast<double>(point - 1) + 1.0e-8);
      radii[point] = sample_radius;
      if (!evaluate_raw(d, d.center_x + d.radi * sample_radius * direction[0],
                        d.center_y + d.radi * sample_radius * direction[1],
                        d.center_z + d.radi * sample_radius * direction[2], true,
                        samples[point]))
        return false;
      Symmetric3 tracefree{};
      double trace = 0.0;
      if (!split_curvature(samples[point], tracefree, trace))
        return false;
      if (point == 1)
        boundary_trace = trace;
      store_tracefree(samples[point], tracefree);
    }
    array<double, bht_fill_points> weights{};
    lagint_4th_weights(radii, radius, weights);
    for (int variable = 0; variable < bht_nvars; ++variable) {
      array<double, bht_fill_points> values{};
      for (int point = 0; point < bht_fill_points; ++point)
        values[point] = samples[point][variable];
      q[variable] = lagint_4th_apply(weights, values);
    }
  } else {
    if (!d.smooth_center_attempted && !compute_smooth_center(d))
      return false;
    if (!d.smooth_center_valid)
      return false;

    array<BhtPoint, 2> samples{};
    for (int point = 0; point < 2; ++point) {
      const double sample_radius =
          d.r_ah + d.extra_points_dr * (static_cast<double>(point) + 1.0e-8);
      if (!evaluate_raw(d, d.center_x + d.radi * sample_radius * direction[0],
                        d.center_y + d.radi * sample_radius * direction[1],
                        d.center_z + d.radi * sample_radius * direction[2], true,
                        samples[point]))
        return false;
      Symmetric3 tracefree{};
      double trace = 0.0;
      if (!split_curvature(samples[point], tracefree, trace))
        return false;
      if (point == 0)
        boundary_trace = trace;
      store_tracefree(samples[point], tracefree);
    }

    const double s = std::clamp(radius / d.r_ah, 0.0, 1.0);
    const double h00 = 2.0 * s * s * s - 3.0 * s * s + 1.0;
    const double h01 = -2.0 * s * s * s + 3.0 * s * s;
    const double h11 = s * s * s - s * s;
    for (int variable = 0; variable < bht_nvars; ++variable) {
      const double derivative = (samples[1][variable] - samples[0][variable]) /
                                d.extra_points_dr;
      q[variable] = h00 * d.smooth_center[variable] + h01 * samples[0][variable] +
                    h11 * d.r_ah * derivative;
    }
  }

  if (!project_tracefree(q))
    return false;
  double interior_trace = 0.0;
  if (d.trace_width > 0.0 && d.r_ah > 0.0) {
    const double inward_distance = std::max(0.0, d.r_ah - radius);
    interior_trace = boundary_trace *
                     std::exp(-std::pow(inward_distance / (d.trace_width * d.r_ah),
                                        d.trace_power));
  }
  add_trace(q, interior_trace);
  zero_hydro(q);
  return validate_state(q);
}

void initialize_bht(BhtData &d, const string &path, bool tabulated_eos,
                    bool use_idpar) {
  d = BhtData{};
  d.grid = read_parameter_file(path + "/rnspar.dat");
  d.parameters = read_bht_parameters(path + "/bhtpar.dat");
  d.grid.r_surf = d.parameters.r_surf;
  if (d.parameters.mass_bh <= 0.0)
    abort_cocal("Invalid BHT mass in bhtpar.dat.");
  const double spin_ratio = std::abs(d.parameters.spin_bh / d.parameters.mass_bh);
  if (spin_ratio > 1.0)
    abort_cocal("Invalid BHT spin in bhtpar.dat.");
  const double sqrt_one_minus_spin2 = std::sqrt(std::max(1.0 - spin_ratio * spin_ratio, 0.0));
  d.grid.rgin = d.parameters.mass_bh * (1.0 + 0.8 * sqrt_one_minus_spin2);

  d.tabulated_eos = tabulated_eos;
  if (tabulated_eos)
    teos_initialize(d.eos, path + "/teos_parameter.dat", path + "/peos_parameter.dat");
  else
    peos_initialize(d.eos, path + "/peos_parameter.dat");

  grid_r_bht(d.grid, d.parameters.mass_bh, d.parameters.spin_bh);
  grid_theta(d.grid);
  grid_phi(d.grid);
  grid_extended(d.grid);

  d.rexc = d.grid.rg(0);
  d.r_eh = d.parameters.mass_bh +
           std::sqrt(d.parameters.mass_bh * d.parameters.mass_bh -
                     d.parameters.spin_bh * d.parameters.spin_bh);
  d.r_ring = std::abs(d.parameters.spin_bh);
  d.r_ah = 0.99 * d.r_eh;
  if (use_idpar)
    read_idpar_center(d, path + "/read_bin_AH.par");
  else
    set_smooth_scales(d);

  allocate_bht_arrays(d);
  read_3d_fields(path + "/rnsgra_3D.las",
                 {&d.psi, &d.alph, &d.bvxd, &d.bvyd, &d.bvzd});
  read_3d_fields(path + "/rnsgra_hij_3D.las",
                 {&d.hxxd, &d.hxyd, &d.hxzd, &d.hyyd, &d.hyzd, &d.hzzd});
  read_3d_fields(path + "/rnsgra_Kij_3D.las",
                 {&d.kxxa, &d.kxya, &d.kxza, &d.kyya, &d.kyza, &d.kzza});
  read_bht_fluid(path + "/rnsflu_3D.las", d.emd, d.omeg, d.ome, d.ber, d.radi);
  invhij_WL_export(d.grid, d.hxxd, d.hxyd, d.hxzd, d.hyyd, d.hyzd, d.hzzd,
                   d.hxxu, d.hxyu, d.hxzu, d.hyyu, d.hyzu, d.hzzu);
  index_vec_down2up_export(d.grid, d.hxxu, d.hxyu, d.hxzu, d.hyyu, d.hyzu, d.hzzu,
                           d.bvxu, d.bvyu, d.bvzu, d.bvxd, d.bvyd, d.bvzd);
  d.have_read_data = true;
}

int gf_index(const cGH *cctkGH, int i, int j, int k) {
  return CCTK_GFINDEX3D(cctkGH, i, j, k);
}

int vector_gf_index(const cGH *cctkGH, int i, int j, int k, int component) {
  return CCTK_VECTGFINDEX3D(cctkGH, i, j, k, component);
}

} // namespace
} // namespace cocal_id_cpp

using namespace cocal_id_cpp;

extern "C" void COCAL_ID_cpp_read_bht_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  std::lock_guard<std::mutex> lock(bht_mutex);
  if (bht.have_read_data)
    return;
  const bool use_idpar = CCTK_Equals(COCAL_ID_bht_fill_method, "IDpar");
  initialize_bht(bht, COCAL_ID_PathToID, COCAL_ID_read_tabulated_eos != 0, use_idpar);
  if (COCAL_ID_verbose) {
    std::ostringstream message;
    message << "BHT radii rexc, r_fill, r_EH, r_ring: " << bht.rexc << " " << bht.r_ah
            << " " << bht.r_eh << " " << bht.r_ring;
    CCTK_INFO(message.str().c_str());
  }
}

extern "C" void COCAL_ID_cpp_deallocate_bht(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  std::lock_guard<std::mutex> lock(bht_mutex);
  bht = BhtData{};
}

extern "C" void COCAL_ID_cpp_bht(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (!bht.have_read_data)
    abort_cocal("COCAL_ID_cpp_bht called before COCAL_ID_cpp_read_bht_data completed.");
  if (CCTK_Equals(initial_Bvec, "Cocal") || CCTK_Equals(initial_Avec, "Cocal") ||
      CCTK_Equals(initial_Aphi, "Cocal"))
    abort_cocal("COCAL BHT is unmagnetized; select zero magnetic initial data.");

  const bool write_metric = CCTK_Equals(initial_data, "CocalBHT");
  const bool write_lapse = CCTK_Equals(initial_lapse, "Cocal");
  const bool write_shift = CCTK_Equals(initial_shift, "Cocal");
  const bool write_hydro = CCTK_Equals(initial_hydro, "Cocal");
  const bool use_idpar = CCTK_Equals(COCAL_ID_bht_fill_method, "IDpar");

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int index = gf_index(cctkGH, i, j, k);
        BhtPoint point{};
        if (!fill_point(bht, x[index], y[index], z[index], use_idpar, point))
          abort_cocal("COCAL BHT interpolation failed.");
        if (write_lapse)
          alp[index] = point[BHT_ALP];
        if (write_shift) {
          betax[index] = point[BHT_BETAX];
          betay[index] = point[BHT_BETAY];
          betaz[index] = point[BHT_BETAZ];
        }
        if (write_hydro) {
          rho[index] = point[BHT_RHO];
          press[index] = point[BHT_PRESS];
          eps[index] = point[BHT_EPS];
          vel[vector_gf_index(cctkGH, i, j, k, 0)] = point[BHT_VX];
          vel[vector_gf_index(cctkGH, i, j, k, 1)] = point[BHT_VY];
          vel[vector_gf_index(cctkGH, i, j, k, 2)] = point[BHT_VZ];
        }
        if (write_metric) {
          gxx[index] = point[BHT_GXX];
          gxy[index] = point[BHT_GXY];
          gxz[index] = point[BHT_GXZ];
          gyy[index] = point[BHT_GYY];
          gyz[index] = point[BHT_GYZ];
          gzz[index] = point[BHT_GZZ];
          kxx[index] = point[BHT_KXX];
          kxy[index] = point[BHT_KXY];
          kxz[index] = point[BHT_KXZ];
          kyy[index] = point[BHT_KYY];
          kyz[index] = point[BHT_KYZ];
          kzz[index] = point[BHT_KZZ];
        }
      }
    }
  }
}
