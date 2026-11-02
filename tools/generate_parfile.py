#!/usr/bin/env python3
"""Generate a Carpet parfile from a COCAL RNS, BNS, BBH or BBHT ID directory."""

import argparse
from array import array as doubles
from collections import namedtuple
import hashlib
import json
import math
from pathlib import Path
import re
import struct

# Shared thorn defaults. A configuration can supply its own active_thorns list.
# Dependencies and the selected evolution thorns are resolved together.
default_thorns = [
    "Carpet",
    "COCAL_ID",
    "Baikal",
    "GRHayLHD",

    "CarpetIOASCII",
    "CarpetIOBasic",
    "CarpetIOHDF5",
    "CarpetIOScalar",
    "TerminationTrigger",

    "ReduceLocTracker",
    "NSSurface",                     # adds SphericalSurface and interpolation
    "QuasiLocalMeasures",            # adds SphericalSurface and TGRtensor
    "VolumeIntegrals_GRMHD",         # always included for moving BNS grids
    "VolumeIntegrals_vacuum",
    "AHFinderDirect",                # adds SphericalSurface
    "Multipole",                     # adds WeylScal4
    "NaNChecker",
    "TimerReport",
    "SystemStatistics",
    "SystemTopology",
    "Trigger",
]

# Vacuum BBH selection. Per-run active_thorns overrides work the same way.
bbh_thorns = [
    "Carpet", "COCAL_ID", "BaikalVacuum",
    "CarpetIOASCII", "CarpetIOBasic", "CarpetIOScalar", "CarpetIOHDF5",
    "AHFinderDirect", "QuasiLocalMeasures", "VolumeIntegrals_vacuum",
    "NaNChecker", "TimerReport", "TerminationTrigger", "Multipole",
]

# BBHT uses the same horizon diagnostics and adds hydrodynamics and rest mass.
bbht_thorns = ["Baikal" if t == "BaikalVacuum" else t for t in bbh_thorns]
bbht_thorns += ["GRHayLHD", "VolumeIntegrals_GRMHD"]



# RUN CONFIGURATIONS: edit here, then run `python tools/generate_parfile.py 0`.
# Add entries for your own runs; keep the committed example generic.
# See tools/README.md for the available grid, run, solver and diagnostic options.
configs = [
    {
        "name": "example",
        "id_dir": "/path/to/ID",
        "profile": "semiprod",       # "local" or "semiprod"
        "grid": {},                  # e.g. {"finest_dx": 0.125, "levels": 8}
        "run": {},                   # e.g. {"iterations": 16, "out_every": 8}
    },
]


# Implementation below; the user settings end above.
# These also determine the conservative nesting allowance in resolve_grid().
GHOST_SIZE = 3
RK_STAGES = 4
G_CGS, C_CGS, MSUN_CGS = 6.67428e-8, 2.99792458e10, 1.98892e33
DENSITY_UNIT = (G_CGS / C_CGS**2)**3 * MSUN_CGS**2  # g/cm^3 -> G = c = Msol = 1
PROLONGATION_ORDER = 5
MIN_DISTANCE = 4

# Each spacetime entry gives its constraint variables, native evolved groups and the name of
# its runtime lapse floor (None: it has none). Adding a solver also requires dependencies,
# parameter settings and validation.
SPACETIME_THORNS = {
    "Baikal": {
        "lapse_floor": "lapse_floor",
        "constraints": ("Baikal::HGF", "Baikal::MU0GF", "Baikal::MU1GF", "Baikal::MU2GF"),
        "evolved": "Baikal::evol_variables",
    },
    "BaikalVacuum": {
        "lapse_floor": "lapse_floor",
        "constraints": ("BaikalVacuum::HGF", "BaikalVacuum::MU0GF", "BaikalVacuum::MU1GF", "BaikalVacuum::MU2GF"),
        "evolved": "BaikalVacuum::evol_variables",
    },
    "ML_BSSN": {
        "lapse_floor": "minimumLapse",
        "constraints": ("ML_BSSN::H", "ML_BSSN::M1", "ML_BSSN::M2", "ML_BSSN::M3"),
        "evolved": " ".join("ML_BSSN::" + group for group in (
            "ML_log_confac", "ML_metric", "ML_Gamma", "ML_trace_curv", "ML_curv",
            "ML_lapse", "ML_dtlapse", "ML_shift", "ML_dtshift")),
    },
    "ML_CCZ4": {
        "lapse_floor": "minimumLapse",
        "constraints": ("ML_CCZ4::H", "ML_CCZ4::M1", "ML_CCZ4::M2", "ML_CCZ4::M3"),
        "evolved": " ".join("ML_CCZ4::" + group for group in (
            "ML_log_confac", "ML_metric", "ML_Gamma", "ML_trace_curv", "ML_curv", "ML_Theta",
            "ML_lapse", "ML_dtlapse", "ML_shift", "ML_dtshift")),
    },
    "LeanBSSNMoL": {
        "lapse_floor": None,
        "constraints": ("LeanBSSNMoL::hc", "LeanBSSNMoL::mcx", "LeanBSSNMoL::mcy", "LeanBSSNMoL::mcz"),
        "evolved": "ADMBase::lapse ADMBase::shift " + " ".join("LeanBSSNMoL::" + group for group in (
            "conf_fac", "hmetric", "hcurv", "trk", "gammat")),
    },
}
# black_hole_interior: density ceiling and in-horizon limits for matter inside a horizon.
HYDRO_THORNS = {
    "GRHayLHD": {"black_hole_interior": True},
    "IllinoisGRMHD": {"black_hole_interior": True},
    "GRHydro": {"black_hole_interior": False},
}
# Gauge settings use the shared vocabulary below. These additional controls
# keep their solver's names and meanings; they cannot override the gauge map.
SPACETIME_OPTIONS = {
    "Baikal": {"diss_strength": float, "newrad_aDD_radpower": int},
    "BaikalVacuum": {"diss_strength": float},
    "ML_BSSN": {"epsDiss": float, "conformalMethod": int, "fixAdvectionTerms": int},
    "ML_CCZ4": {"epsDiss": float, "conformalMethod": int, "fixAdvectionTerms": int,
                "dampk1": float, "dampk2": float, "GammaShift": float},
    "LeanBSSNMoL": {"conf_fac_floor": float, "chi_gamma": float,
                    "make_aa_tracefree": bool, "reset_dethh": bool},
}


# Dependencies for this generator's Carpet evolution and diagnostic recipe.
# This is not a replacement for the build's thornlist or Cactus configuration.ccl.
THORN_DEPENDENCIES = {
    "Carpet": ("CarpetLib", "CarpetRegrid2", "CartGrid3D", "CoordBase", "CoordGauge",
               "Time", "SymBase", "LoopControl", "CycleClock", "Timers"),
    "COCAL_ID": ("ADMBase", "HydroBase", "StaticConformal", "InitBase"),
    "Baikal": ("ADMBase", "ADMCoupling", "ADMMacros", "Boundary", "MoL", "NewRad",
               "TmunuBase", "Vectors", "SummationByParts"),
    "BaikalVacuum": ("ADMBase", "ADMCoupling", "ADMMacros", "Boundary", "MoL",
                     "NewRad", "Vectors", "SummationByParts"),
    "GRHayLHD": ("GRHayLib", "HydroBase", "ADMBase", "ADMMacros", "TmunuBase",
                 "Boundary", "MoL", "SymBase", "SpaceMask"),
    "GRHayLID": ("GRHayLib", "HydroBase"),
    "IllinoisGRMHD": ("GRHayLib", "HydroBase", "ADMBase", "ADMMacros", "TmunuBase",
                      "Boundary", "MoL", "SymBase", "SpaceMask"),
    "GRHydro": ("EOS_Omni", "HydroBase", "ADMBase", "ADMMacros", "TmunuBase",
                "Boundary", "MoL", "SpaceMask", "CartGrid3D"),
    "ML_BSSN": ("ML_BSSN_Helper", "ADMBase", "ADMCoupling", "Boundary", "MoL",
                "NewRad", "TmunuBase", "GenericFD", "LoopControl", "Vectors"),
    "ML_BSSN_Helper": ("ML_BSSN", "CartGrid3D", "CoordGauge"),
    "ML_CCZ4": ("ML_CCZ4_Helper", "ADMBase", "ADMCoupling", "Boundary", "MoL",
                "NewRad", "TmunuBase", "GenericFD", "LoopControl", "Vectors"),
    "ML_CCZ4_Helper": ("ML_CCZ4", "CartGrid3D", "CoordGauge"),
    "LeanBSSNMoL": ("ADMBase", "TmunuBase", "Boundary", "CoordGauge", "MoL", "NewRad", "Dissipation"),
    "Dissipation": ("MoL", "SpaceMask", "SphericalSurface"),
    "CarpetLib": ("Vectors",),
    "CarpetRegrid2": ("Carpet",),
    "CarpetTracker": ("CarpetRegrid2", "SphericalSurface"),
    "CarpetInterp": ("Carpet", "AEILocalInterp"),
    "CarpetReduce": ("Carpet", "LocalReduce"),
    "CarpetSlab": ("Carpet", "Slab"),
    "CarpetIOASCII": ("Carpet", "IOUtil", "CarpetSlab"),
    "CarpetIOBasic": ("Carpet", "IOUtil", "CarpetReduce"),
    "CarpetIOHDF5": ("Carpet", "IOUtil", "CarpetSlab"),
    "CarpetIOScalar": ("Carpet", "IOUtil", "CarpetReduce"),
    "NSSurface": ("SphericalSurface", "CarpetInterp", "HydroBase"),
    "QuasiLocalMeasures": ("SphericalSurface", "TGRtensor", "CarpetInterp", "ADMMacros", "TmunuBase"),
    "AHFinderDirect": ("SphericalSurface", "CarpetInterp", "ADMMacros", "SpaceMask"),
    "Multipole": ("WeylScal4", "CarpetInterp"),
    "WeylScal4": ("ADMBase", "ADMMacros", "SummationByParts", "GenericFD"),
    "ReduceLocTracker": ("CarpetReduce", "HydroBase"),
    "VolumeIntegrals_GRMHD": ("CarpetReduce", "CarpetRegrid2", "HydroBase"),
    "VolumeIntegrals_vacuum": ("CarpetReduce", "ADMBase", "ADMMacros"),
    "Trigger": ("CarpetReduce",),
    "TimerReport": ("Timers",),
    "TerminationTrigger": ("IOUtil",),
    "ADMBase": (), "ADMCoupling": (), "ADMMacros": (), "AEILocalInterp": (),
    "Boundary": (), "CartGrid3D": (), "CoordBase": (), "CoordGauge": (),
    "CycleClock": (), "EOS_Omni": (), "GenericFD": (), "GRHayLib": (), "HydroBase": (), "InitBase": (), "IOUtil": (),
    "LocalReduce": (), "LoopControl": (), "MoL": (), "NaNChecker": (), "NewRad": (),
    "Slab": (), "SpaceMask": (), "StaticConformal": (), "SummationByParts": (),
    "SymBase": (), "SphericalSurface": (), "SystemStatistics": (), "SystemTopology": (),
    "TGRtensor": (), "Time": (), "Timers": (), "TmunuBase": (), "Vectors": (),
}


def resolve_thorns(requested, bns, tabulated=False, vacuum=False):
    """Keep the requested order, then append each dependency once."""
    if not isinstance(requested, list) or not all(isinstance(name, str) for name in requested):
        raise ValueError("active_thorns must be a list of thorn names")
    if len(requested) != len(set(requested)):
        raise ValueError("active_thorns contains duplicate names")

    active = list(requested)
    if bns and "VolumeIntegrals_GRMHD" not in active:
        active.append("VolumeIntegrals_GRMHD")
    if tabulated:
        for thorn in ("GRHayLID", "EOS_Omni"):
            if thorn not in active:
                active.append(thorn)

    # Appending while walking the list also resolves transitive dependencies.
    for thorn in active:
        if thorn not in THORN_DEPENDENCIES:
            raise ValueError(f"Unknown thorn: {thorn}")
        for dependency in THORN_DEPENDENCIES[thorn]:
            if dependency not in active:
                active.append(dependency)

    spacetime = [name for name in SPACETIME_THORNS if name in active]
    hydro = [name for name in HYDRO_THORNS if name in active]
    if len(spacetime) != 1:
        raise ValueError("Select exactly one spacetime thorn: " + ", ".join(SPACETIME_THORNS))
    if vacuum and hydro:
        raise ValueError("Vacuum BBH must not activate a hydro evolution thorn")
    if not vacuum and len(hydro) != 1:
        raise ValueError("Matter needs exactly one hydro thorn: " + ", ".join(HYDRO_THORNS))
    if not vacuum and spacetime == ["BaikalVacuum"]:
        raise ValueError("BaikalVacuum has no matter source terms; select a matter-capable spacetime thorn")
    if tabulated and hydro == ["GRHydro"]:
        raise ValueError("GRHydro tabulated EOS needs a separately validated EOS_Omni table/initialization recipe; use GRHayLHD or IllinoisGRMHD here")
    for thorn in ("Carpet", "COCAL_ID"):
        if thorn not in active:
            raise ValueError(f"This evolution recipe needs {thorn}; keep it in active_thorns")
    return active


def resolve_evolution(cfg, requested, data):
    """Replace explicitly selected solvers; otherwise respect active_thorns."""
    if not isinstance(requested, list) or not all(isinstance(t, str) for t in requested):
        raise ValueError("active_thorns must be a list of thorn names")
    if len(requested) != len(set(requested)):
        raise ValueError("active_thorns contains duplicate names")
    requested = list(requested)
    for key, choices in (("spacetime_thorn", SPACETIME_THORNS), ("hydro_thorn", HYDRO_THORNS)):
        selected = cfg[key]
        if selected is not None:
            position = next((i for i, name in enumerate(requested) if name in choices), len(requested))
            requested = [name for name in requested if name not in choices]
            if key == "spacetime_thorn":
                old_helpers = {name + "_Helper" for name in SPACETIME_THORNS if name != selected}
                requested = [name for name in requested if name not in old_helpers]
            requested.insert(position, selected)
    active = resolve_thorns(requested, data["kind"] == "bns",
                            tabulated=data["EOS_type"] == "Tabulated", vacuum=data["kind"] == "bbh")
    spacetime = cfg["spacetime_thorn"] = next(t for t in SPACETIME_THORNS if t in active)
    hydro = cfg["hydro_thorn"] = next((t for t in HYDRO_THORNS if t in active), None)
    if hydro and data["kind"] == "bbht":
        # Matter around a black hole: the puncture lapse collapses, and con2prim (u^0 ~ W/alpha)
        # needs a lapse floor and in-horizon limits to survive it.
        if not SPACETIME_THORNS[spacetime]["lapse_floor"]:
            raise ValueError(f"{spacetime} has no lapse floor, which matter around a black hole needs")
        if not HYDRO_THORNS[hydro]["black_hole_interior"]:
            raise ValueError(f"{hydro} has no black-hole interior limits here; use GRHayLHD or IllinoisGRMHD")
    return requested, active


def resolve_gauge(cfg, data):
    """Check runtime gauge capabilities against the selected implementation."""
    thorn = cfg["spacetime_thorn"]
    # Lapse floor 1e-8 for every ID (the BHT value). A puncture lapse collapses toward 0;
    # vacuum survives that, GRHayL's u^0 ~ W/alpha does not (BBHT 21641319, 21693998).
    gauge = dict(lapse="1+log", shift="gamma-driver", advect_lapse=True, advect_shift=True,
                 driver_order=1 if thorn == "LeanBSSNMoL" else 2,
                 eta=1 / data["source_adm_mass"],
                 lapse_floor=1e-8 if SPACETIME_THORNS[thorn]["lapse_floor"] else None)
    unknown = cfg["gauge"].keys() - gauge.keys()
    if unknown:
        raise ValueError("Unknown gauge settings: " + ", ".join(sorted(unknown)))
    gauge.update(cfg["gauge"])
    if gauge["lapse"] not in ("1+log", "harmonic") or gauge["shift"] != "gamma-driver":
        raise ValueError('gauge supports lapse="1+log"/"harmonic" and shift="gamma-driver"')
    for name in ("advect_lapse", "advect_shift"):
        if type(gauge[name]) is not bool:
            raise ValueError(f"gauge.{name} must be True or False")
    if type(gauge["driver_order"]) is not int or gauge["driver_order"] not in (1, 2):
        raise ValueError("gauge.driver_order must be 1 or 2")
    for name in ("eta", "lapse_floor"):
        value = gauge[name]
        if value is None and name == "lapse_floor" and not SPACETIME_THORNS[thorn]["lapse_floor"]:
            continue
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError(f"gauge.{name} must be a finite nonnegative number")
    if thorn in ("Baikal", "BaikalVacuum") and (gauge["lapse"] != "1+log" or gauge["driver_order"] != 2
                                                  or not gauge["advect_lapse"] or not gauge["advect_shift"]):
        raise ValueError(f"{thorn}'s generated equations fix advective 1+log and the second-order Gamma driver; changing them requires regenerating the thorn")
    if gauge["lapse_floor"] is not None and not SPACETIME_THORNS[thorn]["lapse_floor"]:
        raise ValueError(f"{thorn} has no runtime lapse floor")
    if thorn == "LeanBSSNMoL" and (gauge["lapse"] != "1+log" or gauge["driver_order"] != 1):
        raise ValueError("LeanBSSNMoL provides 1+log and a first-order Gamma driver, without a runtime lapse floor")
    return gauge


def spacetime_parameters(thorn, gauge, options, *, black_holes=False):
    """Fourth-order moving-puncture settings in each solver's own convention."""
    eta = gauge["eta"]
    if thorn in ("Baikal", "BaikalVacuum"):
        parameters = {**({"FD_order": 4} if black_holes else {}),
                      "diss_strength": "0.2", "eta": real(eta)}
        if gauge["lapse_floor"] != 1e-15:
            parameters["lapse_floor"] = real(gauge["lapse_floor"])
    elif thorn == "LeanBSSNMoL":
        parameters = dict(derivs_order=4, calculate_constraints="yes", calculate_constraints_every=1,
                          precollapsed_lapse="no", eta_beta=real(eta), beta_Gamma="0.75",
                          zeta_alpha=int(gauge["advect_lapse"]), zeta_beta=int(gauge["advect_shift"]))
    else:
        parameters = dict(fdOrder=4, conformalMethod=1,
                          harmonicN=1 if gauge["lapse"] == "1+log" else 2,
                          harmonicF=2 if gauge["lapse"] == "1+log" else 1,
                          evolveA=0, evolveB=gauge["driver_order"] - 1,
                          advectLapse=int(gauge["advect_lapse"]), advectShift=int(gauge["advect_shift"]),
                          shiftGammaCoeff=0.75,
                          # With evolveB=0, McLachlan multiplies the entire
                          # (Gamma - betaDriver*beta) term by shiftGammaCoeff.
                          betaDriver=real(eta if gauge["driver_order"] == 2 else eta / 0.75),
                          minimumLapse=real(gauge["lapse_floor"]), epsDiss="0.2",
                          initial_boundary_condition='"extrapolate-gammas"', rhs_boundary_condition='"NewRad"')
    allowed = SPACETIME_OPTIONS[thorn]
    for name, value in options.items():
        if name not in allowed:
            raise ValueError(f"Unsupported {thorn} spacetime_options key: {name}; gauge controls belong in gauge")
        kind = allowed[name]
        if (kind in (int, bool) and type(value) is not kind) or (kind is float and type(value) not in (int, float)):
            raise ValueError(f"{thorn}.{name} needs a {kind.__name__} value")
        parameters[name] = ("yes" if value else "no") if kind is bool else real(value)
    return parameters


def hydro_variables(thorn):
    """Native fields for finite checks and output; HydroBase stays common."""
    if thorn == "GRHydro":
        return "GRHydro::dens GRHydro::scon GRHydro::tau", ""
    magnetic = (" IllinoisGRMHD::Ax IllinoisGRMHD::Ay IllinoisGRMHD::Az"
                " IllinoisGRMHD::phitilde IllinoisGRMHD::grmhd_B_center") if thorn == "IllinoisGRMHD" else ""
    return (f"{thorn}::grmhd_velocities {thorn}::grmhd_conservatives" + magnetic,
            f"{thorn}::failure_checker")


# Matter limits for every ID, scaled by the ID's peak rest-mass density: the atmosphere
# far below any fluid, the ceiling well above a merger remnant or collapsing core
# (BNS SP_K123.6: 7.3e-13 and 0.022, the validated 1e-12 and 0.02). Psi6threshold marks
# black-hole interiors for the con2prim fixes.
ATMOSPHERE_FRACTION = 1e-9
CEILING_FACTOR = 30
PSI6_THRESHOLD = 60


def matter_limits(data):
    peak = data["peak_density"]
    return dict(atmosphere=ATMOSPHERE_FRACTION * peak, ceiling=CEILING_FACTOR * peak)


def grhayl_limits(data, table=None):
    """GRHayLib floors and ceiling, shared by GRHayLHD and IllinoisGRMHD. GRHayL clamps
    rho_min/rho_max to a table itself, but evaluates the atmosphere unchecked, so keep it
    inside the table."""
    limits = matter_limits(data)
    atmosphere, ceiling = limits["atmosphere"], limits["ceiling"]
    if table:
        atmosphere, ceiling = max(atmosphere, table["rho"][0]), min(ceiling, table["rho"][1])
    return dict(rho_b_atm=real(atmosphere), rho_b_min=real(atmosphere),
                rho_b_max=real(ceiling), Psi6threshold=PSI6_THRESHOLD)


def table_bounds(path):
    """rho (code units), T (MeV) and Y_e ranges of a GRHayL table; needs h5py, like
    make_grhayl_teos.py."""
    try:
        import h5py
    except ImportError:
        raise ValueError("Tabulated IDs need h5py to read the EOS table bounds")
    with h5py.File(path, "r") as table:
        logrho, logtemp, ye = (table[name][()] for name in ("logrho", "logtemp", "ye"))
    return dict(rho=(10**logrho.min() * DENSITY_UNIT, 10**logrho.max() * DENSITY_UNIT),
                T=(10**logtemp.min(), 10**logtemp.max()), Y_e=(ye.min(), ye.max()))


def grhydro_parameters(data):
    atmosphere = real(matter_limits(data)["atmosphere"])
    return dict(GRHydro_eos_type='"General"', GRHydro_eos_table='"Hybrid"',
                recon_method='"ppm"', riemann_solver='"HLLE"', GRHydro_stencil=3,
                sources_spatial_order=4,
                bound='"flat"', rho_abs_min=atmosphere, initial_rho_abs_min=atmosphere,
                GRHydro_enable_internal_excision="no")


def horizon_shape_every(run, find_every):
    """Horizon shapes are 2D-type output, written on iterations where a search runs."""
    every = run["out2d_every"]
    return every * find_every // math.gcd(every, find_every) if every > 0 else 0


def thermal_gamma(data):
    """BBHT's disk is one polytrope and keeps its Gamma; stars use 5/3."""
    return data["EOS_gamma"][0] if data["kind"] == "bbht" else 5 / 3


def hybrid_eos_parameters(data):
    """Map the same cold EOS and thermal Gamma to EOS_Omni's Hybrid handle."""
    parameters = dict(n_pieces=len(data["EOS_gamma"]), hybrid_k0=real(data["EOS_k0"]),
                      hybrid_gamma_th=real(thermal_gamma(data)),
                      poly_k=real(data["EOS_k0"]))
    parameters.update({f"hybrid_gamma[{i}]": real(gamma) for i, gamma in enumerate(data["EOS_gamma"])})
    parameters.update({f"hybrid_rho[{i}]": real(rho) for i, rho in enumerate(data["EOS_boundaries"])})
    return parameters


def quoted(value):
    value = str(value)
    if any(char in value for char in ('"', "\n", "\r", "$", "@")):
        raise ValueError(f"Unsupported character in parfile string: {value!r}")
    return f'"{value}"'


def real(value):
    if not math.isfinite(value):
        raise ValueError("Numeric settings must be finite")
    return f"{value:.16g}"


def array(values):
    return "[" + ", ".join(real(value) for value in values) + "]"


def variables(*rows):
    return '"\n  ' + '\n  '.join(rows) + '\n"'


def resolve_grid(data, overrides, local):
    """Build surface-sized stellar stacks and shared parents, then fit the domain.

    All widths are coordinate half-widths. Grid spacings double outward; box
    widths approximately double, with extra room where numerical nesting needs it.
    The outer domain never sets an inner box's size.
    """
    bns = data["kind"] == "bns"
    centers = data["centers"]
    radii = data["max_coordinate_star_radii"]

    # 1. Apply user overrides to the profile defaults.
    default_levels = 4 if local else (8 if bns else 6)
    grid = {
        "finest_dx": 0.5 if local else 0.125,
        "levels": default_levels,
        "domain": None,
        "star_margin": 0.20,
        "star_padding": 0.0,
        "common_levels": None,
        "movement_cells": 2.0,
    }
    unknown = overrides.keys() - grid.keys()
    if unknown:
        raise ValueError("Unknown grid settings: " + ", ".join(sorted(unknown)))
    grid.update(overrides)

    dx = grid["finest_dx"]
    levels = grid["levels"]
    if type(levels) is not int or not 2 <= levels <= 10:
        raise ValueError("grid levels must be an integer in 2..10, including the base level")
    for name in ("finest_dx", "movement_cells", "star_margin", "star_padding"):
        value = grid[name]
        positive = name in ("finest_dx", "movement_cells")
        if (type(value) not in (int, float) or not math.isfinite(value)
                or value < 0 or (positive and value == 0)):
            bound = "positive" if positive else "nonnegative"
            raise ValueError(f"grid {name} must be a {bound} finite number")
    coarse_dx = dx * 2**(levels - 1)

    def aligned(width, spacing):
        return math.ceil(width / (2 * spacing)) * (2 * spacing)

    def nesting_margin(parent_dx):
        # CarpetRegrid2 nests the buffered fine region plus its ghosts, then
        # adds min_distance and the prolongation stencil on the parent grid.
        # Two further parent cells allow for coordinate rounding and snapping.
        fine_support = RK_STAGES * GHOST_SIZE * (parent_dx / 2)
        coarse_support = (MIN_DISTANCE + PROLONGATION_ORDER // 2 + 2) * parent_dx
        return fine_support + coarse_support

    # 2. Cover the surface allowance after tracking lag and centre rounding.
    # The local profile keeps a minimum half-width for the small fixtures.
    tracking_margin = grid["movement_cells"] * dx if bns else 0.0
    rounding_margin = dx / 2 if bns else 0.0
    fine = []
    for radius in radii:
        allowance = max(grid["star_margin"] * radius, grid["star_padding"])
        half_width = aligned(radius + allowance + tracking_margin + rounding_margin, dx)
        fine.append(max(12 if local else 0, half_width))

    star_boxes = []
    for radius in fine:
        boxes = [radius]
        for level in range(levels - 2, 0, -1):
            parent_dx = dx * 2**(levels - 1 - level)
            parent = max(2 * boxes[0], boxes[0] + nesting_margin(parent_dx))
            boxes.insert(0, aligned(parent, parent_dx))
        star_boxes.append(boxes)

    # 3. Use shared outer boxes once the stellar boxes substantially overlap.
    # Keep at least the finest level star-centred. The overlap must cover half
    # the smaller box's x-width; a mere first touch is a poor handoff to a cube.
    separation = abs(centers[1] - centers[0]) if bns else 0
    # Only for the diagnostic spheres, not the grid.
    wide = bns and separation > 8 * max(radii)
    common_levels = grid["common_levels"]
    if common_levels is not None:
        if not bns or local:
            raise ValueError("common_levels applies only to semiprod BNS grids")
        if type(common_levels) is not int or not 1 <= common_levels <= levels:
            raise ValueError("common_levels must be in 1..levels, including the base level")
    elif bns and not local:
        common_levels = 1  # base level only, if no shared refined parent fits
        for level in range(levels - 2, 0, -1):
            widths = [boxes[level - 1] for boxes in star_boxes]
            lower = max(x - width for x, width in zip(centers, widths))
            upper = min(x + width for x, width in zip(centers, widths))
            if upper - lower >= min(widths):
                common_levels = level + 1
                break
    else:
        common_levels = 0

    common_boxes = []
    if bns and not local:
        # Inactive central radii follow the stellar size, for a later collapse.
        common_boxes = [max(boxes[level - 1] for boxes in star_boxes)
                        for level in range(1, levels)]
        for level in range(levels - 2, 0, -1):
            parent_dx = dx * 2**(levels - 1 - level)
            child_extent = common_boxes[level]
            if level < common_levels:
                # Support the next-finer union. Same-level stellar requests may
                # extend beyond this central box; Carpet takes their union.
                child_extent = max(child_extent, max(
                    abs(x) + boxes[level] for x, boxes in zip(centers, star_boxes)))
            required = child_extent + nesting_margin(parent_dx)
            common_boxes[level - 1] = aligned(
                max(common_boxes[level - 1], required), parent_dx)

    # 4. Fit the stacks inside the domain, leaving room for boundary stencils.
    boundary_clearance = nesting_margin(coarse_dx)
    stack_extent = [
        max(abs(x) + boxes[0] for x, boxes in zip(centers, star_boxes)),
        max(boxes[0] for boxes in star_boxes),
        max(boxes[0] for boxes in star_boxes),
    ]
    if common_levels > 1:
        stack_extent = [max(extent, common_boxes[0]) for extent in stack_extent]
    domain = grid["domain"]
    if domain is None:
        if local:
            domain = [max(64, extent + boundary_clearance) for extent in stack_extent]
        else:
            stellar_extent = max(abs(x) + r for x, r in zip(centers, radii))
            half_width = max(
                1536 if bns else 384,
                4 * stellar_extent,
                max(stack_extent) + boundary_clearance,
            )
            domain = [half_width] * 3

        # Eight coarse cells give the domain convenient divisibility.
        domain_unit = 8 * coarse_dx
        domain = [math.ceil(width / domain_unit) * domain_unit for width in domain]

    if not isinstance(domain, (list, tuple)) or len(domain) != 3:
        raise ValueError("domain needs three positive half-widths aligned to coarse cells")
    for width in domain:
        if type(width) not in (int, float) or not math.isfinite(width) or width <= 0:
            raise ValueError("domain needs three positive half-widths aligned to coarse cells")
        coarse_cells = width / coarse_dx
        if not math.isclose(coarse_cells, round(coarse_cells)):
            raise ValueError("domain half-widths must align to coarse cells")

    if any(width - extent < boundary_clearance
           for width, extent in zip(domain, stack_extent)):
        raise ValueError(
            f"Refinement stacks need {boundary_clearance / coarse_dx:g} coarse cells of boundary clearance; "
            "enlarge domain or reduce levels/surface padding"
        )

    # 5. CoM and density searches must contain their star, but not its companion.
    region = []
    for half_width, radius in zip(fine, radii):
        search_radius = half_width + 1
        if bns:
            search_radius = min(search_radius, separation / 2 - dx)
        if search_radius <= radius:
            raise ValueError("Stellar diagnostic spheres must enclose their star without reaching the companion")
        region.append(search_radius)

    grid.update(
        domain=list(domain),
        common_levels=common_levels,
        coarse_dx=coarse_dx,
        wide=wide,
        fine_box_radii=fine,
        star_boxes=star_boxes,
        common_boxes=common_boxes,
        tracking_radii=region,
    )
    return grid


def settings(id_dir, overrides=None):
    """Apply one run's overrides; grid and run defaults come from its profile."""
    values = {
        "name": None,
        "profile": "semiprod",
        "compact": False,
        "xunit": False,
        "ecc_cor_velx": 0.0,
        "ecc_cor_velz": 0.0,
        "eos_table": None,
        "surface_density": 1e-7,
        "grid": {},
        "run": {},
        "active_thorns": default_thorns,
        "spacetime_thorn": None,   # BBH: automatic BaikalVacuum; stellar: Baikal
        "hydro_thorn": None,       # matter: GRHayLHD; vacuum: no hydro solver
        "gauge": {},
        "spacetime_options": {},
    }
    if overrides is not None:
        if not isinstance(overrides, dict):
            raise ValueError("Each configuration must be a dictionary")
        unknown = overrides.keys() - values.keys()
        if unknown:
            raise ValueError("Unknown configuration settings: " + ", ".join(sorted(unknown)))
        values.update(overrides)

    if values["profile"] not in ("local", "semiprod"):
        raise ValueError('profile must be "local" or "semiprod"')
    if values['spacetime_thorn'] is not None and (not isinstance(values['spacetime_thorn'], str) or values['spacetime_thorn'] not in SPACETIME_THORNS):
        raise ValueError('spacetime_thorn must be ' + ', '.join(SPACETIME_THORNS))
    if values['hydro_thorn'] is not None and values['hydro_thorn'] not in HYDRO_THORNS:
        raise ValueError('hydro_thorn must be ' + ', '.join(HYDRO_THORNS))
    for name in ("compact", "xunit"):
        if type(values[name]) is not bool:
            raise ValueError(f"{name} must be True or False")
    for name in ("ecc_cor_velx", "ecc_cor_velz"):
        velocity_pair(values[name], name)
    for name in ("grid", "run", "gauge", "spacetime_options"):
        if not isinstance(values[name], dict):
            raise ValueError(f"{name} overrides must be a dictionary")
        values[name] = dict(values[name])
    output_name = values["name"]
    if output_name is not None:
        if (not isinstance(output_name, str) or not output_name.strip()
                or output_name in (".", "..") or Path(output_name).name != output_name):
            raise ValueError("name must be a filename stem, not a directory path")
    if not isinstance(id_dir, (str, Path)) or not str(id_dir).strip():
        raise ValueError("Set id_dir to the COCAL ID directory")
    values["id_dir"] = Path(id_dir).expanduser().resolve()
    if values["eos_table"] is not None:
        table = Path(values["eos_table"]).expanduser().resolve()
        if not table.is_file():
            raise ValueError(f"EOS table does not exist: {table}")
        values["eos_table"] = table
    return values


def velocity_pair(value, name):
    """Validate a shared component or a pair in COCAL patch-1/patch-2 order."""
    pair = value if isinstance(value, (list, tuple)) else (value, value)
    if len(pair) != 2 or any(type(v) not in (int, float) or not math.isfinite(v) for v in pair):
        raise ValueError(f"{name} must be a finite number or two finite numbers")
    return pair


def has_bns_velocity(cfg):
    return any(v != 0 for key in ("ecc_cor_velx", "ecc_cor_velz")
               for v in velocity_pair(cfg[key], key))


def select_config(target):
    """A bare integer selects configs[index]; any other argument is an ID path."""
    if target.isdecimal():
        index = int(target)
        if not isinstance(configs, list) or index >= len(configs):
            raise ValueError(f"No configuration {index}; add it to configs near the top of this file")
        if not isinstance(configs[index], dict):
            raise ValueError(f"Configuration {index} must be a dictionary")
        overrides = dict(configs[index])
        if "id_dir" not in overrides:
            raise ValueError(f"Configuration {index} needs id_dir")
        id_dir = overrides.pop("id_dir")
        return id_dir, overrides, index
    return target, {}, None


def resolve_run(cfg, data, grid, active):
    """Resolve cadence and runtime controls, then check numerical limits."""
    local = cfg["profile"] == "local"
    # Every stellar ID can collapse: the trigger adds two finer levels at its centre.
    collapse = data["kind"] in ("rns", "bns") and not local and "Trigger" in active
    levels = grid["levels"]
    max_levels = min(10, levels + 2) if collapse else levels
    # Carpet iterations include levels reserved for a later collapse trigger.
    sync_every = 2**(max_levels - 1)
    run = {
        "iterations": 8 if local else 100000000,
        "final_time": 1000,
        "dtfac": 0.25 if local else 0.4,
        # Output cadences in iterations; 1D/2D default to multiples of out_every, 3D is off.
        "out_every": sync_every,
        "out1d_every": None,
        "out2d_every": None,
        "out3d_every": -1,
        "nan_check_every": sync_every,
        "poison": local,
        "walltime_hours": 1 if local else 12,
        "termination_margin_minutes": 5 if local or data["kind"] in ("bbh", "bbht") else 45,
        "checkpoint_hours": -1 if local or data["kind"] in ("bbh", "bbht") else 3,
        "checkpoint_on_terminate": not local and data["kind"] not in ("bbh", "bbht"),
    }
    if cfg["run"].keys() - run.keys():
        raise ValueError("run accepts " + ", ".join(run))
    run.update(cfg["run"])
    every = run["out_every"]
    for name in ("out1d_every", "out2d_every"):
        if run[name] is None:
            run[name] = every if local else 4 * every
    for name in ("out_every", "out1d_every", "out2d_every", "out3d_every", "nan_check_every"):
        value = run[name]
        if value == -1 and name in ("out1d_every", "out2d_every", "out3d_every"):
            continue
        # Divisors of the coarse step repeat on the same levels; multiples include every level.
        if type(value) is not int or value < 1 or (sync_every % value and value % sync_every):
            raise ValueError(f"{name} must divide or be a multiple of the coarse step ({sync_every} iterations)"
                             + (", or -1 to disable" if name != "out_every" and name != "nan_check_every" else ""))
    if type(run["iterations"]) is not int or run["iterations"] < 0:
        raise ValueError("iterations must be a nonnegative integer")
    if type(run["poison"]) is not bool:
        raise ValueError("poison must be True or False")
    for name in ("final_time", "walltime_hours", "dtfac"):
        if not math.isfinite(run[name]) or run[name] <= 0:
            raise ValueError(f"{name} must be positive")
    margin = run["termination_margin_minutes"]
    if not math.isfinite(margin) or not 0 < margin < 60 * run["walltime_hours"]:
        raise ValueError("termination_margin_minutes must be positive and less than the allocated walltime")
    checkpoint_hours = run["checkpoint_hours"]
    if not math.isfinite(checkpoint_hours) or (checkpoint_hours != -1 and checkpoint_hours <= 0):
        raise ValueError("checkpoint_hours must be positive, or -1 to disable periodic checkpoints")
    if type(run["checkpoint_on_terminate"]) is not bool:
        raise ValueError("checkpoint_on_terminate must be True or False")
    cfg["gauge"] = resolve_gauge(cfg, data)
    eta = cfg["gauge"]["eta"]
    if eta * run["dtfac"] * grid["coarse_dx"] >= 2.785293563:
        raise ValueError("Coarse-grid eta*dt exceeds the isolated RK4 damping limit; reduce spacing or levels")

    # HydroBase is refreshed as often as anything reads it: outputs and the per-step trackers.
    used = [run[name] for name in ("out_every", "out1d_every", "out2d_every", "out3d_every") if run[name] > 0]
    run.update(track_every=sync_every, max_levels=max_levels, collapse=collapse, eta=eta,
               hydro_every=math.gcd(sync_every, *used))
    return run


def resolve_stellar(id_dir, overrides):
    """Settings, metadata, thorns, grid and run for an RNS or BNS ID."""
    cfg = settings(id_dir, overrides)
    data = CocalMetadata(cfg["id_dir"]).read()
    if data["EOS_type"] == "Tabulated" and cfg["eos_table"] is None:
        raise ValueError("This ID has a tabulated EOS; set eos_table to the matching GRHayL HDF5 table")
    if data["EOS_type"] == "Hybrid" and cfg["eos_table"] is not None:
        raise ValueError("eos_table applies only to tabulated EOS IDs")
    if data["kind"] != "bns" and (cfg["compact"] or cfg["xunit"] or has_bns_velocity(cfg)):
        raise ValueError("compact, xunit and eccentricity velocities apply to BNS only")
    requested, active = resolve_evolution(cfg, cfg["active_thorns"], data)
    grid = resolve_grid(data, cfg["grid"], cfg["profile"] == "local")
    run = resolve_run(cfg, data, grid, active)
    return cfg, data, grid, run, requested, active


def resolve_black_holes(id_dir, overrides, disk):
    """Settings, metadata, thorns, grid and run for a BBH or BBHT ID; local by default."""
    cfg = settings(id_dir, {"profile": "local", **(overrides or {})})
    if any(cfg[key] for key in ("compact", "xunit", "eos_table")) or has_bns_velocity(cfg):
        raise ValueError("Stellar EOS/export options do not apply to BBH/BBHT")
    data = read_bbht_metadata(cfg["id_dir"]) if disk else read_bbh_metadata(cfg["id_dir"])
    grid = (resolve_bbht_grid if disk else resolve_bbh_grid)(data, cfg["grid"], cfg["profile"] == "local")
    requested = list((overrides or {}).get("active_thorns", bbht_thorns if disk else bbh_thorns))
    if grid["track_horizons"]:
        requested += [thorn for thorn in ("AHFinderDirect", "CarpetTracker") if thorn not in requested]
    elif "CarpetTracker" in requested:
        raise ValueError("Enable grid.track_horizons to configure CarpetTracker for BBH")
    requested, active = resolve_evolution(cfg, requested, data)
    levels = grid["levels"]
    cfg["run"].setdefault("iterations", (1 if disk else 4) * 2**(levels - 1))
    run = resolve_run(cfg, data, grid, active)
    run.update(finest_dt=run["dtfac"] * grid["finest_dx"], coarse_dt=run["dtfac"] * grid["coarse_dx"],
               synchronized_steps=run["iterations"] / 2**(levels - 1))
    return cfg, data, grid, run, requested, active


# Tabulated EOS atmosphere temperature (MeV) and electron fraction.
TABULATED_ATMOSPHERE = dict(T=0.1, Y_e=0.05)


class Setup:
    """The resolved inputs, and the layout several thorn blocks share: the
    spherical surfaces, the integration spheres and, for black holes, the
    horizon searches."""

    def __init__(self, cfg, data, grid, run, active):
        self.cfg, self.data, self.grid, self.run, self.active = cfg, data, grid, run, active
        self.kind = data["kind"]
        self.bns = self.kind == "bns"
        self.holes = self.kind in ("bbh", "bbht")
        self.disk = self.kind == "bbht"
        self.local = cfg["profile"] == "local"
        self.evolution, self.hydro = cfg["spacetime_thorn"], cfg["hydro_thorn"]
        self.constraint_vars = SPACETIME_THORNS[self.evolution]["constraints"]
        self.constraints = " ".join(self.constraint_vars)
        self.centers = data["centers"]
        self.table = None
        if self.holes:
            self.hole_layout()
        else:
            self.star_layout()

    def star_layout(self):
        grid, data = self.grid, self.data
        self.radii = data["max_coordinate_star_radii"]
        self.region = grid["tracking_radii"]
        cutoff = self.cfg["surface_density"]
        atmosphere = matter_limits(data)["atmosphere"]
        if not math.isfinite(cutoff) or cutoff <= atmosphere:
            raise ValueError(f"surface_density must exceed the atmosphere density {atmosphere:g}")
        self.native, self.failure = hydro_variables(self.hydro)
        if data["EOS_type"] == "Tabulated":
            self.table = table_bounds(self.cfg["eos_table"])
            for name, value in TABULATED_ATMOSPHERE.items():
                low, high = self.table[name]
                if not low <= value <= high:
                    raise ValueError(f"{name}_atm = {value} is outside the EOS table ({low:g}..{high:g})")

        coarse_dx = grid["coarse_dx"]
        self.outer = outer = min(grid["domain"])
        self.stellar_extent = max(abs(x) + r for x, r in zip(self.centers, self.radii))
        self.global_radii = []
        for numerator in (6, 7, 8, 9):
            radius = outer * numerator / 10
            if not self.stellar_extent + 4 * coarse_dx < radius < outer - 4 * coarse_dx:
                radius = None
            self.global_radii.append(radius)
        # BNS semiproduction adds shared outer levels to the two moving stacks.
        # Local runs omit the common stack to keep their memory footprint bounded.
        self.offset = 1 if self.bns and not self.local else 0

        # Stable surface slots: stellar exterior spheres, outer spheres, horizon, contours.
        nstar = len(self.centers)
        self.horizon_slot = nstar + 4
        self.contour_slots = list(range(self.horizon_slot + 1, self.horizon_slot + 1 + nstar))
        if "NSSurface" in self.active:
            count = self.contour_slots[-1] + 1
        elif "AHFinderDirect" in self.active:
            count = self.horizon_slot + 1
        else:
            count = self.horizon_slot
        self.surfaces = []
        for slot in range(count):
            if slot < nstar:
                name = ("Lhs" if slot == 0 else "Rhs") if self.bns else "star"
                spec = dict(name=name, x=self.centers[slot], radius=1.05 * self.radii[slot], enabled=True)
            elif slot < self.horizon_slot:
                fraction = (0.6, 0.7, 0.8, 0.9)[slot - nstar]
                radius = self.global_radii[slot - nstar]
                spec = dict(name=f"outer {fraction:g}", x=0, radius=radius, enabled=radius is not None)
            elif slot == self.horizon_slot:
                spec = dict(name="apparent horizon", x=0, radius=None, enabled="AHFinderDirect" in self.active)
            else:
                star = slot - self.contour_slots[0]
                name = ("Lhs contour" if star == 0 else "Rhs contour") if self.bns else "stellar contour"
                spec = dict(name=name, x=self.centers[star], radius=None, enabled=True)
            self.surfaces.append(spec)

        # Vacuum integrals, in slot order: whole-domain totals, BNS angular momentum
        # and constraints about each star, then the finite ADM spheres.
        self.vacuum_integrals, self.finite_adm = [], []

        def integral(quantity, radius=-1, star=None):
            slot = len(self.vacuum_integrals) + 1
            self.vacuum_integrals.append(dict(slot=slot, quantity=quantity, radius=radius, star=star))
            return slot

        integral("ADM_Mass")
        integral("ADM_Angular_Momentum")
        if self.bns:
            for radius in ([256, 384, 512] if grid["wide"] else [40, 60, 80]):
                if self.stellar_extent < radius < outer - 4 * coarse_dx:
                    integral("ADM_Angular_Momentum", radius)
        integral("H_M_CnstraintsL2")
        if self.bns:
            for star, radius in enumerate(self.region):
                integral("H_M_CnstraintsL2", radius, star)
        integral("ADM_Momentum")
        for radius in self.global_radii[:3]:
            if radius is not None:
                self.finite_adm.append(dict(radius=radius, mass_slot=integral("ADM_Mass", radius),
                                            angular_momentum_slot=integral("ADM_Angular_Momentum", radius)))

    def hole_layout(self):
        grid, data, run = self.grid, self.data, self.run
        self.tracking = grid["track_horizons"]
        # Search and regrid every 1/8 coarse step: a hole moves ~0.08 per coarse step, beyond the
        # tracking clearance (probe 21639370 lost both horizons; 16 at 8 levels found them all).
        self.horizon_every = max(1, run["track_every"] // 8)
        self.find_every = self.horizon_every if self.tracking else run["out_every"]
        # Stable slots: two actual horizons, followed by two finite-radius spheres.
        domain = grid["domain"]
        outer = [.5 * min(domain), .7 * min(domain)]
        if self.disk:
            # Extraction spheres must enclose the torus, including on small grids.
            outer[0] = max(outer[0], data["disk_outer_radius"] + 4 * grid["coarse_dx"])
            outer[1] = max(outer[1], outer[0] + 2 * grid["coarse_dx"])
            if outer[1] >= min(domain) - 4 * grid["coarse_dx"]:
                raise ValueError("BBHT extraction spheres need a larger outer domain")
        self.outer_spheres = outer
        self.surfaces = ([dict(name=f"BH{m + 1}", x=x, radius=None, enabled=True) for m, x in enumerate(self.centers)]
                         + [dict(name=f"outer{s}", x=0, radius=radius, enabled=True) for s, radius in enumerate(outer)])
        # Common horizon after a merger, about the origin (centre of mass). Coordinate-separation
        # Peters time underestimates the merger (qc0: ~2 vs ~18), so searches start at a quarter
        # of it: early enough never to miss the merger, cheap until then.
        separation = abs(self.centers[1] - self.centers[0])
        self.common_after = 0.25 * 5 / 64 * separation**4 / data["source_adm_mass"]**3
        # Output groups. GRHydro declares its three conservative groups separately.
        self.stress = (["TmunuBase::stress_energy_scalar TmunuBase::stress_energy_vector TmunuBase::stress_energy_tensor"]
                       if self.evolution != "BaikalVacuum" else [])
        self.matter_outputs = []
        if self.disk:
            native, failure = hydro_variables(self.hydro)
            conserved = native.partition(" ")[2] if self.hydro != "GRHydro" else native
            self.matter_outputs = ["HydroBase::rho HydroBase::press HydroBase::eps HydroBase::vel HydroBase::w_lorentz",
                                   " ".join(filter(None, (conserved, failure)))]


# One function per thorn returns that thorn's parameters (unqualified names),
# or a Block when it also sets the section title or integral labels. Stellar
# (RNS/BNS) and black-hole (BBH/BBHT) IDs share each function; where their
# settings differ, both are written there, the black-hole case first. Only
# active thorns are called, and an empty result writes no block.
Block = namedtuple("Block", "parameters title labels", defaults=(None, None))


def yes_no(flag):
    return "yes" if flag else "no"


# Execution and resource controls

def cactus(s):
    return {"terminate": '"either"', "cctk_itlast": s.run["iterations"],
            "cctk_final_time": real(s.run["final_time"])}


def termination_trigger(s):
    run = s.run
    limits = {"max_walltime": real(run["walltime_hours"]),
              "on_remaining_walltime": real(run["termination_margin_minutes"])}
    if s.holes:
        return limits
    return {**limits, "check_file_every": 256, "termination_file": '"TerminationTrigger.txt"',
            "termination_from_file": '"yes"', "create_termination_file": '"yes"'}


def timer_report(s):
    if s.holes:
        return {"out_every": s.run["out_every"], "out_filename": '"TimerReport"'}
    return {"out_every": "$out_every", "out_filename": '"TimerReport"', "output_all_timers": '"yes"',
            "output_all_timers_together": '"yes"', "output_all_timers_readable": '"yes"', "n_top_timers": 40}


# Grid, refinement and timestep

def cart_grid_3d(s):
    return {"type": '"coordbase"', "domain": '"full"', "avoid_origin": '"no"'}


def carpet(s):
    poison = yes_no(s.run["poison"])
    if s.holes:
        return {"domain_from_coordbase": '"yes"', "ghost_size": GHOST_SIZE, "max_refinement_levels": s.grid["levels"],
                "refinement_centering": '"vertex"', "use_buffer_zones": '"yes"',
                "prolongation_order_space": PROLONGATION_ORDER, "prolongation_order_time": 2,
                "init_fill_timelevels": '"yes"', "poison_new_timelevels": poison, "check_for_poison": '"no"'}
    return {"ghost_size": GHOST_SIZE, "domain_from_coordbase": '"yes"', "refinement_centering": '"vertex"',
            "enable_all_storage": "no", "use_buffer_zones": '"yes"', "poison_new_timelevels": poison,
            "check_for_poison": '"no"', "init_3_timelevels": "no", "init_fill_timelevels": '"yes"',
            "output_timers_every": "$out_every", "max_refinement_levels": "$max_levels",
            "prolongation_order_space": PROLONGATION_ORDER, "prolongation_order_time": 2}


def carpet_lib(s):
    poison = {"poison_new_memory": yes_no(s.run["poison"]), "poison_value": 114}
    if s.holes:
        return poison
    return {**poison, "support_staggered_operators": '"yes"', "print_memstats_every": "$out_every"}


def carpet_interp(s):
    return {"check_tree_search": '"no"', "tree_search": '"yes"'}


def loop_control(s):
    # LoopControl's initial random search can request an empty aligned tile range.
    return Block({"random_jump_probability": 0}, "Fixed loop tiling")


def coord_base(s):
    coord = {}
    for axis, extent in zip("xyz", s.grid["domain"]):
        coord.update({f"{axis}min": real(-extent), f"{axis}max": real(extent), f"d{axis}": real(s.grid["coarse_dx"])})
        for side in ("lower", "upper"):
            coord[f"boundary_size_{axis}_{side}"] = GHOST_SIZE
            if not s.holes:
                coord[f"boundary_shiftout_{axis}_{side}"] = 0
    return coord


def carpet_regrid2(s):
    grid = s.grid
    nesting = {"ensure_proper_nesting": "yes", "min_distance": MIN_DISTANCE,
               "freeze_unaligned_levels": "yes", "freeze_unaligned_parent_levels": "yes"}
    if s.holes:
        amr = {"regrid_every": s.horizon_every if s.tracking else 0, "num_centres": 3 if s.disk else 2, **nesting}
        if s.tracking:
            amr["snap_to_coarse"] = "yes"
        for m, (x, boxes) in enumerate(zip(s.centers, grid["hole_boxes"]), 1):
            amr.update({f"num_levels_{m}": grid["levels"], f"position_x_{m}": real(x)})
            amr.update({f"radius_{m}[{level}]": real(radius) for level, radius in enumerate(boxes, 1)})
            if s.tracking:
                amr[f"movement_threshold_{m}"] = real(grid["movement_cells"] * grid["finest_dx"])
        if s.disk:
            amr["num_levels_3"] = grid["disk_levels"]
            amr.update({f"radius_3[{level}]": real(r) for level, (r, z) in enumerate(grid["disk_boxes"], 1)})
            amr.update({f"radius_z_3[{level}]": real(z) for level, (r, z) in enumerate(grid["disk_boxes"], 1)})
        return amr
    max_levels, collapse = s.run["max_levels"], s.run["collapse"]
    amr = {"regrid_every": "$track_every" if s.bns or collapse else 0, "num_centres": len(s.centers) + s.offset,
           "snap_to_coarse": "yes", **nesting}
    if s.offset:
        amr.update(num_levels_1=grid["common_levels"], position_x_1=0)
        common_boxes = list(grid["common_boxes"])
        while len(common_boxes) < max_levels - 1:
            common_boxes.append(common_boxes[-1] / 2)
        amr.update({f"radius_1[{level}]": real(radius) for level, radius in enumerate(common_boxes, 1)})
    for c, (x, boxes) in enumerate(zip(s.centers, grid["star_boxes"]), s.offset + 1):
        if collapse and c == 1:
            # The trigger steers centre 1, here the star itself: give it the collapse radii.
            boxes = list(boxes)
            while len(boxes) < max_levels - 1:
                boxes.append(boxes[-1] / 2)
        amr.update({f"num_levels_{c}": "$num_levels", f"position_x_{c}": real(x),
                    f"movement_threshold_{c}": real(grid["movement_cells"] * grid["finest_dx"])})
        amr.update({f"radius_{c}[{level}]": real(radius) for level, radius in enumerate(boxes, 1)})
    return amr


def carpet_tracker(s):
    # The two hole boxes follow the horizons in surfaces 0 and 1.
    return {"surface[0]": 0, "surface[1]": 1} if s.holes and s.tracking else {}


def time(s):
    return {"dtfac": real(s.run["dtfac"]), "timestep_method": '"courant_static"'}


def mol(s):
    return {"ODE_Method": '"rk4"', "MoL_Intermediate_Steps": RK_STAGES, "MoL_Num_Scratch_Levels": 1}


def init_base(s):
    return {"initial_data_setup_method": '"init_all_levels"'}


# Initial data and evolution

def cocal_id(s):
    cfg, data = s.cfg, s.data
    reader = {"COCAL_ID_PathToID": quoted(data["id"]), "COCAL_ID_verbose": "no"}
    if s.holes:
        reader["COCAL_ID_bbh_fill_method"] = '"smooth"'
    elif s.bns:
        if cfg["compact"]:
            reader["COCAL_ID_bns_compact"] = "yes"
        if cfg["xunit"]:
            reader["COCAL_ID_bns_xunit"] = "yes"
        if isinstance(cfg["ecc_cor_velx"], (list, tuple)):
            for star, value in enumerate(cfg["ecc_cor_velx"]):
                reader[f"COCAL_ID_ecc_cor_velx_star[{star}]"] = real(value)
        else:
            # A scalar is the radial/xunit coefficient for both stars.
            reader["COCAL_ID_ecc_cor_velx"] = real(cfg["ecc_cor_velx"])
        velz = velocity_pair(cfg["ecc_cor_velz"], "ecc_cor_velz")
        if any(velz):
            reader.update({f"COCAL_ID_ecc_cor_velz_star[{star}]": real(value) for star, value in enumerate(velz)})
        reader.update({f"COCAL_ID_cen_x{star}": real(x) for star, x in enumerate(s.centers, 1)})
    else:
        reader["COCAL_ID_rnstype"] = quoted(data["rnstype"])
    if data["EOS_type"] == "Tabulated":
        reader["COCAL_ID_read_tabulated_eos"] = "yes"
    return Block(reader, "COCAL initial data")


def adm_base(s):
    evolution = quoted(s.evolution)
    initial = {"rns": "CocalRNS", "bns": "CocalBNS", "bbh": "CocalBBH", "bbht": "CocalBBHT"}[s.kind]
    return {"metric_type": '"physical"', "evolution_method": evolution, "lapse_evolution_method": evolution,
            "shift_evolution_method": evolution, "dtlapse_evolution_method": evolution,
            "dtshift_evolution_method": evolution, "initial_data": quoted(initial), "initial_lapse": '"Cocal"',
            "initial_shift": '"zero"' if s.bns else '"Cocal"', "initial_dtlapse": '"zero"',
            "initial_dtshift": '"zero"', "metric_timelevels": 3, "lapse_timelevels": 3, "shift_timelevels": 3}


def hydro_base(s):
    grhydro = {"evolution_method": '"GRHydro"', "prolongation_type": '"ENO"'} if s.hydro == "GRHydro" else {}
    if s.holes:
        if not s.disk:
            return {"initial_hydro": '"zero"', **grhydro}
        return {"initial_hydro": '"Cocal"', "timelevels": 3, "initial_Bvec": '"zero"', "initial_Avec": '"zero"',
                "initial_Aphi": '"zero"', **grhydro}
    thermal = '"GRHayLID"' if s.table else '"none"'
    return {"timelevels": 3, "initial_hydro": '"Cocal"', "initial_Y_e": thermal, "initial_temperature": thermal,
            "initial_entropy": thermal, "initial_Avec": '"zero"', "initial_Aphi": '"zero"',
            "initial_Bvec": '"none"' if s.hydro == "GRHayLHD" else '"zero"', **grhydro}


def adm_macros(s):
    return {"spatial_order": 4}


def space_mask(s):
    return {"use_mask": '"yes"'}


def tmunu_base(s):
    storage = yes_no(s.evolution != "BaikalVacuum")
    tmunu = {"stress_energy_storage": storage, "stress_energy_at_RHS": storage, "prolongation_type": '"none"'}
    if s.kind != "bbh":
        tmunu["timelevels"] = 3
    return tmunu


def spacetime(s):
    return spacetime_parameters(s.evolution, s.cfg["gauge"], s.cfg["spacetime_options"], black_holes=s.holes)


def dissipation(s):
    if s.evolution != "LeanBSSNMoL":
        return {}
    return {"order": 5, "epsdis": "0.2", "vars": quoted(SPACETIME_THORNS[s.evolution]["evolved"])}


def hydro_solver(s):
    if s.hydro == "GRHydro":
        return grhydro_parameters(s.data)
    return {"Symmetry": '"none"', "Convert_to_HydroBase_every": s.run["hydro_every"] if s.holes else "$hydro_every"}


def grhayl_lib(s):
    """The EOS, con2prim and limits shared by GRHayLHD and IllinoisGRMHD."""
    data, table = s.data, s.table
    con2prim = {"con2prim_routine": '"Palenzuela1D"', "con2prim_backup_routines[0]": '"Newman1D"',
                "calc_primitive_guess": "yes"}
    if table:
        atmosphere = TABULATED_ATMOSPHERE
        return Block({
            "EOS_type": '"Tabulated"', "EOS_tablepath": quoted(s.cfg["eos_table"]),
            "evolve_temperature": "yes", "evolve_entropy": "yes", **con2prim,
            "con2prim_backup_routines[1]": '"Palenzuela1D_entropy"',
            "con2prim_backup_routines[2]": '"Newman1D_entropy"',
            **grhayl_limits(data, table),
            "Y_e_atm": real(atmosphere["Y_e"]), "Y_e_min": real(table["Y_e"][0]), "Y_e_max": real(table["Y_e"][1]),
            "T_atm": real(atmosphere["T"]), "T_min": real(table["T"][0]), "T_max": real(table["T"][1]),
        }, "Tabulated EOS")
    return Block({
        "EOS_type": '"Hybrid"', "neos": len(data["EOS_gamma"]), "gamma_th": real(thermal_gamma(data)),
        "k_ppoly0": real(data["EOS_k0"]), "gamma_ppoly_in": array(data["EOS_gamma"]),
        "rho_ppoly_in": array(data["EOS_boundaries"] or [0]), "evolve_temperature": "no",
        **con2prim, **grhayl_limits(data),
    }, "EOS from the native COCAL input")


def grhayl_id(s):
    if not s.table:
        return {}
    return {"impose_beta_equilibrium": "yes", "beq_temperature": real(TABULATED_ATMOSPHERE["T"]),
            "initialize_magnetic_quantities": "no"}


def eos_omni(s):
    # GRHydro's EOS handle; a tabulated ID activates EOS_Omni without settings.
    return hybrid_eos_parameters(s.data) if s.hydro == "GRHydro" else {}


# Surfaces and tracking

def spherical_surface(s):
    """The angular grid and prescribed sphere geometry of every surface slot."""
    surfaces = {"nsurfaces": len(s.surfaces), "maxntheta": 39, "maxnphi": 76}
    for slot, spec in enumerate(s.surfaces):
        surfaces.update({f"ntheta[{slot}]": 39, f"nphi[{slot}]": 76, f"nghoststheta[{slot}]": 2,
                         f"nghostsphi[{slot}]": 2, f"name[{slot}]": quoted(spec["name"]),
                         f"origin_x[{slot}]": real(spec["x"]),
                         f"set_spherical[{slot}]": yes_no(spec["radius"] is not None)})
        if spec["radius"] is not None:
            surfaces[f"radius[{slot}]"] = real(spec["radius"])
    return surfaces


def quasi_local_measures(s):
    # The same slots as SphericalSurface; unused slots stay disabled.
    qlm = {"num_surfaces": len(s.surfaces), "spatial_order": 4, "interpolator_options": '"order=4"'}
    if s.holes:
        qlm["killing_vector_method"] = '"axial"'
    else:
        qlm["interpolator"] = '"Lagrange polynomial interpolation"'
    qlm.update({f"surface_index[{slot}]": slot if spec["enabled"] else -1 for slot, spec in enumerate(s.surfaces)})
    return qlm


def ns_surface(s):
    if s.holes:
        return {}
    surface = {"compute_every": "$out_every"}
    for slot, radius in zip(s.contour_slots, s.region):
        surface.update({f"find_surface[{slot}]": "yes", f"rho_surface[{slot}]": real(s.cfg["surface_density"]),
                        f"max_radius[{slot}]": real(radius)})
    return surface


def ahfinder_direct(s):
    interpolators = {"surface_interpolator_name": '"Lagrange polynomial interpolation"',
                     "surface_interpolator_pars": '"order=4"'}
    common = {"output_BH_diagnostics": "yes", "output_initial_guess": "no", "output_Theta_every": 0,
              "output_mean_curvature_every": 0, "run_at_CCTK_POST_RECOVER_VARIABLES": "no", **interpolators}
    if s.holes:
        data, run = s.data, s.run
        ah = {"N_horizons": 3, "find_every": s.find_every, "output_h_every": horizon_shape_every(run, s.find_every),
              "run_at_CCTK_ANALYSIS": "yes", "run_at_CCTK_POSTSTEP": "no", **common,
              # AHFinderDirect requires C1 geometry interpolation for Newton convergence.
              "geometry_interpolator_name": '"Hermite polynomial interpolation"',
              "geometry_interpolator_pars": '"order=2"',
              "max_Newton_iterations__initial": 40, "verbose_level": '"algorithm highlights"'}
        if s.tracking:
            ah.update(move_origins="yes", reshape_while_moving="yes", predict_origin_movement="yes")
        radii = data["bh_excision_radii"]
        for m, (x, r) in enumerate(zip(s.centers, radii), 1):
            ah.update({f"origin_x[{m}]": real(x), f"initial_guess_method[{m}]": '"coordinate sphere"',
                       f"initial_guess__coord_sphere__x_center[{m}]": real(x),
                       f"initial_guess__coord_sphere__radius[{m}]": real(1.1 * r),
                       f"max_allowable_horizon_radius[{m}]": real(4 * r),
                       f"which_surface_to_store_info[{m}]": m - 1})
        # The common horizon: guess sphere twice the summed individual guesses, as in the
        # ET qc0 parfiles, searched once per coarse step after s.common_after.
        separation = abs(s.centers[1] - s.centers[0])
        ah.update({"initial_guess_method[3]": '"coordinate sphere"',
                   "find_after_individual_time[3]": real(s.common_after),
                   "find_every_individual[3]": run["track_every"] if s.tracking else run["out_every"],
                   "initial_guess__coord_sphere__radius[3]": real(2 * 1.1 * sum(radii)),
                   "max_allowable_horizon_radius[3]": real(2 * (sum(radii) + separation))})
        return ah
    return {"N_horizons": 1, "find_every": "$track_every*4",  # a forming horizon, every 4 coarse steps
            "output_h_every": horizon_shape_every(s.run, 4 * s.run["track_every"]), **common,
            "geometry_interpolator_name": '"Lagrange polynomial interpolation"', "geometry_interpolator_pars": '"order=4"',
            "initial_guess_method[1]": '"coordinate sphere"', "origin_x[1]": "0.0", "origin_y[1]": "0.0",
            "origin_z[1]": "0.0", "initial_guess__coord_sphere__x_center[1]": "0.0",
            "initial_guess__coord_sphere__y_center[1]": "0.0", "initial_guess__coord_sphere__z_center[1]": "0.0",
            "initial_guess__coord_sphere__radius[1]": "2.0", "max_allowable_horizon_radius[1]": "20.0",
            "which_surface_to_store_info_by_name[1]": '"apparent horizon"'}


def reduce_loc_tracker(s):
    if s.holes:
        return {}
    tracker = {"num_tracked_vars": len(s.centers), "compute_every": "$track_every"}
    for star, (x, radius) in enumerate(zip(s.centers, s.region)):
        tracker.update({f"tracked_var[{star}]": '"HydroBase::rho"', f"selection_region[{star}]": '"follow_previous"',
                        f"target_x0[{star}]": real(x), f"target_radius[{star}]": real(radius)})
    return tracker


# Volume integrals

def volume_integrals_grmhd(s):
    output = {"enable_file_output": 1, "verbose": 0}
    if s.holes:
        return {"NumIntegrals": 1, "VolIntegral_out_every": s.run["out_every"], **output,
                "Integration_quantity_keyword[1]": '"restmass"'}
    if not s.bns:
        return Block({"NumIntegrals": 1, "VolIntegral_out_every": "$out_every", **output,
                      "Integration_quantity_keyword[1]": '"restmass"'}, labels={1: "Whole-domain rest mass"})
    # Fixed slots, so the integral columns are the same for every run.
    matter = {"NumIntegrals": 8, "VolIntegral_out_every": "$track_every", **output,
              "Integration_quantity_keyword[1]": '"centerofmass"', "Integration_quantity_keyword[6]": '"restmass"'}
    labels = {1: "Whole-domain CoM", 6: "Whole-domain rest mass"}
    for star, (x, radius, com, mass) in enumerate(zip(s.centers, s.region, (2, 4), (7, 8))):
        name, centre = ("Lhs" if star == 0 else "Rhs"), star + s.offset
        labels[com] = f"{name} CoM"
        labels[com + 1] = f"{name} contour origin (reuse CoM)" if "NSSurface" in s.active else "Reserved contour slot"
        labels[mass] = f"{name} rest mass"
        matter.update({f"Integration_quantity_keyword[{com}]": '"centerofmass"',
                       f"volintegral_sphere__center_x_initial[{com}]": real(x),
                       f"volintegral_inside_sphere__radius[{com}]": real(radius),
                       f"volintegral_sphere__tracks__amr_centre[{com}]": centre,
                       f"amr_centre__tracks__volintegral_inside_sphere[{com}]": centre,
                       f"Integration_quantity_keyword[{com + 1}]": '"one"'})
        if "SphericalSurface" in s.active:
            matter[f"sphericalsurface__tracks__volintegral_inside_sphere[{com}]"] = star
        if "NSSurface" in s.active:
            # Reuse the same masked CoM integrands, without applying a second mask.
            matter.update({f"Integration_quantity_keyword[{com + 1}]": '"usepreviousintegrands"',
                           f"volintegral_usepreviousintegrands_num_integrands[{com + 1}]": 4,
                           f"sphericalsurface__tracks__volintegral_inside_sphere[{com + 1}]": s.contour_slots[star]})
        matter.update({f"Integration_quantity_keyword[{mass}]": '"restmass"',
                       f"volintegral_sphere__center_x_initial[{mass}]": real(x),
                       f"volintegral_inside_sphere__radius[{mass}]": real(radius),
                       f"volintegral_sphere__tracks__amr_centre[{mass}]": centre})
    return Block(matter, labels=labels)


def volume_integrals_vacuum(s):
    names = ("HamiltonianVarString", "Momentum0VarString", "Momentum1VarString", "Momentum2VarString")
    vacuum = {"enable_file_output": 1, "verbose": 0,
              **{name: quoted(var) for name, var in zip(names, s.constraint_vars)}}
    if s.holes:
        data, outer = s.data, s.outer_spheres[0]
        # Constraints and volume inside outer sphere 0, then outside each hole's
        # fill (reusing the integrands), then ADM mass and J and the constraints.
        vacuum.update(NumIntegrals=9, VolIntegral_out_every=s.run["out_every"])
        for slot, kind in [(1, "H_M_CnstraintsL2"), (2, "usepreviousintegrands"), (3, "usepreviousintegrands"),
                           (4, "one"), (5, "usepreviousintegrands"), (6, "usepreviousintegrands"),
                           (7, "ADM_Mass"), (8, "ADM_Angular_Momentum"), (9, "H_M_CnstraintsL2")]:
            vacuum[f"Integration_quantity_keyword[{slot}]"] = quoted(kind)
            if slot in (1, 4, 7, 8, 9):
                vacuum[f"volintegral_inside_sphere__radius[{slot}]"] = real(outer)
            if slot in (2, 3, 5, 6):
                hole = 0 if slot in (2, 5) else 1
                vacuum[f"volintegral_usepreviousintegrands_num_integrands[{slot}]"] = 4 if slot < 4 else 1
                vacuum[f"volintegral_sphere__center_x_initial[{slot}]"] = real(s.centers[hole])
                vacuum[f"volintegral_outside_sphere__radius[{slot}]"] = real(max(.45, 3 * data["bh_excision_radii"][hole]))
                if s.tracking:
                    vacuum[f"volintegral_sphere__tracks__amr_centre[{slot}]"] = hole
        return vacuum
    vacuum.update(NumIntegrals=len(s.vacuum_integrals), VolIntegral_out_every="$out_every")
    labels = {}
    for integral in s.vacuum_integrals:
        slot, radius, star = integral["slot"], integral["radius"], integral["star"]
        vacuum[f"Integration_quantity_keyword[{slot}]"] = quoted(integral["quantity"])
        vacuum[f"volintegral_inside_sphere__radius[{slot}]"] = real(radius)
        region = "Whole domain" if radius < 0 else f"r={radius:g}"
        if star is not None:
            vacuum[f"volintegral_sphere__center_x_initial[{slot}]"] = real(s.centers[star])
            vacuum[f"volintegral_sphere__tracks__amr_centre[{slot}]"] = star + s.offset
            region = ("Lhs" if star == 0 else "Rhs") + " " + region
        labels[slot] = f"{integral['quantity']}, {region}"
    return Block(vacuum, labels=labels)


# Wave extraction and field checks

PSI4 = '"WeylScal4::Psi4r{sw=-2 cmplx=\'WeylScal4::Psi4i\' name=\'psi4\'}"'


def multipole(s):
    if s.holes:
        return {"nradii": 1, "radius[0]": real(s.outer_spheres[0]), "out_every": s.run["out_every"], "l_max": 4,
                "ntheta": 32, "nphi": 64, "output_ascii": "yes", "output_hdf5": "no", "variables": PSI4}
    radii = [radius for radius in s.global_radii if radius is not None]
    return {"out_dir": '"./data/Multipole"', "l_max": 4, "nphi": 240, "ntheta": 120, "out_every": "$out_every",
            "output_ascii": "yes", "output_hdf5": "yes", "variables": PSI4,
            **{f"radius[{slot}]": real(radius) for slot, radius in enumerate(radii)}, "nradii": len(radii)}


def nan_checker(s):
    evolved = SPACETIME_THORNS[s.evolution]["evolved"]
    if s.holes:
        return {"check_every": s.run["nan_check_every"], "action_if_found": '"abort"',
                "check_vars": variables("ADMBase::metric ADMBase::curv ADMBase::lapse ADMBase::shift", evolved,
                                        *s.stress, *s.matter_outputs,
                                        *([s.hydro + "::grmhd_velocities"] if s.hydro == "IllinoisGRMHD" else []))}
    return Block({
        "check_every": "$nan_check_every", "action_if_found": '"abort"',
        "check_vars": variables("ADMBase::metric ADMBase::lapse ADMBase::shift",
                                "HydroBase::rho HydroBase::eps HydroBase::press", s.native, s.constraints,
                                *([evolved] if s.evolution != "Baikal" else [])),
    }, "Check native velocities before the ANALYSIS HydroBase conversion" if s.hydro != "GRHydro"
       else "Check native evolved fields")


# Output and checkpointing

ADM_OUTPUT = "ADMBase::metric ADMBase::curv ADMBase::lapse ADMBase::shift"


def carpet_io_basic(s):
    if s.holes:
        return {"outInfo_every": s.run["out_every"], "outInfo_reductions": '"minimum maximum"',
                "outInfo_vars": '"ADMBase::alp"'}
    return {"outInfo_every": "$out_every",
            "outInfo_vars": ('"Carpet::physical_time_per_hour HydroBase::rho{reductions=\'maximum\'} '
                             + s.constraint_vars[0] + '{reductions=\'maximum\'}"')}


def carpet_io_scalar(s):
    files = {"one_file_per_group": '"yes"', "all_reductions_in_one_file": '"yes"'}
    if s.holes:
        return {"outScalar_every": s.run["out_every"], **files, "outScalar_reductions": '"minimum maximum norm2"',
                "outScalar_vars": variables(ADM_OUTPUT, s.constraints, *s.stress, *s.matter_outputs)}
    return {"outScalar_every": "$out_every", **files,
            "outScalar_reductions": '"norm1 norm2 norm_inf sum maximum minimum"',
            "outScalar_vars": variables("HydroBase::rho HydroBase::press HydroBase::vel HydroBase::eps",
                                        "ADMBase::lapse ADMBase::shift", s.constraints,
                                        " ".join(filter(None, (s.failure, "Carpet::timing{reductions='average'}"))))}


def carpet_io_ascii(s):
    active = s.active
    if s.holes:
        scalars = []
        if "QuasiLocalMeasures" in active:
            scalars.append("QuasiLocalMeasures::qlm_state QuasiLocalMeasures::qlm_scalars")
        if "SphericalSurface" in active:
            scalars.append("SphericalSurface::sf_active SphericalSurface::sf_valid SphericalSurface::sf_origin")
        if s.tracking:
            scalars += ["SphericalSurface::sf_info", "CarpetRegrid2::positions CarpetRegrid2::active"]
        return {"out0D_every": s.run["out_every"], "out0D_vars": variables(*scalars),
                "out1D_every": s.run["out1d_every"], "one_file_per_group": "yes",
                "out1D_vars": variables(ADM_OUTPUT, s.constraints, *s.matter_outputs)}
    scalars = ["Carpet::timing"]
    if "ReduceLocTracker" in active:
        scalars.append("ReduceLocTracker::reduce_loc_reals")
    if "QuasiLocalMeasures" in active:
        scalars += ["QuasiLocalMeasures::qlm_state", "QuasiLocalMeasures::qlm_scalars",
                    "QuasiLocalMeasures::qlm_multipole_moments", "QuasiLocalMeasures::qlm_grid_real"]
    if "SphericalSurface" in active:
        scalars += ["SphericalSurface::sf_active", "SphericalSurface::sf_valid", "SphericalSurface::sf_origin",
                    "SphericalSurface::sf_info", "SphericalSurface::sf_coordinate_descriptors"]
    if "NSSurface" in active:
        scalars.append("NSSurface::info_valid")
    return {"out0D_every": "$out_every", "out0D_vars": variables(*scalars), "out1D_every": "$out_every1d",
            "one_file_per_group": "yes", "output_symmetry_points": "no",
            "out1D_vars": '"\n' + "".join(f" {name}\n" for name in (
                "HydroBase::rho", "HydroBase::vel", "HydroBase::press", "ADMBase::lapse", "ADMBase::metric",
                "ADMBase::curv", *s.constraint_vars[:2])) + '"',
            "out2D_every": "$out_every" if "NSSurface" in active else -1, "out2D_vars": shape_output(s)}


def shape_output(s):
    return '"SphericalSurface::sf_radius"' if "NSSurface" in s.active else '""'


def carpet_io_hdf5(s):
    run = s.run
    checkpoint = yes_no(run["checkpoint_hours"] > 0 or run["checkpoint_on_terminate"])
    if s.holes:
        matter = s.matter_outputs
        return {"out2D_every": run["out2d_every"], "out2D_vars": variables("ADMBase::lapse", s.constraints, *matter),
                "out3D_every": run["out3d_every"], "out3D_vars": variables("ADMBase::lapse ADMBase::metric", *matter[:1]),
                "checkpoint": checkpoint}
    return {"checkpoint": checkpoint, "out2D_every": "$out_every2d",
            "out2D_vars": variables("HydroBase::rho HydroBase::vel HydroBase::press",
                                    "ADMBase::lapse ADMBase::shift ADMBase::metric", s.constraints),
            "out3D_every": "$out_every3d", "out3D_vars": '"\n HydroBase::rho\n ADMBase::lapse\n "',
            "out_vars": shape_output(s), "out_every": "$out_every", "output_buffer_points": "no"}


def io_util(s):
    run = s.run
    io = {"out_dir": "${parfile}", "checkpoint_dir": '"checkpoints"', "checkpoint_ID": "no", "checkpoint_every": -1,
          "checkpoint_every_walltime_hours": run["checkpoint_hours"],
          "checkpoint_on_terminate": yes_no(run["checkpoint_on_terminate"]), "checkpoint_keep": 2,
          "recover_dir": '"checkpoints"'}
    if s.holes:
        return {**io, "parfile_write": '"no"', "recover": '"autoprobe"'}
    return {**io, "out_mode": '"onefile"', "recover": '"no"' if s.local else '"autoprobe"'}


# Collapse trigger and thread topology

def trigger(s):
    if s.holes:
        return {}
    return {"Trigger_Number": 1 if s.run["collapse"] else 0, "Trigger_Checked_Variable[0]": '"ADMBase::alp"',
            "Trigger_Reduction[0]": '"minimum"', "Trigger_Relation[0]": '"<"', "Trigger_Checked_Value[0]": "0.1",
            "Trigger_Reaction[0]": '"steerscalar"', "Trigger_Steered_Scalar[0]": '"CarpetRegrid2::num_levels[0]"',
            "Trigger_Steered_Scalar_Value[0]": '"$max_levels"', "Trigger_Once[0]": 1}


def system_topology(s):
    return {} if s.holes else {"set_thread_bindings": '"no"'}


# Every thorn with parameters, in output order.
BLOCKS = {
    "Cactus": cactus, "TerminationTrigger": termination_trigger, "TimerReport": timer_report,
    "CartGrid3D": cart_grid_3d, "Carpet": carpet, "CarpetLib": carpet_lib, "CarpetInterp": carpet_interp,
    "LoopControl": loop_control, "CoordBase": coord_base, "CarpetRegrid2": carpet_regrid2,
    "CarpetTracker": carpet_tracker, "Time": time, "MoL": mol, "InitBase": init_base,
    "COCAL_ID": cocal_id, "ADMBase": adm_base, "HydroBase": hydro_base, "ADMMacros": adm_macros,
    "SpaceMask": space_mask, "TmunuBase": tmunu_base,
    **{thorn: spacetime for thorn in SPACETIME_THORNS}, "Dissipation": dissipation,
    **{thorn: hydro_solver for thorn in HYDRO_THORNS}, "GRHayLib": grhayl_lib, "GRHayLID": grhayl_id,
    "EOS_Omni": eos_omni,
    "SphericalSurface": spherical_surface, "QuasiLocalMeasures": quasi_local_measures, "NSSurface": ns_surface,
    "AHFinderDirect": ahfinder_direct, "ReduceLocTracker": reduce_loc_tracker,
    "VolumeIntegrals_GRMHD": volume_integrals_grmhd, "VolumeIntegrals_vacuum": volume_integrals_vacuum,
    "Multipole": multipole, "NaNChecker": nan_checker,
    "CarpetIOBasic": carpet_io_basic, "CarpetIOScalar": carpet_io_scalar, "CarpetIOASCII": carpet_io_ascii,
    "CarpetIOHDF5": carpet_io_hdf5, "IOUtil": io_util, "Trigger": trigger, "SystemTopology": system_topology,
}


def thorn_blocks(s):
    """Each active thorn's block, from its function in BLOCKS."""
    blocks = {}
    for thorn, function in BLOCKS.items():
        if thorn == "Cactus" or thorn in s.active:
            block = function(s)
            block = block if isinstance(block, Block) else Block(block)
            blocks[thorn] = dict(parameters=block.parameters, title=block.title or thorn, labels=block.labels or {})
    return blocks


def render_parfile(blocks, active, macros, title, profile):
    """Assemble selected blocks and check their structure, not Cactus expressions.

    Cactus -P handles expression syntax and compiled parameter names/types/ranges.
    The generator only checks errors it can identify without parsing that language.
    """
    sections = [
        f"# {title}\n"
        f"# Generated by generate_parfile.py; settings: {profile}.\n",
        'ActiveThorns = "\n  ' + '\n  '.join(active) + '\n"',
    ]
    if macros:
        sections.append("# Resolution and output cadence\n"
                        + "\n".join(f"${name} = {value}" for name, value in macros.items()))
    emitted = set()
    # These Cactus implementations are provided by the named Carpet IO thorns.
    namespaces = {"CarpetIOBasic": "IOBasic", "CarpetIOScalar": "IOScalar"}

    for thorn, block in blocks.items():
        namespace = namespaces.get(thorn, thorn)
        parameters = {key: str(value) for key, value in block["parameters"].items()
                      if value is not None}
        if not parameters:
            continue

        labels = block["labels"]
        if labels:
            def slot_number(item):
                match = re.search(r"\[(\d+)\]$", item[0])
                return int(match[1]) if match else 0
            parameters = dict(sorted(parameters.items(), key=slot_number))

        lines = []
        width = max(len(key) for key in parameters)
        previous_slot = None
        for key, value in parameters.items():
            # An index may itself be a Cactus expression. Check only the base name.
            base_name = key.partition("[")[0]
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", base_name):
                raise ValueError(f"Invalid parameter name in {thorn} block: {key}")
            full_name = f"{namespace}::{key}"
            if full_name.lower() in emitted:
                raise ValueError(f"Duplicate generated parameter: {full_name}")
            emitted.add(full_name.lower())

            if not value.strip():
                raise ValueError(f"Empty value for {full_name}")

            if labels or thorn in ("SphericalSurface", "ReduceLocTracker", "NSSurface"):
                match = re.search(r"\[(\d+)\]$", key)
                slot = match[1] if match else None
                if slot != previous_slot and lines:
                    lines.append("")
                    if labels and slot is not None:
                        lines.append(f"# Integral {slot}: {labels[int(slot)]}")
                previous_slot = slot
            elif thorn == "CarpetRegrid2" and key.startswith("num_levels_"):
                lines.append("")
            lines.append(f"{namespace}::{key:<{width}} = {value}")
        sections.append(f"# {block['title']}\n" + "\n".join(lines))

    return "\n\n".join(sections) + "\n"


def generate(id_dir, overrides=None):
    """Read the ID, resolve the setup, assemble thorn blocks, and record provenance."""
    directory = Path(id_dir).expanduser()
    if (directory / "bhtpar.dat").is_file():
        raise ValueError("BHT uses the dedicated par/bht recipes; this generator supports RNS/BNS/BBH/BBHT")
    if (directory / "bbhgrids_3D_mpt1.las").is_file():
        resolved = resolve_black_holes(id_dir, overrides, disk=(directory / "disk_state.las").is_file())
    else:
        resolved = resolve_stellar(id_dir, overrides)
    cfg, data, grid, run, requested, active = resolved
    s = Setup(cfg, data, grid, run, active)
    # Stellar parfiles name their resolution and cadences; black-hole parfiles use values.
    macros = {} if s.holes else {
        "num_levels": grid["levels"], "max_levels": run["max_levels"], "track_every": run["track_every"],
        "out_every": run["out_every"], "out_every1d": run["out1d_every"], "out_every2d": run["out2d_every"],
        "out_every3d": run["out3d_every"], "nan_check_every": run["nan_check_every"],
        "hydro_every": run["hydro_every"],
    }
    hydro = f"{cfg['hydro_thorn']} + " if cfg["hydro_thorn"] else ""
    title = f"Carpet + {hydro}{cfg['spacetime_thorn']}: COCAL {data['kind'].upper()}"
    text = render_parfile(thorn_blocks(s), active, macros, title, cfg["profile"])

    metadata = dict(
        data, profile=cfg["profile"], config_name=cfg["name"],
        spacetime_thorn=cfg["spacetime_thorn"], hydro_thorn=cfg["hydro_thorn"],
        gauge=cfg["gauge"], spacetime_options=cfg["spacetime_options"],
        requested_thorns=requested, active_thorns=active,
        added_thorns=[name for name in active if name not in requested],
        grid=grid, run=run,
    )
    if s.holes:
        tracking = s.tracking
        metadata.update(
            common_horizon_after=s.common_after,
            tracking=dict(enabled=tracking, every=s.find_every,
                          movement_threshold=grid["movement_cells"] * grid["finest_dx"] if tracking else None),
            outer_qlm_radii=s.outer_spheres, exterior_constraint_slots=[1, 2, 3], exterior_volume_slots=[4, 5, 6])
    else:
        metadata.update(
            surface_density=cfg["surface_density"], compact=cfg["compact"], xunit=cfg["xunit"],
            ecc_cor_velx=cfg["ecc_cor_velx"], ecc_cor_velz=cfg["ecc_cor_velz"],
            eos_table=str(cfg["eos_table"]) if cfg["eos_table"] else None,
            fine_box_radii=grid["fine_box_radii"],
            tracking=dict(every=run["track_every"],
                          movement_threshold=grid["movement_cells"] * grid["finest_dx"] if s.bns else None),
            stellar_surface_slots=list(range(len(s.centers))) if "SphericalSurface" in active else [],
            contour_slots=s.contour_slots if "NSSurface" in active else [],
            finite_radius_ADM_volumes=s.finite_adm if "VolumeIntegrals_vacuum" in active else [],
            vacuum_integrals=s.vacuum_integrals if "VolumeIntegrals_vacuum" in active else [],
            outer_qlm_radii=s.global_radii if "QuasiLocalMeasures" in active else [])
    metadata.update(generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    parfile_sha256=hashlib.sha256(text.encode()).hexdigest())
    return text, metadata


class CocalMetadata:
    """Read native ID metadata without loading the large fluid grids."""

    def __init__(self, directory):
        self.directory = Path(directory).resolve()

    def _numbers(self, line):
        return [float(word.replace("D", "e").replace("d", "e"))
                for word in line.split(":", 1)[0].split("#", 1)[0].split()]

    def _read_eos(self, path):
        """Use the same units and phase ordering as COCAL's peos_initialize."""
        lines = path.read_text().splitlines()
        phases = int(self._numbers(lines[0])[0])
        if not 1 <= phases <= 10:
            raise ValueError(f"{path}: expected a piecewise-polytropic EOS")
        rho0, pressure0 = self._numbers(lines[1])[:2]
        table = [self._numbers(line)[:2] for line in lines[2:phases + 3]][::-1]
        if len(table) != phases + 1 or any(len(row) != 2 for row in table):
            raise ValueError(f"{path}: incomplete EOS table")
        rho, gamma = zip(*table)
        if (not all(math.isfinite(x) and x > 0 for x in (*rho, rho0, pressure0))
                or not all(math.isfinite(x) and x > 1 for x in gamma)
                or any(a >= b for a, b in zip(rho, rho[1:]))):
            raise ValueError(f"{path}: invalid density boundaries or adiabatic indices")
        phase = next((i for i in range(1, phases + 1)
                      if rho[i - 1] <= rho0 <= rho[i]), None)
        if phase is None:
            raise ValueError(f"{path}: EOS reference density is outside the table")
        g, c, solar_mass = G_CGS, C_CGS, MSUN_CGS
        density_unit = DENSITY_UNIT
        pressure_unit = g**3 * solar_mass**2 / c**8
        k0 = pressure0 / rho0**gamma[phase] * pressure_unit / density_unit**gamma[phase]
        for i in range(phase - 1, 0, -1):
            k0 *= (rho[i] * density_unit)**(gamma[i + 1] - gamma[i])
        return dict(EOS_gamma=list(gamma[1:]), EOS_k0=k0,
                    EOS_boundaries=[x * density_unit for x in rho[1:-1]])

    def _fluid_constants(self, path):
        # Only the trailer is needed; do not load the large fluid grid.
        with path.open("rb") as stream:
            shape = [int(word) for word in stream.readline().split()]
            if len(shape) != 3 or any(size < 1 for size in shape):
                raise ValueError(f"{path}: expected three positive grid dimensions")
            start = stream.tell()
            record_width = len(stream.readline())
            stream.seek(start + math.prod(size + 1 for size in shape) * record_width)
            values = self._numbers(stream.readline().decode("ascii"))
        if len(values) < 3 or not all(math.isfinite(value) for value in values[:3]) or values[2] <= 0:
            raise ValueError(f"{path}: missing fluid constants after the grid records")
        return values[:3]

    def _surface_radius(self, path):
        with path.open() as stream:
            nt, np = map(int, stream.readline().split())
            values = [self._numbers(line)[0] for line in stream if line.strip()]
        if len(values) != (nt + 1) * (np + 1):
            raise ValueError(f"{path}: surface dimensions do not match its records")
        if not all(math.isfinite(x) and x > 0 for x in values):
            raise ValueError(f"{path}: surface radii must be finite and positive")
        return max(values)

    def _is_table(self, path):
        """A COCAL EOS table starts with its row count (teos_initialize reads i5, then the
        rows); polytrope exports may carry a placeholder file without one."""
        if not path.is_file():
            return False
        lines = [line for line in path.read_text().splitlines() if line.strip()]
        words = lines[0].split() if lines else []
        return len(words) == 1 and words[0].isdigit() and int(words[0]) == len(lines) - 1 > 1

    def _rns_peak_density(self, directory):
        """rho_max (off-centre peaks) or else rho_c, in G = c = Msol = 1 units."""
        path = next((directory / name for name in ("rnsphyseq.dat", "rnsphyseq_v1.dat")
                     if (directory / name).is_file()), None)
        if path is None:
            raise ValueError(f"{directory}: no rnsphyseq.dat or rnsphyseq_v1.dat")
        values = {}
        for line in path.read_text().splitlines():
            key, _, value = line.partition("=")
            if key.strip() in ("rho_max", "rho_c"):
                values.setdefault(key.strip(), self._numbers(value)[0])
        if not values:
            raise ValueError(f"{path}: no rho_max or rho_c")
        return values.get("rho_max", values.get("rho_c"))

    def read(self):
        """Read unmagnetized RNS or BNS data; refuse ambiguous metadata."""
        directory = self.directory
        if (directory / "bnsflu_3D_mpt1.las").exists():
            teos1 = directory / "teos_parameter_mpt1.dat"
            teos2 = directory / "teos_parameter_mpt2.dat"
            if self._is_table(teos1) or self._is_table(teos2):
                if not self._is_table(teos1) or not self._is_table(teos2) or teos1.read_bytes() != teos2.read_bytes():
                    raise ValueError("Both BNS stars need identical tabulated EOS inputs")
                eos = {"EOS_type": "Tabulated"}
            else:
                eos = self._read_eos(directory / "peos_parameter_mpt1.dat")
                other = self._read_eos(directory / "peos_parameter_mpt2.dat")
                for key in eos:
                    a = eos[key] if isinstance(eos[key], list) else [eos[key]]
                    b = other[key] if isinstance(other[key], list) else [other[key]]
                    if len(a) != len(b) or not all(math.isclose(x, y, rel_tol=1e-12)
                                                  for x, y in zip(a, b)):
                        raise ValueError("GRHayLHD needs the same EOS for both stars")
                eos["EOS_type"] = "Hybrid"
            phy = {}
            for line in (directory / "bnsphyseq.dat").read_text().splitlines():
                if "=" in line:
                    key, value = line.split("=", 1)
                    if key.strip() in ("R0", "d*R0", "restmass", "admmass_asymp", "rho_max"):
                        phy[key.strip()] = self._numbers(value.split("nc:", 1)[0])
            radii, scales, distances = [], [], []
            for patch in (1, 2):
                par = (directory / f"rnspar_mpt{patch}.dat").read_text().splitlines()
                kind = par[2].split()[2]
                if kind not in ("CO", "IR", "SP"):
                    raise ValueError(f"Unsupported BNS fluid type {kind}")
                omega, _, scale = self._fluid_constants(
                    directory / f"bnsflu_3D_mpt{patch}.las")
                if not math.isclose(scale, phy["R0"][0], rel_tol=1e-10):
                    raise ValueError("bnsphyseq.dat and the fluid export disagree on R0")
                # COCAL derives dis from its excision-grid radius, not rgmid.
                grid_data = (directory / f"grid_data_mpt{patch}.dat").read_text().splitlines()
                distance_line = next(i + 1 for i, line in enumerate(grid_data) if "dis<" in line)
                distance = self._numbers(grid_data[distance_line])[1]
                if not math.isclose(distance * scale, phy["d*R0"][patch - 1], rel_tol=1e-5):
                    raise ValueError("BNS patch centres disagree with bnsphyseq.dat")
                rsurf = self._numbers((directory / f"r_surf_mpt{patch}.dat").read_text())[0]
                radii.append(self._surface_radius(directory / f"bnssur_3D_mpt{patch}.las") * rsurf * scale)
                scales.append(scale)
                distances.append(distance)
            if not math.isclose(scales[0], scales[1], rel_tol=1e-10):
                raise ValueError("BNS patches use different coordinate scales")
            if not math.isclose(distances[0], distances[1], rel_tol=1e-8):
                raise ValueError("BNS patches use different midpoint distances")
            centers = [-distances[0] * scales[0], distances[1] * scales[1]]
            metadata = dict(kind="bns", centers=centers, max_coordinate_star_radii=radii,
                            source_adm_mass=phy["admmass_asymp"][0],
                            source_rest_masses=phy["restmass"], peak_density=max(phy["rho_max"]))
        elif (directory / "rnsflu_3D.las").exists():
            if any(directory.glob("rnsEMF*.las")) or (directory / "rns4ve_3D.las").exists():
                raise ValueError("This generator supports unmagnetized RNS/BNS only")
            if self._is_table(directory / "teos_parameter.dat"):
                eos = {"EOS_type": "Tabulated"}
            else:
                eos = self._read_eos(directory / "peos_parameter.dat")
                eos["EOS_type"] = "Hybrid"
            omega, _, scale = self._fluid_constants(directory / "rnsflu_3D.las")
            rows = [self._numbers(line) for line in (directory / "rnsphyplot_all.dat").read_text().splitlines()
                    if line.strip() and not line.lstrip().startswith("#")]
            matches = [row for row in rows if len(row) > 23
                       and math.isclose(row[2], scale, rel_tol=1e-9)
                       and math.isclose(row[18], omega / scale, rel_tol=1e-9, abs_tol=1e-14)]
            if len(matches) != 1:
                raise ValueError("Expected one RNS sequence matching the fluid radius and rotation rate")
            row = matches[0]
            metadata = dict(kind="rns", centers=[0.0],
                            rnstype="RNS_WL" if (directory / "rnsgra_hij_3D.las").exists() else "RNS_CF",
                            max_coordinate_star_radii=[self._surface_radius(directory / "rnssur_3D.las") * scale],
                            source_adm_mass=row[19], source_rest_masses=[row[21]], sequence=int(row[0]),
                            peak_density=self._rns_peak_density(directory))
        else:
            raise ValueError(f"{directory}: expected a COCAL RNS or BNS fluid export")
        if not all(math.isfinite(x) and x > 0 for x in
                   [metadata["source_adm_mass"], *metadata["source_rest_masses"],
                    *metadata["max_coordinate_star_radii"], metadata["peak_density"]]):
            raise ValueError("ID masses, stellar radii and peak density must be finite and positive")
        return dict(metadata, id=str(directory), **eos)



def read_bbh_metadata(directory, *, allow_disk=False):
    """Read saved BBH geometry; generator parameters never rebuild its mesh."""
    directory = Path(directory).expanduser().resolve()
    if (directory / 'disk_state.las').exists() and not allow_disk:
        raise ValueError('BBHT needs its own matter reader; this case is vacuum BBH')
    centers, radii, constants = [], [], []
    for m in (1, 2, 3):
        suffix = f'_3D_mpt{m}.las'
        for stem in ('bbhgrids', 'bbhgra', 'bnsgra_Kij'):
            if not (directory / (stem + suffix)).is_file():
                raise ValueError(f'Missing BBH export: {stem + suffix}')
        if not (directory / f'bin_ex_par_mpt{m}.dat').is_file():
            raise ValueError(f'Missing companion-excision metadata for patch {m}')
        with (directory / ('bbhgra' + suffix)).open('rb') as f:
            shape = list(map(int, f.readline().split()))
            start = f.tell()
            record = f.readline()
            width = 23 if len(record.rstrip()) > 100 else 20
            f.seek(start + math.prod(n + 1 for n in shape) * len(record))
            footer = f.readline().decode().rstrip()
        count = 6 if len(footer) > 5 * width else 4
        values = [float(footer[i*width:(i+1)*width].replace('D','E')) for i in range(count)]
        constants.append(values[-3:])  # dis,rgin,sepa in both writer formats
        if m < 3:
            centers.append((-1 if m == 1 else 1) * values[-3])
            radii.append(values[-2])
    if not math.isclose(centers[1] - centers[0], constants[0][2], abs_tol=1e-10):
        raise ValueError('Saved BBH centers disagree with the saved separation')
    source = {}
    for line in (directory / 'bbhphyseq.dat').read_text().splitlines():
        if '=' in line:
            key, value = line.split('=', 1)
            try:
                source[key.strip()] = [float(x.replace('D','E')) for x in value.split()]
            except ValueError:
                continue
    mass = source['M_ADM (asymptotic)'][0]
    if not all(math.isfinite(x) and x > 0 for x in [mass, *radii]):
        raise ValueError('BBH mass and inner radii must be finite and positive')
    return dict(kind='bbh', id=str(directory), EOS_type='Vacuum', centers=centers,
                bh_excision_radii=radii, source_adm_mass=mass,
                source_ah_areas=[source[f'AH area   {s} BH'][0] for s in ('1st','2nd')],
                source_ah_masses=[source[f'AH mass   {s} BH'][0] for s in ('1st','2nd')])


def read_bbht_metadata(directory):
    """Read the saved disk EOS, rather than the solver's old input K."""
    data = read_bbh_metadata(directory, allow_disk=True)
    directory = Path(data['id'])
    with (directory / 'disk_state.las').open('rb') as f:
        def record(size):
            marker = f.read(4)
            if len(marker) != 4 or struct.unpack('<i', marker)[0] != size:
                raise ValueError('BBHT requires little-endian, 4-byte Fortran records')
            value = f.read(size)
            if len(value) != size or f.read(4) != marker:
                raise ValueError('Truncated BBHT disk record')
            return value
        header = record(29)
        if header[:17] != b'BBHT_BHT_FLUID_V1':
            raise ValueError('Unsupported BBHT disk format')
        shape = struct.unpack('<3i', header[17:])
        with (directory / 'bbhgrids_3D_mpt3.las').open() as grid:
            if tuple(map(int, grid.readline().split())) != shape:
                raise ValueError('BBHT disk and patch-3 dimensions differ')
        if any(n < 1 for n in shape):
            raise ValueError('Invalid BBHT disk dimensions')
        nr, nt, nphi = shape
        coords = struct.unpack(f'<{sum(shape) + 3}d', record(8 * (sum(shape) + 3)))
        k, gamma, _, _ = struct.unpack('<4d', record(32))
        # rho then omega, Fortran order (r fastest); the disk height sizes the static boxes.
        plane = (nr + 1) * (nt + 1)
        rho = doubles("d")
        rho.frombytes(record(16 * plane * (nphi + 1))[:8 * plane * (nphi + 1)])
        occupied = bytearray(plane)
        for index, value in enumerate(rho):
            if value > 0:
                occupied[index % plane] = 1
        height = max((abs(coords[i] * math.cos(coords[nr + 1 + j]))
                      for j in range(nt + 1) for i in range(nr + 1) if occupied[i + (nr + 1)*j]),
                     default=0.)
    if not math.isfinite(k) or k <= 0 or not math.isfinite(gamma) or gamma <= 1:
        raise ValueError('Invalid saved BBHT EOS')
    edges = []
    namelist = re.sub(r'!.*', '', (directory / 'diskpar.dat').read_text())
    for key in ('disk_rin', 'disk_rout', 'disk_rho_target'):
        match = re.search(r'\b' + key + r'\s*=\s*([+\-\d.eEdD]+)', namelist, re.I)
        if not match:
            raise ValueError(f'BBHT diskpar.dat must specify {key}')
        edges.append(float(match[1].lower().replace('d', 'e')))
    peak = edges.pop()
    if not all(math.isfinite(r) for r in edges) or not 0 < edges[0] < edges[1]:
        raise ValueError('Invalid BBHT disk edges')
    if not math.isfinite(peak) or peak <= 0:
        raise ValueError('Invalid BBHT disk_rho_target')
    if not 0 < height < edges[1]:
        raise ValueError('Invalid BBHT disk height')
    data.update(kind='bbht', EOS_type='Hybrid', EOS_k0=k, EOS_gamma=[gamma],
                EOS_boundaries=[], disk_inner_radius=edges[0], disk_outer_radius=edges[1], disk_height=height, peak_density=peak,
                disk_format='BBHT_BHT_FLUID_V1', magnetic_field=False,
                eos_constant_source='disk_state.las')
    return data


def resolve_bbh_grid(data, overrides, local):
    """Two equal-size stacks in the saved physical coordinates.

    Profiles refine each hole while retaining the same coarse spacing. Carpet
    constructs the union and adds the RK buffers and proper-nesting support.
    """
    grid = dict(finest_dx=min(data['bh_excision_radii'])/(6 if local else 12),
                levels=5 if local else 6, domain=[8.,8.,8.], fine_box_radius=None,
                track_horizons=False, movement_cells=2.0)
    if overrides.keys() - grid.keys():
        raise ValueError('BBH grid accepts finest_dx, levels, domain, fine_box_radius, track_horizons and movement_cells')
    grid.update(overrides)
    if type(grid['track_horizons']) is not bool:
        raise ValueError('track_horizons must be True or False')
    if not math.isfinite(grid['movement_cells']) or grid['movement_cells'] <= 0:
        raise ValueError('movement_cells must be positive')
    dx, levels = grid['finest_dx'], grid['levels']
    if type(levels) is not int or not 2 <= levels <= 8 or not math.isfinite(dx) or dx <= 0:
        raise ValueError('BBH needs positive finest_dx and 2..8 levels')
    grid['coarse_dx'] = dx * 2**(levels-1)
    domain = grid['domain']
    if len(domain) != 3 or any(not math.isfinite(x) or x <= 0 for x in domain):
        raise ValueError('BBH domain needs three positive half-widths')
    # Round the outer domain outward to whole coarse cells, recording the result.
    domain = [math.ceil(x/grid['coarse_dx']) * grid['coarse_dx'] for x in domain]
    grid['domain'] = domain
    fine = [math.ceil(max(1.8*r, 12*dx)/(2*dx)) * 2*dx for r in data['bh_excision_radii']]
    clearance = (4 + grid['movement_cells']) * dx if grid['track_horizons'] else 4*dx
    if grid['track_horizons']:
        fine = [max(width, math.ceil((r+clearance)/(2*dx))*2*dx)
                for width,r in zip(fine,data['bh_excision_radii'])]
    if grid['fine_box_radius'] is not None:
        radius = grid['fine_box_radius']
        if not math.isfinite(radius) or radius < max(data['bh_excision_radii']) + clearance:
            raise ValueError('fine_box_radius must enclose both inner spheres plus stencil and movement clearance')
        fine = [math.ceil(radius/(2*dx)) * 2*dx] * 2
    # Equal-sized stacks avoid thin slivers in Carpet's box union/granulation
    # for unequal holes. Physical centres and source radii remain independent.
    fine = [max(fine)] * 2
    boxes = [[r*2**(levels-2-l) for l in range(levels-1)] for r in fine]
    extent = [max(abs(x)+b[0] for x,b in zip(data['centers'],boxes)),
              max(b[0] for b in boxes), max(b[0] for b in boxes)]
    if any(e+4*grid['coarse_dx'] >= width for e,width in zip(extent,domain)):
        raise ValueError('BBH boxes need more outer-domain clearance')
    grid.update(hole_boxes=boxes, fine_box_radii=fine)
    return grid


def resolve_bbht_grid(data, overrides, local):
    """Two BH stacks and a static central stack covering the entire disk.

    The static stack reaches the first level whose moving hole boxes stay inside
    the disk's inner radius at every orbital phase (box corners included), so no
    hole-box regrid or buffer zone touches the disk: every coarser hole box
    regridding through the disk lost rest mass (jobs 21641319-21709407). That
    level's boxes are shrunk to fit when the next finer level still nests."""
    options = dict(overrides)
    disk_levels = options.pop('disk_levels', None)
    disk_radius = options.pop('disk_radius', 1.15 * data['disk_outer_radius'])
    disk_height = options.pop('disk_height', 1.15 * data['disk_height'])
    defaults = dict(finest_dx=min(data['bh_excision_radii'])/(4 if local else 8),
                    levels=7 if local else 8, domain=[64.,64.,64.])
    grid = resolve_bbh_grid(data, {**defaults, **options}, local)
    if disk_levels is None:
        # Carpet's buffers and nesting extend each box 10-13 cells past its radius (21709407).
        overhang = lambda level: 14*grid['coarse_dx']/2**level
        offset = max(abs(x) for x in data['centers'])
        boxes = grid['hole_boxes']
        for level in range(1, grid['levels']):
            dx = grid['coarse_dx']/2**level
            # In-plane box corner after a quarter orbit, two parent cells clear of the edge.
            fits = (data['disk_inner_radius'] - 4*dx - offset)/math.sqrt(2) - overhang(level)
            fits = math.floor(fits/(2*dx)) * 2*dx
            finer = (boxes[0][level] + overhang(level+1) + MIN_DISTANCE*dx
                     if level < grid['levels']-1 else 0.)
            if boxes[0][level-1] <= fits or finer <= fits:
                for stack in boxes:
                    stack[level-1] = min(stack[level-1], fits)
                disk_levels = level
                break
        else:
            raise ValueError('No hole level fits inside the disk inner radius; refine the hole boxes')
    if type(disk_levels) is not int or not 1 <= disk_levels <= grid['levels']:
        raise ValueError('disk_levels must be between 1 and grid.levels, including the base')
    if not math.isfinite(disk_radius) or disk_radius < data['disk_outer_radius']:
        raise ValueError('disk_radius must enclose the saved disk outer edge')
    if not math.isfinite(disk_height) or disk_height < data['disk_height']:
        raise ValueError('disk_height must enclose the saved disk height')
    # Static disk boxes, (in-plane radius, z-radius) per level: the disk level is flattened
    # to the disk; each coarser level is a cube of twice the radius, the empty outer region
    # adding one level per factor two in distance, and the domain grows to fit.
    boxes = []
    if disk_levels > 1:
        dx = grid['coarse_dx'] / 2**(disk_levels-1)
        boxes.append((math.ceil(disk_radius/(2*dx)) * 2*dx, math.ceil(disk_height/(2*dx)) * 2*dx))
        for level in range(disk_levels-2, 0, -1):
            dx = grid['coarse_dx'] / 2**level
            radius = math.ceil(2*boxes[0][0]/(2*dx)) * 2*dx
            boxes.insert(0, (radius, radius))
        need = 2 * boxes[0][0]
        grid['domain'] = [max(width, math.ceil(need/grid['coarse_dx']) * grid['coarse_dx'])
                          for width in grid['domain']]
    if boxes and boxes[-1][1] < data['disk_height'] + 4*grid['coarse_dx']/2**(disk_levels-1):
        raise ValueError('disk_height needs four disk-level cells beyond the saved disk height')
    extent = boxes[0][0] if boxes else disk_radius
    if extent + 4*grid['coarse_dx'] >= min(grid['domain']):
        raise ValueError('BBHT disk boxes need more outer-domain clearance')
    if boxes and boxes[-1][0] < data['disk_outer_radius'] + 4*grid['coarse_dx']/2**(disk_levels-1):
        raise ValueError('disk_radius needs four disk-level cells beyond the saved disk outer edge')
    grid.update(disk_levels=disk_levels, disk_radius=disk_radius, disk_height=disk_height, disk_boxes=boxes)
    return grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", metavar="INDEX_OR_ID", help="configuration index or COCAL ID directory")
    parser.add_argument("--force", action="store_true", help="replace existing .par and .json outputs")
    args = parser.parse_args()
    try:
        id_dir, overrides, index = select_config(args.target)
        text, metadata = generate(id_dir, overrides)
        metadata["config_index"] = index
        output_name = metadata["config_name"] or Path(metadata["id"]).name
        output = Path.cwd() / (output_name + ".par")
        for path in (output, output.with_suffix(".json")):
            if path.exists() and not args.force:
                raise ValueError(f"{path.name} exists; use another directory or --force to replace it")
        output.write_text(text)
        output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise SystemExit(f"Cannot generate parfile: {error}") from None
    print(output)
    if index is not None:
        print(f'Configuration {index}: {metadata["config_name"] or Path(metadata["id"]).name}')
    if metadata["added_thorns"]:
        print(f'Thorns: {len(metadata["active_thorns"])} active, '
              f'{len(metadata["added_thorns"])} dependencies added (listed in JSON)')
    eos_description = ("vacuum" if metadata["EOS_type"] == "Vacuum" else "tabulated EOS" if metadata["EOS_type"] == "Tabulated"
                       else f'{len(metadata["EOS_gamma"])} EOS segment(s)')
    print(f'{metadata["kind"].upper()}, {eos_description}, '
          f'{metadata["profile"]}, dx={metadata["grid"]["finest_dx"]:g}, '
          f'{metadata["grid"]["levels"]} levels, dtfac={metadata["run"]["dtfac"]:g}')
    grid = metadata["grid"]
    print(f'Domain half-widths={array(grid["domain"])}, coarse dx={grid["coarse_dx"]:g}, '
          f'finest box half-widths={array([boxes[-1] for boxes in grid.get("hole_boxes", grid.get("star_boxes", []))])}')
    cadence = metadata['tracking']['every']
    activity = f'tracking every {cadence} iterations'
    if metadata['kind'] in ('bbh', 'bbht'):
        activity = (f'horizon/AMR tracking every {cadence} iterations'
                    if metadata['tracking']['enabled']
                    else f'fixed AMR, horizon checks every {cadence} iterations')
    print(f'Object centres x={array(metadata["centers"])}, '
          f'{activity}, output every {metadata["run"]["out_every"]}')
    if metadata["kind"] == "bns":
        print(f'Common refined levels={max(0, grid["common_levels"] - 1)}, '
              f'box movement threshold={metadata["tracking"]["movement_threshold"]:g}')
        print(f'BNS export settings: compact={metadata["compact"]}, '
              f'xunit={metadata["xunit"]}, ecc_cor_velx={metadata["ecc_cor_velx"]}, '
              f'ecc_cor_velz={metadata["ecc_cor_velz"]}')
    if any(radius is None for radius in metadata["outer_qlm_radii"]):
        print("Outer spheres that cannot enclose the stars with stencil clearance are omitted.")


if __name__ == "__main__":
    main()
