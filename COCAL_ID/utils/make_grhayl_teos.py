#!/usr/bin/env python3
"""Prepare a hot tabulated EOS file for GRHayL/GRHayLET.

COCAL TEOS files are cold beta-equilibrium slices and are not sufficient for
GRHayLHD evolution.  This utility validates and optionally downsamples an
existing NRPyEOS / NuLib-style HDF5 table, or converts a raw CompOSE
eos.zip/eos directory, into the schema read by GRHayLib.
"""

from __future__ import annotations

import argparse
import json
import zipfile
from contextlib import contextmanager
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:  # pragma: no cover - exercised by user environment
    raise SystemExit(
        "make_grhayl_teos.py requires h5py and numpy. Run it in a Python "
        "environment that provides them, e.g. `uv run python ...`."
    ) from exc


TABLE_VARIABLES = (
    "logpress",
    "logenergy",
    "entropy",
    "munu",
    "cs2",
    "dedt",
    "dpdrhoe",
    "dpderho",
    "muhat",
    "mu_e",
    "mu_p",
    "mu_n",
    "Xa",
    "Xh",
    "Xn",
    "Xp",
    "Abar",
    "Zbar",
    "gamma",
)

# GRHayL's tabulated EOS reader expects cgs units in the file, then converts
# them internally to code units.  CompOSE gives n_b in fm^-3 and thermo columns
# per baryon, so keep the conversion factors explicit here.
BARYON_MASS_G = 1.66053906660e-24
AMU_MEV = 931.494061
MEV_TO_ERG = 1.602176634e-6
C_CGS = 2.99792458e10
C2_CGS = C_CGS * C_CGS
FM3_TO_CM3_INV = 1.0e39
MEVFM3_TO_DYNCM2 = MEV_TO_ERG * FM3_TO_CM3_INV
MEV_OVER_C2_TO_G = MEV_TO_ERG / C2_CGS

# Normal COCAL_ID workflow defaults.  The matching CLI switches are kept for
# reproducibility/debugging, but are intentionally hidden from --help.
DEFAULT_COMPOSE_ENERGY_ZERO = "amu"
DEFAULT_COMPOSE_DENSITY_MASS = "amu"
DEFAULT_COMPOSE_CHEMICAL_POTENTIAL_ZERO = "rest-mass-subtracted"
DEFAULT_COMPOSE_DERIVATIVE_METHOD = "pchip-log"
DEFAULT_DEPSDT_RELATIVE_FLOOR = 1.0e-14


def dataset_array(handle: h5py.File, name: str) -> np.ndarray:
    if name not in handle:
        raise KeyError(name)
    return np.asarray(handle[name][()])


def scalar_int(handle: h5py.File, name: str, fallback: int | None = None) -> int:
    if name not in handle:
        if fallback is None:
            raise KeyError(name)
        return fallback
    value = np.asarray(handle[name][()]).reshape(-1)
    if value.size != 1:
        raise ValueError(f"{name} should contain one integer, found shape {handle[name].shape}")
    return int(value[0])


def read_axes(handle: h5py.File) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if "logrho" in handle:
        logrho = dataset_array(handle, "logrho").astype(np.float64).reshape(-1)
    elif "nb" in handle:
        nb = dataset_array(handle, "nb").astype(np.float64).reshape(-1)
        if np.any(nb <= 0.0):
            raise ValueError("Cannot build logrho from non-positive Compose nb values")
        logrho = np.log10(nb * FM3_TO_CM3_INV * BARYON_MASS_G)
    else:
        raise KeyError("logrho")

    if "logtemp" in handle:
        logtemp = dataset_array(handle, "logtemp").astype(np.float64).reshape(-1)
    elif "t" in handle:
        temp = dataset_array(handle, "t").astype(np.float64).reshape(-1)
        if np.any(temp <= 0.0):
            raise ValueError("Cannot build logtemp from non-positive Compose t values")
        logtemp = np.log10(temp)
    else:
        raise KeyError("logtemp")

    if "ye" in handle:
        ye = dataset_array(handle, "ye").astype(np.float64).reshape(-1)
    elif "yq" in handle:
        ye = dataset_array(handle, "yq").astype(np.float64).reshape(-1)
    else:
        raise KeyError("ye")

    return logrho, logtemp, ye


def check_axis(name: str, values: np.ndarray, atol: float = 1.0e-12) -> None:
    if values.ndim != 1 or values.size < 2:
        raise ValueError(f"{name} must be a 1D axis with at least two points")
    diffs = np.diff(values)
    if np.any(diffs <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    step = diffs[0]
    err = float(np.max(np.abs(diffs - step)))
    scale = max(1.0, abs(float(step)), abs(float(values[-1] - values[0])))
    if err > atol * scale:
        raise ValueError(
            f"{name} is not uniformly spaced enough for GRHayL indexing: "
            f"max step error {err:.6e}"
        )


def variable_as_3d(
    handle: h5py.File,
    name: str,
    nrho: int,
    ntemp: int,
    nye: int,
) -> np.ndarray:
    data = dataset_array(handle, name).astype(np.float64)
    if data.shape == (nye, ntemp, nrho):
        return data
    if data.size == nrho * ntemp * nye:
        return np.reshape(data, (nye, ntemp, nrho), order="C")
    raise ValueError(
        f"{name} has shape {data.shape}; expected {(nye, ntemp, nrho)} "
        f"or {nrho * ntemp * nye} flat values"
    )


def stride_indices(n: int, stride: int, axis_name: str) -> np.ndarray:
    if stride < 1:
        raise ValueError(f"{axis_name} stride must be positive")
    idx = np.arange(0, n, stride, dtype=np.int64)
    if idx.size < 2:
        raise ValueError(f"{axis_name} stride leaves fewer than two points")
    return idx


def validate_input(handle: h5py.File) -> dict[str, object]:
    logrho, logtemp, ye = read_axes(handle)
    nrho = scalar_int(handle, "pointsrho", len(logrho))
    ntemp = scalar_int(handle, "pointstemp", len(logtemp))
    nye = scalar_int(handle, "pointsye", len(ye))

    if (nrho, ntemp, nye) != (len(logrho), len(logtemp), len(ye)):
        raise ValueError(
            "Axis lengths do not match pointsrho/pointstemp/pointsye: "
            f"{(len(logrho), len(logtemp), len(ye))} vs {(nrho, ntemp, nye)}"
        )
    check_axis("logrho", logrho)
    check_axis("logtemp", logtemp)
    check_axis("ye", ye)

    missing = [name for name in TABLE_VARIABLES if name not in handle]
    if missing:
        raw_compose_hint = any(name in handle for name in ("Thermo_qty", "eos.thermo"))
        hint = ""
        if raw_compose_hint:
            hint = (
                " This looks like a raw or partially processed Compose table; "
                "derive the full NRPyEOS thermodynamic fields first."
            )
        raise KeyError(f"Missing GRHayL table datasets: {', '.join(missing)}.{hint}")
    if "energy_shift" not in handle:
        raise KeyError("Missing GRHayL table dataset: energy_shift")

    for name in TABLE_VARIABLES:
        variable_as_3d(handle, name, nrho, ntemp, nye)

    have_rel_cs2 = None
    if "have_rel_cs2" in handle:
        value = np.asarray(handle["have_rel_cs2"][()]).reshape(-1)
        if value.size != 1:
            raise ValueError(f"have_rel_cs2 should contain one integer, found shape {handle['have_rel_cs2'].shape}")
        have_rel_cs2 = int(value[0])

    return {
        "pointsrho": nrho,
        "pointstemp": ntemp,
        "pointsye": nye,
        "rho_cgs_min": float(10.0 ** logrho[0]),
        "rho_cgs_max": float(10.0 ** logrho[-1]),
        "temp_min_MeV": float(10.0 ** logtemp[0]),
        "temp_max_MeV": float(10.0 ** logtemp[-1]),
        "ye_min": float(ye[0]),
        "ye_max": float(ye[-1]),
        "have_rel_cs2": have_rel_cs2,
    }


def copy_attrs(src: h5py.AttributeManager, dst: h5py.AttributeManager) -> None:
    for key, value in src.items():
        dst[key] = value


def write_table(
    input_path: Path,
    output_path: Path,
    rho_stride: int,
    temp_stride: int,
    ye_stride: int,
    force: bool,
) -> dict[str, object]:
    if output_path.exists() and not force:
        raise FileExistsError(f"{output_path} exists; pass --force to overwrite")

    with h5py.File(input_path, "r") as src:
        before = validate_input(src)
        logrho, logtemp, ye = read_axes(src)
        nrho = int(before["pointsrho"])
        ntemp = int(before["pointstemp"])
        nye = int(before["pointsye"])

        irho = stride_indices(nrho, rho_stride, "rho")
        itemp = stride_indices(ntemp, temp_stride, "temperature")
        iye = stride_indices(nye, ye_stride, "Ye")

        if output_path.exists():
            output_path.unlink()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(output_path, "w") as dst:
            copy_attrs(src.attrs, dst.attrs)
            dst.attrs["source_table"] = str(input_path)
            dst.attrs["generated_by"] = "COCAL_ID/utils/make_grhayl_teos.py"
            dst.attrs["rho_stride"] = rho_stride
            dst.attrs["temp_stride"] = temp_stride
            dst.attrs["ye_stride"] = ye_stride

            dst.create_dataset("pointsrho", data=np.asarray([irho.size], dtype=np.int32))
            dst.create_dataset("pointstemp", data=np.asarray([itemp.size], dtype=np.int32))
            dst.create_dataset("pointsye", data=np.asarray([iye.size], dtype=np.int32))
            dst.create_dataset("logrho", data=logrho[irho])
            dst.create_dataset("logtemp", data=logtemp[itemp])
            dst.create_dataset("ye", data=ye[iye])
            if "have_rel_cs2" in src:
                dst.create_dataset(
                    "have_rel_cs2",
                    data=np.asarray(src["have_rel_cs2"][()]).astype(np.int32).reshape(-1),
                )
            dst.create_dataset(
                "energy_shift",
                data=np.asarray(src["energy_shift"][()]).astype(np.float64).reshape(-1),
            )

            for name in TABLE_VARIABLES:
                table = variable_as_3d(src, name, nrho, ntemp, nye)
                reduced = table[np.ix_(iye, itemp, irho)]
                dst.create_dataset(name, data=reduced)

    with h5py.File(output_path, "r") as check:
        after = validate_input(check)
    after["source_table"] = str(input_path)
    after["output_table"] = str(output_path)
    after["strides"] = {"rho": rho_stride, "temperature": temp_stride, "Ye": ye_stride}
    return after


def read_compose_axis(lines: list[str], name: str) -> np.ndarray:
    if len(lines) < 2:
        raise ValueError(f"{name} is too short")
    n = int(lines[1].split()[0])
    values = np.asarray([float(line.split()[0]) for line in lines[2 : 2 + n]], dtype=np.float64)
    if values.size != n:
        raise ValueError(f"{name} declares {n} values but contains {values.size}")
    return values


@contextmanager
def open_compose_member(input_path: Path, member: str):
    if input_path.is_dir():
        with (input_path / member).open("rb") as f:
            yield f
    else:
        with zipfile.ZipFile(input_path) as zf:
            with zf.open(member) as f:
                yield f


def read_compose_text(input_path: Path, member: str) -> list[str]:
    with open_compose_member(input_path, member) as f:
        return [line.decode("ascii").strip() for line in f if line.strip()]


def is_raw_compose_input(path: Path) -> bool:
    if path.is_dir():
        return (path / "eos.thermo").exists() and (path / "eos.nb").exists()
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = set(zf.namelist())
        return {"eos.thermo", "eos.nb", "eos.t", "eos.yq"}.issubset(names)
    return False


def finite_gradient(values: np.ndarray, coords: np.ndarray, axis: int, name: str) -> np.ndarray:
    grad = np.gradient(values, coords, axis=axis, edge_order=2)
    if not np.all(np.isfinite(grad)):
        raise ValueError(f"Non-finite derivative while computing {name}")
    return grad


def pchip_slopes(values: np.ndarray, coords: np.ndarray, axis: int, name: str) -> np.ndarray:
    """Shape-preserving derivatives along one monotone coordinate axis."""

    values_axis = np.moveaxis(values, axis, 0)
    original_shape = values_axis.shape
    y = values_axis.reshape(original_shape[0], -1)
    x = coords.astype(np.float64)
    if x.ndim != 1 or x.size != y.shape[0]:
        raise ValueError(f"{name} coordinate axis has incompatible shape")
    if x.size < 3:
        return finite_gradient(values, coords, axis=axis, name=name)

    h = np.diff(x)
    if np.any(h <= 0.0):
        raise ValueError(f"{name} coordinate axis must be strictly increasing")
    delta = np.diff(y, axis=0) / h[:, None]
    slopes = np.zeros_like(y)

    same_sign = (delta[:-1] * delta[1:]) > 0.0
    w1 = (2.0 * h[1:] + h[:-1])[:, None]
    w2 = (h[1:] + 2.0 * h[:-1])[:, None]
    slopes_mid = slopes[1:-1]
    denom = np.zeros_like(slopes_mid)
    tmp = np.zeros_like(slopes_mid)
    np.divide(w1, delta[:-1], out=denom, where=same_sign)
    np.divide(w2, delta[1:], out=tmp, where=same_sign)
    denom += tmp
    numerator = np.broadcast_to(w1 + w2, slopes_mid.shape)
    np.divide(numerator, denom, out=slopes_mid, where=same_sign)

    first = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1])
    first = np.where(first * delta[0] <= 0.0, 0.0, first)
    first = np.where(
        ((delta[0] * delta[1]) < 0.0) & (np.abs(first) > 3.0 * np.abs(delta[0])),
        3.0 * delta[0],
        first,
    )
    slopes[0] = first

    last = ((2.0 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]) / (h[-1] + h[-2])
    last = np.where(last * delta[-1] <= 0.0, 0.0, last)
    last = np.where(
        ((delta[-1] * delta[-2]) < 0.0) & (np.abs(last) > 3.0 * np.abs(delta[-1])),
        3.0 * delta[-1],
        last,
    )
    slopes[-1] = last

    result = slopes.reshape(original_shape)
    result = np.moveaxis(result, 0, axis)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"Non-finite derivative while computing {name}")
    return result


def derivative_summary(values: np.ndarray, prefix: str) -> dict[str, float | int]:
    finite = np.isfinite(values)
    summary: dict[str, float | int] = {
        f"{prefix}_nonfinite_points": int(values.size - np.count_nonzero(finite)),
    }
    if np.any(finite):
        finite_values = values[finite]
        summary.update(
            {
                f"{prefix}_min": float(np.min(finite_values)),
                f"{prefix}_max": float(np.max(finite_values)),
                f"{prefix}_negative_points": int(np.count_nonzero(finite_values < 0.0)),
                f"{prefix}_zero_points": int(np.count_nonzero(finite_values == 0.0)),
            }
        )
    return summary


def sound_speed_summary(cs2_dimensionless: np.ndarray) -> dict[str, float | int]:
    finite = np.isfinite(cs2_dimensionless)
    summary: dict[str, float | int] = {
        "cs2_nonfinite_points": int(cs2_dimensionless.size - np.count_nonzero(finite)),
    }
    if np.any(finite):
        finite_values = cs2_dimensionless[finite]
        summary.update(
            {
                "cs2_min_dimensionless": float(np.min(finite_values)),
                "cs2_max_dimensionless": float(np.max(finite_values)),
                "cs2_negative_points": int(np.count_nonzero(finite_values < 0.0)),
                "cs2_superluminal_points": int(np.count_nonzero(finite_values > 1.0)),
            }
        )
    return summary


def uniform_log_axis(values: np.ndarray, name: str) -> np.ndarray:
    log_values = np.log10(values)
    check_axis(name, log_values, atol=1.0e-6)
    return 10.0 ** np.linspace(log_values[0], log_values[-1], values.size)


def uniform_linear_axis(values: np.ndarray, name: str) -> np.ndarray:
    check_axis(name, values, atol=1.0e-6)
    return np.linspace(values[0], values[-1], values.size)


def read_compose_thermo(
    input_path: Path,
    nb: np.ndarray,
    temp: np.ndarray,
    ye: np.ndarray,
    energy_zero: str,
    depsdt_relative_floor: float,
    derivative_method: str,
    density_mass: str,
    chemical_potential_zero: str,
    clean_sound_speed: bool,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    nrho, ntemp, nye = nb.size, temp.size, ye.size
    dims = (nye, ntemp, nrho)
    q1 = np.empty(dims, dtype=np.float64)
    entropy = np.empty(dims, dtype=np.float64)
    q3 = np.empty(dims, dtype=np.float64)
    q4 = np.empty(dims, dtype=np.float64)
    q5 = np.empty(dims, dtype=np.float64)
    q7 = np.empty(dims, dtype=np.float64)
    filled = np.zeros(dims, dtype=bool)

    with open_compose_member(input_path, "eos.thermo") as f:
        header = f.readline().split()
        if len(header) < 3:
            raise ValueError("eos.thermo header should contain neutron/proton masses and lepton flag")
        neutron_mass_mev = float(header[0])
        proton_mass_mev = float(header[1])
        lepton_flag = int(header[2])

        for raw in f:
            parts = raw.split()
            if not parts:
                continue
            it = int(parts[0]) - 1
            ir = int(parts[1]) - 1
            iy = int(parts[2]) - 1
            if not (0 <= it < ntemp and 0 <= ir < nrho and 0 <= iy < nye):
                raise ValueError(f"eos.thermo index out of range: {(it + 1, ir + 1, iy + 1)}")
            q1[iy, it, ir] = float(parts[3])
            entropy[iy, it, ir] = float(parts[4])
            q3[iy, it, ir] = float(parts[5])
            q4[iy, it, ir] = float(parts[6])
            q5[iy, it, ir] = float(parts[7])
            q7[iy, it, ir] = float(parts[9])
            filled[iy, it, ir] = True

    if not np.all(filled):
        raise ValueError(f"eos.thermo is missing {np.size(filled) - int(np.count_nonzero(filled))} grid points")

    density_mass_g = BARYON_MASS_G
    if density_mass == "neutron":
        density_mass_g = neutron_mass_mev * MEV_OVER_C2_TO_G
    elif density_mass != "amu":
        raise ValueError(f"Unknown density_mass={density_mass!r}")

    nb3 = nb.reshape(1, 1, nrho)
    rho_cgs = nb * FM3_TO_CM3_INV * density_mass_g
    rho3 = rho_cgs.reshape(1, 1, nrho)
    press_cgs = q1 * nb3 * MEVFM3_TO_DYNCM2

    total_energy_mev_per_baryon = (1.0 + q7) * neutron_mass_mev
    if energy_zero == "amu":
        eps_mev_per_baryon = total_energy_mev_per_baryon - AMU_MEV
    elif energy_zero == "neutron":
        eps_mev_per_baryon = q7 * neutron_mass_mev
    elif energy_zero == "compose":
        eps_mev_per_baryon = q7 * density_mass_g / MEV_OVER_C2_TO_G
    else:
        raise ValueError(f"Unknown energy_zero={energy_zero!r}")
    eps_cgs = eps_mev_per_baryon * MEV_TO_ERG / density_mass_g
    eps_min = float(np.min(eps_cgs))
    eps_abs_max = float(np.max(np.abs(eps_cgs)))
    energy_shift = 0.0
    if eps_min <= 0.0:
        energy_shift = -eps_min + max(1.0e-12 * C2_CGS, 1.0e-12 * eps_abs_max)

    if derivative_method == "finite":
        temp_axis = temp.astype(np.float64)
        press_t = finite_gradient(press_cgs, temp_axis, axis=1, name="dP/dT")
        eps_t = finite_gradient(eps_cgs, temp_axis, axis=1, name="deps/dT")
        press_rho = finite_gradient(press_cgs, rho_cgs, axis=2, name="dP/drho")
        eps_rho = finite_gradient(eps_cgs, rho_cgs, axis=2, name="deps/drho")
    elif derivative_method == "pchip-log":
        log_temp = np.log(temp.astype(np.float64))
        log_rho = np.log(rho_cgs.astype(np.float64))
        eps_shifted = eps_cgs + energy_shift
        log_press = np.log(press_cgs)
        log_eps_shifted = np.log(eps_shifted)

        press_t = (
            press_cgs
            / temp.reshape(1, ntemp, 1)
            * pchip_slopes(log_press, log_temp, axis=1, name="dlnP/dlnT")
        )
        eps_t = (
            eps_shifted
            / temp.reshape(1, ntemp, 1)
            * pchip_slopes(log_eps_shifted, log_temp, axis=1, name="dlneps/dlnT")
        )
        press_rho = (
            press_cgs
            / rho3
            * pchip_slopes(log_press, log_rho, axis=2, name="dlnP/dlnrho")
        )
        eps_rho = (
            eps_shifted
            / rho3
            * pchip_slopes(log_eps_shifted, log_rho, axis=2, name="dlneps/dlnrho")
        )
    else:
        raise ValueError(f"Unknown derivative_method={derivative_method!r}")

    depsdt_floor = depsdt_relative_floor * float(np.max(np.abs(eps_t)))
    depsdt_floor = max(depsdt_floor, np.finfo(np.float64).tiny)
    eps_t_regularized_count = int(np.count_nonzero(np.abs(eps_t) < depsdt_floor))
    eps_t_floor = np.maximum(np.abs(eps_t), depsdt_floor)
    eps_t_safe = np.where(eps_t >= 0.0, eps_t_floor, -eps_t_floor)
    dpderho = press_t / eps_t_safe
    dpdrhoe = press_rho - press_t * eps_rho / eps_t_safe
    sound_numerator = dpdrhoe + press_cgs * dpderho / (rho3 * rho3)
    enthalpy_cgs = C2_CGS + eps_cgs + press_cgs / rho3
    cs2_dimensionless = sound_numerator / enthalpy_cgs
    cs2_before_cleaning = cs2_dimensionless.copy()
    if clean_sound_speed:
        cs2_dimensionless = np.nan_to_num(cs2_dimensionless, nan=1.0 - 1.0e-15, posinf=1.0, neginf=0.0)
        cs2_dimensionless = np.clip(cs2_dimensionless, 0.0, 1.0)
        sound_numerator = cs2_dimensionless * enthalpy_cgs
    gamma = sound_numerator * rho3 / press_cgs

    if not np.all(np.isfinite(dpderho)):
        raise ValueError("Non-finite dpderho generated")
    if not np.all(np.isfinite(dpdrhoe)):
        raise ValueError("Non-finite dpdrhoe generated")
    if not np.all(np.isfinite(cs2_dimensionless)):
        raise ValueError("Non-finite cs2 generated")
    if not np.all(np.isfinite(gamma)):
        raise ValueError("Non-finite gamma generated")
    if np.any(enthalpy_cgs <= 0.0):
        raise ValueError("Non-positive enthalpy generated")
    if np.any(press_cgs <= 0.0):
        raise ValueError("Non-positive pressure generated")
    if np.any(eps_cgs + energy_shift <= 0.0):
        raise ValueError("Energy shift did not make logenergy positive")

    if chemical_potential_zero == "rest-mass-subtracted":
        mu_b = neutron_mass_mev * q3
    elif chemical_potential_zero == "full":
        mu_b = neutron_mass_mev * (1.0 + q3)
    else:
        raise ValueError(f"Unknown chemical_potential_zero={chemical_potential_zero!r}")
    mu_q = neutron_mass_mev * q4
    mu_l = neutron_mass_mev * q5
    mu_n = mu_b
    mu_p = mu_b + mu_q
    mu_e = mu_l - mu_q

    tables = {
        "logpress": np.log10(press_cgs),
        "logenergy": np.log10(eps_cgs + energy_shift),
        "entropy": entropy,
        "munu": mu_l,
        "cs2": sound_numerator,
        "dedt": eps_t,
        "dpdrhoe": dpdrhoe,
        "dpderho": dpderho,
        "muhat": mu_n - mu_p,
        "mu_e": mu_e,
        "mu_p": mu_p,
        "mu_n": mu_n,
        "gamma": gamma,
    }
    metadata = {
        "neutron_mass_mev": neutron_mass_mev,
        "proton_mass_mev": proton_mass_mev,
        "lepton_flag": lepton_flag,
        "density_mass": density_mass,
        "density_mass_g": float(density_mass_g),
        "energy_zero": energy_zero,
        "chemical_potential_zero": chemical_potential_zero,
        "derivative_method": derivative_method,
        "clean_sound_speed": bool(clean_sound_speed),
        "energy_shift_cgs": energy_shift,
        "eps_cgs_min_unshifted": eps_min,
        "eps_cgs_max_unshifted": float(np.max(eps_cgs)),
        "cs2_table_cgs_min": float(np.min(sound_numerator)),
        "cs2_table_cgs_max": float(np.max(sound_numerator)),
        "gamma_min": float(np.min(gamma)),
        "gamma_max": float(np.max(gamma)),
        "depsdt_relative_floor": depsdt_relative_floor,
        "depsdt_floor_cgs": depsdt_floor,
        "depsdt_regularized_points": eps_t_regularized_count,
    }
    metadata.update(sound_speed_summary(cs2_before_cleaning))
    if clean_sound_speed:
        metadata.update({f"cleaned_{key}": value for key, value in sound_speed_summary(cs2_dimensionless).items()})
    metadata.update(derivative_summary(eps_t, "dedt"))
    metadata.update(derivative_summary(dpdrhoe, "dpdrhoe"))
    metadata.update(derivative_summary(dpderho, "dpderho"))
    return tables, metadata


def read_compose_composition(
    input_path: Path,
    nrho: int,
    ntemp: int,
    nye: int,
) -> dict[str, np.ndarray]:
    dims = (nye, ntemp, nrho)
    x_n = np.zeros(dims, dtype=np.float64)
    x_p = np.zeros(dims, dtype=np.float64)
    x_a = np.zeros(dims, dtype=np.float64)
    x_h = np.zeros(dims, dtype=np.float64)
    abar = np.zeros(dims, dtype=np.float64)
    zbar = np.zeros(dims, dtype=np.float64)

    with open_compose_member(input_path, "eos.compo") as f:
        for raw in f:
            parts = raw.split()
            if not parts:
                continue
            it = int(parts[0]) - 1
            ir = int(parts[1]) - 1
            iy = int(parts[2]) - 1
            pos = 4
            npairs = int(parts[pos])
            pos += 1
            for _ in range(npairs):
                idx = int(parts[pos])
                frac = float(parts[pos + 1])
                pos += 2
                if idx == 10:
                    x_n[iy, it, ir] = frac
                elif idx == 11:
                    x_p[iy, it, ir] = frac
                elif idx == 4002:
                    x_a[iy, it, ir] = 4.0 * frac
            nquad = int(parts[pos])
            pos += 1
            for _ in range(nquad):
                idx = int(parts[pos])
                a_val = float(parts[pos + 1])
                z_val = float(parts[pos + 2])
                frac = float(parts[pos + 3])
                pos += 4
                if idx == 999:
                    abar[iy, it, ir] = a_val
                    zbar[iy, it, ir] = z_val
                    x_h[iy, it, ir] = a_val * frac

    return {"Xn": x_n, "Xp": x_p, "Xa": x_a, "Xh": x_h, "Abar": abar, "Zbar": zbar}


def write_compose_table(
    input_path: Path,
    output_path: Path,
    rho_stride: int,
    temp_stride: int,
    ye_stride: int,
    force: bool,
    energy_zero: str,
    depsdt_relative_floor: float,
    derivative_method: str,
    density_mass: str,
    chemical_potential_zero: str,
    clean_sound_speed: bool,
) -> dict[str, object]:
    if output_path.exists() and not force:
        raise FileExistsError(f"{output_path} exists; pass --force to overwrite")

    nb_raw = read_compose_axis(read_compose_text(input_path, "eos.nb"), "eos.nb")
    temp_raw = read_compose_axis(read_compose_text(input_path, "eos.t"), "eos.t")
    ye_raw = read_compose_axis(read_compose_text(input_path, "eos.yq"), "eos.yq")
    nb = uniform_log_axis(nb_raw, "logrho")
    temp = uniform_log_axis(temp_raw, "logtemp")
    ye = uniform_linear_axis(ye_raw, "ye")
    irho = stride_indices(nb.size, rho_stride, "rho")
    itemp = stride_indices(temp.size, temp_stride, "temperature")
    iye = stride_indices(ye.size, ye_stride, "Ye")

    thermo, meta = read_compose_thermo(
        input_path,
        nb,
        temp,
        ye,
        energy_zero,
        depsdt_relative_floor,
        derivative_method,
        density_mass,
        chemical_potential_zero,
        clean_sound_speed,
    )
    comp = read_compose_composition(input_path, nb.size, temp.size, ye.size)
    tables = {**thermo, **comp}

    if output_path.exists():
        output_path.unlink()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, "w") as dst:
        dst.attrs["source_table"] = str(input_path)
        dst.attrs["source_format"] = "raw_compose"
        dst.attrs["generated_by"] = "COCAL_ID/utils/make_grhayl_teos.py"
        dst.attrs["energy_zero"] = energy_zero
        dst.attrs["density_mass"] = density_mass
        dst.attrs["chemical_potential_zero"] = chemical_potential_zero
        dst.attrs["derivative_method"] = derivative_method
        dst.attrs["clean_sound_speed"] = clean_sound_speed
        dst.attrs["compose_lepton_flag"] = meta["lepton_flag"]
        if meta["lepton_flag"] == 0:
            dst.attrs["warning"] = (
                "For no-electron CompOSE tables, munu is the tabulated lepton "
                "chemical potential and is zero by construction."
            )
        dst.create_dataset("pointsrho", data=np.asarray([irho.size], dtype=np.int32))
        dst.create_dataset("pointstemp", data=np.asarray([itemp.size], dtype=np.int32))
        dst.create_dataset("pointsye", data=np.asarray([iye.size], dtype=np.int32))
        dst.create_dataset("have_rel_cs2", data=np.asarray([0], dtype=np.int32))
        dst.create_dataset("logrho", data=np.log10(nb[irho] * FM3_TO_CM3_INV * float(meta["density_mass_g"])))
        dst.create_dataset("logtemp", data=np.log10(temp[itemp]))
        dst.create_dataset("ye", data=ye[iye])
        dst.create_dataset("energy_shift", data=np.asarray([meta["energy_shift_cgs"]], dtype=np.float64))
        for name in TABLE_VARIABLES:
            dst.create_dataset(name, data=tables[name][np.ix_(iye, itemp, irho)])

    with h5py.File(output_path, "r") as check:
        after = validate_input(check)
    after.update(meta)
    after["source_table"] = str(input_path)
    after["source_format"] = "raw_compose"
    after["output_table"] = str(output_path)
    after["energy_zero"] = energy_zero
    after["density_mass"] = density_mass
    after["chemical_potential_zero"] = chemical_potential_zero
    after["derivative_method"] = derivative_method
    after["clean_sound_speed"] = clean_sound_speed
    after["strides"] = {"rho": rho_stride, "temperature": temp_stride, "Ye": ye_stride}
    after["composition_mass_fraction_sum_min"] = float(np.min(comp["Xn"] + comp["Xp"] + comp["Xa"] + comp["Xh"]))
    after["composition_mass_fraction_sum_max"] = float(np.max(comp["Xn"] + comp["Xp"] + comp["Xa"] + comp["Xh"]))
    return after


def reject_cocal_cold_table(path: Path) -> None:
    if path.is_file() and path.suffix.lower() in {".dat", ".txt"}:
        with path.open(errors="ignore") as f:
            first = f.readline().split()
        if first and first[0].isdigit():
            raise SystemExit(
                f"{path} looks like a COCAL cold TEOS slice. Use the full hot "
                "parent Compose/NRPyEOS HDF5 table for GRHayL evolution."
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "Normal workflow:\n"
            "  make_grhayl_teos.py eos.zip DD2_GRHayL_NRPyEOS.h5 --force\n\n"
            "The default raw-CompOSE conversion preserves thermodynamic data and "
            "uses shape-preserving derivatives only for GRHayL auxiliary arrays."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", type=Path, help="Input hot EOS HDF5 table or raw CompOSE zip/directory")
    parser.add_argument("output", type=Path, nargs="?", help="Output GRHayL HDF5 table")
    parser.add_argument("--rho-stride", type=int, default=1)
    parser.add_argument("--temp-stride", type=int, default=1)
    parser.add_argument("--ye-stride", type=int, default=1)
    parser.add_argument(
        "--energy-zero",
        choices=("amu", "neutron", "compose"),
        default=DEFAULT_COMPOSE_ENERGY_ZERO,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--density-mass",
        choices=("amu", "neutron"),
        default=DEFAULT_COMPOSE_DENSITY_MASS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--chemical-potential-zero",
        choices=("rest-mass-subtracted", "full"),
        default=DEFAULT_COMPOSE_CHEMICAL_POTENTIAL_ZERO,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--derivative-method",
        choices=("pchip-log", "finite"),
        default=DEFAULT_COMPOSE_DERIVATIVE_METHOD,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--depsdt-relative-floor",
        type=float,
        default=DEFAULT_DEPSDT_RELATIVE_FLOOR,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--clean-sound-speed",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--json", action="store_true", help="Print summary as JSON")
    args = parser.parse_args()

    reject_cocal_cold_table(args.input)
    if args.validate_only:
        with h5py.File(args.input, "r") as handle:
            summary = validate_input(handle)
            summary["input_table"] = str(args.input)
    else:
        if args.output is None:
            raise SystemExit("output is required unless --validate-only is passed")
        if is_raw_compose_input(args.input):
            summary = write_compose_table(
                args.input,
                args.output,
                args.rho_stride,
                args.temp_stride,
                args.ye_stride,
                args.force,
                args.energy_zero,
                args.depsdt_relative_floor,
                args.derivative_method,
                args.density_mass,
                args.chemical_potential_zero,
                args.clean_sound_speed,
            )
        else:
            summary = write_table(
                args.input,
                args.output,
                args.rho_stride,
                args.temp_stride,
                args.ye_stride,
                args.force,
            )

    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        for key, value in summary.items():
            print(f"{key}: {value}")


if __name__ == "__main__":
    main()
