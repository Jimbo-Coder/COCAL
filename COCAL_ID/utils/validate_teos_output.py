#!/usr/bin/env python3
"""Check COCAL_ID TEOS hydro output against a COCAL teos_parameter.dat table."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def read_ascii_1d(path: Path, iteration: int) -> dict[tuple[float, float, float], float]:
    values: dict[tuple[float, float, float], float] = {}
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            if int(cols[0]) != iteration:
                continue
            coord = (float(cols[9]), float(cols[10]), float(cols[11]))
            values[coord] = float(cols[12])
    return values


def lagrange4(xs: list[float], fs: list[float], x: float) -> float:
    total = 0.0
    for i in range(4):
        term = fs[i]
        for j in range(4):
            if i != j:
                term *= (x - xs[j]) / (xs[i] - xs[j])
        total += term
    return total


def read_teos(path: Path) -> tuple[list[float], list[float], list[float]]:
    with path.open() as f:
        npts = int(f.readline().split()[0])
        rho: list[float] = []
        press: list[float] = []
        q: list[float] = []
        for _ in range(npts):
            ene, pre, _h, rrho = (float(x) for x in f.readline().split()[:4])
            del ene
            rho.append(rrho)
            press.append(pre)
            q.append(pre / rrho)
    return q, rho, press


def read_peos_crust(path: Path) -> tuple[float, float]:
    with path.open() as f:
        f.readline()
        kappa_crust, gamma_crust = (float(x) for x in f.readline().split()[:2])
    return kappa_crust, gamma_crust


def crust_from_q(q: float, kappa: float, gamma: float) -> tuple[float, float]:
    fac1 = 1.0 / (gamma - 1.0)
    fac2 = gamma / (gamma - 1.0)
    fack = kappa**(-fac1)
    qloc = max(q, 1.0e-60)
    rho = fack * qloc**fac1
    press = fack * qloc**fac2
    return rho, press


def interpolate(table_x: list[float], table_f: list[float], x: float) -> float:
    if x <= table_x[0]:
        return table_f[0]
    idx = 1
    while idx < len(table_x) and x > table_x[idx]:
        idx += 1
    idx = min(max(idx, 2), len(table_x) - 2)
    start = idx - 2
    return lagrange4(table_x[start : start + 4], table_f[start : start + 4], x)


def expected_rho_press_from_q(
    q: float,
    q_table: list[float],
    rho_table: list[float],
    press_table: list[float],
    kappa_crust: float,
    gamma_crust: float,
) -> tuple[float, float]:
    if q <= q_table[0]:
        return crust_from_q(q, kappa_crust, gamma_crust)
    return interpolate(q_table, rho_table, q), interpolate(q_table, press_table, q)


def expected_from_rho(
    rho: float,
    q_table: list[float],
    rho_table: list[float],
    press_table: list[float],
    kappa_crust: float,
    gamma_crust: float,
) -> tuple[float, float, float]:
    if rho <= rho_table[0]:
        q = kappa_crust * rho ** (gamma_crust - 1.0)
        expected_rho, expected_press = expected_rho_press_from_q(
            q, q_table, rho_table, press_table, kappa_crust, gamma_crust
        )
        return q, expected_rho, expected_press

    iphase = 1
    for i in range(1, len(rho_table)):
        if (rho - rho_table[i]) * (rho - rho_table[i - 1]) <= 0.0:
            iphase = i
            break

    lo = q_table[iphase - 1]
    hi = q_table[iphase]

    def rho_residual(q: float) -> float:
        expected_rho, _expected_press = expected_rho_press_from_q(
            q, q_table, rho_table, press_table, kappa_crust, gamma_crust
        )
        return expected_rho - rho

    flo = rho_residual(lo)
    fhi = rho_residual(hi)
    if flo * fhi > 0.0:
        for i in range(max(0, iphase - 4), min(len(q_table) - 1, iphase + 4)):
            flo = rho_residual(q_table[i])
            fhi = rho_residual(q_table[i + 1])
            if flo * fhi <= 0.0:
                lo = q_table[i]
                hi = q_table[i + 1]
                break

    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fmid = rho_residual(mid)
        if flo * fmid <= 0.0:
            hi = mid
            fhi = fmid
        else:
            lo = mid
            flo = fmid

    q = 0.5 * (lo + hi)
    expected_rho, expected_press = expected_rho_press_from_q(
        q, q_table, rho_table, press_table, kappa_crust, gamma_crust
    )
    return q, expected_rho, expected_press


def scale(value: float, src_min: float, src_max: float, dst_min: float, dst_max: float) -> float:
    if src_max == src_min:
        return 0.5 * (dst_min + dst_max)
    return dst_min + (value - src_min) * (dst_max - dst_min) / (src_max - src_min)


def polyline(points: list[tuple[float, float]], color: str) -> str:
    coords = " ".join(f"{x:.3f},{y:.3f}" for x, y in points)
    return f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="1.5"/>'


def nice_ticks(low: float, high: float, max_ticks: int = 7) -> list[float]:
    if high == low:
        return [low]
    raw_step = abs(high - low) / max(1, max_ticks - 1)
    exponent = math.floor(math.log10(raw_step))
    scale_factor = 10.0**exponent
    for multiplier in (1.0, 2.0, 5.0, 10.0):
        step = multiplier * scale_factor
        if abs(high - low) / step <= max_ticks:
            break
    start = math.ceil(low / step) * step
    ticks: list[float] = []
    tick = start
    while tick <= high + 0.5 * step:
        if low - 1.0e-12 <= tick <= high + 1.0e-12:
            ticks.append(0.0 if abs(tick) < 1.0e-14 else tick)
        tick += step
    return ticks


def format_tick(value: float) -> str:
    if value == 0.0:
        return "0"
    if 1.0e-3 <= abs(value) < 1.0e4:
        return f"{value:g}"
    return f"{value:.1e}"


def write_svg_plot(
    path: Path,
    x_values: list[float],
    rho_values: list[float],
    abs_residuals: list[float],
) -> None:
    width = 900
    height = 420
    margin_left = 70
    margin_right = 30
    margin_top = 35
    margin_bottom = 55
    panel_gap = 65
    panel_width = (width - margin_left - margin_right - panel_gap) / 2
    panel_height = height - margin_top - margin_bottom

    safe_residuals = [max(v, 1.0e-300) for v in abs_residuals]
    log_residuals = [math.log10(v) for v in safe_residuals]
    y_min = math.floor(min(log_residuals))
    y_max = math.ceil(max(log_residuals))
    if y_max == y_min:
        y_max += 1

    rho_log = [math.log10(v) for v in rho_values]

    def panel(
        left: float,
        title: str,
        xlabel: str,
        xs: list[float],
        ys: list[float],
        x_tick_values: list[float] | None = None,
        x_tick_labels_override: list[str] | None = None,
    ) -> list[str]:
        x_min = min(xs)
        x_max = max(xs)
        x0 = left
        x1 = left + panel_width
        y0 = margin_top
        y1 = margin_top + panel_height
        pairs = sorted(
            (x, math.log10(max(y, 1.0e-300))) for x, y in zip(xs, ys)
        )
        pts = [
            (scale(x, x_min, x_max, x0, x1), scale(y, y_min, y_max, y1, y0))
            for x, y in pairs
        ]
        y_ticks = list(range(int(y_min), int(y_max) + 1))
        x_ticks = x_tick_values if x_tick_values is not None else nice_ticks(x_min, x_max)
        y_tick_lines = [
            f'<line x1="{x0:.1f}" y1="{scale(tick, y_min, y_max, y1, y0):.3f}" '
            f'x2="{x1:.1f}" y2="{scale(tick, y_min, y_max, y1, y0):.3f}" '
            'stroke="#ddd" stroke-width="1"/>'
            for tick in y_ticks
        ]
        x_tick_lines = [
            f'<line x1="{scale(tick, x_min, x_max, x0, x1):.3f}" y1="{y0:.1f}" '
            f'x2="{scale(tick, x_min, x_max, x0, x1):.3f}" y2="{y1:.1f}" '
            'stroke="#eee" stroke-width="1"/>'
            for tick in x_ticks
        ]
        y_tick_labels = [
            f'<text x="{x0 - 8:.1f}" y="{scale(tick, y_min, y_max, y1, y0) + 4:.1f}" '
            f'text-anchor="end" font-size="11">1e{tick}</text>'
            for tick in y_ticks
        ]
        x_tick_text = (
            x_tick_labels_override
            if x_tick_labels_override is not None
            else [format_tick(tick) for tick in x_ticks]
        )
        x_tick_labels = [
            f'<line x1="{scale(tick, x_min, x_max, x0, x1):.3f}" y1="{y1:.1f}" '
            f'x2="{scale(tick, x_min, x_max, x0, x1):.3f}" y2="{y1 + 5:.1f}" '
            'stroke="black" stroke-width="1"/>'
            f'<text x="{scale(tick, x_min, x_max, x0, x1):.3f}" y="{y1 + 18:.1f}" '
            f'text-anchor="middle" font-size="11">{label}</text>'
            for tick, label in zip(x_ticks, x_tick_text)
        ]
        return [
            f'<text x="{0.5 * (x0 + x1):.1f}" y="22" text-anchor="middle" font-size="15">{title}</text>',
            f'<rect x="{x0:.1f}" y="{y0:.1f}" width="{panel_width:.1f}" height="{panel_height:.1f}" fill="white" stroke="black"/>',
            *y_tick_lines,
            *x_tick_lines,
            polyline(pts, "#1f77b4"),
            *[
                f'<circle cx="{x:.3f}" cy="{y:.3f}" r="2.2" fill="#1f77b4"/>'
                for x, y in pts
            ],
            f'<text x="{0.5 * (x0 + x1):.1f}" y="{height - 16}" text-anchor="middle" font-size="13">{xlabel}</text>',
            *y_tick_labels,
            *x_tick_labels,
        ]

    left_panel = panel(margin_left, "TEOS error vs x", "x", x_values, abs_residuals)
    rho_display = rho_log
    rho_tick_values = nice_ticks(min(rho_log), max(rho_log))
    rho_tick_labels = [format_tick(tick) for tick in rho_tick_values]
    rho_xlabel = "log10(rho)"
    rho_sorted_values = sorted(set(rho_log))
    if len(rho_sorted_values) > 1:
        gaps = [
            (rho_sorted_values[i + 1] - rho_sorted_values[i], i)
            for i in range(len(rho_sorted_values) - 1)
        ]
        gap, gap_index = max(gaps)
        if gap > 20.0:
            low_min = rho_sorted_values[0]
            low_max = rho_sorted_values[gap_index]
            high_min = rho_sorted_values[gap_index + 1]
            high_max = rho_sorted_values[-1]
            low_span = max(low_max - low_min, 1.0)
            high_span = max(high_max - high_min, 1.0)
            high_offset = 1.0

            def broken_log_rho(value: float) -> float:
                if value <= low_max:
                    return 0.4 * (value - low_min) / low_span
                return high_offset + (value - high_min) / high_span

            rho_display = [broken_log_rho(value) for value in rho_log]
            high_tick_logs = nice_ticks(high_min, high_max)
            rho_tick_logs = [low_min, *high_tick_logs]
            rho_tick_values = [broken_log_rho(value) for value in rho_tick_logs]
            rho_tick_labels = [format_tick(value) for value in rho_tick_logs]
            rho_xlabel = "log10(rho)"
    right_panel = panel(
        margin_left + panel_width + panel_gap,
        "TEOS error vs rho",
        rho_xlabel,
        rho_display,
        abs_residuals,
        rho_tick_values,
        rho_tick_labels,
    )

    svg = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="18" y="{0.5 * height:.1f}" transform="rotate(-90 18 {0.5 * height:.1f})" text-anchor="middle" font-size="13">|P / P_TEOS(rho) - 1|, log scale</text>',
        *left_panel,
        *right_panel,
        "</svg>",
    ]
    path.write_text("\n".join(svg) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--press", type=Path, required=True)
    parser.add_argument("--teos", type=Path, required=True)
    parser.add_argument("--peos", type=Path)
    parser.add_argument("--iteration", type=int, default=0)
    parser.add_argument("--rho-cut", type=float, default=0.0)
    parser.add_argument("--rtol", type=float)
    parser.add_argument("--plot", type=Path, help="Write an SVG residual plot")
    args = parser.parse_args()

    rho_values = read_ascii_1d(args.rho, args.iteration)
    press_values = read_ascii_1d(args.press, args.iteration)
    q_table, rho_table, press_table = read_teos(args.teos)
    peos_path = args.peos if args.peos is not None else args.teos.with_name("peos_parameter.dat")
    kappa_crust, gamma_crust = read_peos_crust(peos_path)

    checked = 0
    worst_rel = -1.0
    worst_abs = -1.0
    worst_coord: tuple[float, float, float] | None = None
    worst_kind = ""
    worst_pair: tuple[float, float] | None = None
    max_press_rho_rel = -1.0
    max_rho_reconstruction_rel = -1.0
    plot_x: list[float] = []
    plot_rho: list[float] = []
    plot_residual: list[float] = []

    for coord, rho in rho_values.items():
        if rho <= args.rho_cut:
            continue
        press = press_values[coord]
        _q, expected_rho, expected_press_of_rho = expected_from_rho(
            rho, q_table, rho_table, press_table, kappa_crust, gamma_crust
        )
        press_rho_residual = press / expected_press_of_rho - 1.0
        press_rho_abs_err = abs(press - expected_press_of_rho)
        rho_abs_err = abs(rho - expected_rho)
        press_rho_rel_err = press_rho_abs_err / max(abs(expected_press_of_rho), 1.0e-300)
        rho_rel_err = rho_abs_err / max(abs(expected_rho), 1.0e-300)
        max_press_rho_rel = max(max_press_rho_rel, press_rho_rel_err)
        max_rho_reconstruction_rel = max(max_rho_reconstruction_rel, rho_rel_err)
        checked += 1
        plot_x.append(coord[0])
        plot_rho.append(rho)
        plot_residual.append(abs(press_rho_residual))
        if press_rho_rel_err > worst_rel:
            worst_rel = press_rho_rel_err
            worst_abs = press_rho_abs_err
            worst_coord = coord
            worst_kind = "press(rho->q->press)"
            worst_pair = (press, expected_press_of_rho)
        if rho_rel_err > worst_rel:
            worst_rel = rho_rel_err
            worst_abs = rho_abs_err
            worst_coord = coord
            worst_kind = "rho(rho->q)"
            worst_pair = (rho, expected_rho)

    if checked == 0:
        raise SystemExit("No non-atmosphere points were checked")

    print(f"checked_points {checked}")
    print(f"max_relative_error {worst_rel:.16e}")
    print(f"max_press_rho_relative_error {max_press_rho_rel:.16e}")
    print(f"max_rho_reconstruction_relative_error {max_rho_reconstruction_rel:.16e}")
    print(f"max_absolute_error {worst_abs:.16e}")
    print(f"worst_quantity {worst_kind}")
    print(f"worst_coord {worst_coord}")
    print(f"worst_output {worst_pair[0]:.16e}")
    print(f"worst_teos {worst_pair[1]:.16e}")

    if args.plot is not None:
        order = sorted(range(len(plot_x)), key=plot_x.__getitem__)
        write_svg_plot(
            args.plot,
            [plot_x[i] for i in order],
            [plot_rho[i] for i in order],
            [plot_residual[i] for i in order],
        )
        print(f"plot {args.plot}")

    if args.rtol is not None and worst_rel > args.rtol:
        raise SystemExit(f"TEOS check failed: relative error {worst_rel:.6e} > {args.rtol:.6e}")


if __name__ == "__main__":
    main()
