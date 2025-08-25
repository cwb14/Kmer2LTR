#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from typing import List, Tuple, Dict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # <-- added

# ----------------------------
# Column map (0-based indices)
# ----------------------------
MODEL_COLS = {
    "p":    (6, 7),    # p-dist, p-time
    "JC69": (8, 9),    # JC69-dist, JC69-time
    "K2P":  (10, 11),  # K2P-dist, K2P-time
}

# ----------------------------
# IO utilities
# ----------------------------
def read_tsv_cols(path: str, dist_idx: int, time_idx: int) -> Tuple[np.ndarray, np.ndarray]:
    """Read two float columns (distance, time) from a TSV with no header."""
    dists = []
    times = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) <= max(dist_idx, time_idx):
                continue
            try:
                d = float(parts[dist_idx])
                t = float(parts[time_idx])
            except ValueError:
                continue
            dists.append(d)
            times.append(t)
    if not dists:
        return np.array([]), np.array([])
    return np.asarray(dists, dtype=float), np.asarray(times, dtype=float)

def file_prefix_for_legend(path: str) -> str:
    """Legend name = substring before first '.' in basename."""
    base = os.path.basename(path)
    return base.split('.', 1)[0] if '.' in base else base

# ----------------------------
# KDE (no external deps)
# ----------------------------
def silverman_bandwidth(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 2:
        return 1.0
    sd = np.std(x, ddof=1)
    if sd == 0:
        return max(1e-6, np.abs(x[0]) * 1e-3 + 1e-6)
    return 1.06 * sd * (n ** (-1 / 5))

def kde_gaussian(x: np.ndarray, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    if x.size == 0:
        return np.zeros_like(grid)
    h = max(bandwidth, 1e-12)
    z = (grid[:, None] - x[None, :]) / h
    K = np.exp(-0.5 * z ** 2) / np.sqrt(2 * np.pi)
    dens = K.sum(axis=1) / (x.size * h)
    return dens

# ----------------------------
# Axis mapping: distance <-> time (MYA)
# ----------------------------
def compute_axis_mapping_mya(all_dists: np.ndarray, all_times_years: np.ndarray):
    """
    Build linear mapping between distance (bottom axis) and time (top axis, MYA):
       time_MYA = (a*dist + b) / 1e6
    using min/max alignment across all inputs.
    """
    if all_dists.size == 0 or all_times_years.size == 0:
        return (lambda x: x, lambda y: y)

    dmin, dmax = float(np.min(all_dists)), float(np.max(all_dists))
    tmin, tmax = float(np.min(all_times_years)), float(np.max(all_times_years))

    if dmax == dmin:
        dmax = dmin + 1e-9
    if tmax == tmin:
        tmax = tmin + 1e-9

    a = (tmax - tmin) / (dmax - dmin)
    b = tmin - a * dmin

    def forward(dist_vals):         # dist -> time_MYA
        return (a * dist_vals + b) / 1e6

    def inverse(time_mya_vals):     # time_MYA -> dist
        return (time_mya_vals * 1e6 - b) / a

    return forward, inverse

# ----------------------------
# Peak finding (simple, no scipy)
# ----------------------------
def find_top_peaks(xgrid: np.ndarray, y: np.ndarray, k: int) -> List[Tuple[float, float, int]]:
    """
    Return up to k peaks as (x_peak, y_peak, idx) sorted by descending height.
    A peak is y[i-1] < y[i] > y[i+1]. Ties resolved by height.
    """
    if y.size < 3:
        return []
    peaks = []
    for i in range(1, len(y) - 1):
        if y[i] > y[i-1] and y[i] > y[i+1]:
            peaks.append((xgrid[i], y[i], i))
    if not peaks:
        j = int(np.argmax(y))
        return [(xgrid[j], y[j], j)]
    peaks.sort(key=lambda tup: tup[1], reverse=True)
    return peaks[:max(1, k)]

# ----------------------------
# Plotting
# ----------------------------
def plot_density(
    inputs: List[str],
    model: str,
    miu_str: str,
    out_pdf: str,
    xmax_override: float = None,
    label_peaks: int = 0,  # 0 = off, >0 = number of peaks to mark with vlines & legend entries
):
    dist_idx, time_idx = MODEL_COLS[model]

    series: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    all_dists = []
    all_times_years = []

    # Load data
    for path in inputs:
        d, t = read_tsv_cols(path, dist_idx, time_idx)
        if d.size == 0:
            print(f"[warn] No usable rows in: {path}", file=sys.stderr)
            continue
        series[path] = (d, t)
        all_dists.append(d)
        all_times_years.append(t)

    if not series:
        raise SystemExit("No valid data found in provided input file(s).")

    all_dists_cat = np.concatenate(all_dists)
    all_times_cat = np.concatenate(all_times_years)

    # Build global grid across all datasets for consistent comparison
    xmin = float(np.min(all_dists_cat))
    xmax = float(np.max(all_dists_cat))

    # Apply optional x-max override (distance axis)
    if xmax_override is not None:
        if xmax_override <= xmin:
            raise SystemExit(f"--xmax ({xmax_override}) must be greater than data min ({xmin}).")
        xmax = float(xmax_override)

    # Small padding
    if xmin == xmax:
        xmax = xmin + 1e-9
    xpad = 0.02 * (xmax - xmin) if xmax > xmin else 0.01
    xmin_plot = xmin - xpad
    xmax_plot = xmax + xpad

    xgrid = np.linspace(xmin_plot, xmax_plot, 800)

    # Compute densities and draw curves
    densities = []
    labels = []
    line_colors = []
    plt.figure(figsize=(8.8, 5.8))
    ax = plt.gca()

    for path, (d, _t) in series.items():
        bw = silverman_bandwidth(d)
        dens = kde_gaussian(d, xgrid, bw)
        line, = ax.plot(xgrid, dens, linewidth=2, label=file_prefix_for_legend(path))
        densities.append(dens)
        labels.append(file_prefix_for_legend(path))
        line_colors.append(line.get_color())  # capture cycle color

    ax.set_xlim(xmin_plot, xmax_plot)
    ax.set_ylabel("Density")
    ax.set_xlabel(f"Genetic distance ({model})")

    # Secondary X-axis (time in MYA)
    forward, inverse = compute_axis_mapping_mya(all_dists_cat, all_times_cat)
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel("Time (MYA)")

    # Title with μ
    plt.title(f"Density of {model} divergence (μ = {miu_str})")

    if label_peaks and label_peaks > 0:
        # Draw vertical dashed lines at peaks, and build legend entries as "<label> - x1, x2, ..."
        legend_items = []  # (primary_peak_x, label_string, color)
        for dens, label, color in zip(densities, labels, line_colors):
            peaks = find_top_peaks(xgrid, dens, k=label_peaks)
            if not peaks:
                # nothing to draw / add
                continue

            # Sort peaks by x (distance) ascending just for a cleaner legend per series
            peaks_sorted_by_x = sorted(peaks, key=lambda p: p[0])

            # Draw vlines for each peak using the same color as the density
            for (xp, _yp, _idx) in peaks_sorted_by_x:
                ax.axvline(x=xp, linestyle="--", linewidth=1, color=color, alpha=0.8)

            # Build legend text: "LABEL - 0.0197, 0.0421"
            peak_strs = [f"{xp:.4f}" for (xp, _yp, _idx) in peaks_sorted_by_x]
            legend_text = f"{label} - {', '.join(peak_strs)}"

            # Use the tallest peak's x as the sorting key across series
            tallest_peak_x = peaks[0][0]  # find_top_peaks returns peaks sorted by height desc
            legend_items.append((tallest_peak_x, legend_text, color))

        # Sort legend by ascending peak value
        legend_items.sort(key=lambda tup: tup[0])  # ascending by primary peak x

        # Build proxy handles so legend color matches the series color
        handles = [Line2D([0], [0], color=col, lw=2) for (_x, _txt, col) in legend_items]
        texts = [txt for (_x, txt, _col) in legend_items]
        ax.legend(handles=handles, labels=texts, loc="upper right", frameon=True)
    else:
        # Original legend when --label-peaks is off
        ax.legend(loc="upper right", frameon=True)

    plt.tight_layout()
    plt.savefig(out_pdf, format="pdf")
    plt.close()
    print(f"[ok] Wrote {out_pdf}")

# ----------------------------
# CLI
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Plot density of LTR-RT divergence (distance bottom X-axis, time top X-axis in MYA)."
    )
    parser.add_argument(
        "-in", dest="inputs", nargs="+", required=True,
        help="One or more TSV files like '<prefix>.LTRs.alns.results'."
    )
    parser.add_argument(
        "-model", choices=["K2P", "JC69", "p"], default="K2P",
        help="Distance/time model to use for axes. Default: K2P."
    )
    parser.add_argument(
        "-miu", default="(not specified)",
        help="Mutation rate (string). Used in the figure title only. Example: 3e-8"
    )
    parser.add_argument(
        "-out", default=None,
        help="Output PDF filename (optional). Default: density_<MODEL>.pdf"
    )
    parser.add_argument(
        "--xmax", type=float, default=None,
        help="Custom max for the distance X-axis (bottom). Must be greater than the observed min."
    )
    parser.add_argument(
        "--label-peaks", nargs="?", const="1", default="0",
        help="Mark top N peaks per curve with vertical dashed lines, and put peak values in the legend "
             "(sorted ascending by the tallest peak). Default N=1 when flag present. "
             "If omitted, no peak markers."
    )

    args = parser.parse_args()
    out_pdf = args.out if args.out else f"density_{args.model}.pdf"

    # parse label-peaks value
    try:
        label_peaks = int(args.label_peaks)
    except ValueError:
        print("[error] --label-peaks must be an integer if provided with a value.", file=sys.stderr)
        sys.exit(1)

    try:
        plot_density(
            inputs=args.inputs,
            model=args.model,
            miu_str=args.miu,
            out_pdf=out_pdf,
            xmax_override=args.xmax,
            label_peaks=label_peaks
        )
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
