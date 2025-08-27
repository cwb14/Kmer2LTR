#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from typing import List, Tuple, Dict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # <-- used for legend proxies

# ----------------------------
# Column map (0-based indices)
# ----------------------------
MODEL_COLS = {
    "p":    (6, 7),    # p-dist, p-time
    "JC69": (8, 9),    # JC69-dist, JC69-time
    "K2P":  (10, 11),  # K2P-dist, K2P-time
}

SUMMARY_KEYS = {
    "p": "raw_d",
    "JC69": "JC69_d",
    "K2P": "K2P_d",
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

def read_summary_value(results_path: str, model: str):
    """
    Read the model-specific distance from '<results_path>.summary'.
    Returns float or None if not found/unreadable.
    """
    sum_path = results_path + ".summary"
    key = SUMMARY_KEYS[model]
    if not os.path.exists(sum_path):
        print(f"[warn] Summary file not found for {results_path}: {sum_path}", file=sys.stderr)
        return None
    try:
        with open(sum_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) != 2:
                    continue
                k, v = parts
                if k == key:
                    try:
                        return float(v)
                    except ValueError:
                        return None
    except Exception as e:
        print(f"[warn] Failed reading {sum_path}: {e}", file=sys.stderr)
    return None

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
    label_peaks: int = 0,          # 0 = off, >0 = number of peaks to mark
    label_summary: bool = False,   # mark vline at model-specific summary value
    no_legend: bool = False,       # suppress legend entirely
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
    plt.figure(figsize=(8.8, 5.8))
    ax = plt.gca()

    # Track series metadata for legend assembly
    per_series = {}  # path -> dict(color, label, peak_xs (list), summary (float|None), sort_key)

    for path, (d, _t) in series.items():
        bw = silverman_bandwidth(d)
        dens = kde_gaussian(d, xgrid, bw)
        line, = ax.plot(xgrid, dens, linewidth=2, label=file_prefix_for_legend(path))
        per_series[path] = {
            "color": line.get_color(),
            "label": file_prefix_for_legend(path),
            "peak_xs": [],
            "summary": None,
            "sort_key": None,
            "dens": dens  # keep to find peaks
        }

    ax.set_xlim(xmin_plot, xmax_plot)
    ax.set_ylabel("Density")
    ax.set_xlabel(f"Genetic distance ({model})")

    # Secondary X-axis (time in MYA)
    forward, inverse = compute_axis_mapping_mya(all_dists_cat, all_times_cat)
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel("Time (MYA)")

    # Title with μ
    plt.title(f"Density of {model} divergence (μ = {miu_str})")

    # Summary-based vlines per series
    if label_summary:
        for path in per_series.keys():
            v = read_summary_value(path, model)
            if v is None:
                continue
            ax.axvline(x=v, linestyle="--", linewidth=1.5, color=per_series[path]["color"], alpha=0.85)
            per_series[path]["summary"] = v

    # Peak-based vlines and capture peak positions
    if label_peaks and label_peaks > 0:
        for path, meta in per_series.items():
            dens = meta["dens"]
            peaks = find_top_peaks(xgrid, dens, k=label_peaks)
            if not peaks:
                continue
            # sort peaks by x for prettier legend text
            peaks_sorted_by_x = sorted(peaks, key=lambda p: p[0])
            xs = [xp for (xp, _yp, _idx) in peaks_sorted_by_x]
            for xp in xs:
                ax.axvline(x=xp, linestyle="--", linewidth=1, color=meta["color"], alpha=0.8)
            meta["peak_xs"] = xs
            # sort key = tallest peak's x (find_top_peaks returns peaks sorted by height desc)
            meta["sort_key"] = peaks[0][0]

    # ----------------------------
    # Legend construction
    # ----------------------------
    if not no_legend:
        any_values = any((meta["peak_xs"] or (meta["summary"] is not None)) for meta in per_series.values())
        if any_values:
            # If any peaks exist, sort by tallest-peak x; else if only summaries, sort by summary.
            if any(meta["peak_xs"] for meta in per_series.values()):
                sort_items = sorted(per_series.values(), key=lambda m: (float('inf') if m["sort_key"] is None else m["sort_key"]))
            else:
                sort_items = sorted(per_series.values(), key=lambda m: (float('inf') if m["summary"] is None else m["summary"]))
            handles = []
            texts = []
            for meta in sort_items:
                color = meta["color"]
                label = meta["label"]
                peaks_text = ""
                summary_text = ""
                if meta["peak_xs"]:
                    peaks_text = ", ".join(f"{x:.4f}" for x in meta["peak_xs"])
                if meta["summary"] is not None:
                    summary_text = f"{meta['summary']:.4f}"

                # Compose legend line based on available parts:
                if meta["peak_xs"] and (meta["summary"] is not None):
                    # both present
                    text = f"{label} - {peaks_text}; summary: {summary_text}"
                elif meta["peak_xs"]:
                    text = f"{label} - {peaks_text}"
                elif meta["summary"] is not None:
                    text = f"{label} - {summary_text}"
                else:
                    text = label  # fallback (shouldn't happen if any_values is True)

                handles.append(Line2D([0], [0], color=color, lw=2))
                texts.append(text)

            ax.legend(handles=handles, labels=texts, loc="upper right", frameon=True)
        else:
            # Default legend = series names
            ax.legend(loc="upper right", frameon=True)

    # Cleanup: drop stored densities to avoid accidental reuse
    for meta in per_series.values():
        meta.pop("dens", None)

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
        help=("Mark top N peaks per curve with vertical dashed lines, and put peak values in the legend "
              "(sorted ascending by the tallest peak). Default N=1 when flag present. "
              "If omitted, no peak markers.")
    )
    parser.add_argument(
        "--label-summary", action="store_true",
        help=("Draw a vertical dashed line per input at the model-specific value read from "
              "'<input>.summary' (keys: raw_d, JC69_d, K2P_d), and include that value in the legend.")
    )
    parser.add_argument(
        "--no-legend", action="store_true",
        help="Suppress the legend entirely (vlines will still be shown if requested)."
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
            label_peaks=label_peaks,
            label_summary=args.label_summary,
            no_legend=args.no_legend,
        )
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
