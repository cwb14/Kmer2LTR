#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot histogram density of LTR-RT divergence.

Bottom X-axis: genetic distance (p, JC69, or K2P)
Top X-axis:    time in millions of years ago (MYA)
Y-axis:        density (histogram area = 1)

Peak detection uses histogram bin heights directly, with optional
moving-average smoothing and prominence filtering to suppress noise
while faithfully representing the data.
"""

import argparse
import os
import sys
from typing import List, Tuple, Dict, Optional, Union
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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


def read_summary_value(results_path: str, model: str) -> Optional[float]:
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
# Auto bin calculation
# ----------------------------
def auto_bin_count(data: np.ndarray, xmin: float, xmax: float) -> int:
    """
    Compute bin count using a tightened Freedman-Diaconis rule, tuned
    for right-skewed genomic distance data.

        bin_width = 0.75 * IQR * n^(-1/3)

    The standard FD multiplier (2.0) produces overly coarse bins for
    skewed distributions where detail near the mode matters.  The 0.75
    multiplier gives roughly 2.7× finer resolution while still being
    data-adaptive (IQR-based, robust to outliers).

    Falls back to Sturges' rule (1 + log2(n)) if IQR is zero (e.g. many
    identical values), and clamps the result to [10, 500].

    Parameters
    ----------
    data : concatenated distance values across all input series
    xmin, xmax : the plot range (used to convert width -> count)

    Returns
    -------
    int : number of bins
    """
    n = data.size
    if n < 2:
        return 10

    q75, q25 = np.percentile(data, [75, 25])
    iqr = q75 - q25
    span = xmax - xmin
    if span <= 0:
        return 10

    if iqr > 0:
        # Tightened Freedman-Diaconis
        bin_width = 1.0 * iqr * (n ** (-1.0 / 3.0))
        # bin_width = 0.75 * iqr * (n ** (-1.0 / 3.0)) I adjust to slightly larger bins with 1.0
        nbins = max(1, int(np.ceil(span / bin_width)))
    else:
        # Fallback: Sturges
        nbins = max(1, int(np.ceil(1.0 + np.log2(n))))

    # Clamp to a sensible range
    nbins = max(10, min(nbins, 500))

    return nbins


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
# Smoothing & peak finding
# ----------------------------
def moving_average(y: np.ndarray, window: int) -> np.ndarray:
    """
    Centred moving average.  Window is forced to an odd number >= 1.
    Edge bins use a smaller (truncated) window so output length = input length.
    """
    w = max(1, window)
    if w % 2 == 0:
        w += 1  # force odd for symmetric centring
    if w == 1 or y.size <= 1:
        return y.copy()
    half = w // 2
    out = np.empty_like(y, dtype=float)
    for i in range(len(y)):
        lo = max(0, i - half)
        hi = min(len(y), i + half + 1)
        out[i] = y[lo:hi].mean()
    return out


def compute_prominence(y: np.ndarray, idx: int) -> float:
    """
    Prominence of a peak at index `idx`:
      Walk left until a higher value or the edge -> track the min along the way.
      Walk right until a higher value or the edge -> track the min along the way.
      prominence = y[idx] - max(left_min, right_min)
    """
    h = y[idx]
    # Walk left
    left_min = h
    for j in range(idx - 1, -1, -1):
        if y[j] > h:
            break
        left_min = min(left_min, y[j])
    else:
        left_min = min(left_min, y[0])
    # Walk right
    right_min = h
    for j in range(idx + 1, len(y)):
        if y[j] > h:
            break
        right_min = min(right_min, y[j])
    else:
        right_min = min(right_min, y[-1])
    return h - max(left_min, right_min)


def find_histogram_peaks(
    bin_centers: np.ndarray,
    bin_heights: np.ndarray,
    k: int,
    smooth_window: int = 3,
    min_prominence: float = 0.0,
) -> List[Tuple[float, float, int]]:
    """
    Find up to `k` peaks from histogram bin heights.

    1. Smooth bin heights with a centred moving average of width `smooth_window`.
    2. Identify local maxima (higher than both neighbours) on the smoothed curve.
    3. Discard peaks whose prominence (on the smoothed curve) is below
       `min_prominence`.  0 means keep all local maxima.
    4. Return top k by prominence as (bin_center, raw_height, index),
       sorted by descending prominence.

    Peak *positions* (bin centres) and reported *heights* always come from the
    original (unsmoothed) histogram so they faithfully represent the data.
    """
    if bin_heights.size < 3:
        if bin_heights.size == 0:
            return []
        j = int(np.argmax(bin_heights))
        return [(bin_centers[j], bin_heights[j], j)]

    sm = moving_average(bin_heights, smooth_window)

    # Local maxima on the smoothed curve
    candidates = []
    for i in range(1, len(sm) - 1):
        if sm[i] > sm[i - 1] and sm[i] > sm[i + 1]:
            prom = compute_prominence(sm, i)
            if prom >= min_prominence:
                candidates.append((bin_centers[i], bin_heights[i], i, prom))

    if not candidates:
        # Fallback: report the global maximum of the raw histogram
        j = int(np.argmax(bin_heights))
        return [(bin_centers[j], bin_heights[j], j)]

    # Sort by prominence (primary) then height (secondary), descending
    candidates.sort(key=lambda c: (c[3], c[1]), reverse=True)

    return [(c[0], c[1], c[2]) for c in candidates[:max(1, k)]]


# ----------------------------
# Plotting
# ----------------------------
def plot_density(
    inputs: List[str],
    model: str,
    miu_str: str,
    out_pdf: str,
    xmax_override: float = None,
    label_peaks: int = 0,
    label_summary: bool = False,
    no_legend: bool = False,
    bins_arg: Union[int, str] = 50,   # int or "auto"
    smooth_window: int = 3,
    min_prominence: float = 0.0,
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

    # X-axis range
    xmin = float(np.min(all_dists_cat))
    xmax = float(np.max(all_dists_cat))

    if xmax_override is not None:
        if xmax_override <= xmin:
            raise SystemExit(f"--xmax ({xmax_override}) must be greater than data min ({xmin}).")
        xmax = float(xmax_override)

    if xmin == xmax:
        xmax = xmin + 1e-9
    xpad = 0.02 * (xmax - xmin) if xmax > xmin else 0.01
    xmin_plot = xmin - xpad
    xmax_plot = xmax + xpad

    # Resolve bin count — "auto" uses Freedman-Diaconis on all data combined
    if isinstance(bins_arg, str) and bins_arg.lower() == "auto":
        bins = auto_bin_count(all_dists_cat, xmin_plot, xmax_plot)
        print(f"[info] Auto bins (tightened Freedman-Diaconis): {bins}", file=sys.stderr)
    else:
        bins = int(bins_arg)

    # Time mapping (must be built before peak labelling)
    forward, inverse = compute_axis_mapping_mya(all_dists_cat, all_times_cat)

    # ------- Figure -------
    plt.figure(figsize=(8.8, 5.8))
    ax = plt.gca()

    per_series: Dict[str, dict] = {}

    for path, (d, _t) in series.items():
        # Draw histogram — all series use the same bins/range
        heights, edges, patches = ax.hist(
            d,
            bins=bins,
            density=True,
            range=(xmin_plot, xmax_plot),
            alpha=0.45,
            edgecolor="white",
            linewidth=0.4,
            label=file_prefix_for_legend(path),
        )
        color = patches[0].get_facecolor()

        bin_centers = 0.5 * (edges[:-1] + edges[1:])

        per_series[path] = {
            "color": color,
            "label": file_prefix_for_legend(path),
            "bin_centers": bin_centers,
            "bin_heights": heights,
            "peak_xs": [],
            "peak_myas": [],
            "summary": None,
            "sort_key": None,
        }

    ax.set_xlim(xmin_plot, xmax_plot)
    ax.set_ylabel("Density")
    ax.set_xlabel(f"Genetic distance ({model})")

    # Secondary X-axis (time in MYA)
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel("Time (MYA)")

    # Title
    plt.title(f"Density of {model} divergence (μ = {miu_str})")

    # Summary vlines
    if label_summary:
        for path in per_series:
            v = read_summary_value(path, model)
            if v is None:
                continue
            ax.axvline(x=v, linestyle="--", linewidth=1.5,
                       color=per_series[path]["color"], alpha=0.85)
            per_series[path]["summary"] = v

    # Peak detection and vlines
    if label_peaks and label_peaks > 0:
        for path, meta in per_series.items():
            peaks = find_histogram_peaks(
                meta["bin_centers"],
                meta["bin_heights"],
                k=label_peaks,
                smooth_window=smooth_window,
                min_prominence=min_prominence,
            )
            if not peaks:
                continue
            # Sort by x for legend readability
            peaks_sorted = sorted(peaks, key=lambda p: p[0])
            xs = [xp for xp, _h, _i in peaks_sorted]
            myas = [float(forward(np.array([xp]))[0]) for xp in xs]

            for xp in xs:
                ax.axvline(x=xp, linestyle="--", linewidth=1.2,
                           color=meta["color"], alpha=0.8)

            meta["peak_xs"] = xs
            meta["peak_myas"] = myas
            tallest_x = peaks[0][0]
            meta["sort_key"] = tallest_x

    # ------- Legend -------
    if not no_legend:
        any_values = any(
            meta["peak_xs"] or (meta["summary"] is not None)
            for meta in per_series.values()
        )
        if any_values:
            if any(meta["peak_xs"] for meta in per_series.values()):
                sort_items = sorted(
                    per_series.values(),
                    key=lambda m: m["sort_key"] if m["sort_key"] is not None else float('inf'),
                )
            else:
                sort_items = sorted(
                    per_series.values(),
                    key=lambda m: m["summary"] if m["summary"] is not None else float('inf'),
                )

            handles = []
            texts = []
            for meta in sort_items:
                color = meta["color"]
                label = meta["label"]
                peaks_str = ""
                summary_str = ""

                if meta["peak_xs"]:
                    parts = []
                    for xp, mya in zip(meta["peak_xs"], meta["peak_myas"]):
                        parts.append(f"{xp:.4f} ({mya:.3f} MYA)")
                    peaks_str = ", ".join(parts)

                if meta["summary"] is not None:
                    summary_mya = float(forward(np.array([meta["summary"]]))[0])
                    summary_str = f"{meta['summary']:.4f} ({summary_mya:.3f} MYA)"

                if meta["peak_xs"] and meta["summary"] is not None:
                    text = f"{label} - {peaks_str}; summary: {summary_str}"
                elif meta["peak_xs"]:
                    text = f"{label} - {peaks_str}"
                elif meta["summary"] is not None:
                    text = f"{label} - {summary_str}"
                else:
                    text = label

                handles.append(Line2D([0], [0], color=color, lw=2))
                texts.append(text)

            ax.legend(handles=handles, labels=texts, loc="upper right", frameon=True)
        else:
            ax.legend(loc="upper right", frameon=True)

    plt.tight_layout()
    plt.savefig(out_pdf, format="pdf")
    plt.close()
    print(f"[ok] Wrote {out_pdf}")


# ----------------------------
# CLI
# ----------------------------
def parse_bins(value: str) -> Union[int, str]:
    """Accept 'auto' or a positive integer for --bins."""
    if value.lower() == "auto":
        return "auto"
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--bins must be 'auto' or a positive integer, got '{value}'"
        )
    if n < 1:
        raise argparse.ArgumentTypeError("--bins must be >= 1")
    return n


def main():
    parser = argparse.ArgumentParser(
        description="Plot histogram density of LTR-RT divergence "
                    "(distance bottom X-axis, time top X-axis in MYA)."
    )
    parser.add_argument(
        "-in", dest="inputs", nargs="+", required=True,
        help="One or more TSV files like '<prefix>.LTRs.alns.results'.",
    )
    parser.add_argument(
        "-model", choices=["K2P", "JC69", "p"], default="K2P",
        help="Distance/time model to use for axes. Default: K2P.",
    )
    parser.add_argument(
        "-miu", default="(not specified)",
        help="Mutation rate (string). Used in the figure title only. Example: 3e-8",
    )
    parser.add_argument(
        "-out", default=None,
        help="Output PDF filename (optional). Default: density_<MODEL>.pdf",
    )
    parser.add_argument(
        "--xmax", type=float, default=None,
        help="Custom max for the distance X-axis (bottom). Must be > observed min.",
    )
    parser.add_argument(
        "--label-peaks", nargs="?", const="1", default="0",
        help=("Mark top N peaks per curve with vertical dashed lines and include "
              "peak distance + MYA in the legend. Default N=1 when flag present. "
              "If omitted, no peak markers."),
    )
    parser.add_argument(
        "--label-summary", action="store_true",
        help=("Draw a vertical dashed line per input at the model-specific value "
              "read from '<input>.summary' (keys: raw_d, JC69_d, K2P_d), and "
              "include that value in the legend."),
    )
    parser.add_argument(
        "--no-legend", action="store_true",
        help="Suppress the legend entirely (vlines still shown if requested).",
    )
    parser.add_argument(
        "--bins", type=parse_bins, default=50,
        help=("Number of histogram bins, or 'auto' to select automatically "
              "using the Freedman-Diaconis rule (robust for skewed data). "
              "When 'auto', bin width is computed from ALL input data combined "
              "so every series shares the same bin edges. Default: 50."),
    )
    parser.add_argument(
        "--smooth", type=int, default=3,
        help=("Moving-average window (in bins) used to smooth histogram heights "
              "ONLY for peak detection. Must be >= 1. Does NOT change the plotted "
              "histogram — only affects which bin is identified as a peak. "
              "Use 1 for no smoothing (raw bin heights). Default: 3."),
    )
    parser.add_argument(
        "--min-prominence", type=float, default=0.0,
        help=("Minimum prominence a peak must have to be reported. "
              "Prominence = how far a peak rises above the highest saddle "
              "connecting it to a taller peak. Filters out noise peaks. "
              "Use 0 to accept all local maxima (after smoothing). Default: 0."),
    )

    args = parser.parse_args()
    out_pdf = args.out if args.out else f"density_{args.model}.pdf"

    try:
        label_peaks = int(args.label_peaks)
    except ValueError:
        print("[error] --label-peaks must be an integer.", file=sys.stderr)
        sys.exit(1)

    if args.smooth < 1:
        print("[error] --smooth must be >= 1.", file=sys.stderr)
        sys.exit(1)

    if args.min_prominence < 0:
        print("[error] --min-prominence must be >= 0.", file=sys.stderr)
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
            bins_arg=args.bins,
            smooth_window=args.smooth,
            min_prominence=args.min_prominence,
        )
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
