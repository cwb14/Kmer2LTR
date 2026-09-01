"""K2P divergence density figure for a finished results TSV.

Reads the TSV back rather than accumulating in memory, so the figure is
independent of thread count and of whether the run was resumed.
"""
from __future__ import annotations

import csv
import math
import sys

import numpy as np

# Fraction of the distribution the x-axis covers. A handful of near-saturated
# elements otherwise compress every real burst into the first centimetre of the
# plot; anything past this is counted in the corner annotation, never dropped
# silently.
_XMAX_Q = 99.5
_FILL = "#0173B2"     # colourblind-safe blue


def read_k2p(tsv, status: str = "pass") -> np.ndarray:
    """K2P values from a results TSV, for rows with the given status."""
    vals: list[float] = []
    with open(tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") != status:
                continue
            try:
                d = float(row["k2p"])
            except (TypeError, ValueError):
                continue          # "NA", or a k2p column that is not a number
            if math.isfinite(d):
                vals.append(d)
    return np.asarray(vals, dtype=float)


def myr_per_unit_k2p(mutation_rate: float) -> float:
    """K2P distance spanned by one million years, at this substitution rate.

    The top axis divides by exactly this, so it is a relabelling of the bottom
    axis rather than a second measurement: `k2p.insertion_time` and this must
    agree, and a test pins that they do.
    """
    return 2.0 * mutation_rate * 1e6


def _kde(x: np.ndarray, grid: np.ndarray, reflect: bool = False) -> np.ndarray | None:
    """Gaussian kernel density on `grid`, Silverman bandwidth, or None.

    `reflect` mirrors the sample about zero and adds the mirrored kernels, which
    is the standard correction for a variable that cannot go negative. Without
    it a plain kernel density puts mass below zero and correspondingly
    understates the density just above it -- which for LTR divergence is the
    recent-insertion peak, the part of the figure anyone is reading it for. The
    estimate still integrates to one over the non-negative half-line.

    Evaluated one grid point at a time so memory stays O(n) rather than
    O(n * grid): at whole-genome scale the outer product is gigabytes.
    """
    n = x.size
    # One value repeated is not a distribution. Testing the bandwidth alone is
    # not enough: the standard deviation of N identical floats is ~1e-17 rather
    # than 0, which passes a `> 0` guard and then collapses the kernel far below
    # the grid spacing, so every evaluated point comes back exactly zero.
    if n < 2 or float(np.max(x) - np.min(x)) < 1e-12:
        return None
    q75, q25 = np.percentile(x, [75, 25])
    sd = float(np.std(x, ddof=1))
    scale = min(sd, (q75 - q25) / 1.349) or sd
    bw = 0.9 * scale * n ** -0.2
    if not bw > 0:
        return None
    dens = np.empty(grid.size)
    for i, g in enumerate(grid):
        total = np.exp(-0.5 * ((g - x) / bw) ** 2).sum()
        if reflect:
            total += np.exp(-0.5 * ((g + x) / bw) ** 2).sum()
        dens[i] = total
    return dens / (n * bw * math.sqrt(2.0 * math.pi))


def density_plot(tsv, out_pdf, mutation_rate: float | None = None) -> int:
    """Write the K2P density figure; return the number of elements plotted.

    Writes nothing and returns 0 when no element passed. With a mutation rate,
    a second axis across the top reads the same data in millions of years --
    an exact linear relabelling of the bottom axis, since t = d / (2*mu).
    """
    x = read_k2p(tsv)
    if x.size == 0:
        print(f"Kmer2LTR: warning: no passing element in {tsv}; no plot written",
              file=sys.stderr)
        return 0

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xmax = float(np.percentile(x, _XMAX_Q)) if x.size > 1 else float(x.max())
    xmax = max(xmax, float(np.min(x)) + 1e-6)
    beyond = int((x > xmax).sum())

    rc = {"font.size": 8, "axes.labelsize": 9, "xtick.labelsize": 8,
          "ytick.labelsize": 8, "axes.linewidth": 0.6, "pdf.fonttype": 42,
          "svg.fonttype": "none"}
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(3.6, 2.7))
        # Freedman-Diaconis, halved and clamped. FD's 2*IQR*n^(-1/3) is derived
        # for roughly symmetric data; an LTR divergence distribution is strongly
        # right-skewed, so its IQR is set by the old tail while the structure
        # worth seeing -- the recent-insertion peak -- sits in the first few
        # percent. At full FD width that peak collapses into one bar that
        # visibly disagrees with the kernel density drawn over it.
        fd = np.histogram_bin_edges(x[x <= xmax], bins="fd", range=(0.0, xmax))
        nbins = min(max(2 * (len(fd) - 1), 10), 500)
        ax.hist(x, bins=np.linspace(0.0, xmax, nbins + 1), density=True,
                color=_FILL, alpha=0.35, edgecolor="white", linewidth=0.3)
        grid = np.linspace(0.0, xmax, 512)
        dens = _kde(x, grid, reflect=True)
        if dens is not None:
            ax.plot(grid, dens, color=_FILL, linewidth=1.2)

        ax.set_xlim(0.0, xmax)
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel("K2P divergence between LTRs")
        ax.set_ylabel("Density")
        ax.spines[["top", "right"]].set_visible(False)

        if mutation_rate:
            ax.spines["top"].set_visible(True)
            scale = myr_per_unit_k2p(mutation_rate)
            sec = ax.secondary_xaxis("top", functions=(lambda d: d / scale,
                                                       lambda t: t * scale))
            sec.set_xlabel("Insertion age (Myr)")
            sec.spines["top"].set_linewidth(0.6)

        note = f"n = {x.size:,}"
        if beyond:
            note += f"\n({beyond:,} beyond axis)"
        ax.text(0.97, 0.95, note, transform=ax.transAxes, ha="right", va="top")

        fig.tight_layout()
        fig.savefig(out_pdf, dpi=300)
        plt.close(fig)
    return int(x.size)
