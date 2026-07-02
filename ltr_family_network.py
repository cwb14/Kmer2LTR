#!/usr/bin/env python3
"""Render Kmer2LTR consensus-LTR clustering as a multipage similarity-network PDF.

Nodes = consensus LTRs; edges = all-vs-all mmseqs identity; families = mmseqs
clusters; node color = insertion age. Bold minimum-spanning-tree backbone per
family = the parsimonious 'pedigree spine'.
"""
# =============================================================================
# DEVELOPER NOTE — biological goal (source-only; not printed or user-facing)
# =============================================================================
# Purpose: visualize INTRA-family (intra-cluster) relationships among LTR
# retrotransposons — within a single LTR-RT family, infer and draw which copy
# most plausibly gave rise to which: a putative *transposition path*
# (a copy-and-paste pedigree of the elements).
#
# Why it works: each element's consensus LTR approximates the template that was
# active when it inserted, so pairwise identity between consensus LTRs tracks how
# recently two copies shared a source. The minimum spanning tree on
# distance = 1 - identity is our hypothesis for that pedigree — the most
# parsimonious "who-was-copied-from-whom" backbone. The aim is to read a history
# like this at a glance:
#
#     RT1 -> RT2, RT3, RT4      # RT1 is the source lineage
#     RT2 -> RT5                # RT2 later spawned RT5
#     RT3 -> (nothing)          # a dead-end copy
#     RT4 -> RT6, RT7           # RT4 spawned two more
#
# Encodings: nodes = elements; a bold MST edge = a hypothesized parent->child
# copy event, its width/opacity = identity = confidence in that link; node
# colour = insertion age (older copies seed younger ones); faint kNN edges =
# near-relatives the tree omits (reticulation a strict tree cannot show:
# multiple active source loci, recombination, deleted/unsampled intermediates).
#
# Caveat kept in view: the MST is a hypothesis, not proof. A long, low-identity
# bridge marks a divergence gap or a missing source copy — NOT a confident
# parent->child event.
# python ltrnet/ltr_family_network.py --clusters Alyr_depth0_ltr.LTRs.alns.consensus_id0.70_cluster.tsv --results Alyr_depth0_ltr.LTRs.alns.results --fasta Alyr_depth0_ltr.LTRs.alns.consensus.fa -o Alyr_ltr_family_network.pdf -v
# =============================================================================

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import warnings
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: set before any pyplot import

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

# ---- Constants ----
MODEL_TIME_COL = {"p": 7, "JC69": 9, "K2P": 11}  # 0-based col of insertion time (years)
M8_FIDENT_COL = 3  # 0-based col of fraction identity in our easy-search --format-output
_SEARCH_FORMAT = "query,target,pident,fident,alnlen,mismatch,qcov,tcov,evalue,bits"


# ---- Parsers ----


def superfamily_of(node_id: str) -> str:
    """Return the classification label after the first '#', or 'unknown'."""
    return node_id.split("#", 1)[1] if "#" in node_id else "unknown"


def parse_fasta_lengths(path) -> dict[str, int]:
    """Map FASTA header (no leading '>') -> sequence length. Streams the file.

    Nameless headers ('>' with no id) are skipped and counted, with one warning.
    """
    lengths: dict[str, int] = {}
    name = None
    n = 0
    skipped = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                parts = line[1:].split()
                if parts:
                    name = parts[0]
                else:
                    name = None
                    skipped += 1
                n = 0
            elif line:
                n += len(line)
    if name is not None:
        lengths[name] = n
    if skipped:
        warnings.warn(f"parse_fasta_lengths: skipped {skipped} nameless header(s) in {path}",
                      stacklevel=2)
    return lengths


def parse_clusters(path) -> dict[str, str]:
    """mmseqs cluster.tsv (representative <TAB> member) -> {member: representative}.

    Blank lines are ignored; non-2-column rows are skipped and counted, with a
    single summary warning reporting the total (avoids per-row spam on big files).
    """
    member_to_rep: dict[str, str] = {}
    skipped = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) != 2:
                skipped += 1
                continue
            rep, member = parts
            member_to_rep[member] = rep
    if skipped:
        warnings.warn(f"parse_clusters: skipped {skipped} malformed row(s) in {path}",
                      stacklevel=2)
    return member_to_rep


def group_families(member_to_rep: dict[str, str]) -> dict[str, list[str]]:
    """Invert {member: rep} into {rep: [members...]} preserving first-seen order.

    The representative appears among its own members because mmseqs cluster.tsv
    includes a rep->rep self-row.
    """
    families: dict[str, list[str]] = {}
    for member, rep in member_to_rep.items():
        families.setdefault(rep, []).append(member)
    return families


def parse_ages(path, model: str = "K2P") -> dict[str, float]:
    """Map node -> insertion age (years) from the model's time column of .results.

    Unknown model raises ValueError. Blank lines are ignored; rows with too few
    columns or a non-numeric age are skipped and counted, with one summary warning.
    """
    if model not in MODEL_TIME_COL:
        raise ValueError(f"unknown age model {model!r}; choose from {list(MODEL_TIME_COL)}")
    col = MODEL_TIME_COL[model]
    ages: dict[str, float] = {}
    skipped = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip().split("\t")
            if len(f) <= col:
                skipped += 1
                continue
            try:
                ages[f[0]] = float(f[col])
            except ValueError:
                skipped += 1
                continue
    if skipped:
        warnings.warn(f"parse_ages: skipped {skipped} malformed row(s) in {path}",
                      stacklevel=2)
    return ages


def parse_m8_edges(path) -> dict[frozenset[str], float]:
    """Undirected edges from an mmseqs m8: frozenset({q,t}) -> max fraction identity.

    Drops self-hits; collapses reciprocal rows keeping the higher fident. Blank
    lines are ignored; rows with too few columns or a non-numeric fident are
    skipped and counted, with one summary warning (consistent with the other parsers).
    """
    edges: dict[frozenset[str], float] = {}
    skipped = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip().split("\t")
            if len(f) <= M8_FIDENT_COL:
                skipped += 1
                continue
            q, t = f[0], f[1]
            if q == t:
                continue
            try:
                fident = float(f[M8_FIDENT_COL])
            except ValueError:
                skipped += 1
                continue
            key = frozenset((q, t))
            if fident > edges.get(key, -1.0):
                edges[key] = fident
    if skipped:
        warnings.warn(f"parse_m8_edges: skipped {skipped} malformed row(s) in {path}",
                      stacklevel=2)
    return edges


# ---- Graph operations ----


def order_families(families: dict[str, list[str]], min_size: int = 2):
    """Families with >= min_size members, sorted by descending size (tie: rep)."""
    kept = [(rep, members) for rep, members in families.items() if len(members) >= min_size]
    kept.sort(key=lambda rm: (-len(rm[1]), rm[0]))
    return kept


def family_subgraph(members: list[str], edges: dict[frozenset, float]) -> "nx.Graph":
    """Subgraph over `members`; edges = within-member m8 pairs (distance = 1 - fident)."""
    G = nx.Graph()
    G.add_nodes_from(members)
    mset = set(members)
    for pair, fident in edges.items():
        a, b = tuple(pair)
        if a in mset and b in mset:
            G.add_edge(a, b, fident=fident, distance=1.0 - fident)
    return G


def mst_forest(G: "nx.Graph") -> tuple[set[frozenset], int]:
    """Minimum spanning forest edges (by 'distance') + connected-component count.

    networkx.minimum_spanning_tree returns a spanning forest for disconnected G.
    """
    mst = nx.minimum_spanning_tree(G, weight="distance")
    edge_set = {frozenset((u, v)) for u, v in mst.edges()}
    return edge_set, nx.number_connected_components(G)


# ---- Layout and visualization ----


def age_normalizer(ages, lo: float = 2, hi: float = 98) -> Normalize:
    """Robust color normalization clipped to the [lo, hi] age percentiles.

    Non-finite/None ages are ignored. Guarantees vmin < vmax even for empty or
    degenerate input. clip=True so outlier ages render at the colormap extremes.
    """
    arr = np.asarray([a for a in ages if a is not None and np.isfinite(a)], dtype=float)
    if arr.size == 0:
        return Normalize(vmin=0.0, vmax=1.0, clip=True)
    vmin, vmax = np.percentile(arr, [lo, hi])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin, vmax = float(arr.min()), float(arr.max())
    if vmin == vmax:
        vmax = vmin + 1.0
    return Normalize(vmin=float(vmin), vmax=float(vmax), clip=True)


def layout_positions(G, engine: str = "neato") -> dict:
    """Force-directed node positions; graphviz if available, else spring_layout.

    Falls back to networkx.spring_layout ONLY when graphviz/pygraphviz is
    unavailable in this environment (module not importable, or the layout
    executable missing). Other errors propagate so real bugs are not masked.
    A warning marks any fallback so the caller knows the requested engine was
    not used.
    """
    try:
        H = G.copy()
        for _, _, d in H.edges(data=True):
            # graphviz neato/fdp read 'len'; more identity -> shorter edge
            d["len"] = max(0.05, float(d.get("distance", 0.5)))
        return nx.nx_agraph.graphviz_layout(H, prog=engine)
    except (ImportError, FileNotFoundError) as exc:
        warnings.warn(
            f"layout_positions: graphviz engine {engine!r} unavailable "
            f"({type(exc).__name__}); falling back to spring_layout",
            stacklevel=2,
        )
        return nx.spring_layout(G, weight="fident", seed=1)


# ---- Rendering ----


def run_mmseqs_search(fasta, out_m8, *, threads=8, sensitivity=7.5, min_seq_id=0.30,
                      coverage=0.5, evalue=1e-3, max_seqs=1000, mmseqs_bin="mmseqs",
                      verbose=False) -> str:
    """Run all-vs-all `mmseqs easy-search` on `fasta`, write a BLAST-tab m8 to
    `out_m8`, PRINT the exact command (transparency), and return `out_m8`.

    Settings match those validated for LTR consensus-LTR networks: nucleotide
    search (--search-type 3), high sensitivity (-s 7.5), no low-complexity mask,
    50% mutual coverage, and a permissive identity floor (0.30) so the m8 is a
    superset that the drawing step (edge_id_floor) thresholds further.
    """
    if shutil.which(mmseqs_bin) is None:
        raise FileNotFoundError(
            f"{mmseqs_bin!r} not found in PATH; install mmseqs2 "
            f"(e.g. 'mamba install -c bioconda mmseqs2') or pass a precomputed --m8")
    out_m8 = str(out_m8)
    tmp_dir = tempfile.mkdtemp(prefix="ltrnet_mmseqs_", dir=str(Path(out_m8).parent))
    cmd = [mmseqs_bin, "easy-search", str(fasta), str(fasta), out_m8, tmp_dir,
           "--min-seq-id", f"{min_seq_id}", "-c", f"{coverage}", "--cov-mode", "0",
           "--mask", "0", "-s", f"{sensitivity}", "--max-seqs", str(max_seqs),
           "--search-type", "3", "-e", f"{evalue}", "--threads", str(threads),
           "--format-output", _SEARCH_FORMAT]
    print("[MMSEQS] " + " ".join(cmd))
    try:
        subprocess.run(cmd, check=True,
                       stdout=(None if verbose else subprocess.DEVNULL),
                       stderr=(None if verbose else subprocess.DEVNULL))
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    return out_m8


def build_model(fasta, clusters, results, m8, model: str = "K2P") -> dict:
    """Parse all four inputs into a single graph model dict."""
    member_to_rep = parse_clusters(clusters)
    return {
        "lengths": parse_fasta_lengths(fasta),
        "member_to_rep": member_to_rep,
        "families": group_families(member_to_rep),
        "ages": parse_ages(results, model=model),
        "edges": parse_m8_edges(m8),
    }


def _node_sizes(members, lengths, smin=40.0, smax=400.0, lo=None, hi=None):
    """Map consensus-LTR length -> marker size on a fixed [smin, smax] scale.

    Pass GLOBAL length bounds via lo/hi (min/max across ALL elements) so a given
    length renders at the same size on every page — consistent with the
    document-wide color scale. Falls back to the local min/max when lo/hi are
    not given. Sizes are clamped into [smin, smax].
    """
    vals = [lengths.get(m, 0) for m in members]
    if lo is None or hi is None:
        lo, hi = (min(vals), max(vals)) if vals else (0, 1)
    if hi == lo:
        return {m: (smin + smax) / 2 for m in members}

    def _s(m):
        x = smin + (smax - smin) * (lengths.get(m, 0) - lo) / (hi - lo)
        return max(smin, min(smax, x))

    return {m: _s(m) for m in members}


def _resolve_cmap(cmap):
    """Return a Colormap with an opaque neutral 'bad' color so missing (NaN) ages
    render as visible gray instead of matplotlib's default transparent."""
    base = matplotlib.colormaps[cmap] if isinstance(cmap, str) else cmap
    return base.with_extremes(bad="0.6")


def _select_faint_edges(G, mst_set, edge_id_floor, knn):
    """Faint reticulation = each node's top-`knn` most-similar NON-MST neighbours
    (edges also passing edge_id_floor). An edge is kept if it is in EITHER
    endpoint's top-knn (union). knn <= 0 disables the faint layer; knn is None
    keeps every non-MST edge >= floor (the old behaviour)."""
    if knn is not None and knn <= 0:
        return set()
    keep = set()
    for node in G.nodes():
        cand = []
        for nbr in G.neighbors(node):
            key = frozenset((node, nbr))
            if key in mst_set:
                continue
            fid = G[node][nbr]["fident"]
            if fid >= edge_id_floor:
                cand.append((fid, key))
        cand.sort(key=lambda t: t[0], reverse=True)
        lim = len(cand) if knn is None else knn
        for _, key in cand[:lim]:
            keep.add(key)
    return keep


_MST_ID_LO, _MST_ID_HI = 0.70, 1.00   # identity range mapped to MST edge weight/opacity


def _mst_edge_style(fident, gray, lw_lo, lw_hi, alpha_lo):
    """Map an MST edge's identity to (linewidth, rgba). Higher-identity (more
    confident) links render thicker and more opaque; weak bridges render thin
    and pale. Identity is clipped to [_MST_ID_LO, _MST_ID_HI] before mapping."""
    t = (fident - _MST_ID_LO) / (_MST_ID_HI - _MST_ID_LO)
    t = 0.0 if t < 0 else 1.0 if t > 1 else t
    lw = lw_lo + t * (lw_hi - lw_lo)
    alpha = alpha_lo + t * (1.0 - alpha_lo)
    return lw, (gray, gray, gray, alpha)


def _draw_edges(ax, G, pos, mst_set, edge_id_floor, knn, *,
                mst_gray, mst_lw_range, mst_alpha_min, faint_style):
    """Draw faint reticulation (each node's top-`knn` non-MST neighbours) plus the
    MST spine, where each MST edge's width and opacity encode its % identity
    (confidence). Batched into LineCollections."""
    faint = _select_faint_edges(G, mst_set, edge_id_floor, knn)
    lw_lo, lw_hi = mst_lw_range
    faint_segs = []
    mst_segs, mst_lws, mst_cols = [], [], []
    for u, v in G.edges():
        key = frozenset((u, v))
        seg = (pos[u], pos[v])
        if key in mst_set:
            lw, col = _mst_edge_style(G[u][v]["fident"], mst_gray, lw_lo, lw_hi, mst_alpha_min)
            mst_segs.append(seg); mst_lws.append(lw); mst_cols.append(col)
        elif key in faint:
            faint_segs.append(seg)
    if faint_segs:
        ax.add_collection(LineCollection(faint_segs, **faint_style))
    if mst_segs:
        ax.add_collection(LineCollection(mst_segs, colors=mst_cols,
                                         linewidths=mst_lws, zorder=1))


def _draw_family(ax, members, edges, ages, lengths, norm, cmap, engine,
                 edge_id_floor, knn, len_lo, len_hi):
    G = family_subgraph(members, edges)
    mst, ncomp = mst_forest(G)
    pos = layout_positions(G, engine=engine)
    sizes = _node_sizes(members, lengths, lo=len_lo, hi=len_hi)
    _draw_edges(ax, G, pos, mst, edge_id_floor, knn,
                mst_gray=0.15, mst_lw_range=(0.8, 2.4), mst_alpha_min=0.45,
                faint_style={"colors": "0.75", "linewidths": 0.5, "alpha": 0.6, "zorder": 0})
    colors = [ages.get(m, np.nan) for m in members]
    ax.scatter([pos[m][0] for m in members], [pos[m][1] for m in members],
               s=[sizes[m] for m in members], c=colors, cmap=cmap, norm=norm,
               edgecolors="white", linewidths=0.4, zorder=2)
    dom = Counter(superfamily_of(m) for m in members).most_common(1)[0][0]
    tag = f"{members[0]}  ·  n={len(members)}  ·  {dom}"
    if ncomp > 1:
        tag += f"  ·  {ncomp} components"
    ax.set_title("")
    ax.text(0.01, 0.99, tag, transform=ax.transAxes, va="top", ha="left",
            fontsize=7, family="monospace")
    ax.set_axis_off()


def _draw_overview(ax, kept, model, norm, cmap, edge_id_floor, knn, engine, len_lo, len_hi):
    nodes = [m for _, members in kept for m in members]
    G = nx.Graph()
    G.add_nodes_from(nodes)
    per_family_mst = set()
    for _, members in kept:
        FG = family_subgraph(members, model["edges"])
        mst, _ = mst_forest(FG)
        per_family_mst |= mst
        for a, b, d in FG.edges(data=True):
            G.add_edge(a, b, fident=d["fident"], distance=d["distance"])
    pos = layout_positions(G, engine=engine)
    _draw_edges(ax, G, pos, per_family_mst, edge_id_floor, knn,
                mst_gray=0.35, mst_lw_range=(0.5, 1.2), mst_alpha_min=0.40,
                faint_style={"colors": "0.85", "linewidths": 0.3, "alpha": 0.5, "zorder": 0})
    sizes = _node_sizes(nodes, model["lengths"], smin=10.0, smax=90.0, lo=len_lo, hi=len_hi)
    colors = [model["ages"].get(m, np.nan) for m in nodes]
    ax.scatter([pos[m][0] for m in nodes], [pos[m][1] for m in nodes],
               s=[sizes[m] for m in nodes], c=colors, cmap=cmap, norm=norm,
               edgecolors="white", linewidths=0.3, zorder=2)
    ax.text(0.01, 0.99, f"overview  ·  {len(kept)} families  ·  {len(nodes)} elements",
            transform=ax.transAxes, va="top", ha="left", fontsize=7, family="monospace")
    ax.set_axis_off()


def _add_colorbar(fig, ax, norm, cmap, age_units, age_model):
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    div = 1e6 if age_units == "Myr" else 1.0
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.02,
                      format=lambda x, _pos: f"{x/div:.2f}")
    cb.set_label(f"Insertion age ({age_units}, {age_model})", fontsize=8)
    cb.ax.tick_params(labelsize=7)


def render_pdf(model, out_pdf, *, min_family_size=2, edge_id_floor=0.70,
               knn=4, cmap="viridis", overview_engine="sfdp", family_engine="neato",
               age_units="Myr", age_model="K2P", verbose=False) -> int:
    """Render overview + per-family pages to a multipage PDF. Returns page count."""
    norm = age_normalizer(list(model["ages"].values()))
    cmap_obj = _resolve_cmap(cmap)
    lvals = list(model["lengths"].values())
    len_lo, len_hi = (min(lvals), max(lvals)) if lvals else (0, 1)
    kept = order_families(model["families"], min_size=min_family_size)
    pages = 0
    with PdfPages(out_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        _draw_overview(ax, kept, model, norm, cmap_obj, edge_id_floor, knn,
                       overview_engine, len_lo, len_hi)
        _add_colorbar(fig, ax, norm, cmap_obj, age_units, age_model)
        pdf.savefig(fig, dpi=300); plt.close(fig); pages += 1
        for rep, members in kept:
            if verbose:
                print(f"[RENDER] family {rep} (n={len(members)})")
            fig, ax = plt.subplots(figsize=(7.5, 7.5))
            _draw_family(ax, members, model["edges"], model["ages"], model["lengths"],
                         norm, cmap_obj, family_engine, edge_id_floor, knn, len_lo, len_hi)
            _add_colorbar(fig, ax, norm, cmap_obj, age_units, age_model)
            pdf.savefig(fig, dpi=300); plt.close(fig); pages += 1
    return pages


# ---- CLI ----


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Render Kmer2LTR consensus-LTR clustering as a multipage "
                    "similarity-network PDF (overview + one page per family, "
                    "MST 'pedigree spine' bold, nodes colored by insertion age).")
    p.add_argument("--m8", help="precomputed all-vs-all mmseqs m8 (edges); if omitted, "
                                "generated from --fasta via mmseqs easy-search")
    p.add_argument("--clusters", required=True, help="mmseqs *_cluster.tsv (families)")
    p.add_argument("--results", required=True, help="Kmer2LTR .LTRs.alns.results (insertion age)")
    p.add_argument("--fasta", required=True, help="consensus LTR FASTA (labels + length)")
    p.add_argument("-o", "--out", required=True, help="output multipage PDF")
    p.add_argument("--m8-out", default=None,
                   help="where to write the generated m8 when --m8 is omitted "
                        "(default: <fasta-stem>.allvall.m8); ignored if --m8 is given")
    p.add_argument("--threads", type=int, default=8,
                   help="threads for the internal mmseqs easy-search (default: 8)")
    p.add_argument("--mmseqs", default="mmseqs",
                   help="mmseqs binary to use for the internal search (default: mmseqs)")
    p.add_argument("--min-family-size", type=int, default=2,
                   help="skip families smaller than this (default: 2, i.e. skip singletons)")
    p.add_argument("--edge-id-floor", type=float, default=0.70,
                   help="draw non-MST edges only at/above this fraction identity (default: 0.70)")
    p.add_argument("--knn", type=int, default=4,
                   help="faint reticulation = each node's K most-similar non-MST "
                        "neighbours (default: 4; 0 = MST spine only; large K = all "
                        "edges >= --edge-id-floor, the old dense behaviour)")
    p.add_argument("--age-model", choices=["K2P", "JC69", "p"], default="K2P",
                   help="divergence model whose insertion-time column to use (default: K2P)")
    p.add_argument("--age-units", choices=["Myr", "years"], default="Myr",
                   help="colorbar units (default: Myr)")
    p.add_argument("--layout", default="neato", choices=["neato", "fdp", "sfdp"],
                   help="per-family graphviz layout engine (default: neato)")
    p.add_argument("--cmap", default="viridis", help="matplotlib colormap (default: viridis)")
    p.add_argument("-v", "--verbose", action="store_true", help="per-family progress")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    # --clusters/--results/--fasta are always required; --m8 is optional (generated if absent)
    for label, path in [("--clusters", args.clusters), ("--results", args.results),
                        ("--fasta", args.fasta)]:
        if not Path(path).is_file():
            print(f"ERROR: {label} file not found: {path}", file=sys.stderr)
            return 2
    if args.m8 is not None and not Path(args.m8).is_file():
        print(f"ERROR: --m8 file not found: {args.m8}", file=sys.stderr)
        return 2
    out_dir = Path(args.out).parent
    if not out_dir.exists():
        print(f"ERROR: output directory does not exist: {out_dir}", file=sys.stderr)
        return 2

    m8_path = args.m8
    if m8_path is None:
        if shutil.which(args.mmseqs) is None:
            print(f"ERROR: {args.mmseqs!r} not found in PATH; install mmseqs2 "
                  f"(e.g. 'mamba install -c bioconda mmseqs2') or pass a precomputed --m8",
                  file=sys.stderr)
            return 2
        m8_path = args.m8_out or f"{Path(args.fasta).with_suffix('')}.allvall.m8"
        run_mmseqs_search(args.fasta, m8_path, threads=args.threads,
                          mmseqs_bin=args.mmseqs, verbose=args.verbose)

    if args.verbose:
        print(f"[LOAD] {args.fasta} / {args.clusters} / {args.results} / {m8_path}")
    model = build_model(args.fasta, args.clusters, args.results, m8_path,
                        model=args.age_model)
    pages = render_pdf(model, args.out, min_family_size=args.min_family_size,
                       edge_id_floor=args.edge_id_floor, knn=args.knn, cmap=args.cmap,
                       family_engine=args.layout, age_units=args.age_units,
                       age_model=args.age_model, verbose=args.verbose)
    print(f"[DONE] wrote {pages} pages -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
