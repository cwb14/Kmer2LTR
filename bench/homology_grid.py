"""Homology-only benchmark: perfect elements, known mutations, known flanks,
known indels -- and an EXACT true alignment.

Why this exists alongside `gold_robustness.py`. That benchmark draws its
elements with `select_gold`, which filters on the tool's own output *and*
requires canonical `TG`..`CA` termini. Both are deliberate there, and both are
disqualifying for the question "how well does the tool recover a boundary from
homology alone?" -- the first is circular, the second smuggles in a structural
prior the tool itself declines to use.

This module's primary source is instead `bench/out/truth.fa`: library
`X-LTR` + `X-I` + `X-LTR` concatenations built by `bench/build_truth.py`. Those
elements are **perfect by construction** -- the two LTR copies are literally the
same string, so their divergence is exactly zero and their boundaries are known
exactly -- and selecting them involves no motif, no TSD and no call by ltrk2p.
Everything the benchmark then measures is a deviation the harness itself
introduced.

Three perturbation axes, applied to those perfect elements:

  * substitutions, parameterised by the TARGET PAIRWISE p-distance
    (5%, 10%, ... 35% of sites differing between the two copies), not by a
    branch length -- so "add 25% mutations" means what it says. The K2P
    pairwise distance that delivers it is solved for numerically and recorded
    as `d_nominal`.
  * flanks of 5, 25, 45 and 65 bp, drawn from a composition-matched shuffle of
    the element's own sequence, so a detector cannot succeed by noticing a
    compositional step.
  * indels, at a rate swept independently of the substitution rate.

**The true alignment is tracked through evolution**, so `realized_p` and
`realized_k2p` are exact rather than approximated. `bench/simulate.py` computes
its `realized_k2p` by truncate-and-compare on unaligned strings, which Task 14
found disagrees materially with a proper alignment even at `indel_frac=0`; with
indels switched on it is simply wrong. Nothing here inherits that.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
import time
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ltrk2p.fasta import read_fasta            # noqa: E402
from ltrk2p.k2p import count_substitutions, k2p_distance   # noqa: E402
from ltrk2p.scoring import k2p_probs           # noqa: E402

_TI = {"A": "G", "G": "A", "C": "T", "T": "C"}
_TV = {"A": "CT", "G": "CT", "C": "AG", "T": "AG"}
_BASES = "ACGT"

# The grid. Substitution levels are TARGET PAIRWISE p-distances.
P_GRID = (0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35)
FLANK_GRID = (0, 5, 25, 45, 65)
INDEL_GRID = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05)
INDEL_P_GRID = (0.10, 0.25)      # substitution levels the indel panel is run at
KAPPA = 2.0
MEAN_INDEL_LEN = 3.0

TRUTH_COLUMNS = [
    "seq_id", "source", "orig_id", "panel", "p_target", "d_nominal", "indel_rate",
    "flank5", "flank3", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end",
    "seq_len", "realized_p", "realized_k2p", "true_n_sites", "true_n_ts",
    "true_n_tv", "true_n_gapcols",
]


# --------------------------------------------------------------------------- #
# Divergence parameterisation
# --------------------------------------------------------------------------- #

def d_for_p(p_target: float, kappa: float = KAPPA) -> float:
    """K2P pairwise distance whose expected p-distance is `p_target`.

    Two branches of d/2 compose to one branch of d (the model is Markov and
    reversible), so the probability the two copies differ at a site is
    `1 - p_same(d, kappa)`. That is strictly decreasing in d, so bisection is
    exact to machine precision and needs no derivative.
    """
    if p_target <= 0.0:
        return 0.0
    if p_target >= 0.75:
        raise ValueError(f"p_target {p_target} is at or beyond K2P saturation (0.75)")
    lo, hi = 0.0, 1.0
    while 1.0 - k2p_probs(hi, kappa)[0] < p_target:
        hi *= 2.0
        if hi > 1e6:
            raise ValueError(f"no K2P distance reaches p_target {p_target}")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 1.0 - k2p_probs(mid, kappa)[0] < p_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------------- #
# Evolution with a tracked true alignment
# --------------------------------------------------------------------------- #

def evolve_tracked(seq: str, branch_d: float, kappa: float, indel_rate: float, rng):
    """Evolve `seq` along one branch and return `(descendant, trace)`.

    `trace` is one entry per emitted base: the ancestral index it descends from,
    or `None` for an inserted base. Deleted ancestral positions simply do not
    appear. Two traces from the same ancestor compose into the exact true
    pairwise alignment (`true_alignment` below) -- which is what makes
    `realized_k2p` here a measurement rather than an approximation.

    Substitutions are drawn from the EXACT K2P transition matrix P(branch_d,
    kappa), not from a single "did it mutate" Bernoulli: a single-draw model
    cannot produce multiple hits at one site, so realised divergence would run
    high relative to nominal exactly in the high-divergence regime this
    benchmark exists to probe.
    """
    p_same, p_ti, p_tv = k2p_probs(branch_d, kappa)
    out: list[str] = []
    trace: list[int | None] = []
    i = 0
    n = len(seq)
    while i < n:
        ch = seq[i]
        if indel_rate and rng.random() < indel_rate:
            length = 1 + int(rng.expovariate(1.0 / MEAN_INDEL_LEN))
            length = min(length, 30)
            if rng.random() < 0.5:                     # deletion of `length` sites
                i += length
                continue
            for _ in range(length):                    # insertion before this site
                out.append(rng.choice(_BASES))
                trace.append(None)
        if ch not in _BASES:
            out.append(ch)
            trace.append(i)
            i += 1
            continue
        u = rng.random()
        if u < p_same:
            out.append(ch)
        elif u < p_same + p_ti:
            out.append(_TI[ch])
        elif u < p_same + p_ti + p_tv:
            out.append(_TV[ch][0])
        else:
            out.append(_TV[ch][1])
        trace.append(i)
        i += 1
    return "".join(out), trace


def true_alignment(a: str, ta, b: str, tb) -> tuple[str, str]:
    """Exact pairwise alignment of two descendants of one ancestor.

    Ancestral positions present in both are aligned to each other; a position
    lost on one branch, and any inserted base, becomes a gap column. Both traces
    are non-decreasing over their non-None entries, so one merge pass suffices.
    """
    qa: list[str] = []
    ra: list[str] = []
    i = j = 0
    while i < len(a) or j < len(b):
        ai = ta[i] if i < len(a) else None
        bj = tb[j] if j < len(b) else None
        if i < len(a) and ai is None:
            qa.append(a[i]); ra.append("-"); i += 1
        elif j < len(b) and bj is None:
            qa.append("-"); ra.append(b[j]); j += 1
        elif i >= len(a):
            qa.append("-"); ra.append(b[j]); j += 1
        elif j >= len(b):
            qa.append(a[i]); ra.append("-"); i += 1
        elif ai == bj:
            qa.append(a[i]); ra.append(b[j]); i += 1; j += 1
        elif ai < bj:
            qa.append(a[i]); ra.append("-"); i += 1
        else:
            qa.append("-"); ra.append(b[j]); j += 1
    return "".join(qa), "".join(ra)


def shuffled_flank(pool: str, n: int, rng) -> str:
    """`n` bases matching `pool`'s composition exactly but none of its order."""
    if n <= 0:
        return ""
    reps = -(-n // max(1, len(pool)))
    chars = list(pool * reps)
    rng.shuffle(chars)
    return "".join(chars[:n])


def make_element(ltr: str, internal: str, p_target: float, flank: int,
                 indel_rate: float, rng, kappa: float = KAPPA):
    """One perturbed element plus its exact truth."""
    d = d_for_p(p_target, kappa)
    l5, t5 = evolve_tracked(ltr, d / 2.0, kappa, indel_rate, rng)
    l3, t3 = evolve_tracked(ltr, d / 2.0, kappa, indel_rate, rng)
    pool = (ltr + internal) or "ACGT"
    f5 = shuffled_flank(pool, flank, rng)
    f3 = shuffled_flank(pool, flank, rng)
    seq = f5 + l5 + internal + l3 + f3
    a, b = true_alignment(l5, t5, l3, t3)
    c = count_substitutions(a, b)
    realized_k2p, _ = k2p_distance(c)
    realized_p = (c.n_ts + c.n_tv) / c.n_sites if c.n_sites else None
    truth = {
        "p_target": p_target, "d_nominal": d, "indel_rate": indel_rate,
        "flank5": flank, "flank3": flank,
        "ltr5_start": len(f5) + 1, "ltr5_end": len(f5) + len(l5),
        "ltr3_start": len(f5) + len(l5) + len(internal) + 1,
        "ltr3_end": len(seq) - len(f3), "seq_len": len(seq),
        "realized_p": realized_p, "realized_k2p": realized_k2p,
        "true_n_sites": c.n_sites, "true_n_ts": c.n_ts, "true_n_tv": c.n_tv,
        "true_n_gapcols": c.n_gapcols,
    }
    return seq, truth


# --------------------------------------------------------------------------- #
# Source elements
# --------------------------------------------------------------------------- #

def load_library_elements(truth_fa, truth_tsv, n: int, seed: int,
                          min_len: int = 500, max_len: int = 12000,
                          min_ltr: int = 100):
    """Perfect elements straight out of `build_truth.py`'s construction.

    Joined by RECORD ORDER, the convention every other module here uses: the
    two files are written together and `elem_id` is not guaranteed unique.
    Selection is by length only -- no motif, no TSD, no call by ltrk2p.
    """
    pool = []
    with open(truth_tsv, newline="") as tf:
        for row, (sid, seq) in zip(csv.DictReader(tf, delimiter="\t"), read_fasta(truth_fa)):
            ltr_len = int(row["ltr_len"])
            l5s, l5e = int(row["ltr5_start"]), int(row["ltr5_end"])
            l3s = int(row["ltr3_start"])
            if l5e - l5s + 1 != ltr_len or ltr_len < min_ltr:
                continue
            if not (min_len <= len(seq) <= max_len) or len(seq) < l3s:
                continue
            pool.append((sid, seq[l5s - 1:l5e], seq[l5e:l3s - 1]))
    rng = random.Random(seed)
    if len(pool) > n:
        pool = [pool[i] for i in sorted(rng.sample(range(len(pool)), n))]
    return pool


def load_real_elements(pred_tsv, fasta, n: int, seed: int, max_len: int = 12000,
                       min_bitscore: float = 200.0, max_k2p: float = 0.05):
    """Real genomic elements, motif-free.

    Same filters as `gold_subset.select_gold` MINUS the `TG`..`CA` requirement,
    and the LTR length is read straight out of `pred_tsv` rather than by
    re-classifying: this panel is then defined entirely by an artifact that
    predates any change under test, so re-running it cannot move the goalposts.

    It is still selected using a prediction, so it is not the non-circular
    panel -- the library panel is. Its value is real composition and real
    internal structure, and it is reported separately for exactly that reason.

    Joined by RECORD ORDER, the convention `select_gold` documents.
    """
    pool = []
    with open(pred_tsv, newline="") as fh:
        for row, (sid, seq) in zip(csv.DictReader(fh, delimiter="\t"), read_fasta(fasta)):
            if row["status"] != "pass":
                continue
            if int(row["ltr5_start"]) != 1 or int(row["ltr3_end"]) != int(row["seq_len"]):
                continue
            k = int(row["ltr5_len"])
            if k != int(row["ltr3_len"]) or float(row["k2p"]) > max_k2p:
                continue
            if float(row["bitscore"]) < min_bitscore or len(seq) > max_len:
                continue
            if 2 * k >= len(seq) or k < 100:
                continue
            pool.append((sid, seq[:k], seq[k:len(seq) - k]))
    rng = random.Random(seed)
    if len(pool) > n:
        pool = [pool[i] for i in sorted(rng.sample(range(len(pool)), n))]
    return pool


# --------------------------------------------------------------------------- #
# Grid construction
# --------------------------------------------------------------------------- #

def cells():
    """(panel, p_target, indel_rate, flank) for the whole grid."""
    for p in P_GRID:
        for f in FLANK_GRID:
            yield "subs", p, 0.0, f
    for ir in INDEL_GRID:
        for p in INDEL_P_GRID:
            for f in FLANK_GRID:
                yield "indel", p, ir, f


def build_grid(sources: dict[str, list], out_fa, out_truth, seed: int,
               verbose: bool = False) -> int:
    """Write the perturbed FASTA and its exact truth TSV.

    One `random.Random` per source, bound outside every loop and threaded
    through the whole grid -- never re-seeded per element or per cell, which
    would correlate cells that are supposed to be independent.
    """
    grid = list(cells())
    n = 0
    with open(out_fa, "w") as fa, open(out_truth, "w") as tsv:
        tsv.write("\t".join(TRUTH_COLUMNS) + "\n")
        for si, (source, elements) in enumerate(sources.items()):
            rng = random.Random(seed + si * 1_000_003)
            t0 = time.time()
            for ei, (orig_id, ltr, internal) in enumerate(elements):
                for panel, p, ir, flank in grid:
                    sid = f"{source}__{orig_id}__{panel}__p{p}__i{ir}__f{flank}"
                    seq, truth = make_element(ltr, internal, p, flank, ir, rng)
                    fa.write(f">{sid}\n{seq}\n")
                    tsv.write("\t".join(_fmt(v) for v in [
                        sid, source, orig_id, panel, truth["p_target"],
                        truth["d_nominal"], truth["indel_rate"], truth["flank5"],
                        truth["flank3"], truth["ltr5_start"], truth["ltr5_end"],
                        truth["ltr3_start"], truth["ltr3_end"], truth["seq_len"],
                        truth["realized_p"], truth["realized_k2p"],
                        truth["true_n_sites"], truth["true_n_ts"],
                        truth["true_n_tv"], truth["true_n_gapcols"],
                    ]) + "\n")
                    n += 1
                if verbose and (ei + 1) % 100 == 0:
                    print(f"  {source}: {ei + 1}/{len(elements)} elements "
                          f"({n} records, {time.time() - t0:.0f}s)", file=sys.stderr)
            if verbose:
                print(f"  {source}: {len(elements)} elements x {len(grid)} cells "
                      f"= {len(elements) * len(grid)} records in "
                      f"{time.time() - t0:.0f}s", file=sys.stderr)
    return n


def _fmt(v) -> str:
    """`+ 0.0` normalises negative zero, matching ltrk2p.runner.format_row."""
    return "NA" if v is None else (f"{v + 0.0:.6g}" if isinstance(v, float) else str(v))


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #

_COORDS = ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end")


# Bin edges are the midpoints of the gold benchmark's D_GRID, so a d_hat bin
# label here is directly comparable to one from `run_bench.dhat_bin`.
_DHAT_LABELS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
_DHAT_EDGES = tuple((_DHAT_LABELS[i] + _DHAT_LABELS[i + 1]) / 2.0
                    for i in range(len(_DHAT_LABELS) - 1))


def dhat_bin(d_hat: float) -> float:
    for edge, label in zip(_DHAT_EDGES, _DHAT_LABELS[:-1]):
        if d_hat < edge:
            return label
    return _DHAT_LABELS[-1]


def score_homology(truth_tsv, pred_tsv, dhat_tsv=None) -> dict:
    """Score a prediction TSV against the exact homology truth.

    Cells are keyed `(source, panel, p_target, indel_rate, flank)`. Every
    boundary statistic is UNCONDITIONAL over the cell's records -- a record the
    tool lost contributes to `n_lost`, never silently to a smaller denominator,
    so a configuration cannot look accurate by answering less often.

    Per cell: `n`, `n_pass`, `n_weak`, `n_lost`; `abs_<coord>` and `exact_<coord>`
    over records with a numeric coordinate; `n_false_flank` (flank==0 only) and
    `n_detected` (flank>0 only, both sides called); `called_abs5/3`, the
    unconditional absolute error of called flank length; and K2P error sums
    against both the exact realized K2P and the nominal target.
    """
    acc: dict = defaultdict(lambda: defaultdict(float))
    by_dhat: dict = defaultdict(lambda: defaultdict(float)) if dhat_tsv else None
    with ExitStack() as stack:
        tf = stack.enter_context(open(truth_tsv, newline=""))
        pf = stack.enter_context(open(pred_tsv, newline=""))
        if dhat_tsv:
            hf = stack.enter_context(open(dhat_tsv, newline=""))
            rows = zip(csv.DictReader(tf, delimiter="\t"),
                       csv.DictReader(pf, delimiter="\t"),
                       csv.DictReader(hf, delimiter="\t"))
        else:
            rows = ((t, p, None) for t, p in zip(csv.DictReader(tf, delimiter="\t"),
                                                 csv.DictReader(pf, delimiter="\t")))
        for t, p, h in rows:
            key = (t["source"], t["panel"], float(t["p_target"]),
                   float(t["indel_rate"]), int(t["flank5"]))
            cells_to_fill = [acc[key]]
            if by_dhat is not None and h is not None and h.get("d_hat", "NA") not in ("NA", ""):
                cells_to_fill.append(by_dhat[(t["source"], t["panel"],
                                              dhat_bin(float(h["d_hat"])),
                                              float(t["indel_rate"]), int(t["flank5"]))])
            c = _Fan(cells_to_fill)
            c["n"] += 1
            status = p["status"]
            located = p["ltr5_start"] not in ("NA", "")
            if status == "pass":
                c["n_pass"] += 1
            elif status == "weak_pair":
                c["n_weak"] += 1
            if not located:
                c["n_lost"] += 1
                # a lost record still counts against flank accuracy: it called
                # no flank at all, so its called length is 0
                if int(t["flank5"]) > 0:
                    c["called_abs5"] += int(t["flank5"])
                    c["called_abs3"] += int(t["flank3"])
                continue
            for f in _COORDS:
                e = abs(int(p[f]) - int(t[f]))
                c[f"abs_{f}"] += e
                c[f"exact_{f}"] += (e == 0)
            c["n_located"] += 1
            f5c, f3c = int(p["flank5_len"]), int(p["flank3_len"])
            ftrue = int(t["flank5"])
            if ftrue == 0:
                if f5c > 0 or f3c > 0:
                    c["n_false_flank"] += 1
            else:
                if f5c > 0 and f3c > 0:
                    c["n_detected"] += 1
                c["called_abs5"] += abs(f5c - ftrue)
                c["called_abs3"] += abs(f3c - int(t["flank3"]))
            if p["k2p"] not in ("NA", ""):
                pk = float(p["k2p"])
                if t["realized_k2p"] not in ("NA", ""):
                    err = pk - float(t["realized_k2p"])
                    c["n_k2p"] += 1
                    c["k2p_err"] += err
                    c["k2p_err2"] += err * err
                errn = pk - float(t["d_nominal"])
                c["n_k2p_nom"] += 1
                c["k2p_err_nom"] += errn
                c["k2p_err2_nom"] += errn * errn
    out = {"by_p": {k: dict(v) for k, v in acc.items()}}
    out["by_dhat"] = ({k: dict(v) for k, v in by_dhat.items()}
                      if by_dhat is not None else None)
    return out


class _Fan:
    """Write-through to several accumulator dicts at once.

    The per-record increments are identical whether a cell is keyed by the known
    perturbation or by the tool's own d_hat; only the key differs. Fanning the
    writes keeps one copy of the accounting logic, so the two views cannot drift.
    """

    __slots__ = ("_targets",)

    def __init__(self, targets):
        self._targets = targets

    def __getitem__(self, k):
        return self._targets[0][k]

    def __setitem__(self, k, v):
        delta = v - self._targets[0][k]
        for t in self._targets:
            t[k] += delta


def cells_to_json(scored: dict, path) -> None:
    def _ser(d):
        return None if d is None else {"|".join(str(x) for x in k): v for k, v in d.items()}
    Path(path).write_text(json.dumps({"by_p": _ser(scored["by_p"]),
                                      "by_dhat": _ser(scored.get("by_dhat"))}, indent=1))


def cells_from_json(path) -> dict:
    raw = json.loads(Path(path).read_text())

    def _de(d):
        if d is None:
            return None
        out = {}
        for k, v in d.items():
            src, panel, p, ir, f = k.split("|")
            out[(src, panel, float(p), float(ir), int(f))] = v
        return out
    return {"by_p": _de(raw["by_p"]), "by_dhat": _de(raw.get("by_dhat"))}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="write the perturbed FASTA and exact truth TSV")
    b.add_argument("--truth-fa", default="bench/out/truth.fa")
    b.add_argument("--truth-tsv", default="bench/out/truth.tsv")
    b.add_argument("--n-lib", type=int, default=500)
    b.add_argument("--real-pred", default=None,
                   help="prediction TSV for the motif-free real panel (optional)")
    b.add_argument("--real-fasta", default=None)
    b.add_argument("--n-real", type=int, default=200)
    b.add_argument("--out-fa", default="bench/out/homology_grid.fa")
    b.add_argument("--out-truth", default="bench/out/homology_truth.tsv")
    b.add_argument("--seed", type=int, default=0)
    b.add_argument("-v", "--verbose", action="store_true")

    s = sub.add_parser("score", help="score a prediction TSV against the truth")
    s.add_argument("--truth", default="bench/out/homology_truth.tsv")
    s.add_argument("--pred", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--dhat", default=None,
                   help="per-record d_hat TSV; adds a view keyed by the tool's "
                        "own estimated divergence instead of the known one")

    args = ap.parse_args(argv)
    if args.cmd == "build":
        sources = {}
        lib = load_library_elements(args.truth_fa, args.truth_tsv, args.n_lib, args.seed)
        print(f"homology_grid: library panel = {len(lib)} perfect elements", file=sys.stderr)
        sources["lib"] = lib
        if args.real_pred and args.real_fasta:
            real = load_real_elements(args.real_pred, args.real_fasta,
                                      args.n_real, args.seed)
            print(f"homology_grid: real panel = {len(real)} elements (motif-free)",
                  file=sys.stderr)
            sources["real"] = real
        n = build_grid(sources, args.out_fa, args.out_truth, args.seed,
                       verbose=args.verbose)
        print(f"homology_grid: wrote {n} records -> {args.out_fa}", file=sys.stderr)
        return 0

    scored = score_homology(args.truth, args.pred, args.dhat)
    cells_to_json(scored, args.out)
    print(f"homology_grid: scored {len(scored['by_p'])} cells -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
