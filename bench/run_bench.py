"""Benchmark harness: synthetic dataset generation, generic scoring, and the
spec section 6.4 ablation grid.

**Context.** The gold-subset
benchmark (`bench/gold_subset.py`, `bench/gold_robustness.py`) found that
`T_BITS = 5.0`, calibrated on SYNTHETIC sequence, badly
under-calls overextension on REAL element sequence (false-flank rate at
d=0.30 measured 26.8% on real gold data vs a synthetic-only estimate of
1.2%, a ~20x gap). This module's single most important job is therefore a
T_BITS sweep measured on the REAL gold-perturbed grid
(`bench/out/gold_perturbed.fa` / `bench/out/gold_truth.tsv`), not a rebuild
of that machinery -- see `score_gold_grid`, `sweep_t_bits`, and
`analyze_dhat_threshold` below. The original interfaces
(`make_dataset`, `score_run`, `ablate`) are implemented in full and are what
drive both the literal spec 6.4 grid on synthetic data AND, reusing the very
same `ablate()`, the gold-data ablation grid.

Everything here is orchestration: every ablation "mechanism" listed
(fixed matrix, no_stage3, no_stage4, wfa_vs_matrix, trim_K) is ALREADY a
keyword-only knob on `kmer2ltr.align.classify` on `kmer2ltr.align.classify`. Nothing in `kmer2ltr/` is modified by this
module or by running it.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
from functools import lru_cache
from itertools import islice
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import parasail  # noqa: E402

from kmer2ltr.align import calibrate, classify, discover  # noqa: E402
from kmer2ltr.fasta import read_fasta  # noqa: E402
from kmer2ltr.runner import COLUMNS, format_row  # noqa: E402
from kmer2ltr.runner import run as run_classify  # noqa: E402
from kmer2ltr.scoring import ALPHABET, GENERIC_MATRIX, SCALE  # noqa: E402
from bench.simulate import simulate_element  # noqa: E402
from bench.gold_robustness import D_GRID, FLANK_GRID, SOURCE_GRID  # noqa: E402,F401


def _fmt(v) -> str:
    """None -> "NA", matching kmer2ltr.runner.format_row's convention everywhere
    a truth/prediction TSV in this module can carry an undefined value."""
    return "NA" if v is None else str(v)


# =========================================================================== #
# make_dataset -- synthetic grid from library-consensus truth
# =========================================================================== #

# Required truth columns for score_run (elem_id, ltr5_start, ltr5_end,
# ltr3_start, ltr3_end, d_nominal, realized_k2p, flank5, flank3) plus two
# informational extras (orig_elem_id, kappa) that score_run never reads.
SIM_TRUTH_COLUMNS = [
    "elem_id", "orig_elem_id", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end",
    "d_nominal", "kappa", "realized_k2p", "flank5", "flank3",
]


def make_dataset(truth_fa, truth_tsv, out_fa, out_truth_tsv, ds, kappas, flanks,
                  n_per_cell, seed) -> int:
    """Build a simulated benchmark grid from the library-consensus ground truth.

    For every truth element, split its sequence into (ltr, internal) using
    the truth TSV's own coordinates, then call `bench.simulate.simulate_element`
    once per (d, kappa, flank, rep) cell -- reusing the evolution model
    rather than writing a second one. Writes a FASTA record and a truth-TSV row per
    generated sequence; a unique id is assigned per cell (`{elem_id}__d{d}
    __k{kappa}__f{flank}__r{rep}`) since one source element produces many
    output records.

    `truth_fa` and `truth_tsv` are joined by RECORD ORDER, not by `elem_id`
    lookup -- the convention used throughout this codebase
    (`bench.gold_subset.select_gold`, `bench.gold_robustness.score_gold`):
    build_truth.py writes the two files together, in the same order, and an
    id-keyed join would be wrong the moment an id is not unique (which
    nothing here guarantees). A single `random.Random(seed)` is created once
    and threaded through the whole grid, never re-seeded per element or per
    cell (same footgun `gold_robustness.py` calls out: many short-lived
    `Random()` instances correlate cells that are supposed to be independent).

    Returns the number of records written (== output FASTA record count ==
    output truth-TSV row count, always).
    """
    rng = random.Random(seed)
    n = 0
    with open(truth_tsv, newline="") as tf, open(out_fa, "w") as fa, \
            open(out_truth_tsv, "w") as tsv:
        treader = csv.DictReader(tf, delimiter="\t")
        tsv.write("\t".join(SIM_TRUTH_COLUMNS) + "\n")
        for row, (elem_id, seq) in zip(treader, read_fasta(truth_fa)):
            ltr_len = int(row["ltr_len"])
            ltr5_start, ltr5_end = int(row["ltr5_start"]), int(row["ltr5_end"])
            ltr3_start = int(row["ltr3_start"])
            if ltr5_end - ltr5_start + 1 != ltr_len or len(seq) < ltr3_start:
                print(f"make_dataset: skipping malformed truth row {row.get('elem_id')!r}",
                      file=sys.stderr)
                continue
            ltr = seq[ltr5_start - 1: ltr5_end]
            internal = seq[ltr5_end: ltr3_start - 1]
            for d in ds:
                for kappa in kappas:
                    for flank in flanks:
                        for rep in range(n_per_cell):
                            cell_id = f"{elem_id}__d{d}__k{kappa}__f{flank}__r{rep}"
                            sseq, truth = simulate_element(ltr, internal, d, kappa,
                                                            flank, flank, rng)
                            fa.write(f">{cell_id}\n{sseq}\n")
                            # realized_k2p can be None (k2p_distance is undefined on a
                            # saturated pair -- a genuine, spec'd outcome, not a bug: short
                            # LTRs at high d can saturate by chance). str(None) == "None"
                            # would silently break score_run's "NA" convention, so format
                            # it the same way runner.format_row does everywhere else.
                            tsv.write("\t".join(_fmt(v) for v in [
                                cell_id, elem_id, truth["ltr5_start"], truth["ltr5_end"],
                                truth["ltr3_start"], truth["ltr3_end"], truth["d_nominal"],
                                truth["kappa"], truth["realized_k2p"], truth["flank5"],
                                truth["flank3"],
                            ]) + "\n")
                            n += 1
    return n


# =========================================================================== #
# score_run -- generic scoring for any (pred_tsv, truth_tsv) pair sharing the
# simulate-style truth schema above
# =========================================================================== #

def score_run(pred_tsv, truth_tsv) -> dict:
    """Score a prediction TSV (full `Kmer2LTR` Result schema, e.g. `ablate`'s or
    `run_classify`'s output) against a truth TSV in the schema `make_dataset`
    writes: `elem_id ltr5_start ltr5_end ltr3_start ltr3_end d_nominal
    realized_k2p flank5 flank3` (extra columns are ignored).

    Joined by RECORD ORDER (see `make_dataset`'s docstring for why: the
    truth's `elem_id` column deliberately repeats the SOURCE element's id
    across every grid cell derived from it, so it is not usable as a unique
    join key -- and the tool's own output-order guarantee, "one row per
    record, in input order, always", is what makes a positional zip correct
    here rather than merely convenient).

    Returns:
      - "n": rows joined.
      - "mae_<coord>" / "bias_<coord>" for coord in
        ltr5_start/ltr5_end/ltr3_start/ltr3_end -- mean absolute / signed
        error, over rows where pred carries a numeric value for that column
        (a no_pair/too_short/... row contributes to neither).
      - "k2p_bias" / "k2p_rmse": mean(pred_k2p - realized_k2p) and its RMSE,
        over rows with a numeric predicted k2p. Compared against REALIZED
        K2P, not nominal d, per spec 6.2: this isolates TOOL error from the
        K2P estimator's own sampling variance.
      - "k2p_bias_by_d" / "k2p_rmse_by_d": the same, binned by `d_nominal`.
      - "flank_roc": {flank_len: {"n", "n_detected"}} over flank5>0 rows
        (the grid is symmetric by construction); "detected" requires BOTH
        sides called nonzero, matching `gold_robustness.score_gold`'s
        convention.
      - "status_counts": dict of predicted status -> count.
    """
    coord_fields = ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end")
    abs_err = {c: 0.0 for c in coord_fields}
    signed_err = {c: 0.0 for c in coord_fields}
    n_coord = {c: 0 for c in coord_fields}
    k2p_err_sum = k2p_err2_sum = 0.0
    n_k2p = 0
    k2p_by_d: dict = defaultdict(lambda: {"n": 0, "sum_err": 0.0, "sum_err2": 0.0})
    flank_roc: dict = defaultdict(lambda: {"n": 0, "n_detected": 0})
    status_counts: Counter = Counter()
    n = 0

    with open(truth_tsv, newline="") as tf, open(pred_tsv, newline="") as pf:
        tr = csv.DictReader(tf, delimiter="\t")
        pr = csv.DictReader(pf, delimiter="\t")
        for t, p in zip(tr, pr):
            n += 1
            status_counts[p.get("status", "NA")] += 1

            for c in coord_fields:
                pv = p.get(c, "NA")
                if pv in ("NA", ""):
                    continue
                err = float(pv) - float(t[c])
                abs_err[c] += abs(err)
                signed_err[c] += err
                n_coord[c] += 1

            pk = p.get("k2p", "NA")
            rk = t.get("realized_k2p", "NA")
            # realized_k2p is "NA" whenever the true evolved pair itself saturated
            # k2p_distance (short LTR + high d; see make_dataset) -- there is no
            # ground truth to score against on those rows, so they are excluded
            # from k2p_bias/rmse entirely rather than raising or silently
            # coercing "NA" to 0.
            if pk not in ("NA", "") and rk not in ("NA", "", "None"):
                err = float(pk) - float(rk)
                k2p_err_sum += err
                k2p_err2_sum += err * err
                n_k2p += 1
                dbin = round(float(t["d_nominal"]), 6)
                row = k2p_by_d[dbin]
                row["n"] += 1
                row["sum_err"] += err
                row["sum_err2"] += err * err

            f5t, f3t = int(float(t["flank5"])), int(float(t["flank3"]))
            if f5t > 0 or f3t > 0:
                flen = f5t or f3t
                row = flank_roc[flen]
                row["n"] += 1
                f5c, f3c = p.get("flank5_len", "NA"), p.get("flank3_len", "NA")
                detected = (f5c not in ("NA", "") and int(f5c) > 0
                            and f3c not in ("NA", "") and int(f3c) > 0)
                if detected:
                    row["n_detected"] += 1

    out: dict = {"n": n, "status_counts": dict(status_counts)}
    for c in coord_fields:
        out[f"mae_{c}"] = abs_err[c] / n_coord[c] if n_coord[c] else None
        out[f"bias_{c}"] = signed_err[c] / n_coord[c] if n_coord[c] else None
    out["k2p_bias"] = k2p_err_sum / n_k2p if n_k2p else None
    out["k2p_rmse"] = (k2p_err2_sum / n_k2p) ** 0.5 if n_k2p else None
    out["k2p_bias_by_d"] = {d: v["sum_err"] / v["n"] for d, v in k2p_by_d.items() if v["n"]}
    out["k2p_rmse_by_d"] = {d: (v["sum_err2"] / v["n"]) ** 0.5
                             for d, v in k2p_by_d.items() if v["n"]}
    out["flank_roc"] = {k: dict(v) for k, v in flank_roc.items()}
    return out


# =========================================================================== #
# Shared bounded-submission parallel driver
# =========================================================================== #

def _parallel_map(fasta_path, work_fn, threads: int, *extra_args):
    """Yield work_fn(seq_id, seq, *extra_args) in input order.

    Bounded ProcessPoolExecutor submission (never Executor.map: it drains the
    whole input generator before dispatching a single task, which would
    materialise all of `fasta_path` in the parent -- see `kmer2ltr/runner.py`'s
    own docstring for the measured RSS cost of that). `work_fn` must be a
    module-level function (picklable by reference); `extra_args` must be
    picklable values only -- this is what keeps a non-picklable
    `parasail.Matrix` from ever being handed to `submit` (see
    `_matrix_for_spec` below for how the fixed-matrix ablations route around
    that instead).
    """
    records = read_fasta(fasta_path)
    if threads <= 1:
        for sid, seq in records:
            yield work_fn(sid, seq, *extra_args)
        return
    max_inflight = threads * 4
    with ProcessPoolExecutor(max_workers=threads) as pool:
        it = iter(records)
        pending = deque(pool.submit(work_fn, sid, seq, *extra_args)
                         for sid, seq in islice(it, max_inflight))
        while pending:
            yield pending.popleft().result()
            nxt = next(it, None)
            if nxt is not None:
                sid, seq = nxt
                pending.append(pool.submit(work_fn, sid, seq, *extra_args))


# =========================================================================== #
# ablate -- run classify() over a FASTA under one named configuration
# =========================================================================== #

# spec 6.4 / the spec 6.4 grid. "calibrated" is classify()'s own
# default in every respect -- listed explicitly (rather than left implicit)
# so this dict is a complete, readable inventory of the grid on its own.
# Every mechanism here is an EXISTING classify() keyword; ablate() adds no
# new tool behaviour, only drives it.
ABLATIONS: dict[str, dict] = {
    "calibrated": {},
    "blastn_1_3": {"_matrix_spec": "blastn_1_3"},
    "fixed_1_1": {"_matrix_spec": "fixed_1_1"},
    "no_stage3": {"use_stage3": False},
    "no_stage4": {"use_stage4": False},
    "period_outermost": {"period_rule": "outermost"},
    "wfa_vs_matrix": {"refine": "matrix"},
    "trim_0": {"trim": 0},
    "trim_3": {"trim": 3},
    "trim_5": {"trim": 5},
    "trim_10": {"trim": 10},
}


def _build_fixed_matrix(match: int, mismatch: int) -> "parasail.Matrix":
    """A flat +match/-mismatch matrix at Kmer2LTR's SCALE, N=0 both directions
    -- the same construction as `kmer2ltr.scoring._generic`, parameterised."""
    m = parasail.matrix_create(ALPHABET, 0, 0)
    for i in range(4):
        for j in range(4):
            m.set_value(i, j, match * SCALE if i == j else mismatch * SCALE)
    n_idx = ALPHABET.index("N")
    for i in range(len(ALPHABET)):
        m.set_value(i, n_idx, 0)
        m.set_value(n_idx, i, 0)
    return m


@lru_cache(maxsize=None)
def _matrix_for_spec(spec: str) -> "parasail.Matrix":
    """Resolve a picklable matrix SPEC (a plain string) to a parasail.Matrix,
    INSIDE whichever process calls it.

    parasail.Matrix wraps a raw ctypes pointer and cannot be pickled
    (confirmed empirically: `pickle.dumps(GENERIC_MATRIX)` raises
    `ValueError: ctypes objects containing pointers cannot be pickled` -- and
    that is true of ANY parasail.Matrix instance, not just custom ones,
    since GENERIC_MATRIX is one too). So a fixed-matrix ablation config never
    puts a Matrix object in the classify_kw dict handed to
    ProcessPoolExecutor.submit; it puts this function's NAME instead (a
    plain string, trivially picklable), and each worker process builds its
    own Matrix locally on first use, cached via lru_cache for the rest of
    that worker's lifetime so the O(1) matrix construction is paid once per
    worker, not once per record.
    """
    if spec == "blastn_1_3":
        return _build_fixed_matrix(1, -3)
    if spec == "fixed_1_1":
        return GENERIC_MATRIX
    raise ValueError(f"unknown matrix spec {spec!r}")


def _ablate_work(sid, seq, kw):
    kw = dict(kw)
    spec = kw.pop("_matrix_spec", None)
    if spec is not None:
        kw["matrix"] = _matrix_for_spec(spec)
    return classify(sid, seq, **kw)


def ablate(name: str, fasta_path, out_tsv, *, threads: int = 1, cs: bool = False,
           **overrides) -> dict:
    """Run classify() over every record of `fasta_path` under ablation
    configuration `name` (see ABLATIONS), writing one prediction row per
    record to `out_tsv` in the SAME column order and input order as
    production `Kmer2LTR` output (`kmer2ltr.runner.COLUMNS`) -- so `out_tsv`
    scores exactly like a normal run, with `score_run` / `score_gold_grid`.

    `overrides` layers extra classify() kwargs on top of `ABLATIONS[name]`
    (overrides win) -- e.g. `ablate("calibrated", fa, out, t_bits=8.0)` runs
    the T_BITS sweep as a variant of the default configuration rather than
    inventing a second mechanism for it.

    Scoring is a separate step -- ablate's only job () is to
    produce the prediction TSV. Returns a small run summary: `{"name",
    "kwargs", "fasta", "out_tsv", "n", "elapsed_s"}`.
    """
    if name not in ABLATIONS:
        raise ValueError(f"unknown ablation {name!r}; choices: {sorted(ABLATIONS)}")
    kw = {**ABLATIONS[name], **overrides}
    needs_matrix = "_matrix_spec" in kw

    t0 = time.time()
    if needs_matrix:
        n = 0
        with open(out_tsv, "w") as out:
            out.write("\t".join(COLUMNS) + "\n")
            for r in _parallel_map(fasta_path, _ablate_work, threads, kw):
                out.write(format_row(r) + "\n")
                n += 1
    else:
        # No non-picklable objects anywhere in kw: production code
        # (already tested, already handles bounded submission + resume)
        # drives this directly rather than duplicating that logic a second
        # time.
        with open(out_tsv, "w") as out:
            n = run_classify(fasta_path, out, threads=threads, cs=cs, **kw)

    elapsed = time.time() - t0
    return {"name": name, "kwargs": {k: v for k, v in kw.items() if not k.startswith("_")},
            "fasta": str(fasta_path), "out_tsv": str(out_tsv), "n": n, "elapsed_s": elapsed}


# =========================================================================== #
# d_hat -- expose the divergence Stage 2 already estimates internally
# =========================================================================== #

def _dhat_work(sid, seq):
    """Stage 1 (GENERIC_MATRIX) + Stage 2 calibration ONLY -- the identical
    pair `classify()` itself always runs first, before Stage 3/4/5. Exposes
    `d_hat`, which `classify()` computes but currently discards, for the
    divergence-aware-threshold analysis; production `classify()` is
    untouched."""
    hit = discover(seq, GENERIC_MATRIX)
    if hit is None:
        return sid, None, None
    _, d_hat, kappa_hat = calibrate(seq, hit)
    return sid, d_hat, kappa_hat


def compute_dhat(fasta_path, out_tsv, threads: int = 1) -> int:
    """Write `seq_id, d_hat, kappa_hat` for every record of `fasta_path`
    (`d_hat`/`kappa_hat` are "NA" when Stage 1 finds no significant hit at
    all). Same input order as `fasta_path` / any matching truth or
    prediction TSV, so it zips positionally with them."""
    n = 0
    with open(out_tsv, "w") as out:
        out.write("seq_id\td_hat\tkappa_hat\n")
        for sid, d_hat, kappa_hat in _parallel_map(fasta_path, _dhat_work, threads):
            dh = "NA" if d_hat is None else f"{d_hat:.6g}"
            kh = "NA" if kappa_hat is None else f"{kappa_hat:.6g}"
            out.write(f"{sid}\t{dh}\t{kh}\n")
            n += 1
    return n


# =========================================================================== #
# score_gold_grid -- score one ablate()-produced prediction TSV against
# the gold-perturbation truth (bench/out/gold_truth.tsv)
# =========================================================================== #

# Bin edges = midpoints between consecutive D_GRID values, so a d_hat bin
# label and a d_nominal bin label are directly comparable (same 7 labels).
_DHAT_EDGES = [(D_GRID[i] + D_GRID[i + 1]) / 2.0 for i in range(len(D_GRID) - 1)]


def dhat_bin(d_hat: float) -> float:
    """Map a continuous calibrated d_hat to the nearest D_GRID label."""
    for edge, label in zip(_DHAT_EDGES, D_GRID[:-1]):
        if d_hat < edge:
            return label
    return D_GRID[-1]


def score_gold_grid(truth_tsv, pred_tsv, dhat_tsv=None) -> dict:
    """Score one `ablate()`-produced prediction TSV against
    `bench/gold_robustness.py`'s perturbation truth.

    Extends `gold_robustness.score_gold` with what it doesn't track: (1)
    per-cell MAE of called flank length, computed UNCONDITIONALLY (a missed
    call counts as a called length of 0, so a config that "detects" fewer
    cells but reports the ones it does find more precisely cannot look
    better on MAE alone than one that detects more), and (2) optionally, the
    same per-record increments ALSO accumulated keyed by the TOOL's own
    calibrated `d_hat` (from `compute_dhat`) instead of the perturbation's
    known `d_nominal` -- for the divergence-aware-threshold analysis, which
    by construction cannot use ground truth at decision time.

    Joins truth / pred / dhat by RECORD ORDER -- same convention as
    `gold_subset.select_gold` / `gold_robustness.score_gold` (record order
    is what the whole pipeline actually guarantees; ids are not unique).

    Returns `{"by_d": {(d_nominal, flank_true): {cell}}, "by_dhat": {...} or
    None}`. Each cell: `n, n_pass, n_lost, n_false_flank` (flank_true==0
    only), `n_detected, sum_abs_err5, sum_abs_err3` (unconditional, over all
    n at that flank_true), `sum_called5_det, sum_called3_det` (conditional
    on detection -- "how accurate once found").
    """
    by_d: dict = defaultdict(lambda: defaultdict(float))
    by_dhat = defaultdict(lambda: defaultdict(float)) if dhat_tsv is not None else None
    # Keyed by (d_nominal, tool_called_a_flank) -- tool_called_a_flank is the
    # PIPELINE's own decision (flank5_len>0 or flank3_len>0), independent of
    # ground truth. K2P bias/RMSE is vs d_nominal (gold_truth.tsv carries no
    # realized-K2P ground truth for the perturbed grid, unlike make_dataset's
    # synthetic truth -- see score_run's docstring), restricted to status=pass
    # rows with a numeric k2p. This is what the trim ablation's "overall vs
    # within the flank-called subset" comparison  reads.
    by_called: dict = defaultdict(lambda: defaultdict(float))

    with ExitStack() as stack:
        tf = stack.enter_context(open(truth_tsv, newline=""))
        pf = stack.enter_context(open(pred_tsv, newline=""))
        tr = csv.DictReader(tf, delimiter="\t")
        pr = csv.DictReader(pf, delimiter="\t")
        if dhat_tsv is not None:
            df = stack.enter_context(open(dhat_tsv, newline=""))
            dr = csv.DictReader(df, delimiter="\t")
            rows = zip(tr, pr, dr)
        else:
            rows = ((t, p, None) for t, p in zip(tr, pr))

        for t, p, dh in rows:
            d_true = float(t["d_nominal"])
            flank_true = int(t["flank5"])          # == flank3 by grid construction
            is_pass = p["status"] == "pass"
            f5c = int(p["flank5_len"]) if is_pass else 0
            f3c = int(p["flank3_len"]) if is_pass else 0

            delta: dict = {"n": 1}
            if not is_pass:
                delta["n_lost"] = 1
            else:
                delta["n_pass"] = 1
                if flank_true == 0:
                    if f5c > 0 or f3c > 0:
                        delta["n_false_flank"] = 1
                elif f5c > 0 and f3c > 0:
                    delta["n_detected"] = 1
                    delta["sum_called5_det"] = f5c
                    delta["sum_called3_det"] = f3c
            if flank_true > 0:
                delta["sum_abs_err5"] = abs(f5c - flank_true)
                delta["sum_abs_err3"] = abs(f3c - flank_true)

            cell = by_d[(d_true, flank_true)]
            for k, v in delta.items():
                cell[k] += v

            if by_dhat is not None and dh is not None and dh.get("d_hat", "NA") not in ("NA", ""):
                dbin = dhat_bin(float(dh["d_hat"]))
                cellh = by_dhat[(dbin, flank_true)]
                for k, v in delta.items():
                    cellh[k] += v

            called = bool(is_pass and (f5c > 0 or f3c > 0))
            ccell = by_called[(d_true, called)]
            ccell["n"] += 1
            pk = p.get("k2p", "NA") if is_pass else "NA"
            if pk not in ("NA", ""):
                err = float(pk) - d_true
                ccell["n_k2p"] += 1
                ccell["sum_k2p_err"] += err
                ccell["sum_k2p_err2"] += err * err

    return {"by_d": {k: dict(v) for k, v in by_d.items()},
            "by_dhat": ({k: dict(v) for k, v in by_dhat.items()}
                        if by_dhat is not None else None),
            "by_called": {k: dict(v) for k, v in by_called.items()}}


def cells_to_json(cells: dict, path) -> None:
    """Persist score_gold_grid's output (tuple-keyed dicts) as JSON."""
    def _ser(d):
        return None if d is None else {f"{k[0]}|{k[1]}": v for k, v in d.items()}
    Path(path).write_text(json.dumps({"by_d": _ser(cells["by_d"]),
                                       "by_dhat": _ser(cells.get("by_dhat")),
                                       "by_called": _ser(cells.get("by_called"))}, indent=1))


def cells_from_json(path) -> dict:
    raw = json.loads(Path(path).read_text())

    def _de(d, second=int):
        if d is None:
            return None
        out = {}
        for k, v in d.items():
            a, b = k.split("|")
            out[(float(a), second(b))] = v
        return out
    return {"by_d": _de(raw["by_d"]), "by_dhat": _de(raw.get("by_dhat")),
            "by_called": _de(raw.get("by_called"), second=lambda s: s == "True")}


# =========================================================================== #
# Stage 4 "does it ever fire" diff on RAW (unperturbed) real data
# =========================================================================== #

def stage4_diff(name: str, fasta_path, with_stage4_pred_tsv, outdir, threads: int = 1) -> dict:
    """Compare use_stage4=True (the `with_stage4_pred_tsv` already on disk --
    reused, not recomputed) against use_stage4=False (`ablate("no_stage4",
    ...)`, computed fresh here) on the SAME raw real FASTA, and count how
    many records the two disagree on.

    This is deliberately run on RAW, unperturbed real data (arabidopsis,
    human, the library-consensus truth.fa) rather than the gold-perturbed
    grid: the question ("does Stage 4 -- outermost-pair recovery for
    retained nested elements -- ever change the answer on real data") is
    about whether nested/nearly-nested structure occurs anywhere in real
    input, which the gold-subset's selection filters (status=pass, no
    overextension on the raw call) would bias against by construction.

    A record counts as "changed" if status, or any of the four boundary
    coordinates, differs between the two runs.
    """
    no4_tsv = Path(outdir) / f"stage4diff_{name}_no_stage4.tsv"
    ablate("no_stage4", fasta_path, no4_tsv, threads=threads)

    coord_fields = ("status", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end")
    n = n_changed = 0
    examples = []
    with open(with_stage4_pred_tsv, newline="") as wf, open(no4_tsv, newline="") as nf:
        wr = csv.DictReader(wf, delimiter="\t")
        nr = csv.DictReader(nf, delimiter="\t")
        for w, r in zip(wr, nr):
            n += 1
            if any(w[c] != r[c] for c in coord_fields):
                n_changed += 1
                if len(examples) < 10:
                    examples.append({
                        "seq_id": w["seq_id"],
                        "with_stage4": {c: w[c] for c in coord_fields},
                        "no_stage4": {c: r[c] for c in coord_fields},
                    })
    return {"dataset": name, "n": n, "n_changed": n_changed,
            "frac_changed": n_changed / n if n else None, "examples": examples}


# =========================================================================== #
# CLI orchestration
# =========================================================================== #

def _milestone(msg: str) -> None:
    print(f"run_bench: {msg}", file=sys.stderr, flush=True)


def _run_synthetic_ablation_grid(truth_fa, truth_tsv, outdir: Path, threads: int,
                                  n_truth_sample: int, seed: int) -> None:
    """The literal spec-6.4 grid: build ONE synthetic dataset with
    `make_dataset` (subsampled to `n_truth_sample` truth elements -- keeps
    this demonstrative run's cost modest; the headline real-data ablation
    numbers come from `_run_gold_analysis`, not from here) and score every
    ABLATIONS config against it with `score_run`."""
    sample_fa = outdir / "synth_truth_sample.fa"
    sample_tsv = outdir / "synth_truth_sample.tsv"
    with open(truth_tsv, newline="") as tf:
        rows = list(islice(csv.DictReader(tf, delimiter="\t"), n_truth_sample))
    keep_ids = {r["elem_id"] for r in rows}
    with open(sample_tsv, "w") as out:
        out.write("elem_id\ttotal_len\tltr_len\tltr5_start\tltr5_end\tltr3_start\tltr3_end\n")
        for r in rows:
            out.write("\t".join(r[c] for c in
                       ["elem_id", "total_len", "ltr_len", "ltr5_start", "ltr5_end",
                        "ltr3_start", "ltr3_end"]) + "\n")
    with open(sample_fa, "w") as out:
        n_written = 0
        for sid, seq in read_fasta(truth_fa):
            if sid in keep_ids:
                out.write(f">{sid}\n{seq}\n")
                n_written += 1
            if n_written >= len(keep_ids):
                break

    grid_fa = outdir / "synth_grid.fa"
    grid_tsv = outdir / "synth_grid_truth.tsv"
    _milestone(f"make_dataset: {len(rows)} truth elements sampled -> building synthetic grid")
    n = make_dataset(sample_fa, sample_tsv, grid_fa, grid_tsv,
                      ds=list(D_GRID), kappas=[2.0], flanks=[0, 10, 20, 50, 100],
                      n_per_cell=1, seed=seed)
    _milestone(f"synthetic grid: {n} records")

    summary_path = outdir / "synth_ablation_summary.json"
    results = {}
    for name in ABLATIONS:
        pred_tsv = outdir / f"synth_pred_{name}.tsv"
        t0 = time.time()
        ablate(name, grid_fa, pred_tsv, threads=threads)
        s = score_run(pred_tsv, grid_tsv)
        results[name] = {"n": s["n"], "mae_ltr5_start": s["mae_ltr5_start"],
                          "mae_ltr3_end": s["mae_ltr3_end"], "k2p_bias": s["k2p_bias"],
                          "k2p_rmse": s["k2p_rmse"], "status_counts": s["status_counts"],
                          "flank_roc": s["flank_roc"], "elapsed_s": time.time() - t0}
        _milestone(f"synthetic ablation {name}: n={s['n']} "
                   f"mae_ltr5_start={s['mae_ltr5_start']} k2p_rmse={s['k2p_rmse']}")
    summary_path.write_text(json.dumps(results, indent=1))
    _milestone(f"wrote {summary_path}")


def _run_gold_analysis(gold_fa: Path, gold_truth: Path, outdir: Path, threads: int,
                        tbits_values: list[float], run_ablations: bool,
                        run_dhat: bool) -> None:
    """The priority-change deliverable: T_BITS sweep + ablation grid measured
    on the real gold-perturbed data, plus (optionally) the d_hat pass
    the divergence-aware-threshold analysis needs. Reuses
    `bench/out/gold_pred.tsv` where it already corresponds to the requested
    configuration, instead of recomputing it.
    """
    from kmer2ltr.align import T_BITS as DEFAULT_T_BITS

    # Resolved relative to this file, not cwd: this is a cross-reference to
    # the one canonical output, independent of --outdir or where the
    # caller's shell happens to be when this script runs.
    default_pred = Path(__file__).resolve().parent / "out" / "gold_pred.tsv"

    if run_dhat:
        dhat_tsv = outdir / "gold_dhat.tsv"
        if not dhat_tsv.exists():
            _milestone("computing d_hat (Stage 1-2 only) over the gold-perturbed grid")
            t0 = time.time()
            n = compute_dhat(gold_fa, dhat_tsv, threads=threads)
            _milestone(f"d_hat: {n} records in {time.time() - t0:.0f}s -> {dhat_tsv}")
    else:
        dhat_tsv = outdir / "gold_dhat.tsv"
        dhat_tsv = dhat_tsv if dhat_tsv.exists() else None

    _milestone(f"T_BITS sweep: {tbits_values}")
    for tb in tbits_values:
        cells_path = outdir / f"cells_tbits_{tb:g}.json"
        if cells_path.exists():
            _milestone(f"t_bits={tb:g}: cells already on disk, skipping")
            continue
        if tb == DEFAULT_T_BITS and default_pred.exists():
            pred_tsv = default_pred
            _milestone(f"t_bits={tb:g}: reusing existing {pred_tsv}")
        else:
            pred_tsv = outdir / f"gold_pred_tbits_{tb:g}.tsv"
            if not pred_tsv.exists():
                _milestone(f"t_bits={tb:g}: running classify over gold-perturbed grid")
                t0 = time.time()
                s = ablate("calibrated", gold_fa, pred_tsv, threads=threads, t_bits=tb)
                _milestone(f"t_bits={tb:g}: {s['n']} records in {s['elapsed_s']:.0f}s")
        cells = score_gold_grid(gold_truth, pred_tsv, dhat_tsv=dhat_tsv)
        cells_to_json(cells, cells_path)
        _milestone(f"t_bits={tb:g}: scored -> {cells_path}")

    if run_ablations:
        _milestone("gold ablation grid")
        for name in ABLATIONS:
            cells_path = outdir / f"gold_cells_{name}.json"
            if cells_path.exists():
                _milestone(f"ablation {name}: cells already on disk, skipping")
                continue
            if name in ("calibrated", "trim_0") and default_pred.exists():
                pred_tsv = default_pred
                _milestone(f"ablation {name}: reusing existing {pred_tsv}")
            else:
                pred_tsv = outdir / f"gold_pred_ablate_{name}.tsv"
                if not pred_tsv.exists():
                    _milestone(f"ablation {name}: running classify over gold-perturbed grid")
                    t0 = time.time()
                    s = ablate(name, gold_fa, pred_tsv, threads=threads)
                    _milestone(f"ablation {name}: {s['n']} records in {s['elapsed_s']:.0f}s")
            cells = score_gold_grid(gold_truth, pred_tsv)
            cells_to_json(cells, cells_path)
            _milestone(f"ablation {name}: scored -> {cells_path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--truth-fa", help="library-consensus truth FASTA")
    ap.add_argument("--truth-tsv", help="matching truth TSV")
    ap.add_argument("--outdir", default="bench/out")
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--all-ablations", action="store_true",
                     help="run the spec 6.4 ablation grid on a synthetic make_dataset grid")
    ap.add_argument("--n-truth-sample", type=int, default=1500,
                     help="cap on truth elements sampled for the synthetic ablation grid")
    ap.add_argument("--gold-fa", default="bench/out/gold_perturbed.fa")
    ap.add_argument("--gold-truth", default="bench/out/gold_truth.tsv")
    ap.add_argument("--tbits-sweep", default="",
                     help="comma-separated T_BITS values to sweep on the gold grid")
    ap.add_argument("--gold-ablations", action="store_true",
                     help="run the ABLATIONS grid on the gold-perturbed FASTA")
    ap.add_argument("--dhat", action="store_true",
                     help="compute per-record calibrated d_hat on the gold-perturbed grid")
    ap.add_argument("--stage4-diff", action="store_true",
                     help="diff use_stage4 True/False on --raw-dataset real FASTAs")
    ap.add_argument("--raw-dataset", action="append", default=[], metavar="NAME=PATH")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.all_ablations:
        if not (args.truth_fa and args.truth_tsv):
            print("run_bench: --all-ablations requires --truth-fa/--truth-tsv", file=sys.stderr)
            return 2
        _run_synthetic_ablation_grid(Path(args.truth_fa), Path(args.truth_tsv), outdir,
                                      args.threads, args.n_truth_sample, args.seed)

    tbits_values = [float(x) for x in args.tbits_sweep.split(",") if x.strip()]
    if tbits_values or args.gold_ablations or args.dhat:
        _run_gold_analysis(Path(args.gold_fa), Path(args.gold_truth), outdir, args.threads,
                            tbits_values, args.gold_ablations, args.dhat)

    if args.stage4_diff:
        for spec in args.raw_dataset:
            name, path = spec.split("=", 1)
            with_pred = outdir / f"gold_raw_pred_{name}.tsv"
            if not with_pred.exists():
                _milestone(f"stage4_diff {name}: no existing with-stage4 prediction at "
                           f"{with_pred}, computing")
                ablate("calibrated", path, with_pred, threads=args.threads)
            _milestone(f"stage4_diff {name}: computing no_stage4 and comparing")
            result = stage4_diff(name, path, with_pred, outdir, threads=args.threads)
            out_path = outdir / f"stage4_diff_{name}.json"
            out_path.write_text(json.dumps(result, indent=1))
            _milestone(f"stage4_diff {name}: {result['n_changed']}/{result['n']} changed "
                       f"-> {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
