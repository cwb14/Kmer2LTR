"""Measurement campaign for the 2026-08-31 improvement round.

Every candidate change is a keyword on `ltrk2p.align.classify`, so this module
is pure orchestration: it names configurations, runs each over both benchmark
grids, and scores them. No tool behaviour lives here.

The two grids answer different questions and are both reported:

  * `bench/out/gold_perturbed.fa` -- 260,876 records, real gold elements under
    known perturbation. Continuity with every prior ltrk2p measurement, but
    motif-filtered and tool-selected (see `gold_subset.py`'s circularity note).
  * `bench/out/homology_grid.fa` -- perfect library-consensus elements under
    known substitutions, flanks and indels, with an EXACT true alignment. No
    motif, no TSD, no tool involvement in selection: homology only.

`BASELINE` reproduces the shipped pre-campaign behaviour exactly, so every
one-at-a-time configuration differs from it in a single knob and the equivalence
of the refactor itself is checkable against `bench/out/gold_pred_tbits_10.tsv`.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ltrk2p.align import (GENERIC_MATRIX, SIG_GAPS, calibrate_full, discover,  # noqa: E402
                          ltr_spans)
from bench.homology_grid import score_homology  # noqa: E402
from bench.homology_grid import cells_to_json as hom_cells_to_json  # noqa: E402
from bench.run_bench import _parallel_map, ablate, cells_to_json, score_gold_grid  # noqa: E402

# The shipped pre-campaign pipeline, stated explicitly rather than left implicit
# as "the defaults" -- the defaults are exactly what this campaign is changing.
BASELINE = dict(t_bits=10.0, snap_mode="binary", inner="none", comp="element",
                gap_scheme="legacy", keep_weak=False, stage4_recal=False)

CONFIGS: dict[str, dict] = {
    "baseline": {},
    "stage4_rerun": dict(stage4_recal=True),      # item 2
    "keep_weak": dict(keep_weak=True),            # item 3
    "graded": dict(snap_mode="graded"),           # item 4
    "inner_joint": dict(inner="joint"),           # item 5
    "comp_core": dict(comp="core"),               # item 6
    "gaps_static": dict(gap_scheme="static"),     # item 7a
    "gaps_adaptive": dict(gap_scheme="adaptive"),  # item 7b
    "schedule": dict(t_bits=None),                # item 1
}


def config_kwargs(name: str, extra: dict | None = None) -> dict:
    if name not in CONFIGS:
        raise ValueError(f"unknown config {name!r}; choices: {sorted(CONFIGS)}")
    # t_bits=None is meaningful (it selects the divergence-aware schedule), so
    # nothing here is filtered out for being None.
    return {**BASELINE, **CONFIGS[name], **(extra or {})}


def _milestone(msg: str) -> None:
    print(f"run_improve: {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# d_hat, computed under a specific configuration
# --------------------------------------------------------------------------- #

def _dhat_work(sid, seq, comp, gap_scheme):
    """Stages 1-2 only -- the identical pair `classify` runs before anything
    else. `comp` and `gap_scheme` are threaded through because they change the
    calibration, hence d_hat, hence which schedule bin a record lands in."""
    hit = discover(seq, GENERIC_MATRIX, SIG_GAPS)
    if hit is None:
        return sid, None
    return sid, calibrate_full(seq, ltr_spans(seq, hit), comp=comp,
                               gap_scheme=gap_scheme).d_hat


def compute_dhat(fasta, out_tsv, threads: int, comp: str, gap_scheme: str) -> int:
    n = 0
    with open(out_tsv, "w") as out:
        out.write("seq_id\td_hat\n")
        for sid, d in _parallel_map(fasta, _dhat_work, threads, comp, gap_scheme):
            out.write(f"{sid}\t{'NA' if d is None else f'{d:.6g}'}\n")
            n += 1
    return n


# --------------------------------------------------------------------------- #
# One configuration over both grids
# --------------------------------------------------------------------------- #

def run_config(name: str, outdir: Path, threads: int, gold_fa, gold_truth,
               hom_fa, hom_truth, extra: dict | None = None, tag: str | None = None,
               gold_dhat=None, hom_dhat=None) -> dict:
    tag = tag or name
    kw = config_kwargs(name, extra)
    summary = {"config": tag, "kwargs": {k: v for k, v in kw.items()}}
    for label, fa, truth, dhat in (("gold", gold_fa, gold_truth, gold_dhat),
                                   ("hom", hom_fa, hom_truth, hom_dhat)):
        if fa is None or not Path(fa).exists():
            continue
        pred = outdir / f"pred_{label}_{tag}.tsv"
        cells = outdir / f"cells_{label}_{tag}.json"
        if not pred.exists():
            t0 = time.time()
            _milestone(f"{tag} [{label}]: classify over {fa}")
            s = ablate("calibrated", fa, pred, threads=threads, **kw)
            _milestone(f"{tag} [{label}]: {s['n']} records in {s['elapsed_s']:.0f}s")
            summary[f"{label}_n"] = s["n"]
            summary[f"{label}_elapsed_s"] = s["elapsed_s"]
        if label == "gold":
            cells_to_json(score_gold_grid(truth, pred, dhat_tsv=dhat), cells)
        else:
            hom_cells_to_json(score_homology(truth, pred, dhat), cells)
        _milestone(f"{tag} [{label}]: scored -> {cells}")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="bench/out/improve")
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--gold-fa", default="bench/out/gold_perturbed.fa")
    ap.add_argument("--gold-truth", default="bench/out/gold_truth.tsv")
    ap.add_argument("--hom-fa", default="bench/out/homology_grid.fa")
    ap.add_argument("--hom-truth", default="bench/out/homology_truth.tsv")
    ap.add_argument("--configs", default="",
                    help="comma-separated CONFIGS names to run")
    ap.add_argument("--dhat-for", default=None, metavar="CONFIG",
                    help="compute per-record d_hat under CONFIG's calibration")
    ap.add_argument("--tbits-sweep", default="",
                    help="comma-separated t_bits values to sweep")
    ap.add_argument("--sweep-base", default="baseline",
                    help="config the t_bits sweep is layered on top of")
    ap.add_argument("--skip-gold", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gold_fa = None if args.skip_gold else args.gold_fa
    summaries = []

    dhat_paths = {"gold": None, "hom": None}
    if args.dhat_for:
        kw = config_kwargs(args.dhat_for)
        for label, fa in (("gold", gold_fa), ("hom", args.hom_fa)):
            if fa is None or not Path(fa).exists():
                continue
            path = outdir / f"dhat_{label}_{args.dhat_for}.tsv"
            if not path.exists():
                t0 = time.time()
                _milestone(f"d_hat [{label}] under {args.dhat_for}")
                n = compute_dhat(fa, path, args.threads,
                                 kw.get("comp", "element"), kw.get("gap_scheme", "legacy"))
                _milestone(f"d_hat [{label}]: {n} records in {time.time() - t0:.0f}s")
            dhat_paths[label] = path

    # One config dying must not take the campaign with it: a crash in an
    # experimental configuration is a result about that configuration, not a
    # reason to lose the eight that already ran.
    for name in [c for c in args.configs.split(",") if c.strip()]:
        try:
            summaries.append(run_config(name, outdir, args.threads, gold_fa,
                                        args.gold_truth, args.hom_fa, args.hom_truth,
                                        gold_dhat=dhat_paths["gold"],
                                        hom_dhat=dhat_paths["hom"]))
        except Exception as exc:                       # noqa: BLE001
            _milestone(f"CONFIG FAILED {name}: {type(exc).__name__}: {exc}")
            summaries.append({"config": name, "error": f"{type(exc).__name__}: {exc}"})

    for tb in [float(x) for x in args.tbits_sweep.split(",") if x.strip()]:
        tag = f"{args.sweep_base}_t{tb:g}"
        try:
            summaries.append(run_config(
                args.sweep_base, outdir, args.threads, gold_fa, args.gold_truth,
                args.hom_fa, args.hom_truth, extra={"t_bits": tb}, tag=tag,
                gold_dhat=dhat_paths["gold"], hom_dhat=dhat_paths["hom"]))
        except Exception as exc:                       # noqa: BLE001
            _milestone(f"SWEEP POINT FAILED {tag}: {type(exc).__name__}: {exc}")
            summaries.append({"config": tag, "error": f"{type(exc).__name__}: {exc}"})

    if summaries:
        path = outdir / f"summary_{int(time.time())}.json"
        path.write_text(json.dumps(summaries, indent=1))
        _milestone(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
