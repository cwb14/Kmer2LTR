"""Grid benchmark: known perturbations of gold real elements, scored against
perturbation-defined truth -- the "does the tool recover the known-correct
answer" benchmark. See `bench/gold_subset.py`'s module docstring and
`bench/out/memo_gold_robustness.md` for the full
circularity argument; the short version: this cannot show the tool calls
boundaries correctly on arbitrary input (selection uses the tool's own
output), but it CAN show whether a known perturbation of an already-clean
call breaks the tool, and whether the known-correct answer comes back.

Pipeline, per dataset:
  1. run ltrk2p on the raw real FASTA -> pred.tsv (for gold selection only).
  2. select_gold(pred.tsv, fasta) -> gold elements; cap to --max-gold.
  3. perturb every gold element over the full grid
     d in {0, .05, .1, .2, .3, .4, .5} x flank in {0,1,2,3,5,10,20,50,100,200,500}
     x flank_source in {other, shuffle}
     -> one combined perturbed FASTA + truth TSV across all datasets.
  4. run ltrk2p on the perturbed FASTA -> pred.tsv.
  5. score_gold joins prediction to truth and produces the three headline
     tables: false-flank rate at flank=0 by divergence, flank-detection rate
     by flank length, K2P bias vs known d.

One `random.Random(seed)` per dataset is created ONCE, outside every loop,
and threaded through the whole grid build for that dataset -- never
re-seeded per element or per cell (the harness footgun this benchmark's brief
explicitly calls out: many short-lived `Random()` instances seeded from a
coarse counter can correlate cells that are supposed to be independent).
"""
from __future__ import annotations

import argparse
import csv
import random
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bench.gold_subset import select_gold, perturb, flank_from, gold_ltr_len  # noqa: E402

D_GRID = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
FLANK_GRID = (0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500)
SOURCE_GRID = ("other", "shuffle")
KAPPA_DEFAULT = 2.0

TRUTH_COLUMNS = ["seq_id", "dataset", "orig_id", "d_nominal", "flank5", "flank3",
                 "flank_source", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end",
                 "seq_len"]


# --------------------------------------------------------------------------- #
# Gold selection
# --------------------------------------------------------------------------- #

def run_ltrk2p(fasta_path, out_tsv, threads: int, python: str = sys.executable,
               verbose: bool = False) -> float:
    """Invoke the production CLI exactly as a user would: `python -m ltrk2p
    INPUT -o OUT -t N`. Subprocess, not an in-process `classify()` call --
    this measures the shipped entry point, process pool included."""
    t0 = time.time()
    cmd = [python, "-m", "ltrk2p", str(fasta_path), "-o", str(out_tsv), "-t", str(threads)]
    if verbose:
        print(f"  $ {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)
    dt = time.time() - t0
    if verbose:
        print(f"  ltrk2p finished in {dt:.1f}s -> {out_tsv}", file=sys.stderr)
    return dt


def load_gold(dataset_fasta, pred_tsv, max_gold: int, seed: int,
              min_bitscore: float = 200.0, max_k2p: float = 0.05,
              require_motif: bool = True, verbose: bool = False) -> list[tuple[str, str]]:
    """Materialise select_gold's stream into a list, capped to `max_gold`.

    Gold subsets are a small, curated fraction of the raw library (tens of
    percent pass the filter, and the filter itself is the point) -- holding
    THIS list in memory is not the same claim as holding the raw multi-GB
    library in memory, which nothing here ever does; `select_gold` and
    `read_fasta` upstream of it stream one record at a time.
    """
    gold = list(select_gold(pred_tsv, dataset_fasta, min_bitscore=min_bitscore,
                             max_k2p=max_k2p, require_motif=require_motif))
    n_total = len(gold)
    if len(gold) > max_gold:
        rng_sel = random.Random(seed)
        idx = sorted(rng_sel.sample(range(len(gold)), max_gold))
        gold = [gold[i] for i in idx]
    if verbose:
        print(f"  gold: {len(gold)} selected (of {n_total} passing filters)", file=sys.stderr)
    return gold


# --------------------------------------------------------------------------- #
# Perturbation grid
# --------------------------------------------------------------------------- #

def build_perturbed_grid(datasets: dict[str, list[tuple[str, str]]], out_fasta: Path,
                          out_truth_tsv: Path, seed: int, kappa: float = KAPPA_DEFAULT,
                          verbose: bool = False) -> int:
    """Write one combined perturbed FASTA + truth TSV spanning every dataset.

    `flank_source="shuffle"` draws from the element's own full sequence.
    `flank_source="other"` draws from a DIFFERENT gold element's internal
    region within the same dataset (its neighbour in list order, so every
    element has a distinct, deterministic donor). Both are real, composition-
    realistic sequence -- shuffle is the harder null (also matches the
    element's own composition), per `flank_from`'s docstring.
    """
    n = 0
    with open(out_fasta, "w") as fa, open(out_truth_tsv, "w") as tsv:
        tsv.write("\t".join(TRUTH_COLUMNS) + "\n")
        for ds_i, (dataset, gold) in enumerate(datasets.items()):
            N = len(gold)
            if N < 2:
                if verbose:
                    print(f"  skip {dataset}: gold subset has {N} element(s), need >= 2 "
                          "for flank_source='other'", file=sys.stderr)
                continue
            # internal regions, computed once per element (gold_ltr_len is memoized
            # in bench.gold_subset, so this and perturb()'s own lookups share one cache)
            internals = [seq[gold_ltr_len(seq): len(seq) - gold_ltr_len(seq)]
                         for _, seq in gold]
            # one Random(seed) for this WHOLE dataset's grid, bound outside every loop
            rng = random.Random(seed + ds_i * 1_000_003)
            t0 = time.time()
            for i, (orig_id, seq) in enumerate(gold):
                for d in D_GRID:
                    for flank in FLANK_GRID:
                        for source in SOURCE_GRID:
                            if source == "shuffle":
                                pool = seq
                            else:
                                pool = internals[(i + 1) % N]
                            fsrc = (lambda nn, rr, _pool=pool, _mode=source:
                                    flank_from(_pool, nn, rr, mode=_mode))
                            pseq, truth = perturb(seq, d, kappa, flank, flank, rng, fsrc)
                            cell_id = f"{dataset}__{orig_id}__d{d}__f{flank}__{source}"
                            fa.write(f">{cell_id}\n{pseq}\n")
                            tsv.write("\t".join(str(v) for v in [
                                cell_id, dataset, orig_id, d, flank, flank, source,
                                truth["ltr5_start"], truth["ltr5_end"],
                                truth["ltr3_start"], truth["ltr3_end"], len(pseq),
                            ]) + "\n")
                            n += 1
                if verbose and (i + 1) % 100 == 0:
                    print(f"  {dataset}: {i + 1}/{N} gold elements perturbed "
                          f"({n} records so far, {time.time() - t0:.0f}s)", file=sys.stderr)
            if verbose:
                print(f"  {dataset}: done, {N} gold elements x {len(D_GRID)} x "
                      f"{len(FLANK_GRID)} x {len(SOURCE_GRID)} = "
                      f"{N * len(D_GRID) * len(FLANK_GRID) * len(SOURCE_GRID)} records "
                      f"in {time.time() - t0:.0f}s", file=sys.stderr)
    return n


# --------------------------------------------------------------------------- #
# Scoring: the three headline tables
# --------------------------------------------------------------------------- #

def score_gold(pred_tsv, truth_tsv) -> dict:
    """Join prediction to truth by RECORD ORDER (same contract as
    `select_gold`: `runner.py` preserves input order, and truth/perturbed
    FASTA were written in that same order), and accumulate the three
    headline tables in one streaming pass -- no per-row list is retained.

    Returns a dict with three sub-dicts, each keyed by `(dataset, x)` (plus
    a pooled `"ALL"` dataset key):
      - `false_flank`: keyed by (dataset, d_nominal); flank=0 cells only.
        `n_pass` excludes cells where the tool lost the pair entirely
        (status != pass) -- that is a different failure mode (see `n_lost`),
        not "wrongly called a flank", and folding it in would inflate the
        headline number with something it does not mean.
      - `flank_detect`: keyed by (dataset, flank_len); flank>0 cells only.
        "detected" requires BOTH sides called (flank5_len>0 AND
        flank3_len>0) -- the grid is symmetric (flank5==flank3==flank by
        construction), so a whole-element pass/fail is the natural unit.
      - `k2p_bias`: keyed by (dataset, d_nominal); restricted to flank=0
        cells that were ALSO correctly bounded (no false flank called) --
        otherwise a boundary miss and a divergence-estimator error would be
        conflated into one number that answers neither question cleanly.
    """
    false_flank: dict = defaultdict(lambda: defaultdict(float))
    flank_detect: dict = defaultdict(lambda: defaultdict(float))
    flank_detect_src: dict = defaultdict(lambda: defaultdict(float))
    k2p_bias: dict = defaultdict(lambda: defaultdict(float))

    with open(truth_tsv, newline="") as tf, open(pred_tsv, newline="") as pf:
        tr = csv.DictReader(tf, delimiter="\t")
        pr = csv.DictReader(pf, delimiter="\t")
        for t, p in zip(tr, pr):
            ds_keys = (t["dataset"], "ALL")
            d = float(t["d_nominal"])
            flank5_true = int(t["flank5"])
            flank3_true = int(t["flank3"])
            source = t["flank_source"]
            is_pass = p["status"] == "pass"

            if flank5_true == 0 and flank3_true == 0:
                for ds in ds_keys:
                    row = false_flank[(ds, d)]
                    row["n"] += 1
                    if not is_pass:
                        row["n_lost"] += 1
                        continue
                    row["n_pass"] += 1
                    f5c, f3c = int(p["flank5_len"]), int(p["flank3_len"])
                    if f5c > 0 or f3c > 0:
                        row["n_false_flank"] += 1
                    elif p["k2p"] not in ("NA", ""):
                        krow = k2p_bias[(ds, d)]
                        err = float(p["k2p"]) - d
                        krow["n"] += 1
                        krow["sum_err"] += err
                        krow["sum_err2"] += err * err

            elif flank5_true > 0:                      # == flank3_true (symmetric grid)
                for ds in ds_keys:
                    row = flank_detect[(ds, flank5_true)]
                    row["n"] += 1
                    row_s = flank_detect_src[(ds, flank5_true, source)]
                    row_s["n"] += 1
                    if is_pass:
                        f5c, f3c = int(p["flank5_len"]), int(p["flank3_len"])
                        if f5c > 0 and f3c > 0:
                            for r in (row, row_s):
                                r["n_detected"] += 1
                                r["sum_called5"] += f5c
                                r["sum_called3"] += f3c

    return {"false_flank": false_flank, "flank_detect": flank_detect,
            "flank_detect_by_source": flank_detect_src, "k2p_bias": k2p_bias}


# --------------------------------------------------------------------------- #
# Markdown rendering
# --------------------------------------------------------------------------- #

def render_tables(scores: dict, datasets: list[str]) -> str:
    out = []
    ds_order = list(datasets) + ["ALL"]

    out.append("### Headline 1 -- false-flank rate at flank=0, by divergence\n")
    out.append("Of gold elements perturbed ONLY by mutation (no flank added), the "
                "fraction the tool wrongly reports as over-extended (rate is of "
                "`n_pass`, i.e. conditional on the tool still finding the pair at "
                "all -- `n_lost` is a separate failure mode, listed alongside).\n")
    out.append("| dataset | d | n | n_pass | n_false_flank | false_flank_rate | n_lost |")
    out.append("|---|---|---|---|---|---|---|")
    for ds in ds_order:
        for d in D_GRID:
            row = scores["false_flank"].get((ds, d))
            if not row:
                continue
            rate = row["n_false_flank"] / row["n_pass"] if row["n_pass"] else float("nan")
            out.append(f"| {ds} | {d:g} | {int(row['n'])} | {int(row['n_pass'])} | "
                        f"{int(row['n_false_flank'])} | {rate:.4f} | {int(row['n_lost'])} |")

    out.append("\n### Headline 2 -- flank-detection rate by flank length\n")
    out.append("flank>0 cells; \"detected\" requires both sides called nonzero. "
                "Called length is the mean over DETECTED elements only.\n")
    out.append("| dataset | flank_true | n | n_detected | detection_rate | "
                "mean_called5 | mean_called3 |")
    out.append("|---|---|---|---|---|---|---|")
    for ds in ds_order:
        for fl in FLANK_GRID:
            if fl == 0:
                continue
            row = scores["flank_detect"].get((ds, fl))
            if not row:
                continue
            rate = row["n_detected"] / row["n"] if row["n"] else float("nan")
            mc5 = row["sum_called5"] / row["n_detected"] if row["n_detected"] else float("nan")
            mc3 = row["sum_called3"] / row["n_detected"] if row["n_detected"] else float("nan")
            out.append(f"| {ds} | {fl} | {int(row['n'])} | {int(row['n_detected'])} | "
                        f"{rate:.4f} | {mc5:.1f} | {mc3:.1f} |")

    out.append("\n### Headline 2b -- flank-detection rate by flank length and source\n")
    out.append("Same as above, split by `flank_source` -- `shuffle` is the harder "
                "null (matches the element's own composition).\n")
    out.append("| dataset | flank_true | source | n | n_detected | detection_rate |")
    out.append("|---|---|---|---|---|---|")
    for ds in ds_order:
        for fl in FLANK_GRID:
            if fl == 0:
                continue
            for source in SOURCE_GRID:
                row = scores["flank_detect_by_source"].get((ds, fl, source))
                if not row:
                    continue
                rate = row["n_detected"] / row["n"] if row["n"] else float("nan")
                out.append(f"| {ds} | {fl} | {source} | {int(row['n'])} | "
                            f"{int(row['n_detected'])} | {rate:.4f} |")

    out.append("\n### Headline 3 -- K2P bias vs known d\n")
    out.append("flank=0 cells that were ALSO correctly bounded (no false flank "
                "called), so this isolates the K2P estimator from boundary-call "
                "error. `bias = mean(k2p_called - d_nominal)`, `rmse` likewise.\n")
    out.append("| dataset | d | n | bias | rmse |")
    out.append("|---|---|---|---|---|")
    for ds in ds_order:
        for d in D_GRID:
            row = scores["k2p_bias"].get((ds, d))
            if not row or not row["n"]:
                continue
            n = row["n"]
            bias = row["sum_err"] / n
            rmse = (row["sum_err2"] / n) ** 0.5
            out.append(f"| {ds} | {d:g} | {int(n)} | {bias:+.4f} | {rmse:.4f} |")

    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dataset", action="append", required=True, metavar="NAME=PATH",
                     help="a real-element FASTA to draw a gold subset from; repeatable")
    ap.add_argument("--outdir", default="bench/out")
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--max-gold", type=int, default=800,
                     help="cap on gold elements per dataset (deterministic sample)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--kappa", type=float, default=KAPPA_DEFAULT)
    ap.add_argument("--min-bitscore", type=float, default=200.0)
    ap.add_argument("--max-k2p", type=float, default=0.05)
    ap.add_argument("--no-motif", action="store_true",
                     help="drop the canonical TG..CA terminus requirement")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ds_paths = dict(a.split("=", 1) for a in args.dataset)

    gold: dict[str, list[tuple[str, str]]] = {}
    for name, path in ds_paths.items():
        print(f"ltrk2p-bench: [{name}] classifying raw input for gold selection", file=sys.stderr)
        raw_pred = outdir / f"gold_raw_pred_{name}.tsv"
        run_ltrk2p(path, raw_pred, args.threads, verbose=args.verbose)
        g = load_gold(path, raw_pred, args.max_gold, args.seed,
                      min_bitscore=args.min_bitscore, max_k2p=args.max_k2p,
                      require_motif=not args.no_motif, verbose=True)
        gold[name] = g
        print(f"ltrk2p-bench: [{name}] gold subset = {len(g)} elements", file=sys.stderr)

    perturbed_fa = outdir / "gold_perturbed.fa"
    truth_tsv = outdir / "gold_truth.tsv"
    print("ltrk2p-bench: building perturbation grid "
          f"({len(D_GRID)} d x {len(FLANK_GRID)} flank x {len(SOURCE_GRID)} source = "
          f"{len(D_GRID) * len(FLANK_GRID) * len(SOURCE_GRID)} cells/element)", file=sys.stderr)
    n_records = build_perturbed_grid(gold, perturbed_fa, truth_tsv, args.seed,
                                      kappa=args.kappa, verbose=True)
    print(f"ltrk2p-bench: wrote {n_records} perturbed records", file=sys.stderr)

    pred_tsv = outdir / "gold_pred.tsv"
    print("ltrk2p-bench: classifying perturbed grid", file=sys.stderr)
    dt = run_ltrk2p(perturbed_fa, pred_tsv, args.threads, verbose=True)
    print(f"ltrk2p-bench: perturbed-grid classification: {n_records} records in {dt:.1f}s "
          f"({n_records / dt:.0f}/s)", file=sys.stderr)

    print("ltrk2p-bench: scoring", file=sys.stderr)
    scores = score_gold(pred_tsv, truth_tsv)
    tables = render_tables(scores, list(gold.keys()))
    (outdir / "gold_headline_tables.md").write_text(tables)
    print(tables)
    return 0


if __name__ == "__main__":
    sys.exit(main())
