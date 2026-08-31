"""Derive the divergence-aware `T_BITS` schedule under a PRE-REGISTERED rule.

Task 14 measured that a divergence-aware threshold is a genuine Pareto
improvement over any fixed one, and then deliberately declined to ship a
schedule, for a specific reason: "the per-bin schedule needs to be chosen with
an explicit detection floor in mind, not just a false-flank target, or it will
overcorrect exactly where overextension is most consequential to still catch."

This file is that explicit floor. It is written and committed BEFORE the sweep
it consumes has been run, so the rule cannot be tuned to its own answer.

--------------------------------------------------------------------------- #
THE RULE
--------------------------------------------------------------------------- #

Reference throughout is the shipped flat `t_bits = 10`.

1. `FF_TARGET` = the false-flank rate the flat reference achieves POOLED over
   all bins. The schedule is therefore calibrated to be no worse overall than
   today; what it changes is the *uniformity* of that rate across divergence.
   A claimed flank should carry the same weight of evidence whether the element
   is 2% or 40% diverged, and at a fixed threshold it does not: Task 14 measured
   the rate swinging 40-80x across the `d_hat` range at constant `t_bits`.

2. A value `t` is ADMISSIBLE in bin `b` only if, in that same bin, relative to
   the flat reference:
       large-flank detection is within 2 points   (gold: 50 and 100 bp;
                                                   homology: 45 and 65 bp)
       mid-flank detection is within 5 points     (gold: 20 bp; homology: 25 bp)
       the pair-loss rate has risen by at most 2 points
   Raising `t_bits` costs pairs outright as well as flank calls -- Task 14 saw
   `n_pass` fall at high divergence -- so pair loss is a constrained quantity,
   not a free one. Admissibility must hold on both grids wherever both have
   evidence (see step 5): the gold grid for continuity with every prior
   measurement, and the homology grid because it is the one that owes nothing to
   a motif prior or to the tool's own selection.

3. `t*(b)` = the SMALLEST admissible `t` whose false-flank rate in bin `b` is at
   or below `FF_TARGET`. If no admissible `t` reaches the target, take the
   LARGEST admissible one. Smallest-that-suffices, because every increment of
   `t_bits` above what the target needs is detection given away for nothing.

4. Bins are then forced monotone non-decreasing by a running maximum: a more
   diverged element may never require LESS evidence to claim a flank than a less
   diverged one. This is a regulariser -- it makes the schedule a shape rather
   than seven independent point estimates, and it cannot be gamed by a noisy bin.

5. A bin holding fewer than `MIN_BIN_N` flank-bearing records on the GOLD grid at
   the reference is not estimated at all; it inherits the previous bin's value.
   Where the HOMOLOGY grid falls below `MIN_BIN_N` in a bin it abstains from
   step 2 rather than vetoing it: the homology grid's substitution axis stops at
   35% mutated, so it has little mass in the top `d_hat` bins by construction,
   and letting a structurally empty bin veto would discard the gold grid's
   evidence exactly where the false-flank problem is worst.

   (Amended 2026-08-31, after the rule was first committed but BEFORE the sweep
   it consumes had produced a single number -- the amendment is forced by the
   grid's known divergence coverage, not by any result.)
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bench.homology_grid import cells_from_json as hom_cells_from_json  # noqa: E402
from bench.run_bench import cells_from_json as gold_cells_from_json  # noqa: E402

REFERENCE_T = 10.0
LARGE_TOL = 0.02
MID_TOL = 0.05
LOST_TOL = 0.02
MIN_BIN_N = 200
DHAT_LABELS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
GOLD_LARGE, GOLD_MID = (50, 100), 20
HOM_LARGE, HOM_MID = (45, 65), 25


def _acc(cells, keep):
    out = defaultdict(float)
    for k, v in cells.items():
        if keep(k):
            for name, x in v.items():
                out[name] += x
    return out


def _rate(num, den):
    return None if not den else num / den


def gold_bin_stats(cells_path) -> dict:
    """{d_hat bin: {ff, lost, det<flank>, n_flank}} from one gold sweep point."""
    by = gold_cells_from_json(cells_path)["by_dhat"]
    if by is None:
        raise SystemExit(f"{cells_path} carries no by_dhat view; rerun with --dhat-for")
    out = {}
    for b in DHAT_LABELS:
        zero = _acc(by, lambda k, b=b: k[0] == b and k[1] == 0)
        allc = _acc(by, lambda k, b=b: k[0] == b)
        row = {"ff": _rate(zero.get("n_false_flank", 0.0), zero.get("n_pass", 0.0)),
               "lost": _rate(allc.get("n_lost", 0.0), allc.get("n", 0.0)),
               "n_flank": sum(v.get("n", 0.0) for k, v in by.items()
                              if k[0] == b and k[1] > 0)}
        for f in GOLD_LARGE + (GOLD_MID,):
            c = _acc(by, lambda k, b=b, f=f: k[0] == b and k[1] == f)
            row[f"det{f}"] = _rate(c.get("n_detected", 0.0), c.get("n", 0.0))
        out[b] = row
    return out


def hom_bin_stats(cells_path, source=None) -> dict:
    by = hom_cells_from_json(cells_path)["by_dhat"]
    if by is None:
        raise SystemExit(f"{cells_path} carries no by_dhat view; rerun with --dhat-for")

    def sel(k, b):
        src, panel, dh, ir, f = k
        return dh == b and panel == "subs" and (source is None or src == source)

    out = {}
    for b in DHAT_LABELS:
        zero = _acc(by, lambda k, b=b: sel(k, b) and k[4] == 0)
        allc = _acc(by, lambda k, b=b: sel(k, b))
        row = {"ff": _rate(zero.get("n_false_flank", 0.0), zero.get("n_located", 0.0)),
               "lost": _rate(allc.get("n_lost", 0.0), allc.get("n", 0.0)),
               "n_flank": sum(v.get("n", 0.0) for k, v in by.items()
                              if sel(k, b) and k[4] > 0)}
        for f in HOM_LARGE + (HOM_MID,):
            c = _acc(by, lambda k, b=b, f=f: sel(k, b) and k[4] == f)
            row[f"det{f}"] = _rate(c.get("n_detected", 0.0), c.get("n", 0.0))
        out[b] = row
    return out


def _admissible(stat, ref, large, mid) -> bool:
    """Rule step 2, on one grid, in one bin."""
    for f in large:
        a, r = stat.get(f"det{f}"), ref.get(f"det{f}")
        if a is not None and r is not None and a < r - LARGE_TOL:
            return False
    a, r = stat.get(f"det{mid}"), ref.get(f"det{mid}")
    if a is not None and r is not None and a < r - MID_TOL:
        return False
    a, r = stat.get("lost"), ref.get("lost")
    if a is not None and r is not None and a > r + LOST_TOL:
        return False
    return True


def derive(gold_by_t: dict, hom_by_t: dict, verbose: bool = True):
    """Apply the rule. `*_by_t` map t_bits -> per-bin stats."""
    ts = sorted(set(gold_by_t) & set(hom_by_t))
    if REFERENCE_T not in ts:
        raise SystemExit(f"reference t={REFERENCE_T} missing from the sweep ({ts})")

    ref_g, ref_h = gold_by_t[REFERENCE_T], hom_by_t[REFERENCE_T]
    # Rule step 1: the target is the reference's OWN pooled false-flank rate.
    num = sum(ref_g[b]["ff"] * ref_g[b]["n_flank"] for b in DHAT_LABELS
              if ref_g[b]["ff"] is not None)
    den = sum(ref_g[b]["n_flank"] for b in DHAT_LABELS if ref_g[b]["ff"] is not None)
    ff_target = num / den if den else 0.0

    chosen: dict[float, float] = {}
    rows = []
    for b in DHAT_LABELS:
        if ref_g[b]["n_flank"] < MIN_BIN_N:
            rows.append({"bin": b, "t": None, "why": "gold bin too small"})
            continue
        hom_votes = ref_h[b]["n_flank"] >= MIN_BIN_N
        adm = [t for t in ts
               if _admissible(gold_by_t[t][b], ref_g[b], GOLD_LARGE, GOLD_MID)
               and (not hom_votes
                    or _admissible(hom_by_t[t][b], ref_h[b], HOM_LARGE, HOM_MID))]
        if not adm:
            adm = [REFERENCE_T]
        meets = [t for t in adm
                 if gold_by_t[t][b]["ff"] is not None
                 and gold_by_t[t][b]["ff"] <= ff_target]
        t = min(meets) if meets else max(adm)
        chosen[b] = t
        rows.append({"bin": b, "t": t, "admissible": adm,
                     "ff_at_t": gold_by_t[t][b]["ff"],
                     "ff_at_ref": ref_g[b]["ff"],
                     "why": "smallest meeting target" if meets else "largest admissible"})

    # Rule step 5 then step 4: fill unestimated bins forward, then running max.
    schedule, last = [], None
    for b in DHAT_LABELS:
        t = chosen.get(b, last)
        if t is None:
            t = REFERENCE_T
        last = t if last is None else max(last, t)
        schedule.append((b, last))

    if verbose:
        print(f"FF_TARGET (pooled false-flank rate of flat t={REFERENCE_T:g}) = "
              f"{ff_target:.4f}\n")
        print("| d_hat bin | chosen t | ff at t | ff at ref | admissible | rule |")
        print("|---|---|---|---|---|---|")
        for r in rows:
            adm = ",".join(f"{x:g}" for x in r.get("admissible", []))
            ff = r.get("ff_at_t"); fr = r.get("ff_at_ref")
            print(f"| {r['bin']:g} | {r['t'] if r['t'] is not None else '--'} | "
                  f"{ff if ff is None else f'{ff:.4f}'} | "
                  f"{fr if fr is None else f'{fr:.4f}'} | {adm} | {r['why']} |")
        print("\nAfter forward-fill and monotone running max:")
        for b, t in schedule:
            print(f"  d_hat < {b:g} ... {t:g}")
    return schedule, ff_target, rows


def as_python(schedule) -> str:
    """The literal to paste into `align.T_BITS_SCHEDULE`.

    Bin LABELS are the D_GRID points; the schedule's cut points are the bin
    upper edges, i.e. the midpoints between consecutive labels, so that
    `t_bits_for(d_hat)` reproduces `dhat_bin(d_hat)` exactly.
    """
    edges = [(DHAT_LABELS[i] + DHAT_LABELS[i + 1]) / 2.0
             for i in range(len(DHAT_LABELS) - 1)] + [float("inf")]
    out, prev = [], None
    for (b, t), e in zip(schedule, edges):
        if prev is not None and t == prev:
            out[-1] = (e, t)               # merge equal-valued adjacent bins
        else:
            out.append((e, t))
        prev = t
    body = ",\n    ".join(
        f"({'float(\"inf\")' if e == float('inf') else f'{e:g}'}, {t:g})" for e, t in out)
    return f"T_BITS_SCHEDULE = (\n    {body},\n)"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="bench/out/improve")
    ap.add_argument("--base", default="baseline")
    ap.add_argument("--tbits", default="2,5,8,10,15,20,30")
    ap.add_argument("--source", default=None, help="homology panel to constrain on")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir)

    gold_by_t, hom_by_t = {}, {}
    for t in [float(x) for x in args.tbits.split(",")]:
        tag = f"{args.base}_t{t:g}"
        gold_by_t[t] = gold_bin_stats(outdir / f"cells_gold_{tag}.json")
        hom_by_t[t] = hom_bin_stats(outdir / f"cells_hom_{tag}.json", args.source)

    schedule, ff_target, rows = derive(gold_by_t, hom_by_t)
    print("\n" + as_python(schedule))
    (outdir / "schedule.json").write_text(json.dumps(
        {"schedule": schedule, "ff_target": ff_target, "rows": rows,
         "reference_t": REFERENCE_T, "large_tol": LARGE_TOL, "mid_tol": MID_TOL,
         "lost_tol": LOST_TOL}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
