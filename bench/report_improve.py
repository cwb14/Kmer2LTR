"""Turn `run_improve.py`'s scored-cell JSON into comparison tables.

Reporting rules, applied uniformly so configurations cannot be flattered by the
choice of denominator:

  * Boundary accuracy is reported over LOCATED records (a record with no
    coordinates has no boundary error to average) but `lost` is always printed
    alongside, because a configuration that answers less often will otherwise
    look more accurate.
  * Called-flank error is UNCONDITIONAL: a missed flank counts as a called
    length of zero, so detecting fewer flanks more precisely cannot win.
  * K2P error on the homology grid is against the EXACT realized K2P of the
    true alignment, not against the nominal target -- that separates tool error
    from the estimator's own sampling variance.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bench.homology_grid import cells_from_json as hom_cells_from_json  # noqa: E402
from bench.run_bench import cells_from_json as gold_cells_from_json  # noqa: E402

_COORDS = ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end")


def _pool(cells, keep):
    """Sum every counter over the cells `keep(key)` selects."""
    acc = defaultdict(float)
    for k, v in cells.items():
        if not keep(k):
            continue
        for name, x in v.items():
            acc[name] += x
    return acc


def _rate(num, den):
    return float("nan") if not den else num / den


def hom_summary(cells, source=None, panel="subs", p=None, indel=None):
    """Pooled homology metrics for one slice of the grid."""
    def keep(k):
        s, pn, pt, ir, f = k
        return ((source is None or s == source) and (panel is None or pn == panel)
                and (p is None or pt == p) and (indel is None or ir == indel))

    all_c = _pool(cells, keep)
    zero = _pool(cells, lambda k: keep(k) and k[4] == 0)
    pos = _pool(cells, lambda k: keep(k) and k[4] > 0)
    loc = all_c.get("n_located", 0.0)
    mae = sum(all_c.get(f"abs_{c}", 0.0) for c in _COORDS) / (4 * loc) if loc else float("nan")
    exact = min(all_c.get(f"exact_{c}", 0.0) for c in _COORDS) / loc if loc else float("nan")
    out = {
        "n": all_c.get("n", 0.0),
        "lost": _rate(all_c.get("n_lost", 0.0), all_c.get("n", 0.0)),
        "weak": _rate(all_c.get("n_weak", 0.0), all_c.get("n", 0.0)),
        "bnd_mae": mae,
        "bnd_exact": exact,
        "ff_rate": _rate(zero.get("n_false_flank", 0.0), zero.get("n_located", 0.0)),
        "flank_mae": _rate(pos.get("called_abs5", 0.0) + pos.get("called_abs3", 0.0),
                           2 * pos.get("n", 0.0)),
        "k2p_bias": _rate(all_c.get("k2p_err", 0.0), all_c.get("n_k2p", 0.0)),
        "k2p_rmse": (all_c.get("k2p_err2", 0.0) / all_c["n_k2p"]) ** 0.5
        if all_c.get("n_k2p") else float("nan"),
    }
    for f in (5, 25, 45, 65):
        c = _pool(cells, lambda k, f=f: keep(k) and k[4] == f)
        out[f"det{f}"] = _rate(c.get("n_detected", 0.0), c.get("n", 0.0))
    return out


def gold_summary(cells, d=None):
    by_d = cells["by_d"]
    if d is not None:
        by_d = {k: v for k, v in by_d.items() if k[0] == d}
    zero = _pool(by_d, lambda k: k[1] == 0)
    out = {
        "n": sum(v.get("n", 0.0) for v in by_d.values()),
        "lost": _rate(sum(v.get("n_lost", 0.0) for v in by_d.values()),
                      sum(v.get("n", 0.0) for v in by_d.values())),
        "ff_rate": _rate(zero.get("n_false_flank", 0.0), zero.get("n_pass", 0.0)),
    }
    for f in (10, 20, 50, 100):
        c = _pool(by_d, lambda k, f=f: k[1] == f)
        out[f"det{f}"] = _rate(c.get("n_detected", 0.0), c.get("n", 0.0))
        out[f"mae{f}"] = _rate(c.get("sum_abs_err5", 0.0) + c.get("sum_abs_err3", 0.0),
                               2 * c.get("n", 0.0))
    return out


_HOM_COLS = [("n", "{:.0f}"), ("lost", "{:.4f}"), ("weak", "{:.4f}"),
             ("bnd_mae", "{:.3f}"), ("bnd_exact", "{:.4f}"), ("ff_rate", "{:.4f}"),
             ("det5", "{:.4f}"), ("det25", "{:.4f}"), ("det45", "{:.4f}"),
             ("det65", "{:.4f}"), ("flank_mae", "{:.2f}"),
             ("k2p_bias", "{:+.5f}"), ("k2p_rmse", "{:.5f}")]
_GOLD_COLS = [("n", "{:.0f}"), ("lost", "{:.4f}"), ("ff_rate", "{:.4f}"),
              ("det10", "{:.4f}"), ("det20", "{:.4f}"), ("det50", "{:.4f}"),
              ("det100", "{:.4f}"), ("mae50", "{:.2f}"), ("mae100", "{:.2f}")]


def _table(title, rows, cols, keycol="config"):
    head = [keycol] + [c for c, _ in cols]
    widths = [max(len(h), 8) for h in head]
    for r in rows:
        widths[0] = max(widths[0], len(str(r[keycol])))
    print(f"\n### {title}\n")
    print("| " + " | ".join(h.ljust(w) for h, w in zip(head, widths)) + " |")
    print("|" + "|".join("-" * (w + 2) for w in widths) + "|")
    for r in rows:
        cells = [str(r[keycol]).ljust(widths[0])]
        for (c, fmt), w in zip(cols, widths[1:]):
            v = r.get(c, float("nan"))
            cells.append((fmt.format(v) if v == v else "--").rjust(w))
        print("| " + " | ".join(cells) + " |")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="bench/out/improve")
    ap.add_argument("--configs", required=True)
    ap.add_argument("--source", default="lib")
    ap.add_argument("--per-d", action="store_true", help="also break the gold grid out by d")
    ap.add_argument("--only", default="", help="restrict to these table names: hom,indel,byp,byindel,gold")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir)
    names = [c for c in args.configs.split(",") if c.strip()]

    hom_rows, gold_rows, sub_rows, ind_rows, gd_rows, ir_rows = [], [], [], [], [], []
    for name in names:
        hp = outdir / f"cells_hom_{name}.json"
        if hp.exists():
            cells = hom_cells_from_json(hp)["by_p"]
            hom_rows.append({"config": name, **hom_summary(cells, args.source, "subs")})
            ind_rows.append({"config": name, **hom_summary(cells, args.source, "indel")})
            for p in (0.0, 0.10, 0.20, 0.30, 0.35):
                sub_rows.append({"config": f"{name} p={p:g}",
                                 **hom_summary(cells, args.source, "subs", p=p)})
            for ir in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05):
                ir_rows.append({"config": f"{name} indel={ir:g}",
                                **hom_summary(cells, args.source, "indel", indel=ir)})
        gp = outdir / f"cells_gold_{name}.json"
        if gp.exists():
            gc = gold_cells_from_json(gp)
            gold_rows.append({"config": name, **gold_summary(gc)})
            for d in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5):
                gd_rows.append({"config": f"{name} d={d:g}", **gold_summary(gc, d)})

    only = {x for x in args.only.split(",") if x.strip()}
    def want(k):
        return not only or k in only
    if hom_rows and want("hom"):
        _table(f"Homology grid, substitution panel, source={args.source} "
               f"(pooled over p and flank)", hom_rows, _HOM_COLS)
    if ind_rows and want("indel"):
        _table(f"Homology grid, INDEL panel, source={args.source}", ind_rows, _HOM_COLS)
    if sub_rows and want("byp"):
        _table("Homology grid, substitution panel, by mutation level",
               sub_rows, _HOM_COLS)
    if ir_rows and want("byindel"):
        _table("Homology grid, indel panel, by indel rate per site per branch",
               ir_rows, _HOM_COLS)
    if gold_rows and want("gold"):
        _table("Gold perturbed grid (pooled over d and flank)", gold_rows, _GOLD_COLS)
    if gd_rows and args.per_d:
        _table("Gold perturbed grid, by divergence", gd_rows, _GOLD_COLS)
    return 0


if __name__ == "__main__":
    sys.exit(main())


def roc_points(outdir, base, ts, grid="gold"):
    """(false-flank rate, detection at each flank size) for one config across a
    t_bits sweep -- the curve a change has to beat at MATCHED false-flank rate
    before it can be called an improvement rather than a re-tuning."""
    out = []
    for t in ts:
        tag = f"{base}_t{t:g}" if t is not None else base
        path = Path(outdir) / f"cells_{'gold' if grid == 'gold' else 'hom'}_{tag}.json"
        if not path.exists():
            continue
        if grid == "gold":
            row = gold_summary(gold_cells_from_json(path))
        else:
            row = hom_summary(hom_cells_from_json(path)["by_p"], "lib", "subs")
        row["t"] = t
        out.append(row)
    return sorted(out, key=lambda r: r["ff_rate"])
