"""Score `--period-rule` against a curated set of real elements.

The fixture directory holds three files: a FASTA of elements whose boundaries
the tool gets wrong, a length-matched FASTA of elements it gets right, and an
`expected.tsv` carrying the true boundaries. Only `expect_ltr5` / `expect_ltr3`
are read as truth. Any `kmer2ltr_current_call` column is IGNORED and recomputed,
because a baseline recorded by hand drifts from the code and then quietly
becomes the thing being measured.

The two sets are scored the same way and reported separately: a rule earns its
keep by fixing the first set without moving anything in the second.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from kmer2ltr.align import PERIOD_RULES, classify   # noqa: E402
from kmer2ltr.fasta import read_fasta               # noqa: E402

FASTAS = ("rebound_failures.fa", "controls_must_not_regress.fa")


def _span(text: str) -> tuple[int, int]:
    lo, hi = text.split("-")
    return int(lo), int(hi)


def load(fixture_dir: str | Path) -> tuple[dict[str, str], list[dict]]:
    d = Path(fixture_dir)
    seqs: dict[str, str] = {}
    for name in FASTAS:
        path = d / name
        if not path.exists():
            raise FileNotFoundError(f"fixture missing: {path}")
        seqs.update(dict(read_fasta(str(path))))
    with open(d / "expected.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    missing = [r["id"] for r in rows if r["id"] not in seqs]
    if missing:
        raise KeyError(f"{len(missing)} id(s) in expected.tsv absent from the "
                       f"FASTAs, first: {missing[0]}")
    return seqs, rows


def call(seq_id: str, seq: str, rule: str) -> tuple[int, int, int, int] | None:
    r = classify(seq_id, seq, period_rule=rule)
    if r.ltr5_start is None:
        return None
    return r.ltr5_start, r.ltr5_end, r.ltr3_start, r.ltr3_end


def score(fixture_dir, rules=PERIOD_RULES, tol: int = 10, verbose: bool = False):
    """Per-record calls under each rule, plus a right/wrong tally per set."""
    seqs, rows = load(fixture_dir)
    out = {"tolerance_bp": tol, "records": [], "summary": {}}
    for r in rows:
        exp = (*_span(r["expect_ltr5"]), *_span(r["expect_ltr3"]))
        rec = {"id": r["id"], "set": r["set"], "expected": exp, "calls": {}}
        for rule in rules:
            got = call(r["id"], seqs[r["id"]], rule)
            rec["calls"][rule] = {
                "call": got,
                "correct": got is not None and all(abs(a - b) <= tol
                                                   for a, b in zip(got, exp)),
            }
        out["records"].append(rec)
        if verbose:
            shown = "  ".join(f"{k}={v['call']}" for k, v in rec["calls"].items())
            print(f"  {r['id']:<26}{shown}", file=sys.stderr)
    for sname in sorted({r["set"] for r in rows}):
        sub = [x for x in out["records"] if x["set"] == sname]
        out["summary"][sname] = {
            "n": len(sub),
            **{rule: sum(x["calls"][rule]["correct"] for x in sub) for rule in rules},
        }
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("fixture_dir", help="directory holding the two FASTAs and expected.tsv")
    ap.add_argument("--tol", type=int, default=10,
                    help="bp a boundary may miss by and still count as correct "
                         "(default: %(default)s)")
    ap.add_argument("--json", help="write the full per-record table here")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-record calls")
    a = ap.parse_args(argv)

    res = score(a.fixture_dir, tol=a.tol, verbose=a.verbose)
    for sname, s in res["summary"].items():
        parts = "  ".join(f"{k}: {v}/{s['n']}" for k, v in s.items() if k != "n")
        print(f"{sname:<28}{parts}")
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1))
        print(f"wrote {a.json}", file=sys.stderr)
    # A rule that fixes nothing, or that moves a control, is a failed experiment
    # and must be visible to the shell, not just readable in the table.
    base, new = "best-score", "outermost"
    ctl = res["summary"].get("controls_must_not_regress")
    return 1 if ctl and ctl[new] < ctl[base] else 0


if __name__ == "__main__":
    raise SystemExit(main())
