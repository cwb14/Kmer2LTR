"""Is the refactored pipeline, run at the shipped settings, the same tool?

Compares a fresh `baseline` prediction against an existing prediction produced
by the pre-campaign code at the same settings, record by record. The refactor
changed two things that were meant to be behaviour-neutral -- Stage 3's
extension partner regions now come from the internal region rather than the
discovery window, and the second discovery pass is seeded at the window the
first one settled on -- so "neutral" is a claim that has to be checked, not
asserted.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter

FIELDS = ("status", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--examples", type=int, default=5)
    args = ap.parse_args(argv)

    n = same = 0
    diff_field = Counter()
    shifts = Counter()
    examples = []
    with open(args.old, newline="") as of, open(args.new, newline="") as nf:
        for o, w in zip(csv.DictReader(of, delimiter="\t"),
                        csv.DictReader(nf, delimiter="\t")):
            n += 1
            if o["seq_id"] != w["seq_id"]:
                print(f"ERROR: record order diverged at row {n}", file=sys.stderr)
                return 1
            d = [f for f in FIELDS if o[f] != w[f]]
            if not d:
                same += 1
                continue
            for f in d:
                diff_field[f] += 1
            for f in d:
                if f != "status" and "NA" not in (o[f], w[f]):
                    shifts[min(abs(int(o[f]) - int(w[f])), 100)] += 1
            if len(examples) < args.examples:
                examples.append((w["seq_id"], {f: (o[f], w[f]) for f in d}))

    print(f"records compared : {n}")
    print(f"identical        : {same} ({same / n:.4%})")
    print(f"differing        : {n - same} ({(n - same) / n:.4%})")
    print(f"fields differing : {dict(diff_field)}")
    if shifts:
        tot = sum(shifts.values())
        cum = 0
        print("coordinate shift distribution (bp, capped at 100):")
        for k in sorted(shifts):
            cum += shifts[k]
            print(f"  {k:>4} bp : {shifts[k]:>7}  (cum {cum / tot:.4%})")
    for sid, fields in examples:
        print(f"  e.g. {sid}: {fields}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
