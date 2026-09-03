"""Sweep the `--genome` parameters, reporting every statistic per detector source.

Why per source. A call set pooled from two structure-based detectors is not one
population. LTRharvest under `-mintsd 0 -maxtsd 0` with no `-motif` places its
boundaries by homology alone and has never looked at a terminal motif or a
target-site duplication; LTR_FINDER snaps its boundaries *onto* both. So a TSD
at an LTR_FINDER boundary restates that tool's own criterion, while a TSD at an
LTRharvest boundary is independent evidence. Pooling them mixes a circular
measurement with an informative one and the aggregate reverses each of them --
Simpson's paradox, and on the arabidopsis set a large one. Every rate here is
therefore reported per source and never summed.

Elements are assigned to a source by EXACT interval match against the stitched
SCN files, which is unambiguous: the FASTA header's `chrom:start-end` is the
detector's own reported interval, carried through the pipeline unchanged.

Controls are shift-budget matched. A wider shift window tries more (k, d5, d3)
combinations and so finds more duplications by chance; a control that tried
fewer would make every widening look like a gain. The control here runs the
identical search with both boundaries displaced outwards by `--control-offset`
bases, so it costs exactly as many trials as the measurement it calibrates.

One caveat on that control: it reaches `--control-offset` bases further into the
pad than the measurement does, so an element close enough to a contig edge that
its pad is shorter than `offset + max(k)` can score a duplication but never a
control, which inflates enrichment slightly. `PAD` is 32 and the offset is 10,
so this needs a pad under 17 bases: 68 of 10,307 arabidopsis records and 37 of
16,336 human ones, 0.7% and 0.2%. Too few to move a rate in the third decimal,
but they are there, and a set of short contigs would have more.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import Counter, defaultdict
from itertools import product

sys.path.insert(0, "src")

from kmer2ltr.fasta import read_fasta
from kmer2ltr.genome import PAD, PROBE, element_context, harvest, locus, orient

K_SETS = ((5,), (6, 5), (6, 5, 4), (7, 6, 5))
SHIFT_SETS = ((0,), (0, 1, -1), (0, 1, -1, 2, -2))
MIN_DISTINCT = (2, 3)
MISMATCHES = (0, 1)
# Which axis is searched first when a longer duplication and a smaller boundary
# shift both fit. Neither is obviously right, and on a set selected through a
# +-1 TSD window the choice is not cosmetic: it decides whether a shifted 6-mer
# or an unshifted 5-mer is the one reported.
ORDERS = ("k_major", "shift_major")
PROBES = (24, 40, 60)
IDENTITIES = (0.85, 0.90, 0.95, 1.0)


def scn_intervals(patterns) -> dict:
    """`{(chrom, start, end): source}` from `LABEL=GLOB` arguments.

    An interval both detectors reported is labelled `both`: those are the
    elements on which two independent methods agree exactly, and they behave
    like neither parent population.
    """
    seen: dict[tuple, set] = {}
    for spec in patterns:
        label, _, pattern = spec.partition("=")
        for path in sorted(glob.glob(pattern)):
            with open(path) as fh:
                for line in fh:
                    f = line.split()
                    if len(f) < 12 or not (f[0].isdigit() and f[1].isdigit()):
                        continue
                    s, e = int(f[0]), int(f[1])
                    if e < s:
                        s, e = e, s
                    # The chromosome is the LAST field: the LTR_FINDER files
                    # carry an extra chunk-name column before it.
                    seen.setdefault((f[-1], s, e), set()).add(label)
    return {k: (next(iter(v)) if len(v) == 1 else "both") for k, v in seen.items()}


def _search(ctx, b5, b3, ks, shifts, min_distinct, max_mismatch, order):
    """`genome.find_tsd` with its three refusals opened up as swept parameters.

    `max_mismatch` allows a decayed duplication to still count. It is swept
    rather than assumed: tolerating a mismatch multiplies the number of k-mer
    pairs that qualify by chance, so it has to pay for itself against the
    matched control before it earns a place.

    `order` decides which of a longer duplication and a smaller boundary shift
    wins when both fit. The first hit is returned either way, so this changes
    what is reported without changing whether anything is found at all.
    """
    if order == "k_major":
        trials = ((k, d5, d3) for k in ks for d5 in shifts for d3 in shifts)
    else:
        trials = ((k, d5, d3) for d5 in shifts for d3 in shifts for k in ks)
    for k, d5, d3 in trials:
        i, j = b5 + d5, b3 - d3
        if i - k < 0 or j + k > len(ctx) or i >= j:
            continue
        left, right = ctx[i - k:i], ctx[j:j + k]
        if "N" in left or "N" in right:
            continue
        if left != right:
            if not max_mismatch:
                continue
            if sum(x != y for x, y in zip(left, right)) > max_mismatch:
                continue
        if len(set(left)) < min_distinct:
            continue
        return left, d5, d3
    return None


def orientation_grid(elements, genome, out):
    """How the orientation call depends on how hard it is asked to agree.

    Reported against the strictest setting that decides at all, because there is
    no external truth here: what matters is that loosening the threshold buys
    decidability without ever flipping a call, and a flip is the failure mode
    that would silently corrupt a `--trim-flanks` header.
    """
    loci = [locus(sid) for sid, _ in read_fasta(elements)]
    rows = []
    for probe in PROBES:
        windows = harvest(genome, [x for x in loci if x], PAD, probe)
        calls: dict[float, dict[int, str]] = {t: {} for t in IDENTITIES}
        for i, (_sid, seq) in enumerate(read_fasta(elements)):
            win = windows.get(loci[i]) if loci[i] else None
            if win is None:
                continue
            for t in IDENTITIES:
                c = orient(seq, win, t)
                if c is not None:
                    calls[t][i] = c.orientation
        strict = calls[max(IDENTITIES)]
        for t in IDENTITIES:
            got = calls[t]
            flips = sum(1 for i, v in got.items() if i in strict and strict[i] != v)
            rows.append({"probe": probe, "min_identity": t, "decided": len(got),
                         "of": len(loci), "flips_vs_strictest": flips,
                         "minus": sum(v == "-" for v in got.values())})
            print(f"  probe={probe} id={t} decided={len(got)} flips={flips}",
                  file=sys.stderr)
    out["orientation"] = rows


def records(elements, genome, tsv, sources):
    """One tuple per record: source, context, sequence, and the called flanks."""
    import csv
    with open(tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    loci = [locus(sid) for sid, _ in read_fasta(elements)]
    windows = harvest(genome, [x for x in loci if x], PAD, PROBE)
    print(f"# {len(windows)} of {len(loci)} loci found in the reference",
          file=sys.stderr)
    for i, (sid, seq) in enumerate(read_fasta(elements)):
        row = rows[i]
        if row["seq_id"] != sid:
            sys.exit(f"row {i}: {row['seq_id']!r} != {sid!r}; the TSV is not "
                     f"this input's")
        key = loci[i]
        win = windows.get(key) if key else None
        ctx = orient(seq, win) if win is not None else None
        if ctx is None or row["status"] != "pass":
            continue
        yield (sources.get(key, "unmatched"), ctx, seq,
               int(row["flank5_len"]), int(row["flank3_len"]), row["motif"])


def report(path) -> None:
    """Print the swept cells as markdown, per source, never pooled.

    Enrichment is against the shift-budget-matched control, so cells with
    different shift windows are comparable to each other and not just to zero.
    """
    d = json.load(open(path))
    rows = d["orientation"]
    print(f"### Orientation (control offset {d['control_offset']} bp)\n")
    print("| probe | min identity | decided | of | reverse | flips vs strictest |")
    print("|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['probe']} | {r['min_identity']} | {r['decided']} | {r['of']} "
              f"| {r['minus']} | {r['flips_vs_strictest']} |")
    groups = sorted({g for c in d["cells"] for g in c["groups"]})
    for g in groups:
        print(f"\n### {g}\n")
        print("| k | shifts | distinct | mismatch | order | n | tsd@call | ctl | "
              "enrich | tsd@input | ctl | enrich |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for c in d["cells"]:
            v = c["groups"].get(g)
            if not v or not v["n"]:
                continue
            n = v["n"]
            def rate(k):
                return v[k] / n
            def enr(a, b):
                return f"{rate(a) / rate(b):.1f}x" if v[b] else "-"
            print(f"| {'.'.join(map(str, c['ks']))} "
                  f"| {','.join(map(str, c['shifts']))} | {c['min_distinct']} "
                  f"| {c['max_mismatch']} | {c['order'][0]} | {n} "
                  f"| {rate('tsd_called'):.4f} | {rate('ctl_called'):.4f} "
                  f"| {enr('tsd_called', 'ctl_called')} "
                  f"| {rate('tsd_input'):.4f} | {rate('ctl_input'):.4f} "
                  f"| {enr('tsd_input', 'ctl_input')} |")


def anchors(specs, sources) -> None:
    """Compare `--tsd-anchor` settings, per source, on a signal the flag cannot see.

    `tsd` at the called boundary is useless here: with the anchor on, a record
    carrying a duplication is snapped onto it *by construction*, so that column
    would report the flag's own action back as a success. `motif` is the check
    that stays honest -- nothing in Kmer2LTR reads a terminal dinucleotide at any
    setting -- so the question is whether suppressing a trim lands the boundary on
    `TG`..`CA` more often, and in which source.
    """
    import csv
    print("| anchor | source | n | flanked | mean flank bp | tg..ca | mean k2p |")
    print("|---|---|---|---|---|---|---|")
    for spec in specs:
        label, _, path = spec.partition("=")
        acc: dict[str, Counter] = defaultdict(Counter)
        k2p: dict[str, list] = defaultdict(list)
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["status"] != "pass" or r["orientation"] == "NA":
                    continue
                key = locus(r["seq_id"])
                src = sources.get(key, "unmatched")
                f5, f3 = int(r["flank5_len"]), int(r["flank3_len"])
                for g in (src, "POOLED"):
                    c = acc[g]
                    c["n"] += 1
                    c["flanked"] += bool(f5 or f3)
                    c["bp"] += f5 + f3
                    c["motif"] += r["motif"] == "tg...ca"
                    if r["k2p"] != "NA":
                        k2p[g].append(float(r["k2p"]))
        for g, c in sorted(acc.items()):
            n = c["n"]
            m = sum(k2p[g]) / len(k2p[g]) if k2p[g] else float("nan")
            print(f"| {label} | {g} | {n} | {c['flanked'] / n:.4f} "
                  f"| {c['bp'] / n:.1f} | {c['motif'] / n:.4f} | {m:.5f} |")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--report", metavar="JSON",
                    help="print an existing sweep as markdown and exit")
    ap.add_argument("--elements")
    ap.add_argument("--genome", nargs="+")
    ap.add_argument("--tsv", help="Kmer2LTR output for --elements")
    ap.add_argument("--scn", nargs="*", default=[], metavar="LABEL=GLOB",
                    help="stitched SCN files, e.g. harvest='*.work/*.ltrharvest*.scn'")
    ap.add_argument("--control-offset", type=int, default=10,
                    help="bases to displace both boundaries outwards for the null")
    ap.add_argument("--anchors", nargs="*", default=None, metavar="LABEL=TSV",
                    help="compare Kmer2LTR runs made at different --tsd-anchor "
                         "settings, per source, and exit")
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.report:
        report(a.report)
        return
    if a.anchors:
        anchors(a.anchors, scn_intervals(a.scn))
        return
    missing = [f for f in ("elements", "genome", "tsv", "out") if not getattr(a, f)]
    if missing:
        ap.error("required without --report: " + ", ".join("--" + m for m in missing))

    sources = scn_intervals(a.scn)
    print(f"# {len(sources)} SCN intervals", file=sys.stderr)
    cached = list(records(a.elements, a.genome, a.tsv, sources))
    print(f"# {len(cached)} classified records with genomic context", file=sys.stderr)

    d = a.control_offset
    out = []
    for ks, shifts, floor, mm, order in product(
            K_SETS, SHIFT_SETS, MIN_DISTINCT, MISMATCHES, ORDERS):
        acc: dict[str, Counter] = defaultdict(Counter)
        for src, context, seq, f5, f3, motif in cached:
            ctx, b5, b3 = element_context(seq, context)
            flanked = bool(f5 or f3)
            for group in (src, f"{src}|{'flanked' if flanked else 'clean'}"):
                c = acc[group]
                c["n"] += 1
                c["tsd_input"] += bool(
                    _search(ctx, b5, b3, ks, shifts, floor, mm, order))
                c["tsd_called"] += bool(
                    _search(ctx, b5 + f5, b3 - f3, ks, shifts, floor, mm, order))
                c["ctl_input"] += bool(
                    _search(ctx, b5 - d, b3 + d, ks, shifts, floor, mm, order))
                c["ctl_called"] += bool(
                    _search(ctx, b5 + f5 - d, b3 - f3 + d, ks, shifts, floor, mm, order))
                c["motif"] += motif == "tg...ca"
        out.append({"ks": list(ks), "shifts": list(shifts), "min_distinct": floor,
                    "max_mismatch": mm, "order": order,
                    "groups": {g: dict(c) for g, c in sorted(acc.items())}})
        print(f"  ks={ks} shifts={shifts} floor={floor} mm={mm} "
              f"order={order} done", file=sys.stderr)

    summary = {"control_offset": d, "cells": out}
    orientation_grid(a.elements, a.genome, summary)
    with open(a.out, "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"# wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
