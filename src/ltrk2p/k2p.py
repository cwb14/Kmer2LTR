"""Kimura 2-parameter distance from a pairwise alignment."""
from __future__ import annotations

import math
from dataclasses import dataclass

_PURINES = frozenset("AG")
_PYRIMIDINES = frozenset("CT")
_ACGT = frozenset("ACGT")


@dataclass(frozen=True)
class SubstCounts:
    n_sites: int      # ungapped, unambiguous columns (= n_match + n_ts + n_tv)
    n_match: int
    n_ts: int         # transitions: A<->G, C<->T
    n_tv: int         # transversions: everything else
    n_gapcols: int
    aln_len: int


def count_substitutions(a: str, b: str) -> SubstCounts:
    """Count substitution classes over two equal-length gapped aligned strings.

    Gap columns are counted only in n_gapcols. Columns containing N (or any
    non-ACGT residue) in either sequence are excluded entirely -- they inform
    neither the numerator nor the denominator of the distance.
    """
    if len(a) != len(b):
        raise ValueError(f"aligned strings differ in length: {len(a)} vs {len(b)}")
    n_match = n_ts = n_tv = n_gap = 0
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            n_gap += 1
            continue
        if x not in _ACGT or y not in _ACGT:
            continue
        if x == y:
            n_match += 1
        elif (x in _PURINES and y in _PURINES) or (x in _PYRIMIDINES and y in _PYRIMIDINES):
            n_ts += 1
        else:
            n_tv += 1
    return SubstCounts(n_sites=n_match + n_ts + n_tv, n_match=n_match, n_ts=n_ts,
                       n_tv=n_tv, n_gapcols=n_gap, aln_len=len(a))


def p_distance(c: SubstCounts) -> float | None:
    if c.n_sites == 0:
        return None
    return (c.n_ts + c.n_tv) / c.n_sites


def k2p_distance(c: SubstCounts) -> tuple[float | None, float | None]:
    """Kimura (1980) distance and its standard error.

    Returns (None, None) when the estimate is undefined -- no sites, or
    saturation making a logarithm argument non-positive. Never clamps.
    """
    if c.n_sites == 0:
        return None, None
    n = c.n_sites
    P = c.n_ts / n
    Q = c.n_tv / n
    w1 = 1.0 - 2.0 * P - Q
    w2 = 1.0 - 2.0 * Q
    if w1 <= 0.0 or w2 <= 0.0:
        return None, None
    d = -0.5 * math.log(w1) - 0.25 * math.log(w2)
    # Kimura (1980) variance: a = 1/w1, b = (1/w1 + 1/w2)/2
    a = 1.0 / w1
    b = 0.5 * (1.0 / w1 + 1.0 / w2)
    var = (a * a * P + b * b * Q - (a * P + b * Q) ** 2) / n
    return d, math.sqrt(var) if var > 0 else 0.0
