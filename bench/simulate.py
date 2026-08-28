"""Simulate diverged LTR pairs with known truth."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ltrk2p.k2p import count_substitutions, k2p_distance   # noqa: E402
from ltrk2p.scoring import k2p_probs                       # noqa: E402

_TI = {"A": "G", "G": "A", "C": "T", "T": "C"}
_TV = {"A": "CT", "G": "CT", "C": "AG", "T": "AG"}
_BASES = "ACGT"


def evolve(seq: str, d: float, kappa: float, indel_frac: float, rng) -> str:
    """Evolve `seq` a branch length d under K2P, plus indels.

    The descendant base is sampled from the EXACT K2P transition-probability
    matrix P(d, kappa), not from a single Bernoulli "did a substitution happen"
    draw. This matters: a single-draw model cannot produce multiple hits at one
    site, so the realised divergence between two such branches runs high
    relative to the nominal d that K2P's correction assumes, and the benchmark
    then reports a large phantom bias exactly in the high-divergence regime the
    tool exists to validate. Measured at nominal d = 0.70, two branches of d/2:
    single-draw model gives realised K2P 0.856; exact-matrix sampling gives
    0.697. Indels are applied at indel_frac of the substitution rate with
    geometric lengths (mean ~3 bp).
    """
    p_same, p_ti, p_tv = k2p_probs(d, kappa)
    p_indel = (1.0 - p_same) * indel_frac
    out: list[str] = []
    for ch in seq:
        if ch not in _BASES:
            out.append(ch)
            continue
        u = rng.random()
        if u < p_indel:
            if rng.random() < 0.5:
                continue                                    # deletion
            k = 1 + min(int(rng.expovariate(1 / 3.0)), 10)   # insertion
            out.append("".join(rng.choice(_BASES) for _ in range(k)))
            out.append(ch)
            continue
        # renormalise the substitution probabilities over the non-indel mass
        v = (u - p_indel) / (1.0 - p_indel) if p_indel < 1.0 else 0.0
        if v < p_same:
            out.append(ch)
        elif v < p_same + p_ti:
            out.append(_TI[ch])
        elif v < p_same + p_ti + p_tv:
            out.append(_TV[ch][0])
        else:
            out.append(_TV[ch][1])
    return "".join(out)


def shuffle_dinuc(seq: str, rng) -> str:
    """Composition-preserving shuffle. Flanks drawn this way are a hard null:
    they match the element's base composition, so flank detection cannot
    succeed merely by spotting a compositional shift."""
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def simulate_element(ltr: str, internal: str, d: float, kappa: float,
                     flank5: int, flank3: int, rng, indel_frac: float = 0.10):
    """Build one simulated element with exact truth.

    Both LTR copies are evolved independently from the ancestral consensus for
    d/2, so their pairwise distance is d.
    """
    l5 = evolve(ltr, d / 2, kappa, indel_frac, rng)
    l3 = evolve(ltr, d / 2, kappa, indel_frac, rng)
    pool = (ltr + internal) * 2
    f5 = shuffle_dinuc(pool[:max(flank5, 1)], rng)[:flank5] if flank5 else ""
    f3 = shuffle_dinuc(pool[-max(flank3, 1):], rng)[:flank3] if flank3 else ""
    seq = f5 + l5 + internal + l3 + f3
    realized, _ = k2p_distance(count_substitutions(
        l5[:min(len(l5), len(l3))], l3[:min(len(l5), len(l3))]))
    truth = {
        "ltr5_start": len(f5) + 1,
        "ltr5_end": len(f5) + len(l5),
        "ltr3_start": len(f5) + len(l5) + len(internal) + 1,
        "ltr3_end": len(f5) + len(l5) + len(internal) + len(l3),
        "ltr_len": len(ltr),
        "d_nominal": d,
        "kappa": kappa,
        "flank5": flank5,
        "flank3": flank3,
        "realized_k2p": realized,
    }
    return seq, truth
