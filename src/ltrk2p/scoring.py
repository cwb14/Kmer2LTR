"""Per-element calibrated scoring: K2P log-odds matrix, bit scores, significance."""
from __future__ import annotations

import math
from collections import Counter

import parasail

from .k2p import SubstCounts, k2p_distance

SCALE = 4          # integer units per bit (0.25-bit resolution)
ALPHABET = "ACGTN"
_BASES = "ACGT"
_PURINES = frozenset("AG")
_PYRIMIDINES = frozenset("CT")

# kappa floor of 1.0, not below: kappa < 1 inverts the ti/tv weighting (transversions
# would score above transitions), which for LTR pairs only ever arises from small-sample
# noise on a short core. kappa = 1 degenerates gracefully to equal weighting. This clamp
# affects ONLY the scoring matrix -- the reported K2P distance is computed from observed
# ts/tv counts and never uses kappa, so no reported distance is biased by it.
KAPPA_MIN, KAPPA_MAX, KAPPA_DEFAULT = 1.0, 10.0, 2.0
D_MIN, D_MAX = 0.01, 1.0   # calibration is clamped: the matrix must stay usable


def k2p_probs(d: float, kappa: float) -> tuple[float, float, float]:
    """(p_same, p_transition, p_each_transversion) under K2P at distance d.

    Parameterisation: with transition rate a and transversion rate b per
    direction, d = a + 2b and kappa = a/b, so b = d/(kappa+2), a = kappa*b.
    """
    b = d / (kappa + 2.0)
    a = kappa * b
    e1 = math.exp(-4.0 * b)
    e2 = math.exp(-2.0 * (a + b))
    p_same = 0.25 + 0.25 * e1 + 0.5 * e2
    p_ti = 0.25 + 0.25 * e1 - 0.5 * e2
    p_tv = 0.25 - 0.25 * e1
    return p_same, p_ti, p_tv


def logodds_bits(d: float, kappa: float, freqs: dict[str, float]) -> dict[tuple[str, str], float]:
    """s(x,y) = log2( P_xy(d, kappa) / (4 * f_x * f_y) ), in bits.

    Log-odds of homology against independence: the numerator is the joint
    probability of the pair under the K2P model (whose ancestral distribution
    is uniform, hence the factor 4), the denominator the joint probability
    under independence at the element's own composition.

    Dividing by f_y alone -- the conditional form -- is asymmetric whenever the
    composition is non-uniform, which makes the alignment score depend on which
    sequence is the query. Both f_x and f_y are needed.
    """
    p_same, p_ti, p_tv = k2p_probs(d, kappa)
    out: dict[tuple[str, str], float] = {}
    for x in _BASES:
        for y in _BASES:
            if x == y:
                p = p_same
            elif (x in _PURINES and y in _PURINES) or (x in _PYRIMIDINES and y in _PYRIMIDINES):
                p = p_ti
            else:
                p = p_tv
            fx = max(freqs.get(x, 0.25), 1e-6)
            fy = max(freqs.get(y, 0.25), 1e-6)
            out[(x, y)] = math.log2(p / (4.0 * fx * fy))
    return out


def _blank_matrix() -> "parasail.Matrix":
    return parasail.matrix_create(ALPHABET, 0, 0)


def parasail_matrix(bits_map: dict[tuple[str, str], float]) -> "parasail.Matrix":
    """Integer-scaled parasail matrix. N scores 0 against everything.

    set_value sets ONE direction only, so the N row and N column are both
    written explicitly -- omitting the second would leave N scoring as a
    mismatch in half of all comparisons.
    """
    m = _blank_matrix()
    for i, x in enumerate(_BASES):
        for j, y in enumerate(_BASES):
            m.set_value(i, j, int(round(bits_map[(x, y)] * SCALE)))
    n_idx = ALPHABET.index("N")
    for i in range(len(ALPHABET)):
        m.set_value(i, n_idx, 0)
        m.set_value(n_idx, i, 0)
    return m


def _generic() -> "parasail.Matrix":
    m = _blank_matrix()
    for i in range(4):
        for j in range(4):
            m.set_value(i, j, SCALE if i == j else -SCALE)
    for i in range(len(ALPHABET)):
        m.set_value(i, 4, 0)
        m.set_value(4, i, 0)
    return m


GENERIC_MATRIX = _generic()


def estimate_params(c: SubstCounts, seq: str) -> tuple[float, float, dict[str, float]]:
    """Estimate (d_hat, kappa_hat, base frequencies) from a core alignment."""
    counts = Counter(ch for ch in seq if ch in _BASES)
    total = sum(counts.values()) or 1
    freqs = {b: counts.get(b, 0) / total for b in _BASES}
    # under K2P there is one transition class and two transversion classes
    kappa = KAPPA_DEFAULT if c.n_tv == 0 else 2.0 * c.n_ts / c.n_tv
    kappa = min(max(kappa, KAPPA_MIN), KAPPA_MAX)
    d, _ = k2p_distance(c)
    if d is None:                       # saturated core: fall back to p-distance
        d = (c.n_ts + c.n_tv) / c.n_sites if c.n_sites else D_MAX
    d = min(max(d, D_MIN), D_MAX)
    return d, kappa, freqs


def bits(raw_score: int) -> float:
    return raw_score / SCALE


def evalue(bitscore: float, m: int, n: int, k_const: float = 0.1) -> float:
    """Karlin-Altschul form. Supplies the correct dependence on window size;
    the constant is calibrated empirically against negative controls.

    An overflowing exponent means the alignment is astronomically
    insignificant, so infinity is the semantically correct answer: every
    caller compares against MAX_EVALUE and rejects, which is the same
    decision the finite arithmetic would have reached. Without this guard a
    single pathological record aborts an entire batch run.
    """
    try:
        return k_const * m * n * (2.0 ** -bitscore)
    except OverflowError:
        return math.inf


def wfa_penalties(match: int, mismatch: int, open_p: int, ext_p: int) -> tuple[int, int, int]:
    """Convert a parasail max-score scheme to equivalent WFA min-penalties.

    Score = m*M - x'*X - gapcost, and M = (n1+n2-G)/2 - X, so maximising the
    score is equivalent to minimising (m+x')*X + gapcost + (m/2)*G. Hence
    x = m + x', e = e_p + m/2, o = o_p - e_p. All doubled to stay integral.
    """
    x = 2 * (match + mismatch)
    e = 2 * ext_p + match
    o = 2 * (open_p - ext_p)
    return x, o, e
