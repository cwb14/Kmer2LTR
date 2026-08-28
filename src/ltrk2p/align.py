"""Stages 1-5: locate the terminal repeat pair and align it."""
from __future__ import annotations

from dataclasses import dataclass

import parasail

from .scoring import SCALE, GENERIC_MATRIX, bits, evalue

GAP_OPEN = 6 * SCALE      # 6 bits; parasail convention: gap of k costs open + (k-1)*ext
GAP_EXTEND = 2 * SCALE    # 2 bits
MIN_LEN = 100             # shorter than this cannot hold two LTRs plus internal
W0 = 1500
EDGE = 20                 # touching within this many bp of the inner edge -> grow
MAX_EVALUE = 1e-3         # significance floor; Task 8 reuses this same constant


@dataclass(frozen=True)
class Hit:
    score: int
    qb: int   # 0-based inclusive, index into prefix window (starts at 0)
    qe: int
    rb: int   # 0-based inclusive, index into suffix window (starts at len(S)-w)
    re: int
    w: int


def discover(S: str, matrix, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTEND,
             w0: int = W0, edge: int = EDGE) -> Hit | None:
    """Locate the terminal repeat pair by adaptive-window local alignment.

    Score-only passes throughout: traceback costs ~7x more and is not needed
    until boundaries are settled. Two passes -- forward for the inner ends,
    then reversed-prefix for the outer ends.
    """
    L = len(S)
    if L < MIN_LEN:
        return None
    w_max = L // 2
    w = min(w_max, w0)
    while True:
        P, Q = S[:w], S[L - w:]
        fwd = parasail.sw_striped_sat(P, Q, gap_open, gap_extend, matrix)
        if fwd.score <= 0:
            if w < w_max:
                w = min(w * 2, w_max)
                continue
            return None
        qe, re_ = fwd.end_query, fwd.end_ref
        rev = parasail.sw_striped_sat(P[:qe + 1][::-1], Q[:re_ + 1][::-1],
                                      gap_open, gap_extend, matrix)
        qb = qe - rev.end_query
        rb = re_ - rev.end_ref
        touches_inner = (qe >= w - edge) or (rb <= edge)
        # A window too small to span the LTR (ltr_len >= 2w) leaves the two windows
        # covering DISJOINT parts of the LTR, so the only hit available is background
        # noise -- which has no reason to touch an edge. Insignificance must therefore
        # trigger growth too, or large-LTR elements are silently mis-called.
        weak = evalue(bits(fwd.score), w, w) > MAX_EVALUE
        if (touches_inner or weak) and w < w_max:
            w = min(w * 2, w_max)
            continue
        if weak:
            # At the ceiling with an insignificant hit: there is no terminal repeat here.
            # Returning it anyway would hand callers noise indistinguishable from a real
            # pair, since Hit carries no significance field.
            return None
        return Hit(score=fwd.score, qb=qb, qe=qe, rb=rb, re=re_, w=w)


from .k2p import count_substitutions
from .scoring import estimate_params, logodds_bits, parasail_matrix

MIN_CALIB_SITES = 50


def ltr_spans(S: str, hit: Hit) -> tuple[int, int, int, int]:
    """Absolute 0-based (l5b, l5e, l3b, l3e); ends inclusive."""
    wstart = len(S) - hit.w
    return hit.qb, hit.qe - 1, wstart + hit.rb, wstart + hit.re - 1


def _parasail_aligned_pair(res, query: str, ref: str) -> tuple[str, str]:
    """Expand a parasail traceback into gapped aligned strings."""
    tb = res.traceback
    return tb.query, tb.ref


def core_alignment(S: str, hit: Hit) -> tuple[str, str]:
    l5b, l5e, l3b, l3e = ltr_spans(S, hit)
    q, r = S[l5b:l5e + 1], S[l3b:l3e + 1]
    res = parasail.nw_trace_striped_sat(q, r, GAP_OPEN, GAP_EXTEND, GENERIC_MATRIX)
    return _parasail_aligned_pair(res, q, r)


def calibrate(S: str, hit: Hit):
    """Derive this element's own scoring matrix from its own divergence."""
    a, b = core_alignment(S, hit)
    counts = count_substitutions(a, b)
    d_hat, kappa_hat, freqs = estimate_params(counts, S)
    if counts.n_sites < MIN_CALIB_SITES:
        return GENERIC_MATRIX, d_hat, kappa_hat
    return parasail_matrix(logodds_bits(d_hat, kappa_hat, freqs)), d_hat, kappa_hat
