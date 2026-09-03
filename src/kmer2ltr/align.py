"""Stages 1-5: locate the terminal repeat pair and align it.

Division of labour, extended from the spec's matrix rule to gap penalties:
**the calibrated model (substitution matrix AND gap penalties) determines the
alignment; the generic model (`GENERIC_MATRIX` + `SIG_GAPS`) determines
significance, everywhere.** Pinning the significance side is what lets the
alignment gap penalties change (`gap_scheme=`) without invalidating
`MAX_EVALUE`, which was calibrated against negative controls under exactly the
generic model and nothing else.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, fields

import numpy as np
import parasail

from .scoring import SCALE, GENERIC_MATRIX, bits, evalue

# --------------------------------------------------------------------------- #
# The significance model: fixed forever, so MAX_EVALUE keeps its meaning.
# --------------------------------------------------------------------------- #
SIG_GAP_OPEN = 6 * SCALE      # 6 bits; parasail convention: gap of k costs open + (k-1)*ext
SIG_GAP_EXTEND = 2 * SCALE    # 2 bits
MIN_LEN = 100                 # shorter than this cannot hold two LTRs plus internal
W0 = 1500
EDGE = 20                     # touching within this many bp of the inner edge -> grow
MAX_EVALUE = 1e-10            # significance floor, shared by every gate in the tool
# Calibrated against 46,823 real non-LTR TEs plus a composition-matched shuffled
# null. The shuffled null is clean at every threshold tested, so the residual
# false-"pass" rate is real terminal-repeat structure inside specific subclasses
# -- satellites definitionally, plus some library-consensus entries that are
# near-exact tandem constructions -- not chance alignment noise. No threshold
# separates those from a true LTR pair without the family classification this
# tool declines to do. 1e-10 is the last point before the true-positive cost
# curve steepens: it nearly halves the false-"pass" rate (9.09% -> 5.23%) while
# real-data pass rate moves 99.80% -> 99.67%. Full sweep: docs/benchmarks.md.


@dataclass(frozen=True)
class Gaps:
    """Affine gap penalties in integer SCALE units."""
    open: int = SIG_GAP_OPEN
    extend: int = SIG_GAP_EXTEND


SIG_GAPS = Gaps(SIG_GAP_OPEN, SIG_GAP_EXTEND)

# Retained as module constants because the spec, the plan and several tests
# refer to them by name; they are the "legacy" alignment scheme's values.
GAP_OPEN = SIG_GAP_OPEN
GAP_EXTEND = SIG_GAP_EXTEND

# Static alignment gap schemes, in BITS. "legacy" is what every published
# Kmer2LTR measurement to date was made under. "static" raises the opening cost
# and drops the extension cost, which is the shape an affine approximation to a
# geometric indel-length model actually takes: opening is rare (~1e-3/site ->
# ~10 bits) while continuing is common (mean length ~3 -> ~0.6 bits).
GAP_SCHEMES: dict[str, tuple[float, float]] = {
    "legacy": (6.0, 2.0),
    "static": (10.0, 1.0),
}
# Clamps for `gap_scheme="adaptive"`. The open range brackets indel rates from
# ~6e-5 to ~6e-2 per site; the extend range brackets mean indel lengths from
# ~1.2 bp to infinite. Outside those the estimate is small-sample noise on a
# short core, not signal.
GAP_OPEN_BITS = (4.0, 16.0)
GAP_EXTEND_BITS = (0.25, 3.0)


@dataclass(frozen=True)
class Hit:
    score: int
    qb: int   # 0-based inclusive, index into prefix window (starts at 0)
    qe: int
    rb: int   # 0-based inclusive, index into suffix window (starts at len(S)-w)
    re: int
    w: int


def discover(S: str, matrix, gap_open: int | Gaps = SIG_GAP_OPEN,
             gap_extend: int = SIG_GAP_EXTEND, w0: int = W0, edge: int = EDGE,
             max_evalue: float = MAX_EVALUE) -> Hit | None:
    """Locate the terminal repeat pair by adaptive-window local alignment.

    Score-only passes throughout: traceback costs ~7x more and is not needed
    until boundaries are settled. Two passes -- forward for the inner ends,
    then reversed-prefix for the outer ends.

    `gap_open` accepts a `Gaps` for callers that carry one; the two-int form is
    kept because the plan, the tests and `bench/` all call it that way.
    """
    if isinstance(gap_open, Gaps):
        gap_open, gap_extend = gap_open.open, gap_open.extend
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
        weak = evalue(bits(fwd.score), w, w) > max_evalue
        if (touches_inner or weak) and w < w_max:
            w = min(w * 2, w_max)
            continue
        if weak:
            # At the ceiling with an insignificant hit: there is no terminal repeat here.
            # Returning it anyway would hand callers noise indistinguishable from a real
            # pair, since Hit carries no significance field.
            return None
        return Hit(score=fwd.score, qb=qb, qe=qe, rb=rb, re=re_, w=w)


from .k2p import (count_gap_runs, count_substitutions, insertion_time,
                  k2p_distance, p_distance)
from .scoring import (estimate_params, expected_random_bits, generic_alpha,
                      logodds_bits, parasail_matrix, wfa_penalties)
from .cigar import aligned_pair_from_cigartuples, cs_string, extended_cigar
from pywfa import WavefrontAligner

MIN_CALIB_SITES = 50


def ltr_spans(S: str, hit: Hit) -> tuple[int, int, int, int]:
    """Absolute 0-based (l5b, l5e, l3b, l3e); ends inclusive.

    parasail's end_query/end_ref are 0-based INCLUSIVE, so they are used as-is.
    Do not subtract 1: that truncates the final base of both LTRs.
    """
    wstart = len(S) - hit.w
    return hit.qb, hit.qe, wstart + hit.rb, wstart + hit.re


def _parasail_aligned_pair(res, query: str, ref: str) -> tuple[str, str]:
    """Expand a parasail traceback into gapped aligned strings."""
    tb = res.traceback
    return tb.query, tb.ref


def core_alignment(S: str, hit: Hit) -> tuple[str, str]:
    l5b, l5e, l3b, l3e = ltr_spans(S, hit)
    return _core_alignment_spans(S, (l5b, l5e, l3b, l3e))


def _core_alignment_spans(S: str, spans) -> tuple[str, str]:
    """Bootstrap alignment of the core pair, always under the generic model.

    This is the alignment Stage 2 estimates FROM, so it cannot use the matrix
    or gaps Stage 2 is about to derive.
    """
    l5b, l5e, l3b, l3e = spans
    q, r = S[l5b:l5e + 1], S[l3b:l3e + 1]
    res = parasail.nw_trace_striped_sat(q, r, SIG_GAP_OPEN, SIG_GAP_EXTEND, GENERIC_MATRIX)
    return _parasail_aligned_pair(res, q, r)


def gaps_from_alignment(a: str, b: str, counts) -> Gaps:
    """Affine penalties implied by the indel statistics of a core alignment.

    Under a geometric indel-length model the log-odds cost of a gap of length k
    is `-log2(mu) - (k-1)*log2(P_continue)`, which is exactly parasail's affine
    form -- so `open = -log2(mu)` and `extend = -log2(P_continue)` with `mu` the
    per-column probability that a gap opens and `P_continue = 1 - 1/mean_len`.
    Both are clamped: a core with no gaps at all would otherwise imply an
    infinite opening cost, and a core with one long gap an infinite mean.

    This is the indel-rate estimate the spec's Stage 2 promises and the shipped
    code never used.
    """
    n_gaps = count_gap_runs(a, b)
    denom = counts.n_sites + counts.n_gapcols
    mu = n_gaps / denom if denom else 0.0
    open_bits = -math.log2(mu) if mu > 0 else GAP_OPEN_BITS[1]
    mean_len = counts.n_gapcols / n_gaps if n_gaps else 1.0
    cont = 1.0 - 1.0 / mean_len if mean_len > 1.0 else 0.0
    ext_bits = -math.log2(cont) if cont > 0 else GAP_EXTEND_BITS[1]
    open_bits = min(max(open_bits, GAP_OPEN_BITS[0]), GAP_OPEN_BITS[1])
    ext_bits = min(max(ext_bits, GAP_EXTEND_BITS[0]), GAP_EXTEND_BITS[1])
    return Gaps(int(round(open_bits * SCALE)), int(round(ext_bits * SCALE)))


def gaps_for_scheme(scheme: str, a: str | None = None, b: str | None = None,
                    counts=None) -> Gaps:
    if scheme == "adaptive":
        if a is None:
            return Gaps(*(int(round(v * SCALE)) for v in GAP_SCHEMES["legacy"]))
        return gaps_from_alignment(a, b, counts)
    if scheme not in GAP_SCHEMES:
        raise ValueError(f"unknown gap_scheme {scheme!r}; "
                         f"choices: {sorted(GAP_SCHEMES) + ['adaptive']}")
    o, e = GAP_SCHEMES[scheme]
    return Gaps(int(round(o * SCALE)), int(round(e * SCALE)))


@dataclass(frozen=True)
class Calib:
    matrix: object
    d_hat: float
    kappa_hat: float
    gaps: Gaps
    alpha: float          # bits lost per non-homologous aligned base


def calibrate_full(S: str, spans, *, comp: str = "element",
                   gap_scheme: str = "legacy") -> Calib:
    """Derive this element's own alignment model from its own core pair.

    `comp` selects the composition the log-odds denominator is built from:
    "element" (the whole record, the shipped behaviour) or "core" (the two LTR
    copies only). The two differ whenever the internal region's composition
    differs from the LTRs', which for a long internal region is most of the
    record.
    """
    a, b = _core_alignment_spans(S, spans)
    counts = count_substitutions(a, b)
    l5b, l5e, l3b, l3e = spans
    src = (S[l5b:l5e + 1] + S[l3b:l3e + 1]) if comp == "core" else S
    d_hat, kappa_hat, freqs = estimate_params(counts, src)
    if counts.n_sites < MIN_CALIB_SITES:
        # Degenerate core: keep the generic model whole rather than manufacturing
        # a degenerate one from it. That includes the gap penalties -- an indel
        # rate estimated from under 50 ungapped columns is noise.
        return Calib(GENERIC_MATRIX, d_hat, kappa_hat, gaps_for_scheme("legacy"),
                     generic_alpha())
    lo = logodds_bits(d_hat, kappa_hat, freqs)
    return Calib(parasail_matrix(lo), d_hat, kappa_hat,
                 gaps_for_scheme(gap_scheme, a, b, counts),
                 expected_random_bits(lo, freqs))


def calibrate(S: str, hit, comp: str = "element"):
    """Backwards-compatible 3-tuple form: (matrix, d_hat, kappa_hat).

    Accepts a `Hit` or an already-resolved `(l5b, l5e, l3b, l3e)` tuple.
    """
    spans = ltr_spans(S, hit) if isinstance(hit, Hit) else tuple(hit)
    c = calibrate_full(S, spans, comp=comp)
    return c.matrix, c.d_hat, c.kappa_hat


# --------------------------------------------------------------------------- #
# Stage 3 -- terminal boundary model selection
# --------------------------------------------------------------------------- #

T_BITS = 10.0    # fixed fallback when no d_hat is available; benchmark-calibrated
# Swept over {2,5,8,10,15,20,30} on 260,876 perturbed real elements. Raising it
# from 5 to 10 cuts the false-flank rate ~3x at d=0.3 while large-flank detection
# barely moves; the cost falls on 10-20 bp flanks, which sit near the detection
# floor `effective_t_bits` describes. Sweep table: docs/benchmarks.md.

# Divergence-aware schedule: (exclusive upper bound on d_hat, t_bits). A single
# constant fits badly -- at fixed t_bits the false-flank rate swings 40-80x
# across the observed d_hat range while the threshold does not move. Derived by
# bench/calibrate_flank_threshold.py under a rule fixed before the sweep it
# consumes is run, against an explicit large-flank detection floor.
#
# It RELAXES the threshold at low divergence rather than tightening it at high,
# which is the opposite of the obvious guess: the floor rules out t=15 and t=20
# in every bin. See docs/design.md, Stage 3, for the trade it buys.
T_BITS_SCHEDULE: tuple[tuple[float, float], ...] = (
    (0.025, 2),
    (0.15, 8),
    (float("inf"), 10),
)


def t_bits_for(d_hat: float | None) -> float:
    """Flank-evidence threshold for an element of estimated divergence `d_hat`."""
    if d_hat is None:
        return T_BITS
    for hi, t in T_BITS_SCHEDULE:
        if d_hat < hi:
            return t
    return T_BITS_SCHEDULE[-1][1]


def effective_t_bits(t_bits: float, k: int, alpha: float | None,
                     beta: float | None) -> float:
    """Cap the demanded evidence at what a k-base flank can physically supply.

    The snap test accepts homology to the terminus when the extension costs less
    than `t_bits`. A non-homologous flank accumulates cost at `alpha` bits per
    base (`scoring.expected_random_bits`), so a flank of length k can never
    present more than `k * alpha` bits of evidence against the snap. Whenever
    `k * alpha < t_bits` the test cannot fire whatever the sequence says: the
    boundary is decided before it is looked at. Since `alpha` falls with
    divergence -- 4.2 bits/base at d=0.01 down to 0.73 at d=0.35 -- that blind
    spot is a few bases wide on a young element and over a dozen on an old one,
    which is exactly the shape of the measured detection curve.

    `beta` sets where inside the available evidence the decision sits. A
    homologous continuation scores about 0; a non-homologous one about
    `-k * alpha`; so `beta = 0.5` is the maximum-likelihood split between the two
    hypotheses for a segment of that length, and `beta = 1.0` tests against the
    flank hypothesis's own expectation. `None` restores a flat threshold.

    Below the floor no threshold recovers certainty -- a 1 bp flank carries one
    base of evidence and is near-chance however it is tested -- so this raises
    detection to the information-theoretic ceiling, not past it.
    """
    if beta is None or alpha is None or alpha <= 0.0 or k <= 0:
        return t_bits
    return min(t_bits, beta * alpha * k)


# How much evidence a flank must present, as a fraction of the most it could.
# Lower is more willing to call a flank. `None` disables the cap entirely, which
# is the behaviour every measurement before 2026-08-31 was made under.
FLANK_SENSITIVITY: dict[str, float | None] = {
    "strict": None,
    "balanced": 1.0,
    "sensitive": 0.5,
}

EXT_CAP = 2000       # max flank bases scored in one graded extension (memory bound)
EXT_REF_SLACK = 200  # ref allowance beyond 2x the query length


@dataclass(frozen=True)
class Bounds:
    l5b: int
    l5e: int
    l3b: int
    l3e: int
    margin_bits: float | None

    @property
    def spans(self) -> tuple[int, int, int, int]:
        return self.l5b, self.l5e, self.l3b, self.l3e


def _extend(outer: str, inner: str, matrix, gaps: Gaps, t_bits: float,
            graded: bool, alpha: float | None = None,
            beta: float | None = None,
            credit: float = 0.0) -> tuple[int, int, float]:
    """One anchored outward extension of an alignment end.

    `outer` is the candidate flank, consumed from the core outward; `inner` is
    its homologous partner region, whose far end is free. Both are oriented so
    that index 0 is adjacent to the core, which is what makes "begins anchored"
    mean "next to the core".

    Returns `(k_outer, k_inner, margin_bits)`: how many bases of each side the
    boundary moves, and how decisively the call was made. `k_outer ==
    len(outer)` means the boundary snapped all the way to the sequence terminus.

    `credit` is evidence from outside the sequence that this terminus is already
    the element's boundary, in bits. It raises what the extension is allowed to
    cost and can therefore only ever move the boundary further out; in binary
    mode "further out" means the sequence terminus and nothing else.

    Binary mode is the shipped rule: consume the whole flank iff doing so costs
    less than `t_bits` (under the penalised objective `score - T*(free ends)`,
    that is `s_full > -T`). Graded mode instead takes the furthest endpoint
    whose running score is still above the `-t_bits` noise floor, so a flank
    call becomes a position rather than a verdict.
    """
    if not outer or not inner:
        return 0, 0, None
    # The cap describes how much evidence the flank itself could ever supply;
    # `credit` is evidence from outside the sequence, so it is added after.
    t_bits = effective_t_bits(t_bits, len(outer), alpha, beta) + credit
    n_ref = min(len(inner), 2 * len(outer) + EXT_REF_SLACK)
    inner = inner[:n_ref]
    res = parasail.sg_de_striped_sat(outer, inner, gaps.open, gaps.extend, matrix)
    s_full = bits(res.score)
    margin = abs(s_full + t_bits)
    if not graded:
        if s_full > -t_bits:
            return len(outer), res.end_ref + 1, margin
        return 0, 0, margin
    if s_full > -t_bits:
        return len(outer), res.end_ref + 1, margin
    # Not homologous all the way out. Find the furthest endpoint still above the
    # noise floor. H(k) = best score consuming exactly k outer bases; H(0) = 0.
    #
    # This is the one call that materialises a full DP table, so both sides are
    # capped: the table is len(q) x len(r) int32, and an uncapped pathological
    # record (a 17 kb flank candidate against a 34 kb partner) would allocate
    # ~2 GB per worker. Capping the partner by the CAPPED query length, not by
    # the original outer length, is what keeps that bound real.
    q = outer[:EXT_CAP]
    r = inner[:2 * len(q) + EXT_REF_SLACK]
    # `res` MUST stay referenced for as long as `tab` is read. parasail's
    # `score_table` is a numpy view onto memory owned by the result object, and
    # it does not keep that object alive -- so
    #     tab = np.asarray(parasail.sg_de_table_...(...).score_table)
    # is a use-after-free. It fails silently for small tables (the freed pages
    # are still mapped, so the read returns plausible-looking garbage) and
    # segfaults once the table is big enough for the allocator to hand the pages
    # back: measured at q=2000 x r=4200 on
    # arabidopsis__LR999453.1:20940665-20952438#LTR/Gypsy/Retand__d0.5__f500.
    res = parasail.sg_de_table_striped_sat(q, r, gaps.open, gaps.extend, matrix)
    tab = np.asarray(res.score_table)
    H = tab.max(axis=1)
    ok = np.nonzero(H > -t_bits * SCALE)[0]
    if not len(ok):
        return 0, 0, margin
    k = int(ok[-1])                      # 0-based row -> k = k+1 outer bases
    return k + 1, int(tab[k].argmax()) + 1, margin


def _joint_inner(S: str, b: Bounds, matrix, gaps: Gaps) -> Bounds:
    """Re-derive both inner boundaries from ONE alignment anchored at the outer ones.

    In the shipped pipeline each inner boundary is a by-product of the opposite
    terminus's snap (5' outer and 3' inner are the same alignment end; so are 3'
    outer and 5' inner), and nothing re-examines it once the outer end is fixed.
    Here the two settled outer boundaries are anchored and the two inner ones
    are the free ends of a single semi-global alignment, so they are chosen
    jointly and optimally rather than inherited.
    """
    L = len(S)
    span = b.l3e + 1 - b.l5b
    if span < 4:
        return b
    cap = 2 * max(b.l5e - b.l5b + 1, b.l3e - b.l3b + 1) + 500
    mid = b.l5b + span // 2
    q = S[b.l5b:min(mid, b.l5b + cap)]
    r = S[max(mid, b.l3e + 1 - cap):b.l3e + 1]
    if len(q) < 2 or len(r) < 2:
        return b
    # query begin anchored (5' outer), ref end anchored (3' outer);
    # query end free and ref begin free are the two inner boundaries.
    fwd = parasail.sg_qe_db_striped_sat(q, r, gaps.open, gaps.extend, matrix)
    if fwd.score <= 0:
        return b
    eq = fwd.end_query
    # Reversed, the same alignment has both begins anchored and only the ref end
    # free -- i.e. sg_de -- so end_ref there locates the ref begin here.
    rev = parasail.sg_de_striped_sat(q[:eq + 1][::-1], r[::-1],
                                     gaps.open, gaps.extend, matrix)
    rb = len(r) - 1 - rev.end_ref
    l5e = b.l5b + eq
    l3b = (b.l3e + 1 - len(r)) + rb
    if not (b.l5b <= l5e < l3b <= b.l3e):
        return b
    return Bounds(b.l5b, l5e, l3b, b.l3e, b.margin_bits)


def snap_bounds(S: str, spans, matrix, gaps: Gaps = SIG_GAPS, t_bits: float = T_BITS,
                *, mode: str = "binary", inner: str = "none",
                alpha: float | None = None, beta: float | None = None,
                credit: float = 0.0) -> Bounds:
    """Decide, per terminus, whether homology reaches the end of the sequence.

    Exact model comparison rather than a greedy endpoint: extending the core to
    the terminus is accepted iff the extension costs less than `t_bits`.

    The two extension partner regions are taken from the INTERNAL region
    (capped), not from the discovery window. That removes `w` from this stage --
    which is what makes it re-runnable after Stage 4, where no window exists --
    and closes two degeneracies of the window-bounded form: an empty partner
    region whenever the hit began exactly at the suffix-window edge, and a
    partner region overlapping the 5' LTR itself whenever the window started
    before it.
    """
    L = len(S)
    l5b, l5e, l3b, l3e = spans
    margins: list[float] = []
    graded = mode == "graded"
    if mode not in ("binary", "graded"):
        raise ValueError(f"unknown snap mode {mode!r}; choices: 'binary', 'graded'")

    # End A: the 5' outer boundary and the 3' inner boundary are the same end of
    # the alignment. Reversed so index 0 is adjacent to the core.
    if l5b > 0 and l3b > l5e + 1:
        k_out, k_in, m = _extend(S[:l5b][::-1], S[l5e + 1:l3b][::-1],
                                 matrix, gaps, t_bits, graded, alpha, beta, credit)
        if m is not None:
            margins.append(m)
        l5b -= k_out
        l3b -= k_in

    # End B: the 3' outer boundary and the 5' inner boundary. Already oriented
    # outward from the core, so no reversal.
    if l3e < L - 1 and l3b > l5e + 1:
        k_out, k_in, m = _extend(S[l3e + 1:], S[l5e + 1:l3b],
                                 matrix, gaps, t_bits, graded, alpha, beta, credit)
        if m is not None:
            margins.append(m)
        l3e += k_out
        l5e += k_in

    l5b = max(0, l5b)
    l3e = min(L - 1, l3e)
    l5e = min(l5e, l3b - 1)
    b = Bounds(l5b, l5e, l3b, l3e, min(margins) if margins else None)
    if inner == "joint":
        b = _joint_inner(S, b, matrix, gaps)
    elif inner != "none":
        raise ValueError(f"unknown inner mode {inner!r}; choices: 'none', 'joint'")
    return b


def terminal_snap(S: str, hit: Hit, matrix, t_bits: float = T_BITS,
                  gaps: Gaps = SIG_GAPS, *, mode: str = "binary",
                  inner: str = "none", alpha: float | None = None,
                  beta: float | None = None, credit: float = 0.0) -> Bounds:
    """`snap_bounds` for a caller holding a `Hit` (the Stage 1 entry point)."""
    return snap_bounds(S, ltr_spans(S, hit), matrix, gaps, t_bits,
                       mode=mode, inner=inner, alpha=alpha, beta=beta,
                       credit=credit)


# --------------------------------------------------------------------------- #
# Stage 4 -- outermost pair
# --------------------------------------------------------------------------- #

MIN_OUTER = 50      # shorter outer segments cannot carry a significant pair


def outermost(S: str, bounds: Bounds, matrix=None, max_evalue: float = MAX_EVALUE) -> Bounds:
    """If a flank was called, prefer a significant repeat pair strictly outside it.

    This is what makes a retained (un-excised) nested element report the OUTER
    element rather than the younger, higher-scoring nested one. Returns the
    input `bounds` unchanged when nothing outside qualifies.

    `matrix` is accepted and ignored: this search is deliberately scored with
    the generic model (see below), and the parameter exists only so callers
    holding the element's calibrated matrix can pass it without special-casing.
    """
    L = len(S)
    if bounds.l5b == 0 and bounds.l3e == L - 1:
        return bounds
    outer5 = S[:bounds.l5b]
    outer3 = S[bounds.l3e + 1:]
    if len(outer5) < MIN_OUTER or len(outer3) < MIN_OUTER:
        return bounds
    # Significance is scored with GENERIC_MATRIX, not the per-element calibrated
    # matrix. MAX_EVALUE was calibrated against GENERIC's score scale, and the
    # calibrated matrix is on a different one -- at low divergence it scores a
    # match +8 where GENERIC scores +4, so the same alignment reports roughly
    # double the bits and an E-value ~2^13 too small. Mixing the two leaked
    # spurious "outer pairs" out of pure chance similarity between unrelated
    # flanks (measured 14/2500 before this change, 0/2500 after).
    res = parasail.sw_striped_sat(outer5, outer3, SIG_GAP_OPEN, SIG_GAP_EXTEND,
                                  GENERIC_MATRIX)
    if res.score <= 0:
        return bounds
    if evalue(bits(res.score), len(outer5), len(outer3)) > max_evalue:
        return bounds
    qe, re_ = res.end_query, res.end_ref
    rev = parasail.sw_striped_sat(outer5[:qe + 1][::-1], outer3[:re_ + 1][::-1],
                                  SIG_GAP_OPEN, SIG_GAP_EXTEND, GENERIC_MATRIX)
    qb = qe - rev.end_query
    rb = re_ - rev.end_ref
    off3 = bounds.l3e + 1
    return Bounds(l5b=qb, l5e=qe, l3b=off3 + rb, l3e=off3 + re_,
                  margin_bits=bounds.margin_bits)


def _rediscover_outer(S: str, lo: int, hi: int, matrix, gaps: Gaps) -> tuple | None:
    """Stage 1, re-run on the outer segments under a matrix calibrated to the
    outer pair itself. `lo` is the length of the 5' outer segment and `hi` the
    start of the 3' outer segment."""
    outer5, outer3 = S[:lo], S[hi:]
    if not outer5 or not outer3:
        return None
    res = parasail.sw_striped_sat(outer5, outer3, gaps.open, gaps.extend, matrix)
    if res.score <= 0:
        return None
    qe, re_ = res.end_query, res.end_ref
    rev = parasail.sw_striped_sat(outer5[:qe + 1][::-1], outer3[:re_ + 1][::-1],
                                  gaps.open, gaps.extend, matrix)
    return (qe - rev.end_query, qe, hi + re_ - rev.end_ref, hi + re_)


# --------------------------------------------------------------------------- #
# Stage 5 and the classify() entry point
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Result:
    seq_id: str
    seq_len: int
    status: str
    ltr5_start: int | None      # 1-based inclusive
    ltr5_end: int | None
    ltr3_start: int | None
    ltr3_end: int | None
    ltr5_len: int | None
    ltr3_len: int | None
    flank5_len: int | None
    flank3_len: int | None
    aln_len: int | None
    n_sites: int | None
    n_ts: int | None
    n_tv: int | None
    n_gapcols: int | None
    identity: float | None
    p_dist: float | None
    k2p: float | None
    k2p_se: float | None
    bitscore: float | None
    flank_margin_bits: float | None
    cigar: str | None
    # Appended after `cigar`, not inserted before it: every column that existed
    # before these two keeps its 1-based index, so a documented recipe like
    # `cut -f1,4-7,19,23` still selects the same fields.
    motif: str | None           # e.g. "tg...ca": the two terminal dinucleotides
    k2p_time: int | None        # years since insertion; needs a mutation rate
    # Appended for the same reason, and filled by `genome.annotate` rather than
    # here: locating an LTR pair needs no reference, and every one of these is
    # `None` unless `--genome` was given.
    orientation: str | None     # '+'/'-' of the record against its header locus
    tsd: str | None             # target-site duplication at the called boundary
    tsd_offset: str | None      # "d5,d3": the shift at which `tsd` was found
    tsd_input: str | None       # the same, at the record's termini as supplied


_N_DATA_FIELDS = len(fields(Result)) - 3   # everything after seq_id/seq_len/status


def _empty(seq_id: str, seq_len: int, status: str) -> Result:
    return Result(seq_id, seq_len, status, *([None] * _N_DATA_FIELDS))


def _refine(q: str, r: str, gaps: Gaps = SIG_GAPS) -> tuple[str, str]:
    """Stage 5: exact global alignment via WFA. Penalties are converted from
    the parasail scheme so WFA optimises the same objective."""
    x, o, e = wfa_penalties(match=SCALE, mismatch=SCALE,
                            open_p=gaps.open, ext_p=gaps.extend)
    al = WavefrontAligner(q, mismatch=x, gap_opening=o, gap_extension=e,
                          scope="full", span="end-to-end")
    al(r)
    return aligned_pair_from_cigartuples(al.cigartuples, q, r)


def _refine_matrix(q: str, r: str, matrix, gaps: Gaps = SIG_GAPS) -> tuple[str, str]:
    """Ablation alternative to _refine: exact global alignment under the
    ti/tv-aware calibrated matrix instead of WFA's uniform mismatch penalty."""
    res = parasail.nw_trace_striped_sat(q, r, gaps.open, gaps.extend, matrix)
    return res.traceback.query, res.traceback.ref


def classify(seq_id: str, S: str, **kw) -> Result:
    """Locate the LTR pair and measure its divergence.

    `t_bits=None` uses the divergence-aware schedule (`t_bits_for(d_hat)`); an
    explicit value pins the threshold, which is what `--flank-bits` does.

    `flank_sensitivity` selects how much of a short flank's available evidence
    must be presented before it is believed (see `effective_t_bits`). It encodes
    a prior about the input, not about the sequence, which is why it is a choice
    and not a calibrated constant: `"strict"` suits tightly-extracted structural
    predictions, the looser settings suit input padded with genomic context.

    `mutation_rate` (substitutions per site per year) turns the divergence into
    a `k2p_time` in years; without it that column is `None`.

    `tsd_credit` is external evidence, in bits, that the sequence's own termini
    are already the element's boundaries -- `--tsd-anchor` turns a genomic
    target-site duplication into it. It is a plain number here on purpose: this
    module locates LTR pairs and has no business knowing what a genome is.

    Keyword knobs after `min_bitscore` drive the benchmark ablations; their
    defaults reproduce production behaviour. See `_classify` for the full
    signature -- this wrapper exists only to keep the public return type a
    plain `Result`.
    """
    return _classify(seq_id, S, **kw)[0]


def _classify(seq_id: str, S: str, *, cs: bool = False, t_bits: float | None = None,
              max_evalue: float = MAX_EVALUE, w0: int = W0,
              min_bitscore: float | None = None, matrix=None,
              use_stage3: bool = True, use_stage4: bool = True,
              trim: int = 0, refine: str = "wfa",
              snap_mode: str = "binary", inner: str = "none",
              comp: str = "element", gap_scheme: str = "adaptive",
              keep_weak: bool = True, stage4_recal: bool = True,
              flank_sensitivity: str = "strict",
              mutation_rate: float | None = None,
              tsd_credit: float = 0.0) -> tuple[Result, tuple[str, str] | None]:
    """`classify` plus the final aligned LTR pair.

    The pair is what `extras.py` builds the IUPAC consensus from. Returning it
    here rather than re-aligning is not just an optimisation: it makes the
    consensus and the reported divergence two readings of one alignment, so
    they cannot disagree. It stays off `Result` because every record would then
    carry two more kilobyte-scale strings back from its worker process, on
    every run, for a payload the TSV never emits.
    """
    L = len(S)
    if L < MIN_LEN:
        return _empty(seq_id, L, "too_short"), None
    if all(ch == "N" for ch in S):
        return _empty(seq_id, L, "all_ambiguous"), None

    hit = discover(S, GENERIC_MATRIX, SIG_GAPS, w0=w0)
    if hit is None:
        return _empty(seq_id, L, "no_pair"), None

    if flank_sensitivity not in FLANK_SENSITIVITY:
        raise ValueError(f"unknown flank_sensitivity {flank_sensitivity!r}; "
                         f"choices: {sorted(FLANK_SENSITIVITY)}")
    beta = FLANK_SENSITIVITY[flank_sensitivity]
    d_hat = None
    alpha = generic_alpha()
    gaps = gaps_for_scheme(gap_scheme if gap_scheme != "adaptive" else "legacy")
    if matrix is None:
        cal = calibrate_full(S, ltr_spans(S, hit), comp=comp, gap_scheme=gap_scheme)
        matrix, d_hat, gaps, alpha = cal.matrix, cal.d_hat, cal.gaps, cal.alpha
        # Seed the second pass at the window the first pass settled on rather
        # than re-growing from w0. This is an efficiency change, not a fix: on a
        # 300-element sweep over LTRs of 1.6-6 kb at p in {0.02,0.1,0.2}, every
        # one of which grew the window, restarting at w0 moved no boundary by
        # more than 10 bp and never returned None. It removes the repeated
        # doubling passes, which for a 6 kb-LTR element is three redundant
        # full-window Smith-Waterman alignments per record.
        hit = discover(S, matrix, gaps, w0=hit.w) or hit
    tb = t_bits if t_bits is not None else t_bits_for(d_hat)

    spans = ltr_spans(S, hit)
    if use_stage3:
        bounds = snap_bounds(S, spans, matrix, gaps, tb, mode=snap_mode,
                             inner=inner, alpha=alpha, beta=beta,
                             credit=tsd_credit)
    else:
        bounds = Bounds(*spans, None)

    if use_stage4:
        outer = outermost(S, bounds, matrix, max_evalue)
        if outer.spans != bounds.spans and stage4_recal:
            # Stage 4 found a different, uncharacterised pair. Re-run Stages 2-3
            # on it, as the spec requires: its divergence is not the inner pair's,
            # so neither the inner pair's matrix nor its boundaries apply, and
            # nothing has yet asked whether the OUTER pair reaches the termini.
            lo, hi = bounds.l5b, bounds.l3e + 1
            cal2 = calibrate_full(S, outer.spans, comp=comp, gap_scheme=gap_scheme)
            re_sp = _rediscover_outer(S, lo, hi, cal2.matrix, cal2.gaps) or outer.spans
            tb2 = t_bits if t_bits is not None else t_bits_for(cal2.d_hat)
            if use_stage3:
                bounds = snap_bounds(S, re_sp, cal2.matrix, cal2.gaps, tb2,
                                     mode=snap_mode, inner=inner,
                                     alpha=cal2.alpha, beta=beta,
                                     credit=tsd_credit)
            else:
                bounds = Bounds(*re_sp, outer.margin_bits)
            matrix, gaps = cal2.matrix, cal2.gaps
        else:
            bounds = outer

    q = S[bounds.l5b:bounds.l5e + 1]
    r = S[bounds.l3b:bounds.l3e + 1]
    if not q or not r:
        return _empty(seq_id, L, "no_pair"), None
    a, b = _refine(q, r, gaps) if refine == "wfa" else _refine_matrix(q, r, matrix, gaps)
    if trim:
        if len(a) <= 2 * trim:
            return _empty(seq_id, L, "no_pair"), None
        a, b = a[trim:-trim], b[trim:-trim]
    counts = count_substitutions(a, b)
    if counts.n_sites == 0:
        return _empty(seq_id, L, "no_pair"), None

    # Significance is ALWAYS scored with the generic model -- GENERIC_MATRIX and
    # SIG_GAPS -- never the calibrated one. MAX_EVALUE is calibrated against that
    # scale; the calibrated matrix is on a different one (at low divergence a
    # match scores +8 vs GENERIC's +4), so the same numeric threshold is not
    # transferable, and mixing them leaked spurious hits in Stage 4 (14/2500 ->
    # 0/2500 once gated on GENERIC). Pinning the gap penalties here too is what
    # lets `gap_scheme` change the ALIGNMENT without moving the significance
    # threshold underneath it.
    score = parasail.nw_striped_sat(q, r, SIG_GAP_OPEN, SIG_GAP_EXTEND,
                                    GENERIC_MATRIX).score
    bitscore = bits(score)
    weak = (evalue(bitscore, hit.w, hit.w) > max_evalue
            or (min_bitscore is not None and bitscore < min_bitscore))
    if weak and not keep_weak:
        return _empty(seq_id, L, "no_pair"), None

    d, se = k2p_distance(counts)
    # Reported straight off the settled boundaries, never searched for: the tool
    # uses no terminal-motif prior anywhere, so this column stays an INDEPENDENT
    # check on the boundary call rather than a restatement of it.
    motif = (f"{q[:2]}...{r[-2:]}".lower()
             if len(q) >= 2 and len(r) >= 2 else None)
    # A located-but-insignificant pair is reported, not deleted: its coordinates,
    # counts and divergence are real measurements, and `status == "pass"` still
    # means exactly what it always did, so the TSV remains a clean filter.
    status = "weak_pair" if weak else ("pass" if d is not None else "k2p_undefined")
    pd = p_distance(counts)
    aln_str = cs_string(a, b) if cs else extended_cigar(a, b)
    return (Result(
        seq_id=seq_id, seq_len=L, status=status,
        ltr5_start=bounds.l5b + 1, ltr5_end=bounds.l5e + 1,
        ltr3_start=bounds.l3b + 1, ltr3_end=bounds.l3e + 1,
        ltr5_len=bounds.l5e - bounds.l5b + 1, ltr3_len=bounds.l3e - bounds.l3b + 1,
        flank5_len=bounds.l5b, flank3_len=L - 1 - bounds.l3e,
        aln_len=counts.aln_len, n_sites=counts.n_sites, n_ts=counts.n_ts,
        n_tv=counts.n_tv, n_gapcols=counts.n_gapcols,
        identity=counts.n_match / counts.n_sites,
        p_dist=pd, k2p=d, k2p_se=se, bitscore=bitscore,
        flank_margin_bits=bounds.margin_bits, cigar=aln_str,
        motif=motif, k2p_time=insertion_time(d, mutation_rate),
        orientation=None, tsd=None, tsd_offset=None, tsd_input=None,
    ), (a, b))
