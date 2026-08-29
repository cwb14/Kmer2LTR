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
MAX_EVALUE = 1e-10        # significance floor; Task 8 reuses this same constant
# Task 16 calibrated this against bench/out/negatives.fa (46,823 real non-LTR
# TEs: DNA transposons, LINEs, SINEs, Helitrons, satellites) and a dinucleotide
# -shuffled null. The shuffled null is clean at any threshold tested (<=0.02%
# throughout), so false "pass" calls on negatives.fa are not chance alignment
# noise -- they come from real, strongly-significant direct terminal repeats
# inside specific TE subclasses (satellites definitionally; a handful of
# library-consensus entries in Helitron/CACTA/hAT/Jockey that are themselves
# near-exact tandem constructions), which no significance threshold can
# separate from a true LTR pair without the family classification this tool
# explicitly declines to do (see "Non-goals"). At the old default (1e-3),
# false-"pass" rate on negatives.fa was 9.09%; tightening alone never reaches
# <1% without unacceptable true-positive cost: 2.86% FP still costs 9.8
# points of arab_ltr_all_clean pass rate at 1e-30, and even 1e-200 (0.41% FP)
# collapses real-data pass rate to 31%. 1e-10 is the last point before that
# cliff: FP nearly halves (9.09%->5.23%) while arab_ltr_all_clean pass rate
# moves only 99.80%->99.67% (13/10307 records). The cost concentrates in the
# tool's already-hardest population (d_nominal in {0.4,0.5}, no added flank,
# on the real gold-perturbed grid: 72.8%->53.8%) -- an amplification of an
# existing high-divergence weakness, not a new failure mode. Full sweep:
# bench/out/memo_real.md.


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


from .k2p import count_substitutions, k2p_distance, p_distance
from .scoring import estimate_params, logodds_bits, parasail_matrix, wfa_penalties
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


T_BITS = 10.0    # evidence required to claim a flank exists; benchmark-calibrated
# Task 14 swept {2,5,8,10,15,20,30} on 260,876 perturbed REAL gold elements
# (bench/out/cells_tbits_*.json, bench/out/memo_bench.md). At the old default
# (5.0), false-flank rate at d=0.3 was 26.8%; at 10.0 it is 8.5% (a 3.2x drop),
# while large-flank detection is nearly unchanged (det@50 93.0%->91.3%,
# det@100 91.7%->91.6%). The cost lands on 10-20bp flanks (det@10 75.5%->
# 39.4%, det@20 88.8%->70.5%), which the spec already documents as sitting
# near the detection floor. Applied in Task 16; see docs/
# 2026-08-28-ltrk2p-design.md and bench/out/
# task-14-report.md for the full sweep and the rejected divergence-aware
# alternative (a real Pareto improvement, deliberately not implemented yet).


@dataclass(frozen=True)
class Bounds:
    l5b: int
    l5e: int
    l3b: int
    l3e: int
    margin_bits: float | None


def terminal_snap(S: str, hit: Hit, matrix, t_bits: float = T_BITS) -> Bounds:
    """Decide, per terminus, whether homology reaches the end of the sequence.

    Exact model comparison rather than a greedy endpoint: extending the core
    to the terminus is accepted iff the extension costs less than t_bits.
    Under the penalised objective score - T*(free ends), that is s_ext > -T.
    """
    L = len(S)
    wstart = L - hit.w
    l5b, l5e, l3b, l3e = ltr_spans(S, hit)
    margins: list[float] = []

    # 5' test: reverse the outer segments so "begins anchored" == "next to the core".
    # sg_de -> query fully consumed (reach S[0]), ref end free.
    if l5b > 0:
        q = S[:l5b][::-1]
        r = S[wstart:l3b][::-1]
        if q and r:
            res = parasail.sg_de_striped_sat(q, r, GAP_OPEN, GAP_EXTEND, matrix)
            s5 = bits(res.score)
            margins.append(abs(s5 + t_bits))
            if s5 > -t_bits:
                l5b = 0
                l3b = l3b - (res.end_ref + 1)

    # 3' test: sg_qe -> ref fully consumed (reach S[L-1]), query end free.
    if l3e < L - 1:
        q = S[l5e + 1:hit.w]
        r = S[l3e + 1:]
        if q and r:
            res = parasail.sg_qe_striped_sat(q, r, GAP_OPEN, GAP_EXTEND, matrix)
            s3 = bits(res.score)
            margins.append(abs(s3 + t_bits))
            if s3 > -t_bits:
                l3e = L - 1
                l5e = l5e + (res.end_query + 1)

    l5b = max(0, l5b)
    l3e = min(L - 1, l3e)
    l5e = min(l5e, l3b - 1)
    return Bounds(l5b, l5e, l3b, l3e, min(margins) if margins else None)


# NOTE: MAX_EVALUE and `evalue` are already defined/imported in align.py by Task 5.
# Do NOT redefine them here.
MIN_OUTER = 50      # shorter outer segments cannot carry a significant pair


def outermost(S: str, bounds: Bounds, matrix, max_evalue: float = MAX_EVALUE) -> Bounds:
    """If a flank was called, prefer a significant repeat pair strictly outside it.

    This is what makes a retained (un-excised) nested element report the OUTER
    element rather than the younger, higher-scoring nested one.
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
    res = parasail.sw_striped_sat(outer5, outer3, GAP_OPEN, GAP_EXTEND, GENERIC_MATRIX)
    if res.score <= 0:
        return bounds
    if evalue(bits(res.score), len(outer5), len(outer3)) > max_evalue:
        return bounds
    qe, re_ = res.end_query, res.end_ref
    rev = parasail.sw_striped_sat(outer5[:qe + 1][::-1], outer3[:re_ + 1][::-1],
                                  GAP_OPEN, GAP_EXTEND, GENERIC_MATRIX)
    qb = qe - rev.end_query
    rb = re_ - rev.end_ref
    off3 = bounds.l3e + 1
    return Bounds(l5b=qb, l5e=qe, l3b=off3 + rb, l3e=off3 + re_,
                  margin_bits=bounds.margin_bits)


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


def _empty(seq_id: str, seq_len: int, status: str) -> Result:
    return Result(seq_id, seq_len, status, *([None] * 20))


def _refine(q: str, r: str) -> tuple[str, str]:
    """Stage 5: exact global alignment via WFA. Penalties are converted from
    the parasail scheme so WFA optimises the same objective."""
    x, o, e = wfa_penalties(match=SCALE, mismatch=SCALE, open_p=GAP_OPEN, ext_p=GAP_EXTEND)
    al = WavefrontAligner(q, mismatch=x, gap_opening=o, gap_extension=e,
                          scope="full", span="end-to-end")
    al(r)
    return aligned_pair_from_cigartuples(al.cigartuples, q, r)


def _refine_matrix(q: str, r: str, matrix) -> tuple[str, str]:
    """Ablation alternative to _refine: exact global alignment under the
    ti/tv-aware calibrated matrix instead of WFA's uniform mismatch penalty."""
    res = parasail.nw_trace_striped_sat(q, r, GAP_OPEN, GAP_EXTEND, matrix)
    return res.traceback.query, res.traceback.ref


def classify(seq_id: str, S: str, *, cs: bool = False, t_bits: float = T_BITS,
             max_evalue: float = MAX_EVALUE, w0: int = W0,
             min_bitscore: float | None = None, matrix=None,
             use_stage3: bool = True, use_stage4: bool = True,
             trim: int = 0, refine: str = "wfa") -> Result:
    """Locate the LTR pair and measure its divergence.

    Keyword knobs after min_bitscore drive Task 14's ablations; their defaults
    reproduce production behaviour.
    """
    L = len(S)
    if L < MIN_LEN:
        return _empty(seq_id, L, "too_short")
    if all(ch == "N" for ch in S):
        return _empty(seq_id, L, "all_ambiguous")

    hit = discover(S, GENERIC_MATRIX, w0=w0)
    if hit is None:
        return _empty(seq_id, L, "no_pair")
    if matrix is None:
        matrix, _, _ = calibrate(S, hit)
        hit = discover(S, matrix, w0=w0) or hit
    if use_stage3:
        bounds = terminal_snap(S, hit, matrix, t_bits)
    else:
        l5b, l5e, l3b, l3e = ltr_spans(S, hit)
        bounds = Bounds(l5b, l5e, l3b, l3e, None)
    if use_stage4:
        bounds = outermost(S, bounds, matrix, max_evalue)

    q = S[bounds.l5b:bounds.l5e + 1]
    r = S[bounds.l3b:bounds.l3e + 1]
    if not q or not r:
        return _empty(seq_id, L, "no_pair")
    a, b = _refine(q, r) if refine == "wfa" else _refine_matrix(q, r, matrix)
    if trim:
        if len(a) <= 2 * trim:
            return _empty(seq_id, L, "no_pair")
        a, b = a[trim:-trim], b[trim:-trim]
    counts = count_substitutions(a, b)

    # Significance is ALWAYS scored with GENERIC_MATRIX, never the calibrated one.
    # MAX_EVALUE is calibrated against GENERIC's score scale; the calibrated matrix
    # is on a different scale (at low divergence a match scores +8 vs GENERIC's +4),
    # so the same numeric threshold is not transferable. Mixing them leaked spurious
    # hits in Stage 4 (14/2500 -> 0/2500 once gated on GENERIC). The calibrated
    # matrix determines the ALIGNMENT; GENERIC determines SIGNIFICANCE.
    score = parasail.nw_striped_sat(q, r, GAP_OPEN, GAP_EXTEND, GENERIC_MATRIX).score
    bitscore = bits(score)
    if (evalue(bitscore, hit.w, hit.w) > max_evalue or counts.n_sites == 0
            or (min_bitscore is not None and bitscore < min_bitscore)):
        return _empty(seq_id, L, "no_pair")

    d, se = k2p_distance(counts)
    status = "pass" if d is not None else "k2p_undefined"
    pd = p_distance(counts)
    aln_str = cs_string(a, b) if cs else extended_cigar(a, b)
    return Result(
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
    )
