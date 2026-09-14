"""Stage 1b: period selection, `--period-rule outermost`.

The failure this rule exists for, as it appears in real elements: the record
carries an internal duplication whose two copies are LONGER than the LTRs, so
the highest-scoring alignment of the two search windows sits on that register
instead of the LTR pair's. The wrong pair is always further in, which is what
makes "outermost" a selection rule rather than a preference.
"""
import random

import numpy as np
import pytest

from kmer2ltr.align import (MIN_INTERNAL, MIN_LEN, SIG_GAPS, Hit, _candidate_offsets,
                            _diagonal_matches, _diagonal_segment, _outermost_period,
                            classify, discover, ltr_spans)
from kmer2ltr.scoring import GENERIC_MATRIX


def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _shifted_register(seed=0):
    """An 1800 bp element carrying a longer duplication that OVERLAPS its LTRs.

    Built to the geometry the real failures have: the competing register is
    longer than the LTR pair, so best-score takes it, and its two copies start
    inside the 5' LTR and end inside the 3' LTR, so what is left outside them is
    two non-corresponding pieces of the same LTR. That is what defeats Stage 4 --
    it searches `S[:ltr5_start]` against `S[ltr3_end:]`, and here those two are
    not homologous to each other. A duplication that sat strictly INSIDE the
    LTRs would be recovered by Stage 4 alone and would test nothing.

        S = A + D + E + D + B,  D = C + M + C
        A = M[400:500], B = M[0:100]      chosen so S[0:300] == S[1500:1800]

    Returns the sequence and the LTR pair's 0-based inclusive spans.
    """
    C, M, E = _rnd(100, seed + 1), _rnd(500, seed + 2), _rnd(200, seed + 3)
    D = C + M + C
    S = M[400:500] + D + E + D + M[0:100]
    assert len(S) == 1800 and S[0:300] == S[1500:1800]
    return S, (0, 299, 1500, 1799)


# --------------------------------------------------------------------------- #
# the diagonal profile
# --------------------------------------------------------------------------- #

def test_diagonal_matches_agrees_with_brute_force():
    P, Q = _rnd(120, 11), _rnd(120, 12)
    counts, n = _diagonal_matches(P, Q)
    for d in range(-(len(P) - 1), len(P)):
        want = sum(1 for i in range(len(P))
                   if 0 <= i + d < len(Q) and P[i] == Q[i + d])
        assert int(counts[d % n]) == want, f"offset {d}"


def test_ambiguous_bases_never_count_as_matches():
    """N scores 0 against everything in GENERIC_MATRIX, including another N."""
    P = Q = "N" * 200
    counts, n = _diagonal_matches(P, Q)
    assert int(counts[0]) == 0
    assert _diagonal_segment(P, Q, 0) is None


def test_candidate_offsets_finds_a_planted_register():
    filler, motif = _rnd(600, 21), _rnd(250, 22)
    P = filler[:200] + motif + filler[200:350]
    Q = filler[350:500] + motif + filler[500:]
    planted = 150 - 200                     # Q index of motif minus P index
    assert any(abs(d - planted) < 32 for d in _candidate_offsets(P, Q))


def test_diagonal_segment_locates_the_conserved_stretch():
    """Windows are always the same length, so the two are built that way here."""
    motif = _rnd(300, 31)
    P = _rnd(100, 32) + motif + _rnd(100, 35)
    Q = _rnd(150, 33) + motif + _rnd(50, 34)
    seg = _diagonal_segment(P, Q, 50)
    assert seg is not None
    i_b, i_e, score = seg
    # A chance match abutting the motif can extend the segment by a base or two.
    assert abs(i_b - 100) <= 5 and abs(i_e - 399) <= 5
    assert score >= 300                     # exact match, +1 bit per base


def test_diagonal_segment_takes_no_break_even_bases_in_front():
    """A match then a mismatch just ahead of the conserved stretch net zero, so
    the maximal segment could begin at either end of them. It has to begin
    where the stretch does: the selection rule reads extra bases as span, and
    would move the 5' boundary out by exactly those two."""
    motif = _rnd(300, 36)
    P = "A" * 98 + "GT" + motif + "A" * 50
    Q = "C" * 98 + "GA" + motif + "C" * 50
    assert _diagonal_segment(P, Q, 0) == (100, 399, 300)


def test_a_diagonal_shorter_than_the_minimum_is_not_a_candidate():
    P, Q = _rnd(400, 41), _rnd(400, 42)
    assert _diagonal_segment(P, Q, 400 - MIN_LEN + 1) is None


# --------------------------------------------------------------------------- #
# the selection rule
# --------------------------------------------------------------------------- #

def test_outermost_recovers_the_ltr_pair_the_best_score_rule_misses():
    S, want = _shifted_register()
    best = discover(S, GENERIC_MATRIX, SIG_GAPS)
    outer = discover(S, GENERIC_MATRIX, SIG_GAPS, rule="outermost")
    assert ltr_spans(S, best) != want, "synthetic does not reproduce the failure"
    assert ltr_spans(S, outer) == want


def test_the_two_rules_agree_on_an_element_with_no_competing_register():
    ltr, internal = _rnd(300, 51), _rnd(1200, 52)
    S = _rnd(80, 53) + ltr + internal + ltr + _rnd(120, 54)
    assert (ltr_spans(S, discover(S, GENERIC_MATRIX, SIG_GAPS))
            == ltr_spans(S, discover(S, GENERIC_MATRIX, SIG_GAPS, rule="outermost")))


def test_a_pair_already_reaching_both_termini_is_returned_untouched():
    """The short-circuit: nothing can out-span a pair that starts at base 0 and
    ends at the last base, so no candidate work is done at all."""
    ltr, internal = _rnd(300, 61), _rnd(1200, 62)
    S = ltr + internal + ltr
    best = discover(S, GENERIC_MATRIX, SIG_GAPS)
    assert (best.qb, best.re) == (0, best.w - 1)
    assert discover(S, GENERIC_MATRIX, SIG_GAPS, rule="outermost") == best


def test_a_wider_candidate_is_refused_when_it_leaves_no_internal_region():
    """Two copies that nearly abut are a tandem duplication, not an LTR pair.

    Same windows and same incumbent both times; only the amount of sequence
    between the two copies differs, so the internal floor is the sole cause of
    the different answer.
    """
    P = Q = _rnd(300, 71)
    inner = Hit(score=200, qb=50, qe=100, rb=150, re=200, w=300)
    tight = _outermost_period(P, Q, inner, wstart=300)      # internal 0 bp
    roomy = _outermost_period(P, Q, inner, wstart=700)      # internal 400 bp
    assert tight == inner, "abutting copies must not win on span alone"
    assert roomy.qb == 0 and roomy.re == 299
    assert 300 + 0 - 299 - 1 < MIN_INTERNAL <= 700 + 0 - 299 - 1


def test_an_insignificant_candidate_never_displaces_the_incumbent():
    P, Q = _rnd(2000, 81), _rnd(2000, 82)          # unrelated: no real register
    inner = Hit(score=200, qb=900, qe=1000, rb=1000, re=1100, w=2000)
    assert _outermost_period(P, Q, inner, wstart=2000) == inner


def test_the_incumbent_holds_a_tie():
    ltr, internal = _rnd(300, 91), _rnd(900, 92)
    S = ltr + internal + ltr
    best = discover(S, GENERIC_MATRIX, SIG_GAPS)
    assert _outermost_period(S[:best.w], S[len(S) - best.w:], best,
                            len(S) - best.w) == best


def test_an_unknown_rule_fails_fast():
    with pytest.raises(ValueError, match="unknown period rule"):
        discover(_rnd(500, 101), GENERIC_MATRIX, SIG_GAPS, rule="outer")
    with pytest.raises(ValueError, match="unknown period rule"):
        discover("ACGT", GENERIC_MATRIX, SIG_GAPS, rule="outer")


def test_degenerate_inputs_do_not_raise_under_either_rule():
    for S in ("", "ACGT", "N" * 500, "A" * 500, _rnd(150, 111), _rnd(3000, 112)):
        discover(S, GENERIC_MATRIX, SIG_GAPS, rule="outermost")


# --------------------------------------------------------------------------- #
# through classify
# --------------------------------------------------------------------------- #

def test_classify_reports_the_outer_pair_under_the_flag():
    S, (b5, e5, b3, e3) = _shifted_register()
    plain = classify("x", S)
    outer = classify("x", S, period_rule="outermost")
    assert (outer.ltr5_start, outer.ltr5_end) == (b5 + 1, e5 + 1)
    assert (outer.ltr3_start, outer.ltr3_end) == (b3 + 1, e3 + 1)
    assert (plain.ltr5_start, plain.ltr3_end) != (outer.ltr5_start, outer.ltr3_end), \
        "Stage 4 already recovers this pair; the synthetic tests nothing"
    assert outer.status == "pass" and outer.k2p == 0.0


def test_classify_rejects_an_unknown_rule():
    with pytest.raises(ValueError, match="unknown period rule"):
        classify("x", _rnd(2000, 121), period_rule="widest")


def test_the_chosen_pair_always_contains_the_one_it_replaces():
    """The flag can only widen a located pair, never trade one side for the other.

    Asserted at Stage 1, where the property is structural. Downstream stages then
    recalibrate and snap from the wider seed, so the guarantee is on discovery.
    """
    r = random.Random(5)
    for i in range(60):
        unit = _rnd(r.randint(30, 200), 900 + i)
        ltr = _rnd(r.randint(80, 300), 950 + i) + unit * r.randint(2, 6)
        S = (_rnd(r.randint(0, 200), 1000 + i) + ltr
             + _rnd(r.randint(200, 1500), 1050 + i) + ltr
             + _rnd(r.randint(0, 200), 1100 + i))
        best = discover(S, GENERIC_MATRIX, SIG_GAPS)
        outer = discover(S, GENERIC_MATRIX, SIG_GAPS, rule="outermost")
        if best is None:
            assert outer is None
            continue
        b5, _, _, b3 = ltr_spans(S, best)
        o5, _, _, o3 = ltr_spans(S, outer)
        assert o5 <= b5 and o3 >= b3, f"record {i}: {(o5, o3)} does not contain {(b5, b3)}"
