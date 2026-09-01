"""Tests for the per-element alignment model: gap schemes, the divergence-aware
flank threshold, composition source, and the non-destructive significance gate."""
import math
import random

import pytest

from ltrk2p.align import (GAP_EXTEND_BITS, GAP_OPEN_BITS, GAP_SCHEMES, Gaps,
                          effective_t_bits,
                          SIG_GAPS, T_BITS, T_BITS_SCHEDULE, calibrate_full,
                          classify, discover, gaps_for_scheme,
                          gaps_from_alignment, ltr_spans, t_bits_for)
from ltrk2p.k2p import count_gap_runs, count_substitutions
from ltrk2p.scoring import GENERIC_MATRIX, SCALE


def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _evolve(s, p, seed):
    r = random.Random(seed)
    ti = {"A": "G", "G": "A", "C": "T", "T": "C"}
    tv = {"A": "CT", "G": "CT", "C": "AG", "T": "AG"}
    out = []
    for c in s:
        if r.random() < p:
            out.append(ti[c] if r.random() < 2 / 3 else r.choice(tv[c]))
        else:
            out.append(c)
    return "".join(out)


# --------------------------------------------------------------------------- #
# count_gap_runs
# --------------------------------------------------------------------------- #

def test_gap_runs_counts_openings_not_columns():
    #      a: one 3-column gap; b: two separate 1-column gaps
    a = "AC---GTAGT"
    b = "ACTTTG-AG-"
    assert count_gap_runs(a, b) == 3
    assert count_substitutions(a, b).n_gapcols == 5


def test_gap_runs_zero_when_ungapped():
    assert count_gap_runs("ACGT", "AGGT") == 0


def test_gap_runs_rejects_ragged_input():
    with pytest.raises(ValueError):
        count_gap_runs("ACGT", "ACG")


# --------------------------------------------------------------------------- #
# Gap schemes
# --------------------------------------------------------------------------- #

def test_named_gap_schemes_convert_bits_to_scale_units():
    assert gaps_for_scheme("legacy") == Gaps(6 * SCALE, 2 * SCALE)
    assert gaps_for_scheme("static") == Gaps(10 * SCALE, 1 * SCALE)
    assert gaps_for_scheme("legacy") == SIG_GAPS


def test_unknown_gap_scheme_fails_loudly():
    with pytest.raises(ValueError):
        gaps_for_scheme("whatever")


def test_adaptive_gaps_track_the_observed_indel_rate():
    """A gap-rich core must imply a CHEAPER opening than a gap-poor one --
    that is the whole point of estimating the rate instead of fixing it."""
    # 1 opening in ~200 columns
    sparse_a = "A" * 100 + "---" + "A" * 100
    sparse_b = "A" * 100 + "CGT" + "A" * 100
    # 10 openings in ~200 columns
    rich_a = ("A" * 17 + "---") * 10
    rich_b = ("A" * 17 + "CGT") * 10
    sparse = gaps_from_alignment(sparse_a, sparse_b, count_substitutions(sparse_a, sparse_b))
    rich = gaps_from_alignment(rich_a, rich_b, count_substitutions(rich_a, rich_b))
    assert rich.open < sparse.open


def test_adaptive_gaps_are_clamped_at_both_ends():
    """No gaps at all implies an infinite opening cost; one endless gap implies
    an infinite mean length. Both must land inside the declared range."""
    a, b = "ACGT" * 50, "ACGA" * 50               # no gaps anywhere
    g = gaps_from_alignment(a, b, count_substitutions(a, b))
    assert g.open == int(round(GAP_OPEN_BITS[1] * SCALE))
    assert GAP_EXTEND_BITS[0] * SCALE <= g.extend <= GAP_EXTEND_BITS[1] * SCALE
    for g2 in (gaps_from_alignment("-" * 40, "ACGT" * 10,
                                   count_substitutions("-" * 40, "ACGT" * 10)),):
        assert GAP_OPEN_BITS[0] * SCALE <= g2.open <= GAP_OPEN_BITS[1] * SCALE
        assert GAP_EXTEND_BITS[0] * SCALE <= g2.extend <= GAP_EXTEND_BITS[1] * SCALE


def test_adaptive_gaps_open_more_than_they_extend():
    """An affine scheme whose extension costs more than its opening is not
    affine in any useful sense; the clamps must make that unreachable."""
    for scheme in list(GAP_SCHEMES) + ["adaptive"]:
        g = gaps_for_scheme(scheme)
        assert g.open > g.extend > 0


def test_every_gap_scheme_runs_end_to_end():
    ltr = _rnd(400, 41)
    S = ltr + _rnd(1000, 42) + _evolve(ltr, 0.08, 43)
    for scheme in list(GAP_SCHEMES) + ["adaptive"]:
        r = classify("x", S, gap_scheme=scheme)
        assert r.status == "pass", scheme
        assert r.ltr5_start == 1 and r.ltr3_end == len(S), scheme


# --------------------------------------------------------------------------- #
# Divergence-aware T_BITS
# --------------------------------------------------------------------------- #

def test_schedule_is_monotone_and_ends_at_infinity():
    """A higher estimated divergence must never require LESS evidence to call a
    flank, and every d_hat must land in some bin."""
    bounds = [hi for hi, _ in T_BITS_SCHEDULE]
    values = [t for _, t in T_BITS_SCHEDULE]
    assert bounds == sorted(bounds)
    assert bounds[-1] == math.inf
    assert values == sorted(values)


def test_t_bits_for_falls_back_to_the_fixed_default_without_an_estimate():
    assert t_bits_for(None) == T_BITS


def test_t_bits_for_picks_the_bin_containing_d_hat():
    for hi, t in T_BITS_SCHEDULE:
        probe = (hi - 1e-9) if hi != math.inf else 10.0
        assert t_bits_for(probe) == t
    assert t_bits_for(0.0) == T_BITS_SCHEDULE[0][1]


def test_explicit_t_bits_overrides_the_schedule():
    """--flank-bits must pin the threshold; that is what makes the old fixed
    behaviour still reachable."""
    ltr = _rnd(400, 51)
    S = _rnd(12, 52) + ltr + _rnd(1000, 53) + _evolve(ltr, 0.30, 54)
    lo = classify("x", S, t_bits=2.0)
    hi = classify("x", S, t_bits=60.0)
    # a HIGHER threshold demands more evidence before believing a flank exists,
    # so it must snap at least as often -> at most as much called flank
    assert hi.flank5_len <= lo.flank5_len
    assert hi.flank5_len == 0 and lo.flank5_len > 0


# --------------------------------------------------------------------------- #
# Composition source
# --------------------------------------------------------------------------- #

def test_core_composition_differs_from_element_composition_when_it_should():
    """A GC-balanced LTR pair inside a strongly AT-rich internal region: the two
    composition sources must produce genuinely different matrices, or the knob
    is measuring nothing."""
    r = random.Random(7)
    ltr = "".join(r.choice("ACGT") for _ in range(400))
    internal = "".join(r.choice("ATATAT") for _ in range(3000))
    S = ltr + internal + ltr
    hit = discover(S, GENERIC_MATRIX, SIG_GAPS)
    spans = ltr_spans(S, hit)
    m_el = calibrate_full(S, spans, comp="element")
    m_core = calibrate_full(S, spans, comp="core")
    import numpy as np
    a = np.asarray(m_el.matrix.matrix)[:4, :4]
    b = np.asarray(m_core.matrix.matrix)[:4, :4]
    assert not np.array_equal(a, b), "composition source made no difference"


def test_both_composition_sources_recover_a_clean_element():
    ltr = _rnd(400, 61)
    S = ltr + _rnd(1500, 62) + _evolve(ltr, 0.10, 63)
    for comp in ("element", "core"):
        res = classify("x", S, comp=comp)
        assert res.status == "pass" and res.ltr5_start == 1 and res.ltr3_end == len(S)


# --------------------------------------------------------------------------- #
# Non-destructive significance gate
# --------------------------------------------------------------------------- #

def test_weak_pair_keeps_every_measurement_instead_of_nulling_the_row():
    """The gate's job is to LABEL an insignificant pair, not to erase a real
    measurement. Coordinates, counts and divergence must all survive."""
    ltr = _rnd(400, 71)
    S = ltr + _rnd(1200, 72) + _evolve(ltr, 0.02, 73)
    strong = classify("x", S)
    assert strong.status == "pass"
    weak = classify("x", S, min_bitscore=strong.bitscore + 100.0)
    assert weak.status == "weak_pair"
    for field in ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end",
                  "n_sites", "k2p", "bitscore", "cigar"):
        assert getattr(weak, field) is not None, field
    assert (weak.ltr5_start, weak.ltr3_end) == (strong.ltr5_start, strong.ltr3_end)
    assert weak.k2p == strong.k2p


def test_keep_weak_false_restores_the_destructive_gate():
    ltr = _rnd(400, 74)
    S = ltr + _rnd(1200, 75) + _evolve(ltr, 0.02, 76)
    strong = classify("x", S)
    gone = classify("x", S, min_bitscore=strong.bitscore + 100.0, keep_weak=False)
    assert gone.status == "no_pair" and gone.ltr5_start is None


def test_pass_still_means_exactly_what_it_meant():
    """The filter contract: a `pass` row is significant, so nothing that would
    previously have been dropped may now be labelled `pass`."""
    ltr = _rnd(400, 77)
    S = ltr + _rnd(1200, 78) + _evolve(ltr, 0.02, 79)
    r = classify("x", S, min_bitscore=1e9)
    assert r.status != "pass"
    assert classify("x", _rnd(3000, 80)).status == "no_pair"


# --------------------------------------------------------------------------- #
# The short-flank detection floor
# --------------------------------------------------------------------------- #

def test_a_flat_threshold_is_unreachable_below_its_own_floor():
    """The motivating arithmetic: a k-base flank can present at most k*alpha
    bits, so a flat t_bits cannot fire below k = t_bits/alpha whatever the
    sequence says. This pins the floor the cap exists to remove."""
    alpha = 1.84                                   # bits/base at d_hat ~ 0.10
    t = 8.0
    floor = t / alpha
    assert floor == pytest.approx(4.35, abs=0.05)
    # below the floor an uncapped threshold demands more than exists
    assert effective_t_bits(t, 3, alpha, None) == t
    assert 3 * alpha < t
    # the cap never demands more than the flank can supply
    for k in range(1, 20):
        assert effective_t_bits(t, k, alpha, 1.0) <= k * alpha + 1e-9


def test_cap_is_inert_above_the_floor_and_binds_below_it():
    alpha, t = 1.84, 8.0
    assert effective_t_bits(t, 100, alpha, 0.5) == t          # far above: inert
    assert effective_t_bits(t, 2, alpha, 0.5) == pytest.approx(1.84)


def test_cap_is_disabled_by_default_and_by_a_missing_alpha():
    assert effective_t_bits(10.0, 1, 2.0, None) == 10.0       # beta None
    assert effective_t_bits(10.0, 1, None, 0.5) == 10.0       # no alpha
    assert effective_t_bits(10.0, 0, 2.0, 0.5) == 10.0        # no flank
    assert effective_t_bits(10.0, 5, 0.0, 0.5) == 10.0        # degenerate alpha


def test_strict_reproduces_the_uncapped_behaviour_exactly():
    """`strict` is the default and must be indistinguishable from having no cap
    at all -- that is what makes every earlier measurement still valid."""
    ltr = _rnd(400, 81)
    for seed, div in ((82, 0.02), (83, 0.15), (84, 0.30)):
        S = _rnd(7, seed) + ltr + _rnd(1100, seed + 1) + _evolve(ltr, div, seed + 2)
        a = classify("x", S, flank_sensitivity="strict")
        b = classify("x", S)
        assert (a.ltr5_start, a.ltr5_end, a.ltr3_start, a.ltr3_end, a.k2p) == \
               (b.ltr5_start, b.ltr5_end, b.ltr3_start, b.ltr3_end, b.k2p)


def test_looser_settings_never_call_a_shorter_flank_than_stricter_ones():
    """Ordering is the contract the flag advertises: strict <= balanced <=
    sensitive in willingness to call a flank."""
    ltr = _rnd(400, 85)
    for seed in range(6):
        S = (_rnd(4, 900 + seed) + ltr + _rnd(1100, 910 + seed)
             + _evolve(ltr, 0.08, 920 + seed) + _rnd(4, 930 + seed))
        got = [classify("x", S, flank_sensitivity=m).flank5_len
               for m in ("strict", "balanced", "sensitive")]
        assert got == sorted(got), (seed, got)


def test_unknown_sensitivity_fails_loudly():
    with pytest.raises(ValueError):
        classify("x", _rnd(400, 86) * 3, flank_sensitivity="aggressive")


def test_alpha_falls_with_divergence_which_is_why_the_floor_widens():
    """alpha is the whole mechanism: it is read off the element's own matrix and
    it shrinks as the pair diverges, so the blind spot widens exactly where
    boundaries are hardest."""
    from ltrk2p.scoring import expected_random_bits, logodds_bits
    f = {b: 0.25 for b in "ACGT"}
    alphas = [expected_random_bits(logodds_bits(d, 2.0, f), f)
              for d in (0.01, 0.05, 0.10, 0.20, 0.35)]
    assert alphas == sorted(alphas, reverse=True)
    assert alphas[0] > 4.0 and alphas[-1] < 1.0
