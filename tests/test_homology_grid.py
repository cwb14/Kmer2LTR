"""The homology grid's truth must be exact, or nothing measured against it means
anything. These tests check the construction, not the tool."""
import random

import pytest

from bench.homology_grid import (KAPPA, cells, d_for_p, evolve_tracked,
                                 make_element, shuffled_flank, true_alignment)
from ltrk2p.k2p import count_substitutions
from ltrk2p.scoring import k2p_probs


def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


# --------------------------------------------------------------------------- #
# d_for_p
# --------------------------------------------------------------------------- #

def test_d_for_p_inverts_the_k2p_difference_probability():
    for p in (0.05, 0.10, 0.20, 0.35, 0.5):
        d = d_for_p(p)
        assert 1.0 - k2p_probs(d, KAPPA)[0] == pytest.approx(p, abs=1e-9)


def test_d_for_p_is_zero_at_zero_and_monotone():
    assert d_for_p(0.0) == 0.0
    ds = [d_for_p(p) for p in (0.05, 0.10, 0.20, 0.35)]
    assert ds == sorted(ds)


def test_d_for_p_rejects_saturation():
    with pytest.raises(ValueError):
        d_for_p(0.80)


def test_target_p_is_realised_in_expectation():
    """The headline promise: 'add 25% mutations' must actually produce two
    copies differing at ~25% of sites."""
    rng = random.Random(0)
    anc = _rnd(4000, 1)
    for p in (0.05, 0.15, 0.25, 0.35):
        d = d_for_p(p)
        a, ta = evolve_tracked(anc, d / 2, KAPPA, 0.0, rng)
        b, tb = evolve_tracked(anc, d / 2, KAPPA, 0.0, rng)
        c = count_substitutions(*true_alignment(a, ta, b, tb))
        realised = (c.n_ts + c.n_tv) / c.n_sites
        assert realised == pytest.approx(p, abs=0.02), p


# --------------------------------------------------------------------------- #
# Tracked evolution and the exact alignment
# --------------------------------------------------------------------------- #

def test_trace_length_matches_the_descendant():
    rng = random.Random(2)
    s = _rnd(500, 3)
    d, tr = evolve_tracked(s, 0.2, KAPPA, 0.01, rng)
    assert len(d) == len(tr)
    kept = [i for i in tr if i is not None]
    assert kept == sorted(kept), "ancestral indices must be non-decreasing"
    assert all(0 <= i < len(s) for i in kept)


def test_true_alignment_is_exact_without_indels():
    """No indels means the true alignment is the ungapped column-by-column one,
    so it must agree with a naive zip exactly."""
    rng = random.Random(4)
    anc = _rnd(600, 5)
    a, ta = evolve_tracked(anc, 0.1, KAPPA, 0.0, rng)
    b, tb = evolve_tracked(anc, 0.1, KAPPA, 0.0, rng)
    qa, ra = true_alignment(a, ta, b, tb)
    assert qa == a and ra == b
    assert "-" not in qa and "-" not in ra


def test_true_alignment_reconstructs_both_sequences_with_indels():
    """With indels the alignment must still be lossless: dropping the gaps has
    to give back exactly the two evolved strings."""
    rng = random.Random(6)
    anc = _rnd(800, 7)
    for rate in (0.005, 0.02, 0.05):
        a, ta = evolve_tracked(anc, 0.2, KAPPA, rate, rng)
        b, tb = evolve_tracked(anc, 0.2, KAPPA, rate, rng)
        qa, ra = true_alignment(a, ta, b, tb)
        assert len(qa) == len(ra)
        assert qa.replace("-", "") == a
        assert ra.replace("-", "") == b
        assert not any(x == "-" and y == "-" for x, y in zip(qa, ra)), \
            "an all-gap column is not an alignment column"


def test_true_alignment_never_aligns_non_homologous_positions():
    """Every ungapped column must join two descendants of the SAME ancestral
    site -- that is the property the naive truncate-and-compare in
    bench/simulate.py does not have."""
    rng = random.Random(8)
    anc = _rnd(400, 9)
    a, ta = evolve_tracked(anc, 0.3, KAPPA, 0.02, rng)
    b, tb = evolve_tracked(anc, 0.3, KAPPA, 0.02, rng)
    qa, ra = true_alignment(a, ta, b, tb)
    i = j = 0
    for x, y in zip(qa, ra):
        if x != "-" and y != "-":
            assert ta[i] is not None and ta[i] == tb[j]
        i += x != "-"
        j += y != "-"


def test_indel_rate_actually_produces_more_gaps():
    rng = random.Random(10)
    anc = _rnd(2000, 11)
    prev = -1
    for rate in (0.0, 0.002, 0.01, 0.05):
        a, ta = evolve_tracked(anc, 0.1, KAPPA, rate, rng)
        b, tb = evolve_tracked(anc, 0.1, KAPPA, rate, rng)
        gaps = count_substitutions(*true_alignment(a, ta, b, tb)).n_gapcols
        assert gaps > prev or rate == 0.0
        prev = gaps
    assert prev > 0


# --------------------------------------------------------------------------- #
# Element assembly
# --------------------------------------------------------------------------- #

def test_truth_coordinates_slice_out_the_actual_ltrs():
    """The one invariant every downstream number depends on."""
    rng = random.Random(12)
    ltr, internal = _rnd(300, 13), _rnd(900, 14)
    for p in (0.0, 0.15, 0.35):
        for flank in (0, 5, 45):
            for ir in (0.0, 0.01):
                seq, t = make_element(ltr, internal, p, flank, ir, rng)
                assert len(seq) == t["seq_len"]
                assert t["flank5"] == flank and t["flank3"] == flank
                assert t["ltr5_start"] - 1 == flank
                assert len(seq) - t["ltr3_end"] == flank
                # the internal region is spliced back untouched
                assert seq[t["ltr5_end"]:t["ltr3_start"] - 1] == internal
                assert t["ltr5_start"] <= t["ltr5_end"] < t["ltr3_start"] <= t["ltr3_end"]


def test_zero_mutation_zero_indel_element_is_literally_perfect():
    rng = random.Random(15)
    ltr, internal = _rnd(300, 16), _rnd(900, 17)
    seq, t = make_element(ltr, internal, 0.0, 0, 0.0, rng)
    assert seq == ltr + internal + ltr
    assert t["realized_k2p"] == 0.0
    assert t["realized_p"] == 0.0
    assert t["true_n_gapcols"] == 0


def test_flank_is_composition_matched_but_reordered():
    rng = random.Random(18)
    pool = "AAAACCGT" * 20
    f = shuffled_flank(pool, 40, rng)
    assert len(f) == 40
    # composition of the shuffle is drawn from the pool's, so an AT-rich pool
    # cannot yield a GC-rich flank
    assert f.count("A") > f.count("G")


def test_grid_shape_is_what_the_brief_asked_for():
    got = list(cells())
    subs = [c for c in got if c[0] == "subs"]
    indel = [c for c in got if c[0] == "indel"]
    assert sorted({c[1] for c in subs}) == [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]
    assert sorted({c[3] for c in subs}) == [0, 5, 25, 45, 65]
    assert all(c[2] == 0.0 for c in subs)
    assert len(subs) == 8 * 5
    assert sorted({c[2] for c in indel}) == [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
    assert len(indel) == 6 * 2 * 5
