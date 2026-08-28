# tests/test_simulate.py
import random
import pytest
from bench.simulate import evolve, shuffle_dinuc, simulate_element
from ltrk2p.k2p import count_substitutions, k2p_distance

def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))

def test_evolve_zero_distance_is_identity():
    s = _rnd(500, 1)
    assert evolve(s, 0.0, 2.0, 0.0, random.Random(2)) == s

def test_evolve_produces_target_divergence():
    """Two independent branches of d/2 each -> pairwise K2P ~= d.

    Tolerance is tight (10%) and extends to d=0.7 deliberately: a single-draw
    substitution model passes a loose low-d check but drifts badly at high d
    (realised 0.856 at nominal 0.70), which would inject a phantom bias into
    the whole benchmark. Sampling from the exact K2P matrix is what keeps this
    test passing at the top of the range.
    """
    anc = _rnd(6000, 3)
    for target in (0.05, 0.2, 0.4, 0.7):
        r = random.Random(int(target * 1000))
        a = evolve(anc, target / 2, 2.0, 0.0, r)
        b = evolve(anc, target / 2, 2.0, 0.0, r)
        d, _ = k2p_distance(count_substitutions(a, b))
        assert d == pytest.approx(target, rel=0.10), f"nominal {target}, realised {d}"


def test_evolve_matches_exact_k2p_transition_probabilities():
    """The per-site outcome distribution must match k2p_probs exactly."""
    from ltrk2p.scoring import k2p_probs
    anc = "A" * 40000
    d, kappa = 0.4, 2.0
    out = evolve(anc, d, kappa, 0.0, random.Random(11))
    p_same, p_ti, p_tv = k2p_probs(d, kappa)
    n = len(out)
    assert out.count("A") / n == pytest.approx(p_same, abs=0.01)
    assert out.count("G") / n == pytest.approx(p_ti, abs=0.01)
    assert (out.count("C") + out.count("T")) / n == pytest.approx(2 * p_tv, abs=0.01)

def test_evolve_respects_kappa():
    anc = _rnd(6000, 4)
    a = evolve(anc, 0.15, 8.0, 0.0, random.Random(5))
    c = count_substitutions(anc, a)
    assert c.n_ts > 2 * c.n_tv       # high kappa -> transition-dominated

def test_evolve_indels_change_length():
    anc = _rnd(2000, 6)
    a = evolve(anc, 0.2, 2.0, 0.15, random.Random(7))
    assert a != anc and abs(len(a) - len(anc)) > 0

def test_shuffle_dinuc_preserves_composition():
    s = _rnd(3000, 8)
    t = shuffle_dinuc(s, random.Random(9))
    assert len(t) == len(s)
    assert sorted(t) == sorted(s)

def test_simulate_element_truth_coordinates_slice_correctly():
    ltr, internal = _rnd(400, 10), _rnd(1200, 11)
    seq, truth = simulate_element(ltr, internal, 0.2, 2.0, 50, 70, random.Random(12))
    assert truth["ltr5_start"] == 51
    assert seq[truth["ltr5_start"]-1:truth["ltr5_end"]] != ""
    assert truth["ltr3_end"] == len(seq) - 70
    assert truth["realized_k2p"] is not None

def test_simulate_element_no_flanks_is_perfectly_bounded():
    ltr, internal = _rnd(400, 13), _rnd(1200, 14)
    seq, truth = simulate_element(ltr, internal, 0.1, 2.0, 0, 0, random.Random(15))
    assert truth["ltr5_start"] == 1
    assert truth["ltr3_end"] == len(seq)
