import math
import pytest
import numpy as np
import parasail
from ltrk2p.k2p import SubstCounts
from ltrk2p.scoring import (SCALE, k2p_probs, logodds_bits, parasail_matrix,
                            GENERIC_MATRIX, estimate_params, bits, evalue,
                            wfa_penalties)

def test_k2p_probs_sum_to_one():
    for d in (0.0, 0.05, 0.3, 0.6, 1.5):
        ps, pti, ptv = k2p_probs(d, 2.0)
        assert ps + pti + 2 * ptv == pytest.approx(1.0)

def test_k2p_probs_zero_distance_is_identity():
    ps, pti, ptv = k2p_probs(0.0, 2.0)
    assert ps == pytest.approx(1.0) and pti == pytest.approx(0.0)

def test_k2p_probs_saturate_to_quarter():
    ps, pti, ptv = k2p_probs(50.0, 2.0)
    assert ps == pytest.approx(0.25, abs=1e-6) and ptv == pytest.approx(0.25, abs=1e-6)

def test_k2p_probs_roundtrip_through_distance():
    """Probabilities generated at d must recover d via the K2P formula."""
    for d in (0.05, 0.15, 0.3, 0.5):
        ps, pti, ptv = k2p_probs(d, 2.0)
        P, Q = pti, 2 * ptv
        back = -0.5 * math.log(1 - 2*P - Q) - 0.25 * math.log(1 - 2*Q)
        assert back == pytest.approx(d, rel=1e-9)

def test_logodds_transitions_cost_less_than_transversions():
    b = logodds_bits(0.3, 2.0, {c: 0.25 for c in "ACGT"})
    assert b[("A", "G")] > b[("A", "C")]      # ts penalised less than tv
    assert b[("A", "A")] > b[("A", "G")]      # match beats ts

def test_logodds_match_positive_mismatch_negative_at_target_divergence():
    b = logodds_bits(0.3, 2.0, {c: 0.25 for c in "ACGT"})
    assert b[("A", "A")] > 0 and b[("A", "C")] < 0

def test_logodds_composition_adjustment_downweights_common_bases():
    at_rich = {"A": 0.35, "T": 0.35, "C": 0.15, "G": 0.15}
    b = logodds_bits(0.3, 2.0, at_rich)
    uniform = logodds_bits(0.3, 2.0, {c: 0.25 for c in "ACGT"})
    # matching a common base is less informative than matching a rare one
    assert b[("A", "A")] < uniform[("A", "A")]
    assert b[("C", "C")] > uniform[("C", "C")]

def test_parasail_matrix_scores_N_as_zero_both_directions():
    m = parasail_matrix(logodds_bits(0.3, 2.0, {c: 0.25 for c in "ACGT"}))
    arr = np.array(m.matrix)
    # alphabet is ACGTN -> index 4 is N; row AND column must both be zero
    assert list(arr[4][:5]) == [0, 0, 0, 0, 0]
    assert [arr[i][4] for i in range(5)] == [0, 0, 0, 0, 0]

def test_generic_matrix_is_plus_one_minus_one_scaled():
    arr = np.array(GENERIC_MATRIX.matrix)
    assert arr[0][0] == SCALE and arr[0][1] == -SCALE

def test_estimate_params_recovers_kappa():
    # 20 transitions, 10 transversions -> kappa_hat = 2*P/Q = 2*20/10 = 4
    c = SubstCounts(n_sites=100, n_match=70, n_ts=20, n_tv=10, n_gapcols=0, aln_len=100)
    d, kappa, freqs = estimate_params(c, "ACGT" * 25)
    assert kappa == pytest.approx(4.0)
    assert freqs["A"] == pytest.approx(0.25)

def test_estimate_params_defaults_kappa_when_no_transversions():
    c = SubstCounts(n_sites=100, n_match=90, n_ts=10, n_tv=0, n_gapcols=0, aln_len=100)
    _, kappa, _ = estimate_params(c, "ACGT" * 25)
    assert kappa == 2.0

def test_estimate_params_clamps_extreme_kappa():
    c = SubstCounts(n_sites=1000, n_match=900, n_ts=99, n_tv=1, n_gapcols=0, aln_len=1000)
    _, kappa, _ = estimate_params(c, "ACGT" * 250)
    assert kappa == 10.0

def test_evalue_decreases_with_bitscore_and_grows_with_window():
    assert evalue(50, 1000, 1000) < evalue(20, 1000, 1000)
    assert evalue(50, 8000, 8000) > evalue(50, 1000, 1000)

@pytest.mark.parametrize("m_,xp,op,ep", [(0, 4, 8, 2), (4, 4, 24, 8), (2, 3, 10, 2)])
def test_wfa_penalties_match_parasail_optimum(m_, xp, op, ep):
    """The conversion must make WFA optimise the SAME objective as parasail.

    Derivation: Score = m*M - x'*X - gapcost with M = (n1+n2-G)/2 - X, so
    maximising Score == minimising (m+x')*X + gapcost + (m/2)*G. Doubling to
    stay integral gives the exact identity, verified for these three schemes
    including alignments containing indels:

        -wfa_score == m*(len(a) + len(b)) - 2*parasail_score
    """
    import random
    from pywfa import WavefrontAligner
    random.seed(5)
    x, o, e = wfa_penalties(m_, xp, op, ep)
    mat = parasail.matrix_create("ACGT", m_, -xp)
    for _ in range(15):
        n = random.randint(100, 300)
        a = "".join(random.choice("ACGT") for _ in range(n))
        b = "".join(random.choice("ACGT") if random.random() < 0.2 else ch for ch in a)
        if random.random() < 0.5:                      # force some indels
            b = b[:len(b) // 2] + b[len(b) // 2 + 3:]
        par = parasail.nw_striped_sat(a, b, op, ep, mat).score
        al = WavefrontAligner(a, mismatch=x, gap_opening=o, gap_extension=e,
                              scope="score", span="end-to-end")
        al(b)
        assert -al.score == m_ * (len(a) + len(b)) - 2 * par
