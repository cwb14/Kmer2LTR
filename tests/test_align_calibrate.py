import random
import numpy as np
import pytest
from kmer2ltr.align import discover, calibrate, ltr_spans, core_alignment, MIN_CALIB_SITES
from kmer2ltr.scoring import GENERIC_MATRIX, SCALE

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
            out.append(ti[c] if r.random() < 2/3 else r.choice(tv[c]))
        else:
            out.append(c)
    return "".join(out)

def test_ltr_spans_are_absolute_and_ordered():
    """Coordinates must be absolute (not window-relative) and correctly ordered.

    Boundaries are asserted with tolerance, not exactly: on a flanked element a
    local alignment legitimately frays a base or two past the true edge where a
    chance match extends it (measured here: +1 on both LTR ends, and confirmed
    not to be a base coincidence). Removing that residual is Stage 3's job, not
    Stage 1's. Exact-boundary correctness is covered by
    test_ltr_spans_recovers_exact_ltr_sequence_unflanked, where an unflanked
    element makes the right answer unambiguous.
    """
    ltr = _rnd(300, 1)
    S = _rnd(50, 2) + ltr + _rnd(1000, 3) + ltr + _rnd(70, 4)
    h = discover(S, GENERIC_MATRIX)
    assert h is not None
    l5b, l5e, l3b, l3e = ltr_spans(S, h)

    # absolute, not window-relative, and correctly ordered
    assert 0 <= l5b <= l5e < l3b <= l3e < len(S)

    # near the true boundaries (true: 50..349 and 1350..1649)
    assert abs(l5b - 50) <= 3
    assert abs(l5e - 349) <= 3
    assert abs(l3b - 1350) <= 3
    assert abs(l3e - 1649) <= 3

    # the two LTRs must come out the same length as each other
    assert abs((l5e - l5b) - (l3e - l3b)) <= 3

def test_core_alignment_reconstructs_both_ltrs():
    ltr = _rnd(300, 5)
    S = ltr + _rnd(900, 6) + _evolve(ltr, 0.10, 7)
    h = discover(S, GENERIC_MATRIX)
    a, b = core_alignment(S, h)
    assert len(a) == len(b)
    l5b, l5e, l3b, l3e = ltr_spans(S, h)
    assert a.replace("-", "") == S[l5b:l5e+1]
    assert b.replace("-", "") == S[l3b:l3e+1]

def test_calibrated_matrix_differs_from_generic_at_high_divergence():
    ltr = _rnd(600, 8)
    S = ltr + _rnd(900, 9) + _evolve(ltr, 0.25, 10)
    h = discover(S, GENERIC_MATRIX)
    m, d_hat, kappa = calibrate(S, h)
    assert 0.1 < d_hat < 0.7
    arr = np.array(m.matrix)
    # transition A>G (index 0,2) must be penalised less than transversion A>C (0,1)
    assert arr[0][2] > arr[0][1]

def test_calibration_falls_back_to_generic_on_tiny_core():
    """Fewer than MIN_CALIB_SITES ungapped sites -> keep the generic matrix
    rather than derive a degenerate one from noise."""
    ltr = _rnd(30, 11)
    S = ltr + _rnd(200, 12) + ltr
    h = discover(S, GENERIC_MATRIX)
    if h is not None:
        m, _, _ = calibrate(S, h)
        n_sites = len(core_alignment(S, h)[0].replace("-", ""))
        if n_sites < MIN_CALIB_SITES:
            assert m is GENERIC_MATRIX

def test_calibrated_matrix_scores_N_zero():
    ltr = _rnd(400, 13)
    S = ltr + _rnd(900, 14) + _evolve(ltr, 0.15, 15)
    h = discover(S, GENERIC_MATRIX)
    m, _, _ = calibrate(S, h)
    arr = np.array(m.matrix)
    assert list(arr[4][:5]) == [0, 0, 0, 0, 0]
    assert [arr[i][4] for i in range(5)] == [0, 0, 0, 0, 0]

def test_ltr_spans_recovers_exact_ltr_sequence_unflanked():
    """Ground truth: the element is built from a known LTR, so the sliced spans
    must reproduce it byte-for-byte. Unflanked, so there is no boundary fraying
    and the correct answer is unambiguous. This is the test that catches an
    off-by-one; the flanked fixture cannot, because drift into the flank can
    coincidentally cancel it."""
    for ltr_len in (200, 400, 1200):
        for seed in range(6):
            ltr = _rnd(ltr_len, seed)
            S = ltr + _rnd(1500, seed + 50) + ltr
            h = discover(S, GENERIC_MATRIX)
            assert h is not None
            l5b, l5e, l3b, l3e = ltr_spans(S, h)
            assert S[l5b:l5e + 1] == ltr, f"5' LTR wrong at len={ltr_len} seed={seed}"
            assert S[l3b:l3e + 1] == ltr, f"3' LTR wrong at len={ltr_len} seed={seed}"
            assert l5e - l5b + 1 == ltr_len and l3e - l3b + 1 == ltr_len
